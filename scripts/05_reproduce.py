"""Arm 0: run PCAn v1.0 on every arm 0 assembly and compare its calls with
Supplementary Data 2; S288C against SGD; the determinism check; the gate.

The rules applied here were committed beforehand in config/arm0_reproduction.md.

Usage:
    python scripts/05_reproduce.py [--jobs N] [--force]

Run inside the tools environment (run_all.sh does this). PCAn itself runs only
through scripts/04_run_pcan.sh.
"""
import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import padlib as P  # noqa: E402

SEED = 20261009
GATE = 0.95
GATE_ALTERNATIVES = [0.90, 0.99]
N_BOOT = 10000
S288C = "GCA_000146045.2"
OUT = P.repo("results", "arm0")


def load_inputs():
    species = pd.read_csv(P.repo("config", "species_arm0.tsv"), sep="\t", dtype=str)
    manifest = pd.read_csv(P.repo("logs", "assemblies_arm0.tsv"), sep="\t", dtype=str)
    stats = pd.read_csv(P.repo("results", "arm0", "assembly_stats.tsv"), sep="\t")
    sd2 = pd.read_csv(P.data_dir("tables", "supp_data2.tsv"), sep="\t", dtype=str)
    return species, manifest, stats, sd2


def run_dir(*parts):
    return P.data_dir("arm0", *parts)


def run_one(row):
    acc = row["accession"]
    out = run_dir("calls", acc + ".tsv")
    meta = P.run_pcan(P.data_dir("assemblies", acc + ".fna.gz"), row["pcan_genus"], out,
                      species=row["pcan_kazachstania_species"] if isinstance(row["pcan_kazachstania_species"], str) else "",
                      label=acc, keep_fimo=run_dir("fimo", acc))
    meta["accession"] = acc
    meta["species"] = row["species"]
    return meta


def read_meta(path):
    meta = {}
    with open(path) as fh:
        for line in fh:
            if "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                meta[k] = v
    return meta


# ----------------------------------------------------------------------------
# locating published calls

def name_maps(acc):
    """refseq name -> genbank name, from the NCBI sequence report."""
    m = {}
    path = P.data_dir("assemblies", acc + ".sequence_report.jsonl")
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("refseqAccession") and r.get("genbankAccession"):
                m[r["refseqAccession"]] = r["genbankAccession"]
    return m


def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i >= 0:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def locate(seq_upper, query):
    """All placements of query on either strand: [(start, end, strand)], 1-based."""
    q = query.upper()
    hits = [(i + 1, i + len(q), "+") for i in find_all(seq_upper, q)]
    rc = P.revcomp(q)
    if rc != q:
        hits += [(i + 1, i + len(q), "-") for i in find_all(seq_upper, rc)]
    return hits


def overlaps(a1, b1, a2, b2):
    return a1 <= b2 and a2 <= b1


# ----------------------------------------------------------------------------
# reasons from the filter trace

def reason_for(pub, cands, pass1, motif_only=False):
    """Why the published call is not reproduced exactly, from this run's
    candidate table and first FIMO pass. With motif_only, the published
    interval is the CDEIII motif itself."""
    if pub["start"] is None:
        return ""
    on = cands[cands["Contig"] == pub["contig"]]
    if not motif_only:
        same = on[(on["Start"] == pub["start"]) & (on["End"] == pub["end"])]
        if len(same):
            step = same.iloc[0]["removed_at"]
            return "published candidate formed (rank %d), removed at %s" % (same.iloc[0]["rank"], step) \
                if step != "called" else "published candidate formed and called"
    # the CDEIII motif of the published call, in forward coordinates
    if motif_only:
        m3 = (pub["start"], pub["end"])
    elif pub["strand"] == "+":
        m3 = (pub["end"] - 25, pub["end"])
    else:
        m3 = (pub["start"], pub["start"] + 25)
    p1 = pass1[(pass1["sequence name"] == pub["contig"]) & (pass1["start"] <= m3[1]) & (pass1["stop"] >= m3[0])]
    exact_hit = p1[(p1["start"] == m3[0]) & (p1["stop"] == m3[1])]
    if not len(p1):
        return "no CDEIII hit at the published motif in FIMO's first pass at this threshold"
    if not len(exact_hit):
        return "first-pass CDEIII hit overlaps but is offset from the published motif"
    near = on[(on["Start"] <= pub["end"]) & (on["End"] >= pub["start"])]
    if not len(near):
        return "CDEIII hit present; no CDEI hit in its window"
    best = near.sort_values("rank").iloc[0]
    return "CDEIII hit present; published CDEI not among candidates; best candidate here CDEII %d, removed at %s" % (
        best["CDEIIlen"], best["removed_at"])


def compare_species(srow, sd2, calls):
    """Status of every published call of one species, and the extra calls."""
    acc = srow["accession"]
    pub = sd2[sd2["Genbank_Assembly"] == srow["sd2_accession"]].copy()
    recs = dict(P.read_fasta(P.data_dir("assemblies", acc + ".fna.gz")))
    upper = {}
    rs2gb = name_maps(acc)
    cand_path = run_dir("calls", acc + ".tsv.candidates.tsv.gz")
    cands = pd.read_csv(cand_path, sep="\t") if os.path.exists(cand_path) else pd.DataFrame(
        columns=["Contig", "Start", "End", "rank", "removed_at", "CDEIIlen"])
    p1_path = run_dir("fimo", acc, "pass1.fimo.txt.gz")
    pass1 = pd.read_csv(p1_path, sep="\t") if os.path.exists(p1_path) else pd.DataFrame(
        columns=["sequence name", "start", "stop"])
    if len(pass1):
        pass1 = pass1.rename(columns={"#pattern name": "pattern"})
    rows = []
    matched_calls = set()
    cdeiii_only = srow["arm0_role"] == "cdeiii_only"
    for _, p in pub.iterrows():
        contig_pub = str(p["Contig/scaffold/chromosome"]).strip()
        contig = rs2gb.get(contig_pub, contig_pub)
        r = {"species": srow["species"], "accession": acc, "published_contig": contig_pub, "contig": contig,
             "chr_no": p["ChrNo"], "published_cdeii_len": p["CDEII_Length"], "start": None, "end": None,
             "strand": "", "status": "", "reproduced_start": "", "reproduced_end": "",
             "reproduced_cdeii_len": "", "reason": ""}
        if contig not in recs:
            r.update(status="unlocated", reason="contig not in the downloaded assembly")
            rows.append(r)
            continue
        if contig not in upper:
            upper[contig] = recs[contig].upper()
        query = str(p["CDEIII"]) if cdeiii_only else str(p["CEN"])
        hits = locate(upper[contig], query)
        if len(hits) != 1:
            r.update(status="unlocated",
                     reason="published sequence found %d times on the contig" % len(hits))
            rows.append(r)
            continue
        s, e, strand = hits[0]
        r.update(start=s, end=e, strand=strand)
        on = calls[calls["contig"] == contig]
        if cdeiii_only:
            same = on[(on["cdeiii_start"] == s) & (on["cdeiii_end"] == e)]
            ov = on[(on["start"] <= e) & (on["end"] >= s)]
            if len(same):
                c = same.iloc[0]
                r.update(status="same_cdeiii")
            elif len(ov):
                c = ov.iloc[0]
                r.update(status="overlapping_other_cdeiii")
            else:
                c = None
                r.update(status="missing")
        else:
            exact = on[(on["start"] == s) & (on["end"] == e) & (on["sequence"].str.upper() == query.upper())]
            ov = on[(on["start"] <= e) & (on["end"] >= s)]
            if len(exact):
                c = exact.iloc[0]
                r.update(status="exact")
            elif len(ov):
                c = ov.iloc[0]
                same_len = str(int(c["cdeii_len"])) == str(int(float(p["CDEII_Length"])))
                r.update(status="same_locus_other_boundaries" if same_len else "different_cdeii_length")
            else:
                c = None
                r.update(status="missing")
        if c is not None:
            matched_calls.add((c["contig"], int(c["start"]), int(c["end"])))
            r.update(reproduced_start=int(c["start"]), reproduced_end=int(c["end"]),
                     reproduced_cdeii_len=int(c["cdeii_len"]))
        if r["status"] not in ("exact", "same_cdeiii"):
            r["reason"] = reason_for(r, cands, pass1, motif_only=cdeiii_only)
        rows.append(r)
    extra = calls[[(c.contig, int(c.start), int(c.end)) not in matched_calls for c in calls.itertuples()]].copy()
    extra.insert(0, "species", srow["species"])
    extra.insert(1, "accession", acc)
    return pd.DataFrame(rows), extra


def failed_species(srow, sd2, error):
    """A species whose PCAn run failed: every published call is missing, and
    the reason is the failure itself."""
    pub = sd2[sd2["Genbank_Assembly"] == srow["sd2_accession"]]
    return pd.DataFrame([{"species": srow["species"], "accession": srow["accession"],
                          "published_contig": p["Contig/scaffold/chromosome"], "contig": "",
                          "chr_no": p["ChrNo"], "published_cdeii_len": p["CDEII_Length"], "start": None,
                          "end": None, "strand": "", "status": "missing", "reproduced_start": "",
                          "reproduced_end": "", "reproduced_cdeii_len": "",
                          "reason": "PCAn run failed: " + str(error)[:300]} for _, p in pub.iterrows()])


# ----------------------------------------------------------------------------
# why a published call is not reproduced exactly (written after the gate
# result, to find its cause; nothing here changes the gate)

def occurrences(recs, query):
    q = query.upper()
    rc = P.revcomp(q)
    hits = []
    for name, s in recs:
        S = s.upper()
        for strand, qq in (("+", q), ("-", rc)):
            i = S.find(qq)
            while i >= 0:
                hits.append((name, i + 1, strand))
                i = S.find(qq, i + 1)
    return hits


def longest_prefix(recs, query):
    q = query.upper()
    best = 0
    for _, s in recs:
        S = s.upper()
        for qq in (q, P.revcomp(q)):
            lo, hi = 0, len(qq)
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if qq[:mid] in S:
                    lo = mid
                else:
                    hi = mid - 1
            best = max(best, lo)
    return best


def cause_of(row, recs, pub_seq, diag_accs, anchor_far):
    st, reason = row["status"], str(row["reason"])
    if st == "unlocated":
        if row["accession"] in diag_accs:
            return "calls published on a replaced assembly version", "Supplementary Data 2 contig names belong to %s" % diag_accs[row["accession"]]
        hits = occurrences(recs, pub_seq)
        if len(hits) == 1:
            return "contig name in Supplementary Data 2 differs from the assembly", "sequence occurs once, on %s" % hits[0][0]
        if len(hits) > 1:
            return "published sequence repeated", "occurs %d times in the assembly" % len(hits)
        return "published sequence differs from the assembly", "longest matching prefix %d of %d bases" % (
            longest_prefix(recs, pub_seq), len(pub_seq))
    if anchor_far and "cdeii_length_outside_anchor_window" in reason:
        return "length anchor set by false candidates", reason
    if reason.startswith("no CDEIII hit"):
        return "CDEIII motif below the released threshold", reason
    if "no CDEI hit in its window" in reason or "published CDEI not among candidates" in reason:
        return "published CDEI not found at the released CDEI threshold", reason
    if "removed at" in reason:
        return "candidate formed, removed by " + reason.rsplit("removed at ", 1)[1], reason
    if reason.startswith("first-pass CDEIII hit overlaps"):
        return "CDEIII motif offset from the published one", reason
    return "other", reason


def explain(per, sd2, species, runs, spc):
    diag = pd.read_csv(P.repo("config", "arm0_diagnostics.tsv"), sep="\t", dtype=str)
    diag_accs = dict(zip(diag["accession_listed"], diag["accession_run"]))
    prim = set(spc.loc[spc["arm0_role"] == "primary", "accession"])
    sub = per[per["accession"].isin(prim) & (per["status"] != "exact")]
    rows = []
    for acc, d in sub.groupby("accession"):
        recs = P.read_fasta(P.data_dir("assemblies", acc + ".fna.gz"))
        s = species[species["accession"] == acc].iloc[0]
        pub = sd2[sd2["Genbank_Assembly"] == s["sd2_accession"]]
        meta = runs[runs["accession"] == acc]
        top5 = float(meta.iloc[0]["median_top5_cdeii_len"]) if len(meta) else float("nan")
        pub_med = pd.to_numeric(pub["CDEII_Length"], errors="coerce").median()
        anchor_far = abs(top5 - pub_med) > 30 if not np.isnan(top5) else False
        for _, r in d.iterrows():
            p = pub[(pub["Contig/scaffold/chromosome"] == r["published_contig"]) & (pub["ChrNo"].astype(str) == str(r["chr_no"]))]
            seq = str(p.iloc[0]["CEN"]) if len(p) else ""
            cause, detail = cause_of(r, recs, seq, diag_accs, anchor_far)
            rows.append({"species": r["species"], "accession": acc, "published_contig": r["published_contig"],
                         "chr_no": r["chr_no"], "status": r["status"], "cause": cause, "detail": detail,
                         "median_top5_cdeii_len": top5, "published_median_cdeii_len": pub_med})
    causes = pd.DataFrame(rows)
    counts = causes.groupby(["cause", "status"]).size().rename("published_calls").reset_index()
    counts = counts.sort_values("published_calls", ascending=False)

    # Against the authors' own expectations. ExpectedFalsePosFalseNegs.txt,
    # shipped with PCAn, gives per species the false negatives and false
    # positives the authors expect from PCAn v1.0 at these settings.
    ex = spc[spc["arm0_role"] == "primary"].merge(
        species[["accession", "authors_false_neg", "authors_false_pos"]], on="accession", how="left")
    ex["authors_false_neg"] = pd.to_numeric(ex["authors_false_neg"], errors="coerce")
    ex["authors_false_pos"] = pd.to_numeric(ex["authors_false_pos"], errors="coerce")
    ex["missing_equals_authors_false_neg"] = ex["missing"] == ex["authors_false_neg"]
    ex["extra_equals_authors_false_pos"] = ex["extra"] == ex["authors_false_pos"]
    ex = ex[["species", "accession", "published_calls", "reproduced_calls", "exact", "different_cdeii_length",
             "same_locus_other_boundaries", "missing", "unlocated", "extra", "authors_false_neg",
             "authors_false_pos", "missing_equals_authors_false_neg", "extra_equals_authors_false_pos"]]
    return causes, counts, ex


def diagnostics(sd2):
    """Run PCAn on assembly versions listed in config/arm0_diagnostics.tsv and
    compare with the published calls. Outside the reproduction fraction."""
    import subprocess
    diag = pd.read_csv(P.repo("config", "arm0_diagnostics.tsv"), sep="\t", dtype=str)
    species = pd.read_csv(P.repo("config", "species_arm0.tsv"), sep="\t", dtype=str)
    rows = []
    for d in diag.itertuples():
        out_dir = P.data_dir("diagnostics", d.accession_run)
        os.makedirs(out_dir, exist_ok=True)
        fna = os.path.join(out_dir, d.accession_run + ".fna.gz")
        if not os.path.exists(fna):
            work = os.path.join(out_dir, "pkg")
            subprocess.run(["datasets", "download", "genome", "accession", d.accession_run, "--include", "genome",
                            "--filename", work + ".zip", "--no-progressbar"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            subprocess.run(["unzip", "-q", "-o", work + ".zip", "-d", work], check=True)
            src = [os.path.join(dp, f) for dp, _, fs in os.walk(work) for f in fs if f.endswith("_genomic.fna")][0]
            subprocess.run("pigz -c '%s' > '%s'" % (src, fna), shell=True, check=True)
            subprocess.run(["rm", "-rf", work, work + ".zip"], check=True)
        s = species[species["species"] == d.species].iloc[0]
        out = os.path.join(out_dir, "calls.tsv")
        P.run_pcan(fna, s["pcan_genus"], out, label=d.accession_run)
        calls = pd.read_csv(out, sep="\t")
        recs = P.read_fasta(fna)
        pub = sd2[sd2["Species"] == d.species]
        exact = 0
        for _, p in pub.iterrows():
            hits = [h for h in occurrences(recs, p["CEN"]) if h[0] == p["Contig/scaffold/chromosome"]]
            if len(hits) == 1:
                name, s0, _ = hits[0]
                e0 = s0 + len(p["CEN"]) - 1
                exact += int(((calls["contig"] == name) & (calls["start"] == s0) & (calls["end"] == e0)).any())
        rows.append({"species": d.species, "accession_listed": d.accession_listed, "accession_run": d.accession_run,
                     "published_calls": len(pub), "reproduced_calls": len(calls), "exact": exact,
                     "extra": len(calls) - exact, "reason": d.reason})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# S288C against SGD

def cdeiii_motif_diagnostic(per, species):
    """The species published with CDEIII loci only. For the three Yueomyces
    species PCAn's own entry uses the Saccharomyces CDEIII motif, while the
    authors' table names a Yueomyces motif that ships with PCAn but no entry
    uses; for Grigorovia jiainica the entry's threshold is 1e-5 and the
    table's 1e-6. Each setting is run as PCAn's first pass runs FIMO (fimo
    --thresh 1.0E-<n> --oc on the whole assembly), and the published loci
    that a hit overlaps are counted. A diagnostic, outside the reproduction
    fraction: no call in any arm depends on it."""
    import gzip
    import shutil
    import subprocess
    fimo = P.data_dir("env", "pcan", "bin", "fimo")
    rows = []
    for s in species[species["arm0_role"] == "cdeiii_only"].to_dict("records"):
        acc = s["accession"]
        pub = per[(per["accession"] == acc)].dropna(subset=["start"])
        work = P.data_dir("work", "cdeiii_motif_" + acc)
        os.makedirs(work, exist_ok=True)
        try:
            genome = os.path.join(work, "genome.fna")
            with gzip.open(P.data_dir("assemblies", acc + ".fna.gz"), "rb") as src, open(genome, "wb") as dst:
                shutil.copyfileobj(src, dst)
            for setting, motif, thresh in (("PCAn entry", s["pcan_cdeiii_motif"], s["pcan_cdeiii_thresh"]),
                                           ("authors' table", s["authors_cdeiii_motif"], s["authors_cdeiii_thresh"])):
                motif_file = os.path.join(P.pcan_dir(), "CDEIII", "CDEIII_%s_MEME.txt" % motif)
                oc = os.path.join(work, "fimo")
                subprocess.run([fimo, "--thresh", "1.0E-%d" % round(-math.log10(float(thresh))), "--oc", oc,
                                motif_file, genome], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                hits = pd.read_csv(os.path.join(oc, "fimo.txt"), sep="\t")
                hit_loci = sum(bool(((hits["sequence name"] == r.contig) & (hits["start"] <= r.end)
                                     & (hits["stop"] >= r.start)).any()) for r in pub.itertuples())
                rows.append({"species": s["species"], "accession": acc, "setting": setting,
                             "cdeiii_motif": os.path.basename(motif_file), "threshold": thresh,
                             "published_loci": len(per[per["accession"] == acc]), "published_loci_located": len(pub),
                             "located_loci_with_a_hit": hit_loci, "fimo_hits": len(hits)})
        finally:
            shutil.rmtree(work, ignore_errors=True)
    return pd.DataFrame(rows)


def sgd_centromeres():
    rows = []
    with open(P.data_dir("tables", "sgd_R64-5-1.gff")) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] in ("centromere", "centromere_DNA_Element_I", "centromere_DNA_Element_II",
                        "centromere_DNA_Element_III"):
                attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                # SGD writes each centromere twice: the CEN feature and a child
                # of the same type and span (ID=CEN1_centromere;Parent=CEN1).
                if f[2] == "centromere" and "Parent" in attrs:
                    continue
                rows.append({"chrom": f[0], "type": f[2], "start": int(f[3]), "end": int(f[4]),
                             "strand": f[6], "name": attrs.get("Name", attrs.get("ID", ""))})
    return pd.DataFrame(rows)


def roman_to_int(s):
    vals = {"I": 1, "V": 5, "X": 10}
    total, prev = 0, 0
    for ch in reversed(s):
        v = vals[ch]
        total = total - v if v < prev else total + v
        prev = max(prev, v)
    return total


def s288c_checks(calls):
    """Sequence identity of NCBI's and SGD's S288C, then the 16 centromeres."""
    ncbi = dict(P.read_fasta(P.data_dir("assemblies", S288C + ".fna.gz")))
    report = {}
    with open(P.data_dir("assemblies", S288C + ".sequence_report.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            report[r["genbankAccession"]] = r
    chrom_to_gb = {}
    for gb, r in report.items():
        if r.get("assignedMoleculeLocationType") == "Chromosome":
            chrom_to_gb["chr" + r["chrName"]] = gb
        elif r.get("assignedMoleculeLocationType") == "Mitochondrion":
            chrom_to_gb["chrmt"] = gb
    # SGD's reference FASTA names chromosomes in its headers as [chromosome=I]
    sgd = {}
    name = None
    with open(P.data_dir("tables", "sgd_R64-5-1_reference.fsa")) as fh:
        chunks = []
        for line in fh:
            if line.startswith(">"):
                if name:
                    sgd[name] = "".join(chunks)
                chunks = []
                if "[chromosome=" in line:
                    name = "chr" + line.split("[chromosome=")[1].split("]")[0]
                elif "[location=mitochondrion]" in line:
                    name = "chrmt"
                else:
                    name = line[1:].split()[0]
            else:
                chunks.append(line.strip())
        if name:
            sgd[name] = "".join(chunks)
    seqcheck = []
    for chrom, gb in sorted(chrom_to_gb.items(), key=lambda kv: kv[0]):
        a, b = ncbi.get(gb, ""), sgd.get(chrom, "")
        seqcheck.append({"chromosome": chrom, "genbank": gb, "ncbi_length": len(a), "sgd_length": len(b),
                         "identical": "yes" if a.upper() == b.upper() and a else "no"})
    seqcheck = pd.DataFrame(seqcheck)

    feats = sgd_centromeres()
    cens = feats[feats["type"] == "centromere"].copy()
    rows = []
    used = set()

    def element(kind, c):
        e = feats[(feats["type"] == kind) & (feats["chrom"] == c["chrom"]) & (feats["start"] >= c["start"])
                  & (feats["end"] <= c["end"])]
        return (int(e.iloc[0]["start"]), int(e.iloc[0]["end"])) if len(e) else None

    for _, c in cens.iterrows():
        gb = chrom_to_gb.get(c["chrom"])
        on = calls[calls["contig"] == gb]
        hit = on[(on["start"] <= c["end"]) & (on["end"] >= c["start"])]
        i1, i2, i3 = (element(k, c) for k in ("centromere_DNA_Element_I", "centromere_DNA_Element_II",
                                                "centromere_DNA_Element_III"))
        sgd_cdeii = i2[1] - i2[0] + 1 if i2 else None
        r = {"sgd_name": c["name"], "chromosome": c["chrom"], "contig": gb, "sgd_start": c["start"],
             "sgd_end": c["end"], "sgd_strand": c["strand"],
             "sgd_cdei_len": i1[1] - i1[0] + 1 if i1 else "", "sgd_cdeii_len": sgd_cdeii,
             "sgd_cdeiii_len": i3[1] - i3[0] + 1 if i3 else "", "called": "yes" if len(hit) else "no"}
        if len(hit):
            h = hit.iloc[0]
            used.add((h["contig"], int(h["start"])))
            # Where CDEII meets CDEIII, in forward coordinates: the first base
            # of CDEIII on the forward strand, the last base on the reverse.
            if h["strand"] == "+":
                boundary = int(h["cdeiii_start"]) - i3[0] if i3 else ""
            else:
                boundary = int(h["cdeiii_end"]) - i3[1] if i3 else ""
            r.update(pcan_start=int(h["start"]), pcan_end=int(h["end"]), pcan_strand=h["strand"],
                     pcan_cdeii_len=int(h["cdeii_len"]), start_offset=int(h["start"]) - c["start"],
                     end_offset=int(h["end"]) - c["end"],
                     cdeii_len_difference=(int(h["cdeii_len"]) - sgd_cdeii) if sgd_cdeii else "",
                     cdeii_cdeiii_boundary_offset=boundary)
        rows.append(r)
    vs = pd.DataFrame(rows)
    vs["sort"] = vs["chromosome"].str.replace("chr", "").map(lambda s: roman_to_int(s) if s != "mt" else 99)
    vs = vs.sort_values("sort").drop(columns="sort")
    extra = calls[[(c.contig, int(c.start)) not in used for c in calls.itertuples()]]
    return seqcheck, vs, extra, len(cens)


# ----------------------------------------------------------------------------
# determinism

def determinism(row):
    outs = []
    for i in range(1, 6):
        out = run_dir("determinism", "run%d.tsv" % i)
        P.run_pcan(P.data_dir("assemblies", S288C + ".fna.gz"), row["pcan_genus"], out, label=S288C,
                   keep_fimo=run_dir("determinism", "fimo%d" % i))
        outs.append(out)
    import gzip
    import hashlib

    def digest(path):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()

    rows = []
    for i, out in enumerate(outs, 1):
        rows.append({"run": i, "calls": sum(1 for _ in open(out)) - 1,
                     "call_table_sha256": digest(out),
                     "fimo_pass1_sha256": digest(run_dir("determinism", "fimo%d" % i, "pass1.fimo.txt.gz")),
                     "fimo_pass2_sha256": digest(run_dir("determinism", "fimo%d" % i, "pass2.fimo.txt.gz")),
                     "candidates_sha256": digest(out + ".candidates.tsv.gz"),
                     "elapsed_s": read_meta(out + ".meta.tsv").get("elapsed_s", "")})
    df = pd.DataFrame(rows)
    identical = all(df[c].nunique() == 1 for c in df.columns if c.endswith("sha256"))
    return df, identical


# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=int(P.conf().get("JOBS", 1)))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if P.is_done("reproduce") and not args.force:
        P.log("arm 0 already complete, skipping (use --force to redo)")
        return
    species, manifest, stats, sd2 = load_inputs()
    have = set(manifest["accession"])
    todo = species[(species["arm0_role"] != "excluded") & species["accession"].isin(have)]
    missing_asm = species[(species["arm0_role"] != "excluded") & ~species["accession"].isin(have)]
    for d in ("calls", "fimo", "determinism"):
        os.makedirs(run_dir(d), exist_ok=True)

    # PCAn on every assembly, measuring after the first ten.
    P.log("arm 0: %d assemblies, %d jobs" % (len(todo), args.jobs))
    metas = []
    t0 = time.time()
    rows = [r for _, r in todo.iterrows()]
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_one, r): r["accession"] for r in rows}
        for f in as_completed(futs):
            m = f.result()
            metas.append(m)
            if len(metas) == 10:
                per = (time.time() - t0) / 10 * args.jobs
                left = (len(rows) - 10) * per / args.jobs / 3600
                P.log("first ten assemblies: %.1f s per assembly per job; about %.2f h left for arm 0" % (per, left))
            if m.get("status") != "ok":
                P.log("%s: PCAn failed: %s" % (m["accession"], m.get("error", m.get("stderr_tail", ""))[:300]))
    runs = pd.DataFrame(metas)
    runs = runs.merge(todo[["accession", "arm0_role"]], on="accession", how="left")
    keep = ["accession", "species", "arm0_role", "genus", "status", "pcan_exit_status", "n_calls", "elapsed_s",
            "max_rss_mb", "fimo_pass1_hits", "fimo_pass2_hits", "fimo_max_stored_scores_warning",
            "trace_matches_pcan", "candidates", "median_top5_cdeii_len", "cdeiii_motif", "cdeiii_threshold",
            "cdei_motif", "cdei_threshold", "error"]
    runs = runs.reindex(columns=keep).sort_values("accession")
    P.atomic_write_tsv(runs, P.repo("logs", "arm0_runs.tsv"))

    # Compare with Supplementary Data 2.
    per, extras, allcalls = [], [], []
    for _, s in todo.iterrows():
        acc = s["accession"]
        path = run_dir("calls", acc + ".tsv")
        if not os.path.exists(path):
            m = runs[runs["accession"] == acc]
            err = m.iloc[0]["error"] if len(m) else "no output"
            if s["arm0_role"] in ("primary", "cdeiii_only"):
                per.append(failed_species(s, sd2, err))
            continue
        calls = pd.read_csv(path, sep="\t")
        c2 = calls.copy()
        c2.insert(0, "species", s["species"])
        allcalls.append(c2)
        if s["arm0_role"] in ("primary", "cdeiii_only"):
            pc, ex = compare_species(s, sd2, calls)
            per.append(pc)
            extras.append(ex)
        else:
            ex = calls.copy()
            ex.insert(0, "species", s["species"])
            ex.insert(1, "accession", acc)
            extras.append(ex)
    per = pd.concat(per, ignore_index=True)
    extras = pd.concat(extras, ignore_index=True)
    allcalls = pd.concat(allcalls, ignore_index=True)
    P.atomic_write_tsv(per, os.path.join(OUT, "per_centromere.tsv"))
    P.atomic_write_tsv(extras.drop(columns=["label"], errors="ignore"), os.path.join(OUT, "extra_calls.tsv"))
    P.atomic_write_tsv(allcalls.drop(columns=["label"], errors="ignore"), os.path.join(OUT, "reproduced_calls.tsv"))

    # Per species.
    spc = []
    for _, s in species[species["group"] == "Saccharomycetaceae"].iterrows():
        acc = s["accession"]
        r = {"species": s["species"], "accession": acc, "arm0_role": s["arm0_role"], "reason": s["reason"],
             "pcan_genus": s["pcan_genus"], "published_calls": s["sd2_calls"]}
        if s["arm0_role"] != "excluded" and acc not in have:
            r["arm0_role"] = "excluded"
            ex = pd.read_csv(P.repo("results", "arm0", "excluded_assemblies.tsv"), sep="\t", dtype=str)
            why = ex.loc[ex["accession"] == acc, "reason"]
            r["reason"] = "assembly not used: " + (why.iloc[0] if len(why) else "not downloaded")
        m = runs[runs["accession"] == acc]
        if len(m):
            r["reproduced_calls"] = m.iloc[0]["n_calls"]
            r["pcan_status"] = m.iloc[0]["status"]
        pc = per[per["accession"] == acc] if len(per) else per
        for st in ["exact", "same_locus_other_boundaries", "different_cdeii_length", "missing", "unlocated",
                   "same_cdeiii", "overlapping_other_cdeiii"]:
            r[st] = int((pc["status"] == st).sum()) if len(pc) else 0
        r["extra"] = int((extras["accession"] == acc).sum()) if len(extras) else 0
        spc.append(r)
    spc = pd.DataFrame(spc)
    P.atomic_write_tsv(spc, os.path.join(OUT, "species_counts.tsv"))

    # The gate: exact fraction of published calls in the reproduction set,
    # bootstrapped over species.
    prim = per[per["accession"].isin(spc.loc[spc["arm0_role"] == "primary", "accession"])]
    units = [g["status"].eq("exact").values for _, g in prim.groupby("accession")]

    def frac(us):
        tot = sum(len(u) for u in us)
        return sum(u.sum() for u in us) / tot if tot else float("nan")

    est, lo, hi = P.cluster_bootstrap(units, frac, n_boot=N_BOOT, seed=SEED)
    n_sp = len(units)
    sp_all_exact = sum(bool(u.all()) for u in units)
    wl, wh = P.wilson(sp_all_exact, n_sp)
    gate = pd.DataFrame([{
        "measure": "fraction of published calls reproduced exactly", "unit_of_replication": "species",
        "n_species": n_sp, "n_calls": int(sum(len(u) for u in units)), "n_exact": int(sum(u.sum() for u in units)),
        "estimate": round(est, 4), "ci_low": round(lo, 4), "ci_high": round(hi, 4),
        "interval": "percentile bootstrap over species, %d resamples, seed %d" % (N_BOOT, SEED),
        "threshold": GATE, "decision": "proceed" if est >= GATE else "stop",
        "species_all_exact": sp_all_exact, "species_all_exact_fraction": round(sp_all_exact / n_sp, 4) if n_sp else "",
        "species_all_exact_wilson_low": round(wl, 4), "species_all_exact_wilson_high": round(wh, 4)}])
    P.atomic_write_tsv(gate, os.path.join(OUT, "gate.tsv"))
    sens = pd.DataFrame([{"threshold": t, "estimate": round(est, 4), "ci_low": round(lo, 4),
                          "decision": "proceed" if est >= t else "stop",
                          "same_as_preregistered_threshold": ("proceed" if est >= t else "stop") == gate.iloc[0]["decision"]}
                         for t in [GATE] + GATE_ALTERNATIVES])
    P.atomic_write_tsv(sens, P.repo("results", "sensitivity", "arm0_gate_threshold.tsv"))
    P.log("gate: %.4f exact (%.4f to %.4f) over %d species; %s" % (est, lo, hi, n_sp, gate.iloc[0]["decision"]))

    # The cause of every non-exact call, the comparison with the authors' own
    # expected false negatives and positives, and the replaced-version
    # diagnostics. Written after the gate result; none of it changes the gate.
    causes, cause_counts, expected = explain(per, sd2, species, runs, spc)
    P.atomic_write_tsv(causes, os.path.join(OUT, "non_exact_causes.tsv"))
    P.atomic_write_tsv(cause_counts, os.path.join(OUT, "non_exact_cause_counts.tsv"))
    P.atomic_write_tsv(expected, os.path.join(OUT, "authors_expected_vs_observed.tsv"))
    diag = diagnostics(sd2)
    P.atomic_write_tsv(diag, os.path.join(OUT, "diagnostic_replaced_versions.tsv"))
    P.atomic_write_tsv(cdeiii_motif_diagnostic(per, species), os.path.join(OUT, "diagnostic_cdeiii_motif.tsv"))
    # Software version: would the AT formula of the published table have
    # changed which calls PCAn makes? Run in PCAn's environment, see the script.
    import subprocess
    prim_per = per[per["accession"].isin(spc.loc[spc["arm0_role"] == "primary", "accession"])]
    tmp = os.path.join(P.data_dir("work"), "per_primary.tsv")
    prim_per.to_csv(tmp, sep="\t", index=False)
    subprocess.run([P.data_dir("env", "pcan", "bin", "python"), P.repo("scripts", "lib", "at_formula_check.py"), tmp,
                    run_dir("calls"), os.path.join(OUT, "diagnostic_at_formula.tsv")], check=True)
    os.remove(tmp)
    P.log("authors' expected false negatives match the missing count in %d of %d species"
          % (expected["missing_equals_authors_false_neg"].sum(), len(expected)))

    # FIMO's stored-score cap.
    cap = runs[["accession", "species", "fimo_pass1_hits", "fimo_pass2_hits", "fimo_max_stored_scores_warning"]].copy()
    cap["max_stored_scores"] = 100000
    P.atomic_write_tsv(cap, os.path.join(OUT, "fimo_cap.tsv"))

    # Positive control.
    s288c_calls = pd.read_csv(run_dir("calls", S288C + ".tsv"), sep="\t")
    seqcheck, vs, s_extra, n_cen = s288c_checks(s288c_calls)
    P.atomic_write_tsv(seqcheck, os.path.join(OUT, "s288c_sequence_check.tsv"))
    P.atomic_write_tsv(vs, os.path.join(OUT, "s288c_vs_sgd.tsv"))
    P.log("S288C: SGD lists %d centromeres; PCAn called %d of them; %d other calls"
          % (n_cen, (vs["called"] == "yes").sum(), len(s_extra)))

    # Negative control.
    det, identical = determinism(todo[todo["accession"] == S288C].iloc[0])
    det["identical_across_runs"] = "yes" if identical else "no"
    P.atomic_write_tsv(det, os.path.join(OUT, "determinism.tsv"))
    P.log("determinism: five runs %s" % ("identical" if identical else "DIFFER"))

    if len(missing_asm):
        P.log("%d arm 0 assemblies were not available; see results/arm0/excluded_assemblies.tsv" % len(missing_asm))
    P.mark_done("reproduce")


if __name__ == "__main__":
    main()
