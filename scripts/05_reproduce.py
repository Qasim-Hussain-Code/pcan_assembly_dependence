"""Arm 0: run PCAn v1.0 on every arm 0 assembly, compare the S288C calls with
SGD's centromere annotation, and run S288C five times to check that the output
is deterministic.

Usage:
    python scripts/05_reproduce.py [--jobs N] [--force]

Run inside the tools environment (run_all.sh does this). PCAn itself runs only
through scripts/04_run_pcan.sh.
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import padlib as P  # noqa: E402

S288C = "GCA_000146045.2"
OUT = P.repo("results", "arm0")


def load_inputs():
    species = pd.read_csv(P.repo("config", "species_arm0.tsv"), sep="\t", dtype=str)
    manifest = pd.read_csv(P.repo("logs", "assemblies_arm0.tsv"), sep="\t", dtype=str)
    return species, manifest


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
# S288C against SGD

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
    species, manifest = load_inputs()
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

    # Every call, and the count per species.
    allcalls = []
    for _, s in todo.iterrows():
        acc = s["accession"]
        path = run_dir("calls", acc + ".tsv")
        if not os.path.exists(path):
            continue
        calls = pd.read_csv(path, sep="\t")
        calls.insert(0, "species", s["species"])
        calls.insert(1, "accession", acc)
        allcalls.append(calls)
    allcalls = pd.concat(allcalls, ignore_index=True)
    P.atomic_write_tsv(allcalls.drop(columns=["label"], errors="ignore"), os.path.join(OUT, "calls.tsv"))
    spc = []
    for _, s in species[species["group"] == "Saccharomycetaceae"].iterrows():
        acc = s["accession"]
        r = {"species": s["species"], "accession": acc, "arm0_role": s["arm0_role"], "reason": s["reason"],
             "pcan_genus": s["pcan_genus"]}
        if s["arm0_role"] != "excluded" and acc not in have:
            r["arm0_role"] = "excluded"
            ex = pd.read_csv(P.repo("results", "arm0", "excluded_assemblies.tsv"), sep="\t", dtype=str)
            why = ex.loc[ex["accession"] == acc, "reason"]
            r["reason"] = "assembly not used: " + (why.iloc[0] if len(why) else "not downloaded")
        m = runs[runs["accession"] == acc]
        if len(m):
            r["calls"] = m.iloc[0]["n_calls"]
            r["pcan_status"] = m.iloc[0]["status"]
        spc.append(r)
    P.atomic_write_tsv(pd.DataFrame(spc), os.path.join(OUT, "species_counts.tsv"))
    P.log("arm 0: %d calls on %d assemblies" % (len(allcalls), allcalls["accession"].nunique()))

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
