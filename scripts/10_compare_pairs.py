"""Arms 2 and 3: compare each strain's long-read assembly with another assembly
of the same strain, by the method registered in config/analysis_plan.md.

  --arm 2   the published short-read assembly of each pair (config/arm2_pairs.tsv)
  --arm 3   the assemblies built by scripts/11_assemble_reads.sh at known depths

For every strain: PCAn on the long-read assembly (once, shared by both arms);
the long-read centromeres and their AT-matched null windows; both flanks of
every element placed in the other assembly with minimap2; the five-way status;
the long-read assembly compared with itself as a positive control for the
liftover; and the contig ends of the other assembly located on the long-read
assembly. Each sensitivity variant repeats the placement with its own flank
length, AT bin or identity threshold.

Usage:
    python scripts/10_compare_pairs.py --arm 2|3 [--jobs N] [--force]
"""
import argparse
import os
import re
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import liftover as L  # noqa: E402
import padlib as P  # noqa: E402
import pairstats as PS  # noqa: E402
import status as S  # noqa: E402

GENUS = "Saccharomyces"
VARIANTS = [
    {"name": "primary", "flank": 2000, "half_bin": 0.025, "identity": 0.95},
    {"name": "flank_1kb", "flank": 1000, "half_bin": 0.025, "identity": 0.95},
    {"name": "flank_5kb", "flank": 5000, "half_bin": 0.025, "identity": 0.95},
    {"name": "at_bin_2.5", "flank": 2000, "half_bin": 0.0125, "identity": 0.95},
    {"name": "at_bin_10", "flank": 2000, "half_bin": 0.05, "identity": 0.95},
    {"name": "identity_0.90", "flank": 2000, "half_bin": 0.025, "identity": 0.90},
    {"name": "identity_0.99", "flank": 2000, "half_bin": 0.025, "identity": 0.99},
]
ELEMENT_SETS = sorted({(v["flank"], v["half_bin"]) for v in VARIANTS})


def threads():
    return max(1, int(P.conf().get("THREADS", 1)) // max(1, int(os.environ.get("PAD_JOBS", 1))))


# ----------------------------------------------------------------------------
# one strain's long-read side, shared by arms 2 and 3

def long_side(code, acc):
    """PCAn calls on the long-read assembly, and the element sets (centromeres
    and nulls) for every flank length and AT bin, cached on disk."""
    d = P.data_dir("pairs", code)
    os.makedirs(d, exist_ok=True)
    fa = P.data_dir("assemblies", acc + ".fna.gz")
    calls_path = os.path.join(d, "long_calls.tsv")
    meta = P.run_pcan(fa, GENUS, calls_path, label=code + "_long")
    if meta.get("status") != "ok":
        raise RuntimeError("PCAn failed on the long-read assembly of %s: %s" % (code, meta.get("error", "")))
    calls = pd.read_csv(calls_path, sep="\t")
    seqs = None
    for flank, half in ELEMENT_SETS:
        tag = "f%d_b%g" % (flank, half)
        cpath, npath = os.path.join(d, "cens_%s.tsv" % tag), os.path.join(d, "nulls_%s.tsv" % tag)
        if os.path.exists(cpath) and os.path.exists(npath):
            continue
        if seqs is None:
            seqs = dict(P.read_fasta(fa))
            idx = L.AtIndex(seqs)
        cens = L.centromere_elements(calls, seqs, flank)
        nulls, short = L.null_windows(code, cens, idx, calls, flank, half)
        cens.to_csv(cpath, sep="\t", index=False)
        nulls.to_csv(npath, sep="\t", index=False)
        short.to_csv(os.path.join(d, "null_shortfall_%s.tsv" % tag), sep="\t", index=False)
    return calls, meta


def read_elements(path):
    """An element or status table. Its kind column holds the word "null",
    which pandas reads as a missing value unless told otherwise. Kind is also
    fixed by the element_id prefix ("cen|" or "null|"), which is what sets it
    here, so a table written with the value already lost reads correctly."""
    t = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])
    t["kind"] = np.where(t["element_id"].astype(str).str.startswith("null|"), "null", "centromere")
    return t


def elements(code, flank, half):
    d = P.data_dir("pairs", code)
    tag = "f%d_b%g" % (flank, half)
    return (read_elements(os.path.join(d, "cens_%s.tsv" % tag)),
            read_elements(os.path.join(d, "nulls_%s.tsv" % tag)))


# ----------------------------------------------------------------------------
# one comparison: the long-read assembly against one other assembly

def compare(code, acc, query_fa, query_calls, workdir, control=False, keep_dir=None):
    """Status of every centromere and null window of `code` in the assembly at
    query_fa, for every variant. Returns {variant name: DataFrame}, plus the
    breakpoints of the query assembly on the long-read assembly. With
    keep_dir, the primary flank sequences are kept there, gzip, for audit."""
    os.makedirs(workdir, exist_ok=True)
    long_fa = P.data_dir("assemblies", acc + ".fna.gz")
    long_seqs = dict(P.read_fasta(long_fa))
    target_seqs = dict(P.read_fasta(query_fa))
    by_seq = {}
    for c in query_calls.itertuples():
        by_seq.setdefault(c.contig, []).append((int(c.start), int(c.end)))
    pafs = {}
    sets = [(VARIANTS[0]["flank"], VARIANTS[0]["half_bin"])] if control else ELEMENT_SETS
    for flank, half in sets:
        cens, nulls = elements(code, flank, half)
        el = pd.concat([cens, nulls], ignore_index=True) if not control else cens
        fq = os.path.join(workdir, "flanks_f%d_b%g.fa" % (flank, half))
        L.flank_fasta(el, long_seqs, fq, flank)
        paf = os.path.join(workdir, "flanks_f%d_b%g.paf" % (flank, half))
        L.run_minimap2(query_fa, fq, paf, "asm10", threads())
        pafs[(flank, half)] = (el, paf)
        if keep_dir and (flank, half) == (VARIANTS[0]["flank"], VARIANTS[0]["half_bin"]):
            import gzip
            with open(fq, "rb") as src, gzip.open(os.path.join(keep_dir, "flanks_primary.fa.gz"), "wb") as dst:
                shutil.copyfileobj(src, dst)
    out = {}
    for v in VARIANTS:
        el, paf = pafs[(v["flank"], v["half_bin"])]
        out[v["name"]] = L.classify_elements(el, paf, target_seqs, by_seq, identity=v["identity"])
        if control:
            break
    if control:
        return out, None, 0
    bpaf = os.path.join(workdir, "contigs_on_long.paf")
    L.run_minimap2(long_fa, query_fa, bpaf, "asm5", threads())
    bps, unlocated = L.breakpoints(bpaf, {n: len(s) for n, s in long_seqs.items()})
    return out, bps, unlocated


def contiguity(records):
    """Contig N50, L50 and count, with sequences split at every run of N (the
    rule of scripts/lib/assembly_stats.py, so that every arm sits on one
    axis), and the sequence count, total length, N runs and N bases."""
    seq_lengths, contig_lengths = [], []
    gap_runs = gap_bases = 0
    for _, seq in records:
        s = seq.upper()
        seq_lengths.append(len(s))
        contig_lengths.extend(P.contig_lengths_split_at_gaps(s))
        gap_bases += s.count("N")
        gap_runs += len(re.findall("N+", s))
    return {"n_sequences": len(seq_lengths), "total_length": sum(seq_lengths), "sequence_n50": P.nx(seq_lengths),
            "contig_n50": P.nx(contig_lengths), "contig_l50": P.lx(contig_lengths), "n_contigs": len(contig_lengths),
            "n_gap_runs": gap_runs, "gap_bases": gap_bases}


def cdeii_differences(st, long_calls, query_calls):
    lc = long_calls.set_index(["contig", "start"])["cdeii_len"].to_dict()
    qc = query_calls.set_index(["contig", "start"])["cdeii_len"].to_dict()
    diffs, qlen = [], []
    for r in st.itertuples():
        if r.kind == "centromere" and r.status == "intact_called":
            q = qc.get((r.target, int(r.call_start)))
            l_ = lc.get((r.contig, int(r.start)))
            diffs.append(q - l_ if q is not None and l_ is not None else np.nan)
            qlen.append(q)
        else:
            diffs.append(np.nan)
            qlen.append(np.nan)
    st = st.copy()
    st["query_cdeii_len"] = qlen
    st["cdeii_difference"] = diffs
    return st


def breakpoint_hits(code, acc, bps):
    """Window table of the long-read assembly and the windows its breakpoints
    fall in, saved for the breakpoint-model fit."""
    wpath = P.data_dir("pairs", code, "windows.tsv.gz")
    if os.path.exists(wpath):
        windows = pd.read_csv(wpath, sep="\t")
    else:
        windows = L.window_table(dict(P.read_fasta(P.data_dir("assemblies", acc + ".fna.gz"))))
        windows.to_csv(wpath, sep="\t", index=False, compression="gzip")
    if bps is None or bps.empty:
        return windows, windows.iloc[0:0]
    key = pd.DataFrame({"contig": bps["target"], "window": bps["position"] // L.WINDOW})
    hits = key.merge(windows, on=["contig", "window"], how="left").dropna(subset=["at"])
    return windows, hits


# ----------------------------------------------------------------------------
# arm 2

def strain_arm2(pair):
    code, acc = pair["code"], pair["long_read_accession"]
    d = P.data_dir("arm2", code)
    done = os.path.join(d, "done")
    if os.path.exists(done):
        return code
    os.makedirs(d, exist_ok=True)
    long_calls, _ = long_side(code, acc)
    short_fa = P.data_dir("assemblies", "peter2018", code + ".fna.gz")
    meta = P.run_pcan(short_fa, GENUS, os.path.join(d, "short_calls.tsv"), label=code + "_short")
    if meta.get("status") != "ok":
        raise RuntimeError("PCAn failed on the short-read assembly of %s" % code)
    short_calls = pd.read_csv(os.path.join(d, "short_calls.tsv"), sep="\t")
    work = P.data_dir("work", "arm2_" + code)
    try:
        res, bps, unlocated = compare(code, acc, short_fa, short_calls, work, keep_dir=d)
        for name, st in res.items():
            st = cdeii_differences(st, long_calls, short_calls)
            st.to_csv(os.path.join(d, "status_%s.tsv" % name), sep="\t", index=False)
        ctrl, _, _ = compare(code, acc, P.data_dir("assemblies", acc + ".fna.gz"), long_calls, work + "_self",
                             control=True)
        ctrl["primary"].to_csv(os.path.join(d, "control_self.tsv"), sep="\t", index=False)
        bps.to_csv(os.path.join(d, "breakpoints.tsv.gz"), sep="\t", index=False, compression="gzip")
        with open(os.path.join(d, "breakpoints_unlocated.txt"), "w") as fh:
            fh.write("%d\n" % unlocated)
        breakpoint_hits(code, acc, bps)
    finally:
        shutil.rmtree(work, ignore_errors=True)
        shutil.rmtree(work + "_self", ignore_errors=True)
    with open(done, "w") as fh:
        fh.write(P.utc_now() + "\n")
    return code


FLANK_REASONS = ["placed", "no full-length match", "two or more full-length matches",
                 "one full-length match, not placed"]


def flank_reasons(d, query_fa, flanks_gz):
    """Exploratory, not in the analysis plan. Why each flank of the primary
    analysis was placed or not: the kept flank sequences are aligned again
    exactly as in the comparison, and every alignment that passes the identity
    and query-coverage thresholds is counted, secondary ones included. A flank
    whose sequence occurs twice in the other assembly, as both haplotypes of a
    heterozygous strain can in a short-read assembly, gets a low mapping
    quality and is left unplaced by the registered rule; this tells that case
    apart from a flank that is missing or broken. Cached in the strain's
    directory."""
    out = os.path.join(d, "flank_reasons.tsv")
    if os.path.exists(out):
        return pd.read_csv(out, sep="\t", keep_default_na=False, na_values=[""])
    work = d + "_flankcheck"
    os.makedirs(work, exist_ok=True)
    try:
        fq = os.path.join(work, "flanks.fa")
        with open(fq, "w") as fh:
            for name, seq in P.read_fasta(flanks_gz):
                fh.write(">%s\n%s\n" % (name, seq))
        paf = os.path.join(work, "flanks.paf")
        L.run_minimap2(query_fa, fq, paf, "asm10", threads())
        with open(paf) as fh:
            alns = S.parse_paf(fh)
        names = [n for n, _ in P.read_fasta(fq)]
    finally:
        shutil.rmtree(work, ignore_errors=True)
    by_q = {}
    for a in alns:
        by_q.setdefault(a.qname, []).append(a)
    rows = []
    for n in names:
        al = by_q.get(n, [])
        full = [a for a in al if a.identity >= L.IDENTITY and a.qcov >= L.QCOV]
        prim = [a for a in full if a.primary]
        if S.place_flank(al, L.IDENTITY, L.QCOV, L.MAPQ) is not None:
            reason = "placed"
        elif not full:
            reason = "no full-length match"
        elif len(full) >= 2:
            reason = "two or more full-length matches"
        else:
            # one match that passes identity and coverage, but its mapping
            # quality is below 20 or it is a secondary alignment
            reason = "one full-length match, not placed"
        element_id, side = n.rsplit("|", 1)
        rows.append({"element_id": element_id, "flank": side, "full_length_matches": len(full),
                     "primary_mapq": max((a.mapq for a in prim), default=""), "reason": reason})
    t = pd.DataFrame(rows)
    P.atomic_write_tsv(t, out)
    return t


def element_reason(left, right, status):
    """One reason per element, from its two flanks: a flank with a second copy
    first, then a flank with no match, then the placement of placed flanks."""
    rs = {left, right}
    if status in ("intact_called", "intact_uncalled", "n_run"):
        return "both flanks placed, element between them"
    if not all(isinstance(x, str) for x in rs):
        return "flanks not checked"
    if "two or more full-length matches" in rs:
        return "a flank matches two or more places"
    if "no full-length match" in rs:
        return "a flank matches nowhere"
    if rs == {"placed"}:
        return "both flanks placed, inconsistent with one intact element"
    return "a flank with one full-length match, not placed"


def possibly_allelic(code, cens_status, long_calls, short_calls):
    """Supplementary Data 4 of the PCAn paper lists, per assembly and per
    chromosome, up to two centromere sequences (columns k_A and k_B). Where
    it lists two different ones for this strain at the chromosome whose
    sequence matches either assembly's call, a disagreement may be allelic."""
    sd4 = pd.read_csv(P.data_dir("tables", "supp_data4.tsv"), sep="\t", dtype=str)
    row = sd4[sd4["Assembly"].astype(str).str.strip() == code.replace("SACE_", "")]
    if row.empty:
        return pd.Series(["no Supplementary Data 4 row"] * len(cens_status), index=cens_status.index)
    row = row.iloc[0]
    pairs = []
    for k in range(1, 17):
        a, b = str(row.get("%d_A" % k, "")), str(row.get("%d_B" % k, ""))
        a = "" if a == "nan" else a.upper()
        b = "" if b == "nan" else b.upper()
        pairs.append((a, b))
    lseq = long_calls.set_index(["contig", "start"])["sequence"].str.upper().to_dict()
    qseq = short_calls.set_index(["contig", "start"])["sequence"].str.upper().to_dict()
    out = []
    for r in cens_status.itertuples():
        seqs = {lseq.get((r.contig, int(r.start)), "")}
        if r.status == "intact_called":
            seqs.add(qseq.get((r.target, int(r.call_start)), ""))
        flag = "no"
        for a, b in pairs:
            if (a and a in seqs) or (b and b in seqs):
                flag = "yes" if a and b and a != b else "no"
                break
        else:
            flag = "chromosome not identified"
        out.append(flag)
    return pd.Series(out, index=cens_status.index)


def aggregate_arm2(pairs):
    out = P.repo("results", "arm2")
    sens = P.repo("results", "sensitivity")
    os.makedirs(out, exist_ok=True)
    runs, counts, cen_all, null_all, ctrl_all, short_all, extra = [], [], {}, {}, [], [], []
    bsets, sstats, explore = {}, [], []
    for pr in pairs.to_dict("records"):
        code, acc = pr["code"], pr["long_read_accession"]
        d = P.data_dir("arm2", code)
        if not os.path.exists(os.path.join(d, "done")):
            continue
        short_fa = P.data_dir("assemblies", "peter2018", code + ".fna.gz")
        recs = P.read_fasta(short_fa)
        sstats.append(dict({"strain": code}, **contiguity(recs)))
        short_len = {n: len(s) for n, s in recs}
        del recs
        explore.append((code, pr["peter_zygosity"], flank_reasons(d, short_fa, os.path.join(d, "flanks_primary.fa.gz")),
                        short_len))
        lc = pd.read_csv(P.data_dir("pairs", code, "long_calls.tsv"), sep="\t")
        sc = pd.read_csv(os.path.join(d, "short_calls.tsv"), sep="\t")
        for side, path in (("long", P.data_dir("pairs", code, "long_calls.tsv")), ("short", os.path.join(d, "short_calls.tsv"))):
            m = pd.read_csv(path + ".meta.tsv", sep="\t", header=None, index_col=0)[1].to_dict()
            runs.append({"strain": code, "assembly": side, "n_calls": m.get("n_calls"), "elapsed_s": m.get("elapsed_s"),
                         "max_rss_mb": m.get("max_rss_mb"), "trace_matches_pcan": m.get("trace_matches_pcan"),
                         "fimo_pass1_hits": m.get("fimo_pass1_hits")})
        counts.append({"strain": code, "isolate_name": pr["isolate_name"], "long_read_accession": acc,
                       "long_calls": len(lc), "short_calls": len(sc), "zygosity": pr["peter_zygosity"],
                       "ploidy": pr["peter_ploidy"]})
        for v in VARIANTS:
            st = read_elements(os.path.join(d, "status_%s.tsv" % v["name"]))
            st.insert(0, "strain", code)
            st["zygosity"] = pr["peter_zygosity"]
            cen_all.setdefault(v["name"], []).append(st[st["kind"] == "centromere"])
            null_all.setdefault(v["name"], []).append(st[st["kind"] == "null"])
        prim = cen_all["primary"][-1]
        prim["possibly_allelic"] = possibly_allelic(code, prim, lc, sc).values
        cen_all["primary"][-1] = prim
        used = set(zip(prim.loc[prim["status"] == "intact_called", "target"],
                       prim.loc[prim["status"] == "intact_called", "call_start"].astype(float).astype(int)))
        extra.append({"strain": code, "short_calls": len(sc),
                      "short_calls_without_long_counterpart": int(sum((c.contig, int(c.start)) not in used
                                                                      for c in sc.itertuples()))})
        ctrl = pd.read_csv(os.path.join(d, "control_self.tsv"), sep="\t")
        ctrl.insert(0, "strain", code)
        ctrl_all.append(ctrl)
        bps = pd.read_csv(os.path.join(d, "breakpoints.tsv.gz"), sep="\t")
        windows, hits = breakpoint_hits(code, acc, bps)
        bsets[code] = (windows, hits)
        with open(os.path.join(d, "breakpoints_unlocated.txt")) as fh:
            unl = int(fh.read().strip())
        short_all.append({"strain": code, "breakpoints": len(hits), "contig_ends_not_located": unl,
                          "mean_at_at_breakpoints": hits["at"].mean() if len(hits) else np.nan,
                          "genome_at": float(np.average(windows["at"], weights=windows["length"]))})
    P.atomic_write_tsv(pd.DataFrame(runs), os.path.join(out, "pcan_runs.tsv"))
    P.atomic_write_tsv(pd.DataFrame(sstats), os.path.join(out, "short_read_assembly_stats.tsv"))
    counts = pd.DataFrame(counts)
    P.atomic_write_tsv(counts, os.path.join(out, "call_counts.tsv"))
    cens = pd.concat(cen_all["primary"], ignore_index=True)
    nulls = pd.concat(null_all["primary"], ignore_index=True)
    P.atomic_write_tsv(cens, os.path.join(out, "centromere_status.tsv"))
    P.atomic_write_tsv(nulls[["strain", "element_id", "for_centromere", "contig", "start", "end", "length", "at",
                              "status", "target"]], os.path.join(out, "null_status.tsv"))
    P.atomic_write_tsv(pd.concat(ctrl_all, ignore_index=True), os.path.join(out, "liftover_control.tsv"))
    P.atomic_write_tsv(pd.DataFrame(extra), os.path.join(out, "short_calls_without_long_counterpart.tsv"))
    P.atomic_write_tsv(pd.DataFrame(short_all), os.path.join(out, "breakpoint_summary.tsv"))

    conf = [PS.c1_counts(counts)] + PS.c2_status(cens) + [PS.c3_cdeii(cens)] + PS.c4_breaks(cens, nulls) \
        + [PS.c5_model(bsets)]
    P.atomic_write_tsv(pd.DataFrame(conf), os.path.join(out, "confirmatory.tsv"))

    srows = []
    for v in VARIANTS:
        c = pd.concat(cen_all[v["name"]], ignore_index=True)
        n = pd.concat(null_all[v["name"]], ignore_index=True)
        for r in PS.c2_status(c) + [PS.c3_cdeii(c)] + PS.c4_breaks(c, n):
            srows.append(dict(r, variant=v["name"]))
    homo = set(counts.loc[counts["zygosity"].str.lower() == "homozygous", "strain"])
    c = cens[cens["strain"].isin(homo)]
    n = nulls[nulls["strain"].isin(homo)]
    for r in [PS.c1_counts(counts[counts["strain"].isin(homo)])] + PS.c2_status(c) + [PS.c3_cdeii(c)] \
            + PS.c4_breaks(c, n):
        srows.append(dict(r, variant="heterozygous_excluded"))
    P.atomic_write_tsv(pd.DataFrame(srows), os.path.join(sens, "arm2_variants.tsv"))
    explore_arm2(explore, cens, nulls, out)
    P.log("arm 2 aggregated: %d strains" % len(counts))


def explore_arm2(explore, cens, nulls, out):
    """Exploratory tables, outside the analysis plan: why each broken element
    is broken, from flank_reasons, and whether each intact_called centromere
    can still be checked for synteny in the short-read assembly (its contig
    extends 10 kb beyond both ends of the call, the arm 1 definition)."""
    el = pd.concat([cens[cens["status"] != "not_liftable"][["strain", "element_id", "kind", "status", "target",
                                                            "call_start", "call_end"]],
                    nulls[["strain", "element_id", "status"]].assign(kind="null")], ignore_index=True)
    reasons, lengths, zyg = [], {}, {}
    for code, z, fr, short_len in explore:
        w = fr.pivot(index="element_id", columns="flank", values="reason").reset_index()
        w["strain"] = code
        reasons.append(w.rename(columns={"L": "left_flank", "R": "right_flank"}))
        lengths[code], zyg[code] = short_len, z
    el = el.merge(pd.concat(reasons, ignore_index=True), on=["strain", "element_id"], how="left")
    el["zygosity"] = el["strain"].map(zyg)
    el["element_reason"] = [element_reason(a, b, s) for a, b, s in zip(el["left_flank"], el["right_flank"], el["status"])]
    syn = []
    for r in el.itertuples():
        if r.kind == "centromere" and r.status == "intact_called":
            n = lengths[r.strain].get(r.target, 0)
            syn.append("yes" if int(r.call_start) - 1 >= 10000 and n - int(r.call_end) >= 10000 else "no")
        else:
            syn.append("")
    el["short_read_synteny_checkable"] = syn
    c = el[el["kind"] == "centromere"]
    P.atomic_write_tsv(c[["strain", "zygosity", "element_id", "status", "left_flank", "right_flank", "element_reason",
                          "short_read_synteny_checkable"]], os.path.join(out, "exploratory_centromere_flanks.tsv"))
    summ = el.groupby(["zygosity", "kind", "status", "element_reason"]).size().rename("elements").reset_index()
    tot = el.groupby(["zygosity", "kind"]).size().rename("of_elements").reset_index()
    summ = summ.merge(tot, on=["zygosity", "kind"])
    summ["fraction"] = summ["elements"] / summ["of_elements"]
    P.atomic_write_tsv(summ, os.path.join(out, "exploratory_flank_summary.tsv"))


# ----------------------------------------------------------------------------
# arm 3

def strain_arm3(row):
    """Every assembly that 11_assemble_reads.sh built for one strain."""
    code, acc = row["code"], row["long_read_accession"]
    d = P.data_dir("arm3", code)
    man_path = os.path.join(d, "manifest.tsv")
    if not os.path.exists(man_path):
        return code
    man = pd.read_csv(man_path, sep="\t")
    long_calls, _ = long_side(code, acc)
    for a in man[man["status"] == "ok"].to_dict("records"):
        tag = "%s_d%s_s%s" % (a["assembler"], a["target_depth"], a["seed"])
        res_dir = os.path.join(d, "compare", tag)
        if os.path.exists(os.path.join(res_dir, "done")):
            continue
        os.makedirs(res_dir, exist_ok=True)
        fa = os.path.join(d, "assemblies", tag + ".fna.gz")
        meta = P.run_pcan(fa, GENUS, os.path.join(res_dir, "calls.tsv"), label=code + "_" + tag)
        if meta.get("status") != "ok":
            with open(os.path.join(res_dir, "pcan_failed"), "w") as fh:
                fh.write(str(meta.get("error", "")) + "\n")
            continue
        calls = pd.read_csv(os.path.join(res_dir, "calls.tsv"), sep="\t")
        work = P.data_dir("work", "arm3_" + code + "_" + tag)
        try:
            res, bps, unlocated = compare(code, acc, fa, calls, work, keep_dir=res_dir)
            for name, st in res.items():
                cdeii_differences(st, long_calls, calls).to_csv(os.path.join(res_dir, "status_%s.tsv" % name),
                                                                sep="\t", index=False)
            bps.to_csv(os.path.join(res_dir, "breakpoints.tsv.gz"), sep="\t", index=False, compression="gzip")
            with open(os.path.join(res_dir, "breakpoints_unlocated.txt"), "w") as fh:
                fh.write("%d\n" % unlocated)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        with open(os.path.join(res_dir, "done"), "w") as fh:
            fh.write(P.utc_now() + "\n")
    return code


def aggregate_arm3(runs3):
    out = P.repo("results", "arm3")
    os.makedirs(out, exist_ok=True)
    asm_rows, cen_rows, null_rows, bsets = [], [], [], {}
    for row in runs3.to_dict("records"):
        code, acc = row["code"], row["long_read_accession"]
        d = P.data_dir("arm3", code)
        if not os.path.exists(os.path.join(d, "manifest.tsv")):
            continue
        man = pd.read_csv(os.path.join(d, "manifest.tsv"), sep="\t")
        lc = pd.read_csv(P.data_dir("pairs", code, "long_calls.tsv"), sep="\t")
        for a in man.to_dict("records"):
            tag = "%s_d%s_s%s" % (a["assembler"], a["target_depth"], a["seed"])
            res_dir = os.path.join(d, "compare", tag)
            r = dict(a, strain=code, long_calls=len(lc))
            if os.path.exists(os.path.join(res_dir, "status_primary.tsv")):
                calls = pd.read_csv(os.path.join(res_dir, "calls.tsv"), sep="\t")
                r["calls"] = len(calls)
                for v in VARIANTS:
                    st = read_elements(os.path.join(res_dir, "status_%s.tsv" % v["name"]))
                    st.insert(0, "strain", code)
                    for k, vv in (("assembler", a["assembler"]), ("target_depth", a["target_depth"]),
                                  ("seed", a["seed"]), ("variant", v["name"])):
                        st[k] = vv
                    cen_rows.append(st[st["kind"] == "centromere"])
                    null_rows.append(st[st["kind"] == "null"])
                if a["assembler"] in ("spades", "megahit") and str(a["target_depth"]) == "full" and int(a["seed"]) == 11:
                    bps = pd.read_csv(os.path.join(res_dir, "breakpoints.tsv.gz"), sep="\t")
                    w, h = breakpoint_hits(code, acc, bps)
                    bsets.setdefault(a["assembler"], {})[code] = (w, h)
            asm_rows.append(r)
    asm = pd.DataFrame(asm_rows)
    P.atomic_write_tsv(asm, os.path.join(out, "assemblies.tsv"))
    if not cen_rows:
        P.log("arm 3: nothing to aggregate yet")
        return
    cens = pd.concat(cen_rows, ignore_index=True)
    nulls = pd.concat(null_rows, ignore_index=True)
    prim = cens[cens["variant"] == "primary"]
    P.atomic_write_tsv(prim, os.path.join(out, "centromere_status.tsv"))

    def per_strain(df, hit):
        """Each strain's fraction, strains in sorted order: a bootstrap over
        these draws exactly the strains one over per-strain tables would."""
        return pd.DataFrame({"strain": df["strain"].values, "hit": np.asarray(hit, dtype=float)}) \
            .groupby("strain")["hit"].mean().tolist()

    def mean_of(us):
        return float(np.mean(us))

    conf, sens = [], []
    for (assembler, depth, seed), c in prim.groupby(["assembler", "target_depth", "seed"]):
        lift = c[c["status"] != "not_liftable"]
        units = per_strain(lift, lift["status"] == "intact_called")
        est, lo, hi = P.cluster_bootstrap(units, mean_of, n_boot=PS.N_BOOT, seed=PS.SEED)
        conf.append({"outcome": "D1", "assembler": assembler, "target_depth": depth, "seed": seed,
                     "measure": "recall relative to the long-read calls", "n_strains": len(units),
                     "estimate": est, "ci_low": lo, "ci_high": hi, "median": float(np.median(units)),
                     "q25": float(np.percentile(units, 25)), "q75": float(np.percentile(units, 75)),
                     "interval": "percentile bootstrap over strains"})
        called = c[c["status"] == "intact_called"]
        cu = per_strain(called, called["cdeii_difference"] == 0)
        if cu:
            e2, l2, h2 = P.cluster_bootstrap(cu, mean_of, n_boot=PS.N_BOOT, seed=PS.SEED)
            conf.append({"outcome": "D2", "assembler": assembler, "target_depth": depth, "seed": seed,
                         "measure": "identical CDEII length among intact_called", "n_strains": len(cu),
                         "estimate": e2, "ci_low": l2, "ci_high": h2, "median": float(np.median(cu)),
                         "q25": float(np.percentile(cu, 25)), "q75": float(np.percentile(cu, 75)),
                         "interval": "percentile bootstrap over strains"})
        if int(seed) == 11:
            n = nulls[(nulls["variant"] == "primary") & (nulls["assembler"] == assembler)
                      & (nulls["target_depth"].astype(str) == str(depth)) & (nulls["seed"] == seed)]
            for r in PS.c4_breaks(c, n, label="D4"):
                conf.append(dict(r, assembler=assembler, target_depth=depth, seed=seed))
    for assembler, sets in bsets.items():
        conf.append(dict(PS.c5_model(sets, label="D5"), assembler=assembler, target_depth="full", seed=11))
    P.atomic_write_tsv(pd.DataFrame(conf), os.path.join(out, "confirmatory.tsv"))
    for (variant, assembler, depth, seed), c in cens.groupby(["variant", "assembler", "target_depth", "seed"]):
        if int(seed) != 11 or variant == "primary":
            continue
        n = nulls[(nulls["variant"] == variant) & (nulls["assembler"] == assembler)
                  & (nulls["target_depth"].astype(str) == str(depth)) & (nulls["seed"] == seed)]
        lift = c[c["status"] != "not_liftable"]
        units = per_strain(lift, lift["status"] == "intact_called")
        est, lo, hi = P.cluster_bootstrap(units, mean_of, n_boot=PS.N_BOOT, seed=PS.SEED)
        sens.append({"variant": variant, "outcome": "D1", "assembler": assembler, "target_depth": depth,
                     "measure": "recall relative to the long-read calls", "estimate": est, "ci_low": lo,
                     "ci_high": hi})
        for r in PS.c4_breaks(c, n, label="D4"):
            sens.append(dict(r, variant=variant, assembler=assembler, target_depth=depth))
    P.atomic_write_tsv(pd.DataFrame(sens), P.repo("results", "sensitivity", "arm3_variants.tsv"))
    P.log("arm 3 aggregated: %d assemblies" % len(asm))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", choices=["2", "3"], required=True)
    ap.add_argument("--jobs", type=int, default=int(P.conf().get("JOBS", 1)))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    os.environ["PAD_JOBS"] = str(args.jobs)
    stage = "compare_arm%s" % args.arm
    if P.is_done(stage) and not args.force:
        P.log("arm %s comparison already complete, skipping (use --force to redo)" % args.arm)
        return
    if args.arm == "2":
        pairs = pd.read_csv(P.repo("config", "arm2_pairs.tsv"), sep="\t", dtype=str)
        rows, worker, agg = pairs.to_dict("records"), strain_arm2, lambda: aggregate_arm2(pairs)
    else:
        runs3 = pd.read_csv(P.repo("config", "arm3_runs.tsv"), sep="\t", dtype=str)
        runs3 = runs3[runs3["arm3_eligible"] == "yes"]
        rows, worker, agg = runs3.to_dict("records"), strain_arm3, lambda: aggregate_arm3(runs3)
    failures = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(worker, r): r["code"] for r in rows}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                f.result()
            except Exception as e:  # one strain failing is recorded, not fatal
                failures.append({"strain": futs[f], "error": str(e)[:500]})
                P.log("%s failed: %s" % (futs[f], e))
            if i % 10 == 0:
                P.log("%d of %d strains compared" % (i, len(rows)))
    if failures:
        P.atomic_write_tsv(pd.DataFrame(failures), P.repo("results", "arm%s" % args.arm, "strain_failures.tsv"))
    agg()
    if not failures:
        P.mark_done(stage)


if __name__ == "__main__":
    main()
