"""Arm 1b, a simulation: plant CDEII length variants in real centromeres of
three arm 1 genomes, run PCAn on each edited genome, record the calls, and
delete the FASTA. Scoring is done by scripts/08_score_perturbations.py.

Edits follow scripts/lib/variants.py: an insertion is a tandem duplication of
adjacent CDEII sequence and a deletion removes a stretch from inside CDEII.
CDEI and the CDEIII motif are never touched. The design was committed in
config/arm1_design.md before this script first ran.

Usage:
    python scripts/07_plant_variants.py [--jobs N] [--force]
"""
import argparse
import hashlib
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import padlib as P  # noqa: E402
import variants as V  # noqa: E402

S288C = "GCA_000146045.2"
SIZES = [-10, -5, 5, 10, 15, 20, 30, 40]
PROGRESSIVE_K = [1, 2, 4, 8, 12, 16]
PROGRESSIVE_SIZE = 10
PROGRESSIVE_REPS = 5
OUT = P.repo("results", "arm1")


def seed_of(*parts):
    return int(hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()[:8], 16)


def truth_calls(acc):
    return pd.read_csv(P.data_dir("arm0", "calls", acc + ".tsv"), sep="\t")


def select_genomes():
    elig = pd.read_csv(os.path.join(OUT, "eligibility.tsv"), sep="\t", dtype={"accession": str})
    elig = elig[elig["eligible"] == "yes"].copy()
    elig["median_cdeii_len"] = [truth_calls(a)["cdeii_len"].median() for a in elig["accession"]]
    others = elig[elig["accession"] != S288C]
    shortest = others.sort_values(["median_cdeii_len", "accession"], ascending=[True, True]).iloc[0]
    longest = others.sort_values(["median_cdeii_len", "accession"], ascending=[False, True]).iloc[0]
    rows = [elig[elig["accession"] == S288C].iloc[0], shortest, longest]
    roles = ["reference", "shortest median CDEII", "longest median CDEII"]
    sel = pd.DataFrame(rows).reset_index(drop=True)
    sel.insert(0, "role", roles)
    P.atomic_write_tsv(sel[["role", "accession", "species", "genus", "kazachstania_species", "calls",
                            "median_cdeii_len"]], os.path.join(OUT, "variant_genomes.tsv"))
    return sel


_GENOMES = {}


def records(acc):
    if acc not in _GENOMES:
        _GENOMES[acc] = P.read_fasta(P.data_dir("assemblies", acc + ".fna.gz"))
    return _GENOMES[acc]


def run_edited(run_id, acc, genus, kaz, edits):
    """edits: [(call dict, size, position)]. Writes the edited genome, runs
    PCAn, deletes the FASTA. Returns the run's metadata."""
    out = P.data_dir("arm1", "variant_calls", run_id + ".tsv")
    if os.path.exists(out) and os.path.exists(out + ".meta.tsv"):
        meta = {}
        with open(out + ".meta.tsv") as fh:
            for line in fh:
                k, v = line.rstrip("\n").split("\t", 1)
                meta[k] = v
        return meta
    recs = dict(records(acc))
    for call, size, pos in edits:
        new, _ = V.apply_edit(recs[call["contig"]], call, size, pos)
        recs[call["contig"]] = new
    work = P.data_dir("work", "arm1b")
    os.makedirs(work, exist_ok=True)
    fasta = os.path.join(work, run_id + ".fna")
    try:
        P.write_fasta(fasta, [(n, recs[n]) for n, _ in records(acc)])
        return P.run_pcan(fasta, genus, out, species=kaz, label=run_id)
    finally:
        if os.path.exists(fasta):
            os.remove(fasta)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=int(P.conf().get("JOBS", 1)))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if P.is_done("plant_variants") and not args.force:
        P.log("arm 1b already complete, skipping (use --force to redo)")
        return
    os.makedirs(P.data_dir("arm1", "variant_calls"), exist_ok=True)
    sel = select_genomes()
    P.log("arm 1b genomes: %s" % ", ".join("%s (%s, median CDEII %.0f)" % (r.accession, r.role, r.median_cdeii_len)
                                         for r in sel.itertuples()))

    jobs, edit_rows = [], []
    for g in sel.itertuples():
        kaz = g.kazachstania_species if isinstance(g.kazachstania_species, str) else ""
        seqs = dict(records(g.accession))
        for c in truth_calls(g.accession).to_dict("records"):
            for size in SIZES:
                rid = "%s_%s_%+d" % (g.accession, c["contig"], size)
                pos = V.draw_position(c, size, np.random.default_rng(seed_of(g.accession, c["contig"], size)))
                row = {"run_id": rid, "design": "single", "accession": g.accession, "contig": c["contig"],
                       "strand": c["strand"], "call_start": c["start"], "call_end": c["end"],
                       "unedited_cdeii_len": c["cdeii_len"], "size": size}
                if pos is None:
                    row.update(position="", status="not possible: CDEII too short for this edit")
                else:
                    _, rec = V.apply_edit(seqs[c["contig"]], c, size, pos)
                    row.update(position=pos, changed_bases=rec["changed_bases"], expected_start=rec["expected_start"],
                               expected_end=rec["expected_end"], expected_cdeii_len=rec["expected_cdeii_len"])
                    jobs.append((rid, g.accession, g.genus, kaz, [(c, size, pos)]))
                edit_rows.append(row)

    s = sel[sel["accession"] == S288C].iloc[0]
    s_calls = truth_calls(S288C).to_dict("records")
    s_seqs = dict(records(S288C))
    for k in PROGRESSIVE_K:
        for rep in range(1, PROGRESSIVE_REPS + 1):
            rng = np.random.default_rng(seed_of("progressive", k, rep))
            chosen = sorted(rng.choice(len(s_calls), size=k, replace=False).tolist())
            edits = []
            rid = "%s_progressive_k%02d_r%d" % (S288C, k, rep)
            for i in chosen:
                c = s_calls[i]
                pos = V.draw_position(c, PROGRESSIVE_SIZE, rng)
                _, rec = V.apply_edit(s_seqs[c["contig"]], c, PROGRESSIVE_SIZE, pos)
                edits.append((c, PROGRESSIVE_SIZE, pos))
                edit_rows.append({"run_id": rid, "design": "progressive", "accession": S288C, "contig": c["contig"],
                                  "strand": c["strand"], "call_start": c["start"], "call_end": c["end"],
                                  "unedited_cdeii_len": c["cdeii_len"], "size": PROGRESSIVE_SIZE, "k": k,
                                  "replicate": rep, "position": pos, "changed_bases": rec["changed_bases"],
                                  "expected_start": rec["expected_start"], "expected_end": rec["expected_end"],
                                  "expected_cdeii_len": rec["expected_cdeii_len"]})
            jobs.append((rid, S288C, s["genus"], "", edits))

    edits_df = pd.DataFrame(edit_rows)
    P.atomic_write_tsv(edits_df, os.path.join(OUT, "variant_edits.tsv"))
    P.log("arm 1b: %d edited genomes to run (%d single, %d progressive)"
          % (len(jobs), (edits_df["design"] == "single").sum(), len(PROGRESSIVE_K) * PROGRESSIVE_REPS))

    t0 = time.time()
    runs = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_edited, *j): j[0] for j in jobs}
        for i, f in enumerate(as_completed(futs), 1):
            m = f.result()
            m["run_id"] = futs[f]
            runs.append(m)
            if i % 100 == 0:
                P.log("%d of %d edited genomes done in %.0f min" % (i, len(jobs), (time.time() - t0) / 60))
    runs = pd.DataFrame(runs)
    keep = ["run_id", "status", "n_calls", "elapsed_s", "max_rss_mb", "trace_matches_pcan", "median_top5_cdeii_len",
            "final_median_cdeii_len", "final_median_cdeii_at", "error"]
    P.atomic_write_tsv(runs.reindex(columns=keep).sort_values("run_id"), os.path.join(OUT, "variant_runs.tsv"))
    P.log("arm 1b finished: %d runs, %d PCAn failures" % (len(runs), (runs["status"] != "ok").sum()))
    P.mark_done("plant_variants")


if __name__ == "__main__":
    main()
