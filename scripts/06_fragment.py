"""Arm 1a, a simulation: fragment real chromosome-level assemblies at chosen
positions, run PCAn on every replicate, record each replicate's contiguity,
and delete the FASTA. Scoring is done by scripts/08_score_perturbations.py.

The design, including the rule that picks the genomes, was committed in
config/arm1_design.md before this script first ran.

Usage:
    python scripts/06_fragment.py [--jobs N] [--force] [--select-only]
"""
import argparse
import gzip
import hashlib
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import fragment as F  # noqa: E402
import padlib as P  # noqa: E402

LEVELS_KB = [1000, 500, 200, 100, 50, 20, 10, 5]
MODELS = {"uniform": None, "at_weighted": math.log(10) / 0.30}
SENS_MODELS = {"at_weighted_g5": math.log(5) / 0.30, "at_weighted_g20": math.log(20) / 0.30}
SENS_LEVELS_KB = [50, 20, 10]
REPS, SENS_REPS, MIN_REPS = 10, 5, 5
OUT = P.repo("results", "arm1")
# A fragmented assembly costs PCAn more than the intact one: its script reparses
# the whole FASTA once per CDEIII hit, and the parse slows with many records.
# Until the first replicates are timed, arm 0's time per assembly is scaled by
# this factor for the projection. The factor is my guess.
SLOWDOWN_PRIOR = 2.0


def seed_of(*parts):
    return int(hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()[:8], 16)


def nuclear_chromosomes(acc):
    import json
    out = set()
    with open(P.data_dir("assemblies", acc + ".sequence_report.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("role") == "assembled-molecule" and r.get("assignedMoleculeLocationType") == "Chromosome":
                out.add(r["genbankAccession"])
    return out


def select_genomes():
    species = pd.read_csv(P.repo("config", "species_arm0.tsv"), sep="\t", dtype=str)
    stats = pd.read_csv(P.repo("results", "arm0", "assembly_stats.tsv"), sep="\t")
    runs = pd.read_csv(P.repo("logs", "arm0_runs.tsv"), sep="\t", dtype=str)
    ran = species[species["arm0_role"] == "run"].merge(stats, on="accession", how="inner")
    rows = []
    for _, s in ran.iterrows():
        if s["assembly_level"] not in ("Chromosome", "Complete Genome"):
            continue
        acc = s["accession"]
        r = {"accession": acc, "species": s["species"], "genus": s["pcan_genus"],
             "kazachstania_species": s["pcan_kazachstania_species"] if isinstance(s["pcan_kazachstania_species"], str) else "",
             "assembly_level": s["assembly_level"], "contig_n50": int(s["contig_n50"]),
             "total_length": int(s["total_length"]), "genome_at": float(s["genome_at"])}
        m = runs[runs["accession"] == acc]
        ok = len(m) and m.iloc[0]["status"] == "ok"
        nuclear = nuclear_chromosomes(acc)
        r["nuclear_chromosome_sequences"] = len(nuclear)
        why = []
        if not ok:
            why.append("arm 0 PCAn run did not finish")
            r["calls"] = ""
        else:
            calls = pd.read_csv(P.data_dir("arm0", "calls", acc + ".tsv"), sep="\t")
            r["calls"] = len(calls)
            if len(calls) != len(nuclear):
                why.append("%d calls for %d nuclear chromosome sequences" % (len(calls), len(nuclear)))
            if calls["contig"].duplicated().any() or not set(calls["contig"]) <= nuclear:
                why.append("calls not one per nuclear chromosome sequence")
        r["eligible"] = "no" if why else "yes"
        r["reason"] = "; ".join(why)
        rows.append(r)
    elig = pd.DataFrame(rows).sort_values(["genus", "accession"])
    chosen = []
    for genus, g in elig[elig["eligible"] == "yes"].groupby("genus"):
        g = g.sort_values(["contig_n50", "accession"], ascending=[False, True])
        chosen.append(g.iloc[0])
    chosen = pd.DataFrame(chosen).reset_index(drop=True)
    elig["chosen"] = elig["accession"].isin(chosen["accession"]).map({True: "yes", False: "no"})
    P.atomic_write_tsv(elig, os.path.join(OUT, "eligibility.tsv"))
    P.atomic_write_tsv(chosen, os.path.join(OUT, "genomes.tsv"))
    P.log("arm 1a: %d chromosome-level assemblies, %d eligible, %d genera chosen"
          % (len(elig), (elig["eligible"] == "yes").sum(), len(chosen)))
    return chosen


# ----------------------------------------------------------------------------
# genome data shared by the workers of one process

_G = {}
_LOCK = __import__("threading").Lock()


def genome(acc, gamma=None):
    """One genome's sequences, gap-free intervals and window weights, loaded
    once per process. A 12 Mb genome takes about 12 MB here."""
    with _LOCK:
        if acc not in _G:
            recs = P.read_fasta(P.data_dir("assemblies", acc + ".fna.gz"))
            _G[acc] = {"recs": recs, "intervals": {n: F.gap_free_intervals(s) for n, s in recs},
                       "at": P.at_fraction("".join(s for _, s in recs)), "weights": {}}
        g = _G[acc]
        if gamma is not None and gamma not in g["weights"]:
            g["weights"][gamma] = {n: F.window_weights(s, g["at"], gamma) for n, s in g["recs"]}
    return g


def model_kind(model):
    return "uniform" if model == "uniform" else "at_weighted"


def calibrate_genome(acc, tasks):
    """All calibrations of one genome in one worker process, which then drops
    the genome so that worker memory does not grow with the number of genomes."""
    out = []
    for _, model, gamma, level_kb in tasks:
        g = genome(acc, gamma)
        rate, n50 = F.calibrate_rate(g["recs"], level_kb * 1000, model_kind(model),
                                     seed=seed_of(acc, model, level_kb, "calibration"), genome_at=g["at"],
                                     gamma=gamma or F.GAMMA_DEFAULT, weights=g["weights"].get(gamma))
        out.append({"accession": acc, "model": model, "gamma": gamma if gamma else "",
                    "target_n50": level_kb * 1000, "rate_per_bp": rate, "calibration_mean_n50": round(n50, 1)})
    _G.pop(acc, None)
    return out


def rid_of(acc, model, level_kb, rep):
    return "%s_%s_%dk_r%02d" % (acc, model, level_kb, rep)


def run_replicate(job):
    acc, genus, kaz, model, gamma, level_kb, rep, rate = job
    rid = rid_of(acc, model, level_kb, rep)
    rec_path = P.data_dir("arm1", "records", rid + ".tsv")
    if os.path.exists(rec_path):
        return pd.read_csv(rec_path, sep="\t").iloc[0].to_dict()
    g = genome(acc, gamma)
    seed = seed_of(acc, model, level_kb, rep)
    rec = {"replicate_id": rid, "accession": acc, "model": model, "gamma": gamma if gamma else "",
           "target_n50": level_kb * 1000, "replicate": rep, "seed": seed, "rate_per_bp": rate}
    if rate == 0:
        lengths = F.genome_contig_lengths(g["intervals"], {})
        rec.update(n_cuts=0, realised_n50=P.nx(lengths), realised_l50=P.lx(lengths), n_contigs=len(lengths),
                   status="unperturbed", n_calls="", elapsed_s="", max_rss_mb="", trace_matches_pcan="")
    else:
        cuts = F.draw_cuts(g["recs"], rate, model_kind(model), np.random.default_rng(seed), g["at"],
                           gamma or F.GAMMA_DEFAULT, g["weights"].get(gamma))
        lengths = F.genome_contig_lengths(g["intervals"], cuts)
        work = P.data_dir("work", "arm1")
        os.makedirs(work, exist_ok=True)
        fasta = os.path.join(work, rid + ".fna")
        try:
            P.write_fasta(fasta, [(n, s) for n, s, _, _ in F.apply_cuts(g["recs"], cuts)])
            out = P.data_dir("arm1", "calls", rid + ".tsv")
            meta = P.run_pcan(fasta, genus, out, species=kaz, label=rid)
        finally:
            if os.path.exists(fasta):
                os.remove(fasta)
        with gzip.open(P.data_dir("arm1", "cuts", rid + ".tsv.gz"), "wt") as fh:
            fh.write("sequence\tcut_after\n")
            for name, cp in cuts.items():
                for p in cp:
                    fh.write("%s\t%d\n" % (name, p))
        rec.update(n_cuts=int(sum(len(c) for c in cuts.values())), realised_n50=P.nx(lengths),
                   realised_l50=P.lx(lengths), n_contigs=len(lengths), status=meta.get("status", "failed"),
                   n_calls=meta.get("n_calls", ""), elapsed_s=meta.get("elapsed_s", ""),
                   max_rss_mb=meta.get("max_rss_mb", ""), trace_matches_pcan=meta.get("trace_matches_pcan", ""))
    pd.DataFrame([rec]).to_csv(rec_path + ".part", sep="\t", index=False)
    os.replace(rec_path + ".part", rec_path)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=int(P.conf().get("JOBS", 1)))
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--select-only", action="store_true")
    args = ap.parse_args()
    if P.is_done("fragment") and not args.force:
        P.log("arm 1a already complete, skipping (use --force to redo)")
        return
    for d in ("calls", "cuts", "records"):
        os.makedirs(P.data_dir("arm1", d), exist_ok=True)
    chosen = select_genomes()
    if args.select_only:
        return

    # Calibrate every rate first, in parallel processes.
    cal_path = os.path.join(OUT, "calibration.tsv")
    tasks = [(a, m, g, lv) for a in chosen["accession"] for m, g in MODELS.items() for lv in LEVELS_KB]
    tasks += [(a, m, g, lv) for a in chosen["accession"] for m, g in SENS_MODELS.items() for lv in SENS_LEVELS_KB]
    cal = pd.read_csv(cal_path, sep="\t") if os.path.exists(cal_path) else pd.DataFrame()
    done = set() if cal.empty else set(zip(cal["accession"], cal["model"], cal["target_n50"] // 1000))
    todo = [t for t in tasks if (t[0], t[1], t[3]) not in done]
    if todo:
        P.log("calibrating %d cut rates" % len(todo))
        by_genome = {}
        for t in todo:
            by_genome.setdefault(t[0], []).append(t)
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(calibrate_genome, a, ts) for a, ts in by_genome.items()]
            new = [row for f in as_completed(futs) for row in f.result()]
        cal = pd.concat([cal, pd.DataFrame(new)], ignore_index=True)
        cal = cal.sort_values(["accession", "model", "target_n50"], ascending=[True, True, False])
        P.atomic_write_tsv(cal, cal_path)
    rates = {(r.accession, r.model, int(r.target_n50) // 1000): r.rate_per_bp for r in cal.itertuples()}

    # Projection from arm 0's measured time per assembly, and cuts if needed.
    runs0 = pd.read_csv(P.repo("logs", "arm0_runs.tsv"), sep="\t")
    t_unfrag = runs0.set_index("accession")["elapsed_s"].astype(float)
    budget_h = float(P.conf().get("HOURS_ARM1", 24))

    def plan(reps):
        jobs = []
        for rep in range(1, reps + 1):
            for _, gr in chosen.iterrows():
                for m, gm in MODELS.items():
                    for lv in LEVELS_KB:
                        jobs.append((gr["accession"], gr["genus"], gr["kazachstania_species"], m, gm, lv, rep,
                                     rates[(gr["accession"], m, lv)]))
        for rep in range(1, min(reps, SENS_REPS) + 1):
            for _, gr in chosen.iterrows():
                for m, gm in SENS_MODELS.items():
                    for lv in SENS_LEVELS_KB:
                        jobs.append((gr["accession"], gr["genus"], gr["kazachstania_species"], m, gm, lv, rep,
                                     rates[(gr["accession"], m, lv)]))
        return jobs

    def projected_hours(jobs):
        secs = sum(t_unfrag.get(j[0], t_unfrag.median()) * SLOWDOWN_PRIOR for j in jobs if j[7] > 0)
        return secs / args.jobs / 3600

    reps = REPS
    cuts = []
    while projected_hours(plan(reps)) > budget_h and reps > MIN_REPS:
        cuts.append({"cut": "replicates per level from %d to %d" % (reps, reps - 1),
                     "projected_hours_before": round(projected_hours(plan(reps)), 2), "budget_hours": budget_h})
        reps -= 1
    jobs = plan(reps)
    proj = projected_hours(jobs)
    if proj > budget_h:
        P.log("refusing: arm 1a projects %.1f h on %d jobs against a %.0f h budget, %.1f h short even at %d replicates"
              % (proj, args.jobs, budget_h, proj - budget_h, MIN_REPS))
        sys.exit(3)
    P.atomic_write_tsv(pd.DataFrame(cuts, columns=["cut", "projected_hours_before", "budget_hours"]),
                       os.path.join(OUT, "cuts.tsv"))
    P.log("arm 1a: %d genomes, %d replicate runs (%d replicates per level), projected %.1f h on %d jobs"
          % (len(chosen), len(jobs), reps, proj, args.jobs))

    t0 = time.time()
    records = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_replicate, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            records.append(f.result())
            if i in (2 * args.jobs, 200) or i % 500 == 0:
                done_t = time.time() - t0
                P.log("%d of %d replicates done in %.0f min; at this pace %.1f h remain"
                      % (i, len(jobs), done_t / 60, done_t / i * (len(jobs) - i) / 3600))
    rec = pd.DataFrame(records).sort_values(["accession", "model", "target_n50", "replicate"],
                                            ascending=[True, True, False, True])
    P.atomic_write_tsv(rec, os.path.join(OUT, "fragmentation_runs.tsv"))
    failed = rec[rec["status"] == "failed"]
    P.log("arm 1a finished: %d replicates, %d unperturbed, %d PCAn failures, %d trace mismatches"
          % (len(rec), (rec["status"] == "unperturbed").sum(), len(failed),
             (rec["trace_matches_pcan"] == "no").sum()))
    P.mark_done("fragment")


if __name__ == "__main__":
    main()
