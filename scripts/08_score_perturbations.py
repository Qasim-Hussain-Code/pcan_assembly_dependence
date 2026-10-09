"""Score arm 1, a simulation. For 1a: recall relative to the unfragmented call,
the geometric expectation, the loss decomposition, the synteny-checkable
fraction and gained calls. For 1b: the detection envelope and the progressive
transitions. Definitions are those committed in config/arm1_design.md.

Usage:
    python scripts/08_score_perturbations.py [--force]
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import fragment as F  # noqa: E402
import padlib as P  # noqa: E402

OUT = P.repo("results", "arm1")
SEED, N_BOOT = 20261009, 10000


def truth(acc):
    return pd.read_csv(P.data_dir("arm0", "calls", acc + ".tsv"), sep="\t")


def seq_lengths(acc):
    return {n: len(s) for n, s in P.read_fasta(P.data_dir("assemblies", acc + ".fna.gz"))}


def read_cuts(rid):
    cuts = {}
    with gzip.open(P.data_dir("arm1", "cuts", rid + ".tsv.gz"), "rt") as fh:
        next(fh)
        for line in fh:
            s, p = line.rstrip("\n").split("\t")
            cuts.setdefault(s, []).append(int(p))
    return {s: np.array(sorted(v), dtype=np.int64) for s, v in cuts.items()}


def to_source(df, contig_col="contig", start_col="start", end_col="end"):
    """Map fragment coordinates back to the unfragmented sequence."""
    if df.empty:
        df = df.copy()
        for c in ("source", "source_start", "source_end", "fragment_length"):
            df[c] = []
        return df
    parsed = df[contig_col].map(F.parse_fragment_name)
    df = df.copy()
    df["source"] = [p[0] for p in parsed]
    off = np.array([p[1] for p in parsed])
    df["source_start"] = df[start_col].values + off - 1
    df["source_end"] = df[end_col].values + off - 1
    df["fragment_length"] = [p[2] - p[1] + 1 for p in parsed]
    return df


def fragment_bounds(cuts, pos, length):
    if cuts is None or len(cuts) == 0:
        return 1, length
    i = np.searchsorted(cuts, pos, side="left")       # first cut >= pos
    start = int(cuts[i - 1]) + 1 if i > 0 else 1
    end = int(cuts[i]) if i < len(cuts) else length
    return start, end


def score_replicate(rec, tr, lengths, unfrag_cands):
    rid, acc = rec["replicate_id"], rec["accession"]
    base = {k: rec[k] for k in ("replicate_id", "accession", "model", "target_n50", "replicate", "realised_n50",
                                "realised_l50", "n_contigs")}
    rows, gained = [], []
    if rec["status"] == "unperturbed":
        for t in tr.to_dict("records"):
            r = dict(base, contig=t["contig"], start=t["start"], end=t["end"], strand=t["strand"],
                     status="retained_exact", window_cut=False, geometric_loss_predicted=False,
                     synteny_checkable=F.synteny_checkable(t, 1, lengths[t["contig"]]),
                     fragment_length=lengths[t["contig"]], trace_step="", reported_cdeii_len=t["cdeii_len"],
                     truth_cdeii_len=t["cdeii_len"])
            rows.append(r)
        return rows, gained
    cuts = read_cuts(rid)
    calls = to_source(pd.read_csv(P.data_dir("arm1", "calls", rid + ".tsv"), sep="\t"))
    cand_path = P.data_dir("arm1", "calls", rid + ".tsv.candidates.tsv.gz")
    cands = to_source(pd.read_csv(cand_path, sep="\t"), "Contig", "Start", "End")
    used = set()
    for t in tr.to_dict("records"):
        cp = cuts.get(t["contig"], np.array([], dtype=np.int64))
        fs, fe = fragment_bounds(cp, t["start"], lengths[t["contig"]])
        r = dict(base, contig=t["contig"], start=t["start"], end=t["end"], strand=t["strand"],
                 window_cut=F.window_cut(t, cp), geometric_loss_predicted=F.geometric_loss(t, cp),
                 synteny_checkable=F.synteny_checkable(t, fs, fe), fragment_length=fe - fs + 1,
                 truth_cdeii_len=t["cdeii_len"], trace_step="", reported_cdeii_len="")
        on = calls[(calls["source"] == t["contig"]) & (calls["source_start"] <= t["end"])
                   & (calls["source_end"] >= t["start"])]
        if len(on):
            ex = on[(on["source_start"] == t["start"]) & (on["source_end"] == t["end"])
                    & (on["cdeii_len"] == t["cdeii_len"])]
            hit = ex.iloc[0] if len(ex) else on.iloc[0]
            used.update(on.index.tolist())
            r["status"] = "retained_exact" if len(ex) else "retained_altered"
            r["reported_cdeii_len"] = int(hit["cdeii_len"])
        else:
            r["status"] = "lost_geometric" if r["geometric_loss_predicted"] else "lost_pipeline"
            c = cands[(cands["source"] == t["contig"]) & (cands["source_start"] == t["start"])
                      & (cands["source_end"] == t["end"])]
            r["trace_step"] = c.iloc[0]["removed_at"] if len(c) else "no candidate formed"
        rows.append(r)
    for i, c in calls.iterrows():
        if i in used:
            continue
        u = unfrag_cands[(unfrag_cands["Contig"] == c["source"]) & (unfrag_cands["Start"] == c["source_start"])
                         & (unfrag_cands["End"] == c["source_end"])]
        gained.append(dict(base, source=c["source"], source_start=int(c["source_start"]),
                           source_end=int(c["source_end"]), strand=c["strand"], cdeii_len=int(c["cdeii_len"]),
                           cdeii_at=c["cdeii_at"], fragment_length=int(c["fragment_length"]),
                           unfragmented_fate=u.iloc[0]["removed_at"] if len(u) else "not a candidate in the intact assembly",
                           unfragmented_rank=int(u.iloc[0]["rank"]) if len(u) else ""))
    return rows, gained


def q(x, p):
    return float(np.percentile(x, p)) if len(x) else float("nan")


def summarise_1a(calls, reps):
    per = calls.groupby("replicate_id").agg(
        accession=("accession", "first"), model=("model", "first"), target_n50=("target_n50", "first"),
        realised_n50=("realised_n50", "first"), n_truth=("status", "size"),
        retained=("status", lambda s: s.str.startswith("retained").sum()),
        retained_exact=("status", lambda s: (s == "retained_exact").sum()),
        lost_geometric=("status", lambda s: (s == "lost_geometric").sum()),
        lost_pipeline=("status", lambda s: (s == "lost_pipeline").sum()),
        predicted_lost=("geometric_loss_predicted", "sum"), window_cut=("window_cut", "sum"),
        synteny_checkable=("synteny_checkable", "sum")).reset_index()
    per["recall"] = per["retained"] / per["n_truth"]
    per["exact_recall"] = per["retained_exact"] / per["n_truth"]
    per["geometric_expectation"] = 1 - per["predicted_lost"] / per["n_truth"]
    per["synteny_checkable_fraction"] = per["synteny_checkable"] / per["n_truth"]
    g = reps.set_index("replicate_id")
    per["gained"] = per["replicate_id"].map(lambda r: g.loc[r, "gained"] if r in g.index else 0)
    per["realised_l50"] = per["replicate_id"].map(g["realised_l50"])
    per["n_contigs"] = per["replicate_id"].map(g["n_contigs"])

    rows = []
    for (acc, model, lv), d in per.groupby(["accession", "model", "target_n50"]):
        rows.append({"accession": acc, "model": model, "target_n50": lv, "replicates": len(d),
                     "realised_n50_median": q(d["realised_n50"], 50), "realised_n50_q25": q(d["realised_n50"], 25),
                     "realised_n50_q75": q(d["realised_n50"], 75), "realised_l50_median": q(d["realised_l50"], 50),
                     "n_contigs_median": q(d["n_contigs"], 50),
                     "recall_median": q(d["recall"], 50), "recall_q25": q(d["recall"], 25),
                     "recall_q75": q(d["recall"], 75), "recall_min": d["recall"].min(), "recall_max": d["recall"].max(),
                     "exact_recall_median": q(d["exact_recall"], 50),
                     "geometric_expectation_median": q(d["geometric_expectation"], 50),
                     "synteny_checkable_median": q(d["synteny_checkable_fraction"], 50),
                     "lost_geometric_total": int(d["lost_geometric"].sum()),
                     "lost_pipeline_total": int(d["lost_pipeline"].sum()),
                     "gained_total": int(d["gained"].sum()), "truth_calls_per_replicate": int(d["n_truth"].iloc[0])})
    by_genome = pd.DataFrame(rows)

    pooled = []
    for (model, lv), d in per.groupby(["model", "target_n50"]):
        units = [x["recall"].values for _, x in d.groupby("accession")]
        est, lo, hi = P.cluster_bootstrap(units, lambda us: float(np.mean([u.mean() for u in us])),
                                          n_boot=N_BOOT, seed=SEED)
        gunits = [x["geometric_expectation"].values for _, x in d.groupby("accession")]
        gest, glo, ghi = P.cluster_bootstrap(gunits, lambda us: float(np.mean([u.mean() for u in us])),
                                             n_boot=N_BOOT, seed=SEED)
        pooled.append({"model": model, "target_n50": lv, "genomes": len(units), "replicates": len(d),
                       "realised_n50_median": q(d["realised_n50"], 50),
                       "recall_median": q(d["recall"], 50), "recall_q25": q(d["recall"], 25),
                       "recall_q75": q(d["recall"], 75), "mean_genome_recall": est, "mean_genome_recall_ci_low": lo,
                       "mean_genome_recall_ci_high": hi, "geometric_expectation_mean": gest,
                       "geometric_expectation_ci_low": glo, "geometric_expectation_ci_high": ghi,
                       "synteny_checkable_median": q(d["synteny_checkable_fraction"], 50),
                       "lost_geometric_total": int(d["lost_geometric"].sum()),
                       "lost_pipeline_total": int(d["lost_pipeline"].sum()), "gained_total": int(d["gained"].sum()),
                       "interval": "percentile bootstrap over genomes, %d resamples, seed %d" % (N_BOOT, SEED)})
    return per, by_genome, pd.DataFrame(pooled)


def score_1a():
    runs = pd.read_csv(os.path.join(OUT, "fragmentation_runs.tsv"), sep="\t")
    all_rows, all_gained, rep_rows = [], [], []
    for acc, rr in runs.groupby("accession"):
        tr = truth(acc)
        lengths = seq_lengths(acc)
        unfrag = pd.read_csv(P.data_dir("arm0", "calls", acc + ".tsv.candidates.tsv.gz"), sep="\t")
        for rec in rr.to_dict("records"):
            if rec["status"] not in ("ok", "unperturbed"):
                rep_rows.append({"replicate_id": rec["replicate_id"], "scored": "no", "gained": 0,
                                 "realised_l50": rec["realised_l50"], "n_contigs": rec["n_contigs"]})
                continue
            rows, gained = score_replicate(rec, tr, lengths, unfrag)
            all_rows.extend(rows)
            all_gained.extend(gained)
            rep_rows.append({"replicate_id": rec["replicate_id"], "scored": "yes", "gained": len(gained),
                             "realised_l50": rec["realised_l50"], "n_contigs": rec["n_contigs"]})
        P.log("scored %s: %d replicates" % (acc, len(rr)))
    calls = pd.DataFrame(all_rows)
    gained = pd.DataFrame(all_gained)
    reps = pd.DataFrame(rep_rows)
    P.atomic_write_tsv(calls, os.path.join(OUT, "fragmentation_calls.tsv"))
    P.atomic_write_tsv(gained, os.path.join(OUT, "fragmentation_gained.tsv"))
    main_models = calls[calls["model"].isin(["uniform", "at_weighted"])]
    per, by_genome, pooled = summarise_1a(main_models, reps)
    P.atomic_write_tsv(per, os.path.join(OUT, "fragmentation_replicates.tsv"))
    P.atomic_write_tsv(by_genome, os.path.join(OUT, "recall_by_genome.tsv"))
    P.atomic_write_tsv(pooled, os.path.join(OUT, "recall_pooled.tsv"))

    # Which filter removed the calls lost with their window intact.
    lp = main_models[main_models["status"] == "lost_pipeline"]
    dec = (lp.groupby(["model", "target_n50", "trace_step"]).size().rename("lost_calls").reset_index()
           if len(lp) else pd.DataFrame(columns=["model", "target_n50", "trace_step", "lost_calls"]))
    P.atomic_write_tsv(dec, os.path.join(OUT, "pipeline_losses_by_filter.tsv"))
    gd = (gained[gained["model"].isin(["uniform", "at_weighted"])].groupby(["model", "target_n50", "unfragmented_fate"])
          .size().rename("gained_calls").reset_index()) if len(gained) else pd.DataFrame(
        columns=["model", "target_n50", "unfragmented_fate", "gained_calls"])
    P.atomic_write_tsv(gd, os.path.join(OUT, "gained_calls_by_origin.tsv"))

    # Sensitivity of the AT-weighted model to gamma, against gamma's main value.
    sens = calls[calls["model"].isin(["at_weighted", "at_weighted_g5", "at_weighted_g20"])
                 & calls["target_n50"].isin([50000, 20000, 10000])]
    if len(sens):
        sp, _, spool = summarise_1a(sens, reps)
        spool = spool.copy()
        spool["gamma"] = spool["model"].map({"at_weighted": "ln(10)/0.30", "at_weighted_g5": "ln(5)/0.30",
                                             "at_weighted_g20": "ln(20)/0.30"})
        P.atomic_write_tsv(spool, P.repo("results", "sensitivity", "arm1_at_weight_gamma.tsv"))
    return per, pooled


def score_1b():
    edits = pd.read_csv(os.path.join(OUT, "variant_edits.tsv"), sep="\t")
    runs = pd.read_csv(os.path.join(OUT, "variant_runs.tsv"), sep="\t").set_index("run_id")
    single = []
    for e in edits[edits["design"] == "single"].to_dict("records"):
        r = dict(e)
        if not isinstance(e.get("status"), str) or not e["status"].startswith("not possible"):
            path = P.data_dir("arm1", "variant_calls", e["run_id"] + ".tsv")
            st = runs.loc[e["run_id"], "status"] if e["run_id"] in runs.index else "missing"
            if st != "ok" or not os.path.exists(path):
                r["status"] = "PCAn failed"
            else:
                c = pd.read_csv(path, sep="\t")
                on = c[(c["contig"] == e["contig"]) & (c["start"] <= e["expected_end"]) & (c["end"] >= e["expected_start"])]
                r["called"] = bool(len(on))
                r["reported_cdeii_len"] = int(on.iloc[0]["cdeii_len"]) if len(on) else ""
                r["correct_length"] = bool(len(on)) and int(on.iloc[0]["cdeii_len"]) == int(e["expected_cdeii_len"])
                r["reported_unedited_length"] = bool(len(on)) and int(on.iloc[0]["cdeii_len"]) == int(e["unedited_cdeii_len"])
                r["calls_in_genome"] = len(c)
                if not len(on):
                    cand = pd.read_csv(path + ".candidates.tsv.gz", sep="\t")
                    m = cand[(cand["Contig"] == e["contig"]) & (cand["Start"] == e["expected_start"])
                             & (cand["End"] == e["expected_end"])]
                    r["trace_step"] = m.iloc[0]["removed_at"] if len(m) else "no candidate at the planted span"
                r["status"] = "scored"
        single.append(r)
    single = pd.DataFrame(single)
    P.atomic_write_tsv(single, os.path.join(OUT, "planted_single.tsv"))

    sc = single[single["status"] == "scored"].copy()
    env = []
    for (acc, size), d in sc.groupby(["accession", "size"]):
        k, kc, n = int(d["called"].sum()), int(d["correct_length"].sum()), len(d)
        lo, hi = P.clopper_pearson(k, n)
        clo, chi = P.clopper_pearson(kc, n)
        env.append({"accession": acc, "size": size, "edits": n, "called": k, "called_fraction": k / n,
                    "called_ci_low": lo, "called_ci_high": hi, "correct_length": kc,
                    "correct_length_fraction": kc / n, "correct_length_ci_low": clo, "correct_length_ci_high": chi,
                    "interval": "Clopper-Pearson over centromeres, within genome"})
    for size, d in sc.groupby("size"):
        units = [x for _, x in d.groupby("accession")]
        est, lo, hi = P.cluster_bootstrap(units, lambda us: sum(u["called"].sum() for u in us) / sum(len(u) for u in us),
                                          n_boot=N_BOOT, seed=SEED)
        cest, clo, chi = P.cluster_bootstrap(units, lambda us: sum(u["correct_length"].sum() for u in us)
                                             / sum(len(u) for u in us), n_boot=N_BOOT, seed=SEED)
        env.append({"accession": "pooled", "size": size, "edits": len(d), "called": int(d["called"].sum()),
                    "called_fraction": est, "called_ci_low": lo, "called_ci_high": hi,
                    "correct_length": int(d["correct_length"].sum()), "correct_length_fraction": cest,
                    "correct_length_ci_low": clo, "correct_length_ci_high": chi,
                    "interval": "percentile bootstrap over %d genomes, %d resamples, seed %d; rough with three units"
                                % (len(units), N_BOOT, SEED)})
    P.atomic_write_tsv(pd.DataFrame(env), os.path.join(OUT, "detection_envelope.tsv"))

    prog = []
    for rid, d in edits[edits["design"] == "progressive"].groupby("run_id"):
        path = P.data_dir("arm1", "variant_calls", rid + ".tsv")
        st = runs.loc[rid, "status"] if rid in runs.index else "missing"
        row = {"run_id": rid, "k": int(d["k"].iloc[0]), "replicate": int(d["replicate"].iloc[0]), "status": st}
        if st == "ok" and os.path.exists(path):
            c = pd.read_csv(path, sep="\t")
            called = correct = 0
            for e in d.to_dict("records"):
                on = c[(c["contig"] == e["contig"]) & (c["start"] <= e["expected_end"]) & (c["end"] >= e["expected_start"])]
                called += bool(len(on))
                correct += bool(len(on)) and int(on.iloc[0]["cdeii_len"]) == int(e["expected_cdeii_len"])
            tr = truth("GCA_000146045.2")
            edited = set(d["contig"])
            unedited = tr[~tr["contig"].isin(edited)]
            kept = sum(((c["contig"] == u.contig) & (c["start"] == u.start) & (c["end"] == u.end)).any()
                       for u in unedited.itertuples())
            row.update(total_calls=len(c), edited_called=called, edited_correct_length=correct,
                       unedited=len(unedited), unedited_unchanged=int(kept),
                       median_top5_cdeii_len=runs.loc[rid, "median_top5_cdeii_len"],
                       final_median_cdeii_len=runs.loc[rid, "final_median_cdeii_len"])
        prog.append(row)
    prog = pd.DataFrame(prog).sort_values(["k", "replicate"])
    P.atomic_write_tsv(prog, os.path.join(OUT, "planted_progressive.tsv"))
    return single, prog


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--part", choices=["1a", "1b", "both"], default="both")
    args = ap.parse_args()
    if P.is_done("score_perturbations") and not args.force:
        P.log("arm 1 scoring already complete, skipping (use --force to redo)")
        return
    if args.part in ("1a", "both"):
        per, pooled = score_1a()
        P.log("arm 1a scored: %d replicates" % len(per))
    if args.part in ("1b", "both"):
        single, prog = score_1b()
        P.log("arm 1b scored: %d single edits, %d progressive runs" % (len(single), len(prog)))
    if args.part == "both":
        P.mark_done("score_perturbations")


if __name__ == "__main__":
    main()
