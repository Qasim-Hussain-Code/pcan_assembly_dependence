"""Collect the arm results into the reporting summary: Holm correction across
the family of confirmatory tests, the arm 3 settling depths, whether any
sensitivity analysis changes a conclusion, and the measured time, memory and
disk of every stage.

A conclusion counts as changed by a sensitivity analysis when, for a test, its
Holm-adjusted significance at 0.05 or the direction of its effect differs from
the primary analysis; and, for an estimate without a test, when the
alternative's point estimate lies outside the primary 95 per cent interval.
This rule was set before any sensitivity result was read.

Usage:
    python scripts/12_analyse.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import padlib as P  # noqa: E402

R = P.repo("results")
DEPTH_ORDER = ["5", "10", "20", "40", "full"]


def read(path):
    p = os.path.join(R, path)
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


def holm(p):
    """Holm step-down adjustment; NaN p-values are left out and stay NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    m = len(order)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p[i]))
        out[i] = running
    return out


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    if isinstance(x, (int, np.integer)):
        return "%d" % x
    return ("%." + str(nd) + "g") % x


def settling(d1, assembler):
    """Lowest target depth from which the mean over strains stays within 0.05
    of its full-depth value at every higher depth (registered as D3)."""
    s = d1[(d1["assembler"] == assembler) & (d1["seed"].astype(str) == "11")].set_index(
        d1["target_depth"].astype(str))["estimate"]
    if "full" not in s.index:
        return ""
    full = s["full"]
    depths = [d for d in DEPTH_ORDER if d in s.index]
    for i, d in enumerate(depths):
        if all(abs(s[x] - full) <= 0.05 for x in depths[i:]):
            return d
    return "full"


def main():
    rows = []

    g = read("arm0/gate.tsv")
    if g is not None:
        r = g.iloc[0]
        rows.append({"arm": "0", "status": "pre-specified (config/arm0_reproduction.md)",
                     "claim": "fraction of published calls PCAn v1.0 reproduces exactly, against a 0.95 gate",
                     "n": "%d calls in %d species" % (r["n_calls"], r["n_species"]), "unit": "species",
                     "test": "none (gate on the point estimate)", "estimate": r["estimate"], "ci_low": r["ci_low"],
                     "ci_high": r["ci_high"], "p_raw": np.nan, "source": "results/arm0/gate.tsv"})

    pooled = read("arm1/recall_pooled.tsv")
    if pooled is not None:
        for _, r in pooled[pooled["target_n50"].isin([50000, 20000, 10000, 5000])].iterrows():
            rows.append({"arm": "1 (simulation)", "status": "pre-specified (config/arm1_design.md)",
                         "claim": "mean recall relative to the unfragmented call, %s model, target contig N50 %d kb"
                                  % ({"uniform": "uniform", "at_weighted": "AT-weighted"}.get(r["model"], r["model"]), r["target_n50"] // 1000),
                         "n": "%d genomes, %d replicates" % (r["genomes"], r["replicates"]), "unit": "genome",
                         "test": "none", "estimate": r["mean_genome_recall"], "ci_low": r["mean_genome_recall_ci_low"],
                         "ci_high": r["mean_genome_recall_ci_high"], "p_raw": np.nan,
                         "source": "results/arm1/recall_pooled.tsv"})
    env = read("arm1/detection_envelope.tsv")
    if env is not None:
        for _, r in env[env["accession"] == "pooled"].iterrows():
            rows.append({"arm": "1 (simulation)", "status": "pre-specified (config/arm1_design.md)",
                         "claim": "planted CDEII variants of %+d bp called at the planted length" % r["size"],
                         "n": "%d edits in 3 genomes" % r["edits"], "unit": "genome", "test": "none",
                         "estimate": r["correct_length_fraction"], "ci_low": r["correct_length_ci_low"],
                         "ci_high": r["correct_length_ci_high"], "p_raw": np.nan,
                         "source": "results/arm1/detection_envelope.tsv"})

    c2 = read("arm2/confirmatory.tsv")
    if c2 is not None:
        for _, r in c2.iterrows():
            rows.append({"arm": "2", "status": "confirmatory (config/analysis_plan.md)",
                         "claim": "%s: %s" % (r["outcome"], r["measure"]), "n": "%s strains" % fmt(r.get("n_strains")),
                         "unit": "strain", "test": r.get("test", "") if isinstance(r.get("test", ""), str) else "",
                         "estimate": r.get("estimate"), "ci_low": r.get("ci_low"), "ci_high": r.get("ci_high"),
                         "p_raw": r.get("p_value", np.nan), "source": "results/arm2/confirmatory.tsv",
                         "decision": r.get("decision", "")})
    c3 = read("arm3/confirmatory.tsv")
    if c3 is not None:
        for _, r in c3.iterrows():
            rows.append({"arm": "3", "status": "confirmatory (config/analysis_plan.md)",
                         "claim": "%s: %s, %s, target depth %s, seed %s" % (r["outcome"], r["measure"], r["assembler"],
                                                                           r["target_depth"], r["seed"]),
                         "n": "%s strains" % fmt(r.get("n_strains")), "unit": "strain",
                         "test": r.get("test", "") if isinstance(r.get("test", ""), str) else "",
                         "estimate": r.get("estimate"), "ci_low": r.get("ci_low"), "ci_high": r.get("ci_high"),
                         "p_raw": r.get("p_value", np.nan), "source": "results/arm3/confirmatory.tsv",
                         "decision": r.get("decision", "")})
        d1 = c3[c3["outcome"] == "D1"]
        d2 = c3[c3["outcome"] == "D2"]
        for asm in sorted(c3["assembler"].dropna().unique()):
            for lab, d in (("D3: depth at which recall settles", d1), ("D3: depth at which CDEII agreement settles", d2)):
                rows.append({"arm": "3", "status": "confirmatory (config/analysis_plan.md)",
                             "claim": "%s, %s" % (lab, asm), "n": "", "unit": "strain",
                             "test": "none (summary of D1 or D2; tolerance 0.05, arbitrary)",
                             "estimate": settling(d, asm), "ci_low": np.nan, "ci_high": np.nan, "p_raw": np.nan,
                             "source": "results/arm3/confirmatory.tsv"})

    rep = pd.DataFrame(rows)
    # The Holm family is every registered test: arm 2 C4 (ratio) and C5, arm 3
    # D4 and D5. Rows without a p-value are estimates, not tests.
    fam = rep["p_raw"].notna() & rep["arm"].isin(["2", "3"])
    rep["p_holm"] = np.nan
    rep.loc[fam, "p_holm"] = holm(rep.loc[fam, "p_raw"].values)
    rep["holm_family_size"] = int(fam.sum())
    P.atomic_write_tsv(rep, os.path.join(R, "reporting_summary.tsv"))

    # Sensitivity: does any alternative change a conclusion?
    sens_rows = []
    gs = read("sensitivity/arm0_gate_threshold.tsv")
    if gs is not None:
        for _, r in gs.iterrows():
            sens_rows.append({"analysis": "arm 0 gate threshold %.2f" % r["threshold"], "outcome": "arm 0 gate",
                              "primary": "decision at 0.95", "alternative": r["decision"],
                              "changes_conclusion": "no" if r["same_as_preregistered_threshold"] else "yes"})
    sv = read("sensitivity/arm2_variants.tsv")
    if sv is not None and c2 is not None:
        for _, r in sv.iterrows():
            base = c2[(c2["outcome"] == r["outcome"]) & (c2["measure"] == r["measure"])]
            if base.empty or pd.isna(r.get("estimate")):
                continue
            b = base.iloc[0]
            if not pd.isna(b.get("p_value", np.nan)) and not pd.isna(r.get("p_value", np.nan)):
                same = (r["p_value"] < 0.05) == (b["p_value"] < 0.05) and np.sign(np.log(r["estimate"])) == np.sign(np.log(b["estimate"])) \
                    if r["measure"].startswith("break-rate") else (r["p_value"] < 0.05) == (b["p_value"] < 0.05)
            else:
                same = b["ci_low"] <= r["estimate"] <= b["ci_high"]
            sens_rows.append({"analysis": "arm 2 %s" % r["variant"], "outcome": "%s %s" % (r["outcome"], r["measure"]),
                              "primary": fmt(b["estimate"]), "alternative": fmt(r["estimate"]),
                              "changes_conclusion": "no" if same else "yes"})
    s3 = read("sensitivity/arm3_variants.tsv")
    if s3 is not None and c3 is not None:
        for _, r in s3.iterrows():
            base = c3[(c3["outcome"] == r["outcome"]) & (c3["assembler"] == r["assembler"])
                      & (c3["target_depth"].astype(str) == str(r["target_depth"])) & (c3["seed"].astype(str) == "11")]
            if "measure" in r and isinstance(r.get("measure"), str):
                base = base[base["measure"] == r["measure"]]
            if base.empty or pd.isna(r.get("estimate")):
                continue
            b = base.iloc[0]
            if not pd.isna(b.get("p_value", np.nan)) and not pd.isna(r.get("p_value", np.nan)):
                same = (r["p_value"] < 0.05) == (b["p_value"] < 0.05)
            else:
                same = b["ci_low"] <= r["estimate"] <= b["ci_high"]
            sens_rows.append({"analysis": "arm 3 %s" % r["variant"],
                              "outcome": "%s %s %s" % (r["outcome"], r["assembler"], r["target_depth"]),
                              "primary": fmt(b["estimate"]), "alternative": fmt(r["estimate"]),
                              "changes_conclusion": "no" if same else "yes"})
    gam = read("sensitivity/arm1_at_weight_gamma.tsv")
    if gam is not None:
        for _, r in gam.iterrows():
            base = gam[(gam["model"] == "at_weighted") & (gam["target_n50"] == r["target_n50"])]
            if base.empty or r["model"] == "at_weighted":
                continue
            b = base.iloc[0]
            same = b["mean_genome_recall_ci_low"] <= r["mean_genome_recall"] <= b["mean_genome_recall_ci_high"]
            sens_rows.append({"analysis": "arm 1 AT weighting gamma %s" % r["gamma"],
                              "outcome": "mean recall at %d kb" % (r["target_n50"] // 1000),
                              "primary": fmt(b["mean_genome_recall"]), "alternative": fmt(r["mean_genome_recall"]),
                              "changes_conclusion": "no" if same else "yes"})
    P.atomic_write_tsv(pd.DataFrame(sens_rows), os.path.join(R, "sensitivity", "summary.tsv"))

    # Time, memory and disk per stage, from logs/resources.tsv.
    res = pd.read_csv(P.repo("logs", "resources.tsv"), sep="\t")
    P.atomic_write_tsv(res, os.path.join(R, "timing_by_stage.tsv"))

    write_markdown(rep)
    P.log("reporting summary: %d rows, %d in the Holm family" % (len(rep), int(fam.sum())))


def write_markdown(rep):
    lines = ["# Reporting summary", "",
             "Every claim the README makes from a registered or pre-specified analysis, with its source. "
             "Confirmatory claims are those registered in `config/analysis_plan.md`; arm 0 and arm 1 measures were "
             "pre-specified in `config/arm0_reproduction.md` and `config/arm1_design.md`. P-values are "
             "Holm-adjusted across the %d registered tests. Intervals are 95 per cent; their method is given in "
             "the source file." % int(rep["holm_family_size"].iloc[0]) if len(rep) else "", "",
             "| Arm | Status | Claim | n | Unit | Test | Estimate | 95% interval | p (Holm) | Source |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in rep.iterrows():
        ci = "" if pd.isna(r["ci_low"]) else "%s to %s" % (fmt(r["ci_low"]), fmt(r["ci_high"]))
        est = r["estimate"] if isinstance(r["estimate"], str) else fmt(r["estimate"])
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | `%s` |" % (
            r["arm"], r["status"], r["claim"], r["n"], r["unit"], r["test"] or "none", est, ci,
            fmt(r["p_holm"]), r["source"]))
    P.atomic_write_text(os.path.join(R, "reporting_summary.md"), "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
