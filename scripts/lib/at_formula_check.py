"""Arm 0 diagnostic, run inside PCAn's own environment (pandas 2.0.3, NumPy
1.24.4) so that tied scores sort exactly as PCAn sorts them.

The CDEII AT percentage that Supplementary Data 2 reports is the AT count over
the CDEII length, CDEII being every base between CDEI and the CDEIII motif.
PCAn v1.0 scores candidates with a different value: the AT count of a stretch
shifted one base towards CDEIII, over the length plus one. If the published
calls had been made with the first formula, rescoring every candidate with it
and replaying PCAn's filters should recover published calls that v1.0 misses.
This replays each arm 0 assembly's candidate table three ways and counts the
published calls left standing.

Usage (called by scripts/05_reproduce.py):
    <pcan env>/bin/python scripts/lib/at_formula_check.py <per_centromere.tsv> <calls dir> <out.tsv>
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pcan_trace as T  # noqa: E402

FORMULAS = {
    "pcan_v1": "as PCAn v1.0 computes it (shifted stretch, count over length plus one)",
    "cdeii_count_over_length": "every base between CDEI and CDEIII, count over length (the published column)",
    "cdeii_count_over_length_plus_one": "every base between CDEI and CDEIII, count over length plus one",
}


def rescore(cands, how):
    c = cands.drop(columns=["rank", "removed_at"]).copy()
    seq = c["Sequence"].astype(str)
    if how != "pcan_v1":
        core = seq.str.slice(8, -26)
        denom = core.str.len() + (1 if how.endswith("plus_one") else 0)
        c["CDEIIAT"] = (core.str.count("[AaTt]") / denom * 100).round(2)
    c["OverallScore"] = (c["FIMOIscore"] * c["FIMOIIscore"] / (1 - c["CDEIIAT"] / 100)).round(0)
    c = c.sort_values(by="OverallScore", ascending=False, ignore_index=True)
    _, survivors, _ = T.replay_filters(c)
    return survivors


def main():
    per_path, calls_dir, out = sys.argv[1:4]
    per = pd.read_csv(per_path, sep="\t")
    rows = []
    for acc, p in per.groupby("accession"):
        path = os.path.join(calls_dir, acc + ".tsv.candidates.tsv.gz")
        if not os.path.exists(path):
            continue
        cands = pd.read_csv(path, sep="\t")
        located = p.dropna(subset=["start"])
        for how in FORMULAS:
            surv = rescore(cands, how)
            n = 0
            for r in located.itertuples():
                n += bool(((surv["Contig"] == r.contig) & (surv["Start"] == int(r.start))
                           & (surv["End"] == int(r.end))).any())
            rows.append({"accession": acc, "formula": how, "published_calls": len(p),
                         "published_calls_standing": n, "calls_standing": len(surv)})
    df = pd.DataFrame(rows)
    tot = df.groupby("formula")[["published_calls", "published_calls_standing", "calls_standing"]].sum().reset_index()
    tot["accession"] = "all"
    df = pd.concat([df, tot], ignore_index=True)
    df["formula_description"] = df["formula"].map(FORMULAS)
    df.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
