"""Replay PCAn v1.0's candidate building and filters from the two FIMO passes of
one run, and record the step at which each candidate is removed.

This does not make calls. PCAn makes the calls; this explains them. Every
statement in build_candidates and replay_filters mirrors a line of
AutomatedCENretrieval.py at commit a5fa46f, cited by line number, including
its quirks: the forward-strand window starts 250 bp before the motif start
(Python slicing makes it 249 bp of upstream sequence), the CDEII used for the
AT percentage skips the first base after CDEI and includes the first base of
the CDEIII motif, and the AT percentage divides by the CDEII length plus one.
It runs inside PCAn's own environment (pandas 2.0.3, NumPy 1.24.4) so that
pandas breaks ties in OverallScore exactly as PCAn's sort does.

scripts/lib/pcan_runner.py runs the replay after every PCAn run and checks
that the candidates it leaves standing are PCAn's calls, row for row. A run
where they differ is flagged in its metadata; none of the attribution in this
repository rests on such a run.
"""
import csv
import re

import pandas as pd

STEPS = ["not_in_top50", "cdeii_length_outside_anchor_window", "cdeii_at_70_or_less",
         "negative_motif_score", "duplicate_cdeiii_hit", "duplicate_sequence",
         "not_best_on_contig", "length_and_at_outlier", "called"]

_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def _rc(s):
    return s.translate(_COMP)[::-1]


def build_candidates(genome, pass1_path, pass2_path, cdeiiiup=250):
    """genome: {sequence id: sequence}. Mirrors FIMO2(), lines 150 to 234."""
    queries = pd.read_csv(pass1_path, skiprows=0, delimiter="\t")
    list_seq = {}
    for index, row in queries.iterrows():
        ident, strand = row.iloc[1], row.iloc[4]
        start, end, qs1 = int(row.iloc[2]), int(row.iloc[3]), row.iloc[5]
        if strand == "+":                       # lines 168 to 171
            new_start, new_end = start - cdeiiiup, end
            position1 = new_end
            sequence = genome[ident][new_start:new_end]
        else:                                   # lines 172 to 175
            new_start, new_end = start - 1, end + cdeiiiup
            position1 = new_start
            sequence = _rc(genome[ident][new_start:new_end])
        name = ident + "_" + str(index + 1)     # line 185
        list_seq[name] = [name, sequence, qs1, position1, strand]

    with open(pass2_path) as fh:                # lines 200 to 203
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        fm_hits = list(reader)
    full = []
    for fm in fm_hits:                          # lines 207 to 228
        start, end = int(fm[2]), int(fm[3])
        hit = list_seq.get(fm[1])
        if hit is None:
            continue
        contig = re.search(".+(?=_.+)", fm[1]).group(0)
        short = hit[1][start - 1:]
        at_seq = hit[1][end + 1:len(hit[1]) - 25]
        at_len = len(at_seq)
        at_per = round((at_seq.count("A") + at_seq.count("a") + at_seq.count("T") + at_seq.count("t"))
                       / (at_len + 1) * 100, 2)
        qs2 = float(fm[5])
        overall = round(hit[2] * qs2 * 1 / (1 - at_per / 100), 0)
        if hit[4] == "+":
            ns, ne = hit[3] - len(short) + 1, hit[3]
        else:
            ns, ne = hit[3] + 1, hit[3] + len(short)
        full.append([contig, hit[0], ns, ne, short, str(at_len), str(at_per), str(hit[2]), str(qs2),
                     str(overall)])
    cols = ["Contig", "ContigHit", "Start", "End", "Sequence", "CDEIIlen", "CDEIIAT", "FIMOIscore",
            "FIMOIIscore", "OverallScore"]
    df = pd.DataFrame(full, columns=cols)       # lines 230 to 233
    df = df.astype({"Start": "int", "End": "int", "CDEIIlen": "int", "CDEIIAT": "float",
                    "FIMOIscore": "float", "FIMOIIscore": "float", "OverallScore": "float"})
    return df.sort_values(by="OverallScore", ascending=False, ignore_index=True)


def replay_filters(full, nrows=50):
    """Mirrors Filtering(), lines 255 to 293, and labels every candidate with
    the first step that removes it. Adding the rank column does not change any
    of PCAn's operations, which all act on named columns."""
    full = full.copy()
    full["rank"] = range(1, len(full) + 1)
    removed = {}

    def drop(before, after, step):
        for r in set(before["rank"]) - set(after["rank"]):
            removed.setdefault(r, step)

    median_top5 = full.head(5)["CDEIIlen"].median()
    head = full.head(nrows)
    drop(full, head, "not_in_top50")
    s1 = head.loc[(head["CDEIIlen"] > median_top5 - 30) & (head["CDEIIlen"] < median_top5 + 30)]
    drop(head, s1, "cdeii_length_outside_anchor_window")
    s2 = s1.loc[s1["CDEIIAT"] > 70]
    drop(s1, s2, "cdeii_at_70_or_less")
    s3 = s2.loc[(s2["FIMOIscore"] > 0) & (s2["FIMOIIscore"] > 0)]
    drop(s2, s3, "negative_motif_score")
    s4 = s3.drop_duplicates(subset="ContigHit", keep="first")
    drop(s3, s4, "duplicate_cdeiii_hit")
    s5 = s4.drop_duplicates(subset="Sequence", keep="first")
    drop(s4, s5, "duplicate_sequence")
    s6 = s5.drop_duplicates(subset="Contig", keep="first")
    drop(s5, s6, "not_best_on_contig")
    medlen = s6["CDEIIlen"].median()
    medat = s6["CDEIIAT"].quantile(0.5)
    s7 = s6.loc[((s6["CDEIIlen"] > medlen - 10) & (s6["CDEIIlen"] < medlen + 10)) | (s6["CDEIIAT"] > medat - 7)]
    drop(s6, s7, "length_and_at_outlier")
    full["removed_at"] = [removed.get(r, "called") for r in full["rank"]]
    anchors = {"median_top5_cdeii_len": median_top5, "final_median_cdeii_len": medlen,
               "final_median_cdeii_at": medat}
    return full, s7, anchors


def trace(genome, pass1_path, pass2_path, pcan_rows):
    """Replay one run and compare the survivors with PCAn's own call rows
    (dicts read from CENsequences.txt). Returns (candidates, anchors, matches)."""
    full = build_candidates(genome, pass1_path, pass2_path)
    labelled, survivors, anchors = replay_filters(full)
    mine = [(r.ContigHit, int(r.Start), int(r.End), int(r.CDEIIlen), r.Sequence) for r in survivors.itertuples()]
    theirs = [(r["ContigHit"], int(r["Start"]), int(r["End"]), int(r["CDEIIlen"]), r["Sequence"]) for r in pcan_rows]
    return labelled, anchors, mine == theirs
