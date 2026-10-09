"""Locate long-read centromeres and null windows in another assembly of the
same strain, and find where that assembly breaks. Used by arms 2 and 3 through
scripts/10_compare_pairs.py, with the rules registered in
config/analysis_plan.md.
"""
from __future__ import annotations

import hashlib
import math
import os
import subprocess

import numpy as np
import pandas as pd

import padlib as P
import status as S

FLANK = 2000
IDENTITY, QCOV, MAPQ = 0.95, 0.90, 20
MAX_EXTRA, MIN_GAP, SLACK = 1000, -50, 10
NULL_PER_CEN, AT_HALF_BIN, NULL_CALL_DISTANCE = 20, 0.025, 10000
BREAK_MIN_ALN, BREAK_END_TOL, BREAK_CONTIG_END = 1000, 100, 1000
WINDOW = 500
GAMMA_ARM1 = math.log(10) / 0.30


def seed_of(*parts):
    return int(hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()[:8], 16)


# ----------------------------------------------------------------------------
# elements: the long-read centromeres, and the null windows drawn for them

def centromere_elements(calls, seqs, flank=FLANK):
    """One row per long-read PCAn call. A call with less than `flank` bases on
    either side on its contig cannot be lifted over and is marked so."""
    rows = []
    for c in calls.to_dict("records"):
        s, e, n = int(c["start"]), int(c["end"]), len(seqs[c["contig"]])
        el = seqs[c["contig"]][s - 1:e]
        rows.append({"element_id": "cen|%s|%d" % (c["contig"], s), "kind": "centromere", "contig": c["contig"],
                     "start": s, "end": e, "length": e - s + 1, "at": P.at_fraction(el),
                     "liftable": s - flank >= 1 and e + flank <= n, "cdeii_len": int(c["cdeii_len"])})
    return pd.DataFrame(rows)


class AtIndex:
    """Prefix sums of AT and N along every contig (32-bit, about 100 MB for a
    12 Mb genome), so that the AT count of any window is one subtraction."""

    def __init__(self, seqs):
        self.seqs = seqs
        self.at, self.nn = {}, {}
        for name, s in seqs.items():
            b = np.frombuffer(s.upper().encode(), dtype=np.uint8)
            self.at[name] = np.concatenate(([0], np.cumsum((b == 65) | (b == 84), dtype=np.int32)))
            self.nn[name] = np.concatenate(([0], np.cumsum(b == 78, dtype=np.int32)))

    def candidates(self, name, length, at_lo, at_hi, flank, exclude):
        """0-based starts of windows of `length` on one contig whose AT count
        lies in [at_lo, at_hi], with no N, at least `flank` from both contig
        ends, and outside every (lo, hi) start range in `exclude`."""
        n = len(self.seqs[name])
        if n < length + 2 * flank:
            return np.array([], dtype=np.int64), np.array([], dtype=np.int32)
        at, nn = self.at[name], self.nn[name]
        a = at[length:] - at[:-length]
        ok = (a >= at_lo) & (a <= at_hi)
        ok &= (nn[length:] - nn[:-length]) == 0
        ok[:flank] = False
        ok[n - length - flank + 1:] = False
        for lo, hi in exclude:
            ok[max(lo + 1, 0):max(hi, 0)] = False
        idx = np.nonzero(ok)[0]
        return idx, a[idx]


def null_windows(strain, cens, idx, calls, flank=FLANK, half_bin=AT_HALF_BIN, per_cen=NULL_PER_CEN):
    """Null windows for each liftable centromere: same length, AT within
    half_bin of the centromere's, no N, at least NULL_CALL_DISTANCE from any
    PCAn call, at least `flank` from a contig end. Drawn uniformly without
    replacement with a seed from the strain and centromere. Works one
    centromere at a time, so memory stays near the size of the prefix sums."""
    rows, short = [], []
    excl = {}
    for k in calls.itertuples():
        excl.setdefault(k.contig, []).append((k.start - 1 - NULL_CALL_DISTANCE - 0, k.end + NULL_CALL_DISTANCE))
    for c in cens[cens["liftable"]].to_dict("records"):
        L = int(c["length"])
        # the AT band as whole counts; the small tolerance keeps windows that
        # sit exactly on the band edge from being lost to rounding
        lo_ct = int(np.ceil((c["at"] - half_bin) * L - 1e-9))
        hi_ct = int(np.floor((c["at"] + half_bin) * L + 1e-9))
        pool_names, pool_starts, pool_at = [], [], []
        for name in idx.seqs:
            ex = [(lo - L, hi) for lo, hi in excl.get(name, [])]
            st, a = idx.candidates(name, L, lo_ct, hi_ct, flank, ex)
            if len(st):
                pool_names.append(np.full(len(st), name, dtype=object))
                pool_starts.append(st)
                pool_at.append(a / L)
        names = np.concatenate(pool_names) if pool_names else np.array([], dtype=object)
        starts = np.concatenate(pool_starts) if pool_starts else np.array([], dtype=np.int64)
        ats = np.concatenate(pool_at) if pool_at else np.array([])
        rng = np.random.default_rng(seed_of(strain, c["contig"], c["start"], "null"))
        k = min(per_cen, len(starts))
        if k < per_cen:
            short.append({"centromere": c["element_id"], "qualifying_windows": len(starts)})
        pick = rng.choice(len(starts), size=k, replace=False) if k else []
        for i in pick:
            s0 = int(starts[i])
            rows.append({"element_id": "null|%s|%s|%d" % (c["element_id"], names[i], s0 + 1), "kind": "null",
                         "for_centromere": c["element_id"], "contig": names[i], "start": s0 + 1, "end": s0 + L,
                         "length": L, "at": float(ats[i]), "liftable": True, "cdeii_len": ""})
    return pd.DataFrame(rows), pd.DataFrame(short)


# ----------------------------------------------------------------------------
# placing flanks and classifying elements

def flank_fasta(elements, seqs, path, flank=FLANK):
    with open(path, "w") as fh:
        for e in elements[elements["liftable"]].to_dict("records"):
            s, t = e["start"], e["end"]
            seq = seqs[e["contig"]]
            fh.write(">%s|L\n%s\n>%s|R\n%s\n" % (e["element_id"], seq[s - 1 - flank:s - 1], e["element_id"], seq[t:t + flank]))


def run_minimap2(target, query, out_paf, preset="asm10", threads=1):
    with open(out_paf, "w") as fh:
        subprocess.run(["minimap2", "-x", preset, "-c", "-t", str(threads), target, query], check=True,
                       stdout=fh, stderr=subprocess.DEVNULL)


def classify_elements(elements, paf_path, target_seqs, calls_by_seq, identity=IDENTITY, qcov=QCOV, mapq=MAPQ):
    with open(paf_path) as fh:
        alns = S.parse_paf(fh)
    by_q = {}
    for a in alns:
        by_q.setdefault(a.qname, []).append(a)
    out = []
    for e in elements.to_dict("records"):
        r = dict(e)
        if not e["liftable"]:
            r.update(status="not_liftable", detail="less than the flank length on the long-read contig")
            out.append(r)
            continue
        left = S.place_flank(by_q.get(e["element_id"] + "|L", []), identity, qcov, mapq)
        right = S.place_flank(by_q.get(e["element_id"] + "|R", []), identity, qcov, mapq)
        st, d = S.classify(left, right, e["length"], lambda n: target_seqs[n], calls_by_seq,
                           MAX_EXTRA, MIN_GAP, SLACK)
        r["status"] = st
        r["detail"] = "; ".join("%s=%s" % kv for kv in d.items())
        r["target"] = d.get("target", "")
        r["region_start"], r["region_end"] = d.get("region_start", ""), d.get("region_end", "")
        r["call_start"], r["call_end"] = d.get("call_start", ""), d.get("call_end", "")
        out.append(r)
    return pd.DataFrame(out)


# ----------------------------------------------------------------------------
# where an assembly breaks, mapped onto the long-read assembly

def breakpoints(paf_path, long_lengths):
    """Each contig end of the query assembly, located on the long-read
    assembly through the query contig's alignments (primary, mapping quality
    at least 20, at least 1 kb). An end is located when the outermost
    alignment reaches within 100 bp of it; ends within 1 kb of a long-read
    contig end are not breaks of the query assembly and are dropped."""
    with open(paf_path) as fh:
        alns = [a for a in S.parse_paf(fh) if a.primary and a.mapq >= MAPQ and a.alen >= BREAK_MIN_ALN]
    by_q = {}
    for a in alns:
        by_q.setdefault(a.qname, []).append(a)
    rows, unlocated = [], 0
    for q, al in by_q.items():
        al.sort(key=lambda a: a.qstart)
        for end, a in (("start", al[0]), ("end", al[-1])):
            if end == "start" and a.qstart > BREAK_END_TOL:
                unlocated += 1
                continue
            if end == "end" and a.qlen - a.qend > BREAK_END_TOL:
                unlocated += 1
                continue
            if end == "start":
                t = a.tstart if a.strand == "+" else a.tend
            else:
                t = a.tend if a.strand == "+" else a.tstart
            if t < BREAK_CONTIG_END or t > long_lengths[a.tname] - BREAK_CONTIG_END:
                continue
            rows.append({"query_contig": q, "contig_end": end, "target": a.tname, "position": int(t)})
    return pd.DataFrame(rows), unlocated


def window_table(seqs):
    """500 bp windows of the long-read assembly with their length and AT."""
    rows = []
    for name, s in seqs.items():
        u = s.upper()
        for a in range(0, len(u), WINDOW):
            w = u[a:a + WINDOW]
            at = w.count("A") + w.count("T")
            acgt = at + w.count("C") + w.count("G")
            rows.append((name, a // WINDOW, len(w), at / acgt if acgt else np.nan))
    df = pd.DataFrame(rows, columns=["contig", "window", "length", "at"])
    df["at"] = df["at"].fillna(df["at"].mean())
    return df


class BreakSummary:
    """Everything the breakpoint likelihood needs from one strain, reduced
    exactly: the windows grouped by AT value (a full 500 bp window can take at
    most 501 values) with their summed lengths, and three sums over the
    strain's breakpoints. Evaluating the likelihood is then a short sum, which
    is what makes 10,000 bootstrap refits affordable."""

    def __init__(self, windows, hits):
        self.abar = float(np.average(windows["at"], weights=windows["length"]))
        g = windows.groupby("at")["length"].sum()
        self.u = g.index.values.astype(float) - self.abar
        self.w = g.values.astype(float)
        self.n = len(hits)
        self.sa = float((hits["at"].values - self.abar).sum()) if self.n else 0.0
        self.sl = float(np.log(hits["length"].values).sum()) if self.n else 0.0

    def loglik(self, gamma):
        z = np.log(self.w) + gamma * self.u
        m = z.max()
        return gamma * self.sa + self.sl - self.n * (m + math.log(np.exp(z - m).sum()))


def break_loglik(gamma, windows, hits):
    """Log-likelihood of one strain's breakpoints under P(window) proportional
    to length x exp(gamma (a - a_genome))."""
    return BreakSummary(windows, hits).loglik(gamma)


def fit_gamma(summaries):
    """Maximum-likelihood gamma shared by all strains, from BreakSummary objects."""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda g: -sum(s.loglik(g) for s in summaries), bounds=(-30, 60), method="bounded",
                          options={"xatol": 1e-4})
    return float(res.x), -float(res.fun)
