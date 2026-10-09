"""Statistics for the registered outcomes of arms 2 and 3
(config/analysis_plan.md). Every interval resamples strains, the independent
unit; nothing here resamples centromeres."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

import liftover as L
import padlib as P
import status as S

SEED, N_BOOT, N_PERM = 20261009, 10000, 10000


def c1_counts(counts):
    """counts: one row per strain with long_calls and short_calls."""
    k = int((counts["long_calls"] == counts["short_calls"]).sum())
    n = len(counts)
    lo, hi = P.wilson(k, n) if n >= 10 and k >= 10 and n - k >= 10 else P.clopper_pearson(k, n)
    method = "Wilson" if n >= 10 and k >= 10 and n - k >= 10 else "Clopper-Pearson"
    return {"outcome": "C1", "measure": "strains with equal call counts", "n_strains": n, "k": k,
            "estimate": k / n if n else float("nan"), "ci_low": lo, "ci_high": hi, "interval": method + " over strains"}


def by_strain(df):
    return [g for _, g in df.groupby("strain")]


def strain_counts(df, hit):
    """Per strain, (elements with hit, elements), strains in the sorted order
    of by_strain. A bootstrap over these tuples draws exactly the strains a
    bootstrap over the per-strain tables would, at a fraction of the cost."""
    g = pd.DataFrame({"strain": df["strain"].values, "hit": np.asarray(hit, dtype=int)}).groupby("strain")["hit"]
    return list(zip(g.sum().tolist(), g.size().tolist()))


def pooled(units):
    return sum(k for k, _ in units) / sum(n for _, n in units)


def c2_status(cens):
    """Fraction of liftable long-read centromeres in each status, pooled."""
    lift = cens[cens["status"] != "not_liftable"]
    out = []
    for st in S.STATUSES:
        units = strain_counts(lift, lift["status"] == st)
        est, lo, hi = P.cluster_bootstrap(units, pooled, n_boot=N_BOOT, seed=SEED)
        out.append({"outcome": "C2", "measure": "fraction %s" % st, "n_strains": len(units),
                    "n_centromeres": len(lift), "k": int((lift["status"] == st).sum()), "estimate": est,
                    "ci_low": lo, "ci_high": hi, "interval": "percentile bootstrap over strains"})
    return out


def c3_cdeii(cens):
    called = cens[cens["status"] == "intact_called"]
    units = strain_counts(called, called["cdeii_difference"] == 0)
    if not units:
        return {"outcome": "C3", "measure": "identical CDEII length among intact_called", "n_strains": 0}
    est, lo, hi = P.cluster_bootstrap(units, pooled, n_boot=N_BOOT, seed=SEED)
    return {"outcome": "C3", "measure": "identical CDEII length among intact_called", "n_strains": len(units),
            "n_centromeres": len(called), "k": int((called["cdeii_difference"] == 0).sum()), "estimate": est,
            "ci_low": lo, "ci_high": hi, "interval": "percentile bootstrap over strains"}


def _sums(units):
    """units: per strain (broken centromeres, centromeres, broken nulls, nulls)."""
    return [sum(u[i] for u in units) for i in range(4)]


def _ratio(units):
    kc, nc, k0, n0 = _sums(units)
    if nc == 0 or n0 == 0 or k0 == 0:
        return float("nan")
    return (kc / nc) / (k0 / n0)


def _diff(units):
    kc, nc, k0, n0 = _sums(units)
    return kc / nc - k0 / n0 if nc and n0 else float("nan")


def abs_log_ratio(kc, nc, k0, n0):
    """|log R| for pooled counts. Infinite when one of the two groups has no
    break, so that such a split counts as at least as extreme as any other;
    None when neither has one, where R is undefined."""
    if kc == 0 and k0 == 0:
        return None
    if kc == 0 or k0 == 0:
        return math.inf
    return abs(math.log((kc / nc) / (k0 / n0)))


def c4_breaks(cens, nulls, label="C4"):
    """Break-rate ratio of centromeres against their AT-matched null windows,
    with a bootstrap over strains and a permutation test within strains."""
    el = pd.concat([cens[cens["status"] != "not_liftable"][["strain", "kind", "status"]],
                    nulls[["strain", "kind", "status"]]], ignore_index=True)
    if not set(el["kind"]) <= {"centromere", "null"}:
        raise ValueError("element kinds other than centromere and null: %s" % sorted(set(el["kind"].astype(str))))
    el["broken"] = el["status"].isin(S.BROKEN).astype(int)
    cen = (el["kind"] == "centromere").values
    b = el["broken"].values
    g = pd.DataFrame({"strain": el["strain"].values, "kc": b * cen, "nc": cen.astype(int), "k0": b * ~cen,
                      "n0": (~cen).astype(int)}).groupby("strain").sum()
    counts = [tuple(int(x) for x in r) for r in g[["kc", "nc", "k0", "n0"]].itertuples(index=False, name=None)]
    est, lo, hi = P.cluster_bootstrap(counts, _ratio, n_boot=N_BOOT, seed=SEED)
    dest, dlo, dhi = P.cluster_bootstrap(counts, _diff, n_boot=N_BOOT, seed=SEED)
    # permutation: within each strain, which elements carry the centromere
    # label is shuffled; the number of centromeres per strain is kept
    units = by_strain(el)
    rng = np.random.default_rng(SEED)
    kc, nc, k0, n0 = _sums(counts)
    obs = abs_log_ratio(kc, nc, k0, n0)
    arrays = [(u["broken"].values, int((u["kind"] == "centromere").sum())) for u in units]
    hits = 0
    valid = 0
    if obs is not None:
        for _ in range(N_PERM):
            pc = p0 = 0
            for bb, m in arrays:
                idx = rng.permutation(len(bb))
                pc += bb[idx[:m]].sum()
                p0 += bb[idx[m:]].sum()
            r = abs_log_ratio(pc, nc, p0, n0)
            if r is None:
                continue
            valid += 1
            if r >= obs - 1e-12:
                hits += 1
    p = (hits + 1) / (valid + 1) if valid else float("nan")
    return [{"outcome": label, "measure": "break-rate ratio, centromere over AT-matched null", "n_strains": len(units),
             "n_centromeres": nc, "k": kc, "n_null": n0, "k_null": k0, "estimate": est, "ci_low": lo, "ci_high": hi,
             "p_value": p, "permutations_used": valid,
             "test": "stratified permutation of log ratio, %d permutations" % N_PERM,
             "interval": "percentile bootstrap over strains"},
            {"outcome": label, "measure": "break proportion, centromere minus null", "n_strains": len(units),
             "n_centromeres": nc, "k": kc, "n_null": n0, "k_null": k0, "estimate": dest, "ci_low": dlo,
             "ci_high": dhi, "interval": "percentile bootstrap over strains"}]


def c5_model(strain_sets, label="C5"):
    """strain_sets: {strain: (windows, hits)}. Registered decision rule: the
    AT-weighted model is favoured if its log-likelihood at arm 1's gamma is
    the higher of the two and the bootstrap interval of the ML gamma excludes
    zero; the uniform model if its log-likelihood is higher and the interval
    includes zero; otherwise neither."""
    from scipy.stats import chi2
    sets = [L.BreakSummary(w, h) for w, h in strain_sets.values() if len(h)]
    if not sets:
        return {"outcome": label, "measure": "breakpoint model", "n_strains": 0}
    ll0 = sum(s.loglik(0.0) for s in sets)
    ll1 = sum(s.loglik(L.GAMMA_ARM1) for s in sets)
    g_hat, ll_hat = L.fit_gamma(sets)
    rng = np.random.default_rng(SEED)
    boots = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(sets), len(sets))
        boots.append(L.fit_gamma([sets[i] for i in idx])[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    lr = 2 * (ll_hat - ll0)
    p = float(chi2.sf(max(lr, 0), 1))
    if ll1 > ll0 and lo > 0:
        decision = "AT-weighted"
    elif ll0 >= ll1 and lo <= 0 <= hi:
        decision = "uniform"
    else:
        decision = "neither"
    return {"outcome": label, "measure": "maximum-likelihood gamma of breakpoint placement", "n_strains": len(sets),
            "n_breakpoints": int(sum(s.n for s in sets)), "estimate": g_hat, "ci_low": float(lo),
            "ci_high": float(hi), "loglik_uniform": ll0, "loglik_arm1_at_weighted": ll1, "loglik_ml": ll_hat,
            "p_value": p, "test": "likelihood ratio, gamma = 0 against free gamma, chi-squared 1 df",
            "interval": "percentile bootstrap over strains, %d resamples" % len(boots), "decision": decision}
