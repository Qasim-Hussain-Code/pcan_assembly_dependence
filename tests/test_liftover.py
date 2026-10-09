"""Known-answer tests for the null-window sampler, the breakpoint locator and
the breakpoint likelihood used in arms 2 and 3."""
import math
import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "lib"))
import liftover as L  # noqa: E402


def contig(n, at_frac, rng):
    p = at_frac / 2
    return "".join(rng.choice(list("ACGT"), size=n, p=[p, (1 - at_frac) / 2, (1 - at_frac) / 2, p]))


class TestNullWindows(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        # contig a: GC-rich except one AT-rich stretch; contig b: AT-rich throughout
        a = contig(30000, 0.3, rng) + contig(5000, 0.9, rng) + contig(30000, 0.3, rng)
        b = contig(40000, 0.9, rng)
        b = b[:20000] + "N" * 50 + b[20050:]
        self.seqs = {"a": a, "b": b}
        self.idx = L.AtIndex(self.seqs)
        # one PCAn call in the AT-rich stretch of contig a
        self.calls = pd.DataFrame([{"contig": "a", "start": 32001, "end": 32120, "cdeii_len": 86}])
        self.cens = L.centromere_elements(self.calls, self.seqs, 2000)

    def test_constraints(self):
        nulls, short = L.null_windows("s", self.cens, self.idx, self.calls, flank=2000, half_bin=0.025, per_cen=20)
        self.assertEqual(len(nulls), 20)
        c = self.cens.iloc[0]
        for r in nulls.itertuples():
            seq = self.seqs[r.contig][r.start - 1:r.end]
            self.assertEqual(len(seq), c["length"])
            self.assertNotIn("N", seq)
            at = (seq.count("A") + seq.count("T")) / len(seq)
            self.assertLessEqual(abs(at - c["at"]), 0.025 + 1e-9)
            self.assertGreaterEqual(r.start - 1, 2000)
            self.assertLessEqual(r.end + 2000, len(self.seqs[r.contig]))
            if r.contig == "a":
                # at least 10 kb from the call on either side
                self.assertTrue(r.end < 32001 - 10000 or r.start > 32120 + 10000)

    def test_seeded(self):
        a, _ = L.null_windows("s", self.cens, self.idx, self.calls, flank=2000)
        b, _ = L.null_windows("s", self.cens, self.idx, self.calls, flank=2000)
        self.assertEqual(list(a["start"]), list(b["start"]))
        c, _ = L.null_windows("other", self.cens, self.idx, self.calls, flank=2000)
        self.assertNotEqual(list(a["start"]), list(c["start"]))

    def test_shortfall_is_reported(self):
        nulls, short = L.null_windows("s", self.cens, self.idx, self.calls, flank=2000, half_bin=0.0001, per_cen=20)
        self.assertLessEqual(len(nulls), 20)
        if len(nulls) < 20:
            self.assertEqual(len(short), 1)


class TestBreakpoints(unittest.TestCase):
    def test_contig_ends_located(self):
        # query contig q1 aligns end to end, forward, to long[10000:30000];
        # q2 aligns reverse to long[40000:41500] but leaves 300 bp unaligned
        # at its start, so that end is not located
        paf = ("q1\t20000\t0\t20000\t+\tlong\t100000\t10000\t30000\t19900\t20000\t60\ttp:A:P\n"
               "q2\t1800\t300\t1800\t-\tlong\t100000\t40000\t41500\t1490\t1500\t60\ttp:A:P\n"
               "q3\t5000\t0\t5000\t+\tlong\t100000\t500\t5500\t4990\t5000\t60\ttp:A:P\n")
        with tempfile.NamedTemporaryFile("w", suffix=".paf", delete=False) as fh:
            fh.write(paf)
        bps, unlocated = L.breakpoints(fh.name, {"long": 100000})
        os.remove(fh.name)
        got = sorted(zip(bps["query_contig"], bps["contig_end"], bps["position"]))
        # q1: start at 10000, end at 30000. q2 (reverse): its end maps to tstart
        # 40000; its start is 300 bp unaligned. q3 starts within 1 kb of the
        # long-read contig end (500), so that end is dropped; its end at 5500 stays.
        self.assertEqual(got, [("q1", "end", 30000), ("q1", "start", 10000), ("q2", "end", 40000),
                               ("q3", "end", 5500)])
        self.assertEqual(unlocated, 1)


class TestLikelihood(unittest.TestCase):
    def test_summary_equals_direct_sum(self):
        rng = np.random.default_rng(4)
        w = pd.DataFrame({"length": rng.choice([500, 500, 500, 230], size=300),
                          "at": rng.integers(150, 451, size=300) / 500})
        hits = w.sample(40, random_state=2)
        for g in (0.0, 2.5, -4.0, L.GAMMA_ARM1):
            abar = np.average(w["at"], weights=w["length"])
            logw = np.log(w["length"]) + g * (w["at"] - abar)
            direct = float((np.log(hits["length"]) + g * (hits["at"] - abar)).sum()
                           - len(hits) * math.log(np.exp(logw).sum()))
            self.assertAlmostEqual(L.BreakSummary(w, hits).loglik(g), direct, places=8)

    def test_gamma_recovered(self):
        # breakpoints drawn from the AT-weighted model with gamma 7.7 should
        # give a maximum-likelihood gamma near 7.7
        rng = np.random.default_rng(8)
        w = pd.DataFrame({"length": np.full(20000, 500), "at": rng.uniform(0.3, 0.9, 20000).round(3)})
        abar = w["at"].mean()
        p = np.exp(7.7 * (w["at"] - abar))
        p /= p.sum()
        hits = w.iloc[rng.choice(len(w), size=4000, p=p)]
        g, _ = L.fit_gamma([L.BreakSummary(w, hits)])
        self.assertAlmostEqual(g, 7.7, delta=0.5)


if __name__ == "__main__":
    unittest.main()
