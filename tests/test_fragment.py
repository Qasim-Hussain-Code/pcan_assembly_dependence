"""Known-answer tests for the fragmenter and for the geometry of a PCAn call."""
import math
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "lib"))
import fragment as F  # noqa: E402
import padlib as P  # noqa: E402


class TestCuts(unittest.TestCase):
    def test_apply_cuts_known_fragments(self):
        recs = [("chr1", "AAAACCCCGGGG")]
        frags = F.apply_cuts(recs, {"chr1": [4, 8]})
        self.assertEqual([f[1] for f in frags], ["AAAA", "CCCC", "GGGG"])
        self.assertEqual([f[0] for f in frags], ["chr1__1-4", "chr1__5-8", "chr1__9-12"])
        self.assertEqual("".join(f[1] for f in frags), recs[0][1])

    def test_no_cuts_returns_the_sequence(self):
        frags = F.apply_cuts([("s", "ACGT")], {})
        self.assertEqual(frags, [("s__1-4", "ACGT", "s", 1)])

    def test_names_map_back(self):
        self.assertEqual(F.parse_fragment_name("NC_001133.9__150001-162345"), ("NC_001133.9", 150001, 162345))
        # underscores inside the source name survive
        self.assertEqual(F.parse_fragment_name("JAJLYE010000002.1__1-500"), ("JAJLYE010000002.1", 1, 500))

    def test_every_base_kept_after_random_cuts(self):
        rng = np.random.default_rng(1)
        seq = "".join(rng.choice(list("ACGT"), size=50000))
        cuts = F.draw_cuts([("x", seq)], 1 / 1000, "uniform", np.random.default_rng(7))
        frags = F.apply_cuts([("x", seq)], cuts)
        self.assertEqual("".join(f[1] for f in frags), seq)
        for name, s, src, start in frags:
            src2, a, b = F.parse_fragment_name(name)
            self.assertEqual(seq[a - 1:b], s)
            self.assertEqual(a, start)

    def test_seed_reproduces_cuts(self):
        seq = "ACGT" * 10000
        a = F.draw_cuts([("x", seq)], 1 / 500, "uniform", np.random.default_rng(42))
        b = F.draw_cuts([("x", seq)], 1 / 500, "uniform", np.random.default_rng(42))
        self.assertTrue(np.array_equal(a["x"], b["x"]))

    def test_uniform_rate(self):
        # 1,600 draws of a 100 kb sequence at one cut per 1 kb: about 100 cuts
        # per draw, so the mean over draws has a standard error of 0.25 per
        # cent and must lie within 1 per cent of 100
        seq = "ACGT" * 25000
        rng = np.random.default_rng(3)
        n = [len(F.draw_cuts([("x", seq)], 1 / 1000, "uniform", rng)["x"]) for _ in range(1600)]
        self.assertAlmostEqual(np.mean(n) / (len(seq) - 1) * 1000, 1.0, delta=0.01)

    def test_at_weighted_prefers_at_rich_sequence(self):
        # half the sequence is 90 per cent AT, half 30 per cent AT; genome AT is
        # 60 per cent. With gamma = ln(10)/0.30 each AT-rich window carries
        # weight exp(gamma * 0.30) = 10 and each GC-rich window exp(-gamma * 0.30)
        # = 0.1, so the expected share of cuts in the AT-rich half is 10/10.1.
        rng = np.random.default_rng(5)
        at_half = "".join(rng.choice(list("ACGT"), size=50000, p=[0.45, 0.05, 0.05, 0.45]))
        gc_half = "".join(rng.choice(list("ACGT"), size=50000, p=[0.15, 0.35, 0.35, 0.15]))
        seq = at_half + gc_half
        genome_at = P.at_fraction(seq)
        cuts = np.concatenate([F.draw_cuts([("x", seq)], 1 / 200, "at_weighted", rng, genome_at)["x"]
                               for _ in range(50)])
        share = (cuts <= 50000).mean()
        expect = math.exp(F.GAMMA_DEFAULT * (P.at_fraction(at_half) - genome_at))
        expect = expect / (expect + math.exp(F.GAMMA_DEFAULT * (P.at_fraction(gc_half) - genome_at)))
        self.assertAlmostEqual(share, expect, delta=0.01)

    def test_at_weighted_keeps_the_genome_wide_rate(self):
        rng = np.random.default_rng(9)
        seq = "".join(rng.choice(list("ACGT"), size=100000, p=[0.35, 0.15, 0.15, 0.35]))
        n = [len(F.draw_cuts([("x", seq)], 1 / 1000, "at_weighted", rng, P.at_fraction(seq))["x"])
             for _ in range(200)]
        self.assertAlmostEqual(np.mean(n) / 100, 1.0, delta=0.02)

    def test_fast_contig_lengths_equal_slow_ones(self):
        # random sequence with N runs, including at both ends and back to back
        rng = np.random.default_rng(13)
        parts = []
        for _ in range(40):
            parts.append("".join(rng.choice(list("ACGT"), size=int(rng.integers(1, 3000)))))
            parts.append("N" * int(rng.integers(1, 50)))
        seq = "NN" + "".join(parts) + "ACGT"
        recs = [("s", seq)]
        for rate in (1 / 200, 1 / 2000, 1 / 20000):
            cuts = F.draw_cuts(recs, rate, "uniform", np.random.default_rng(int(1 / rate)))
            slow = sorted(F.fragment_lengths_split_at_gaps(F.apply_cuts(recs, cuts)))
            fast = sorted(F.genome_contig_lengths({"s": F.gap_free_intervals(seq)}, cuts))
            self.assertEqual(slow, fast)

    def test_calibration_reaches_target(self):
        rng = np.random.default_rng(11)
        recs = [("c%d" % i, "".join(rng.choice(list("ACGT"), size=400000))) for i in range(4)]
        rate, n50 = F.calibrate_rate(recs, 20000, "uniform", seed=123)
        self.assertGreater(rate, 0)
        self.assertLess(abs(n50 - 20000) / 20000, 0.03)

    def test_calibration_returns_zero_above_the_original(self):
        recs = [("c", "A" * 1000)]
        rate, n50 = F.calibrate_rate(recs, 5000, "uniform", seed=1)
        self.assertEqual(rate, 0.0)
        self.assertEqual(n50, 1000)


class TestCallGeometry(unittest.TestCase):
    # A forward-strand call as PCAn reports it on S288C chromosome IX:
    # CDEI 355629-355636, CDEIII motif 355721-355746.
    fwd = {"strand": "+", "start": 355629, "end": 355746, "cdei_start": 355629, "cdei_end": 355636,
           "cdeiii_start": 355721, "cdeiii_end": 355746}
    # A reverse-strand call on chromosome X: motif 436306-436331, CDEI 436418-436425.
    rev = {"strand": "-", "start": 436306, "end": 436425, "cdei_start": 436418, "cdei_end": 436425,
           "cdeiii_start": 436306, "cdeiii_end": 436331}

    def test_window_coordinates(self):
        self.assertEqual(F.extraction_window(self.fwd), (355721 - 249, 355746))
        self.assertEqual(F.extraction_window(self.fwd)[1] - F.extraction_window(self.fwd)[0] + 1, 275)
        self.assertEqual(F.extraction_window(self.rev), (436306, 436331 + 250))
        self.assertEqual(F.extraction_window(self.rev)[1] - F.extraction_window(self.rev)[0] + 1, 276)

    def test_forward_cut_upstream_of_cdei_inside_window_loses_the_call(self):
        self.assertTrue(F.geometric_loss(self.fwd, [355721 - 249]))
        self.assertTrue(F.geometric_loss(self.fwd, [355600]))

    def test_forward_cut_just_outside_window_keeps_the_call(self):
        self.assertFalse(F.geometric_loss(self.fwd, [355721 - 250]))
        self.assertFalse(F.geometric_loss(self.fwd, [355746]))

    def test_reverse_cut_beyond_cdei_only_shortens_the_window(self):
        self.assertTrue(F.window_cut(self.rev, [436500]))
        self.assertFalse(F.geometric_loss(self.rev, [436500]))
        self.assertFalse(F.geometric_loss(self.rev, [436425]))

    def test_reverse_cut_inside_the_element_loses_the_call(self):
        self.assertTrue(F.geometric_loss(self.rev, [436306]))
        self.assertTrue(F.geometric_loss(self.rev, [436424]))

    def test_synteny_flanks(self):
        self.assertTrue(F.synteny_checkable(self.fwd, 345629, 365746))
        self.assertFalse(F.synteny_checkable(self.fwd, 345630, 365746))
        self.assertFalse(F.synteny_checkable(self.fwd, 345629, 365745))


if __name__ == "__main__":
    unittest.main()
