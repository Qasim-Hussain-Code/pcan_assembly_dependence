"""Known-answer tests for the interval, correction and break-ratio code that
the confirmatory outcomes rest on. Each expected value is worked out by hand
or is a published textbook value, not taken from the code under test."""
import importlib.util
import math
import os
import sys
import unittest

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts", "lib"))
import padlib as P  # noqa: E402
import pairstats as PS  # noqa: E402

spec = importlib.util.spec_from_file_location("analyse", os.path.join(HERE, "..", "scripts", "12_analyse.py"))
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


class TestIntervals(unittest.TestCase):
    def test_wilson_half(self):
        # 5 of 10: the Wilson 95 per cent interval is 0.2366 to 0.7634
        lo, hi = P.wilson(5, 10)
        self.assertAlmostEqual(lo, 0.2366, places=4)
        self.assertAlmostEqual(hi, 0.7634, places=4)

    def test_clopper_pearson_zero(self):
        # 0 of 10: the exact upper bound solves (1 - p)^10 = 0.025, p = 0.3085
        lo, hi = P.clopper_pearson(0, 10)
        self.assertEqual(lo, 0.0)
        self.assertAlmostEqual(hi, 1 - 0.025 ** (1 / 10), places=4)


class TestHolm(unittest.TestCase):
    def test_three_tests(self):
        # sorted 0.01, 0.03, 0.04: 3 x 0.01 = 0.03; 2 x 0.03 = 0.06; 1 x 0.04,
        # raised to the running maximum 0.06; returned in the input order
        self.assertTrue(np.allclose(A.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06]))

    def test_fixed_family_with_tests_missing(self):
        # two of fourteen tests run: the smallest is multiplied by 14 and the
        # next by 13, as if the twelve missing tests had p = 1
        out = A.holm([0.001, 0.002], 14)
        self.assertTrue(np.allclose(out, [0.014, 0.026]))

    def test_nan_left_out(self):
        out = A.holm([0.01, float("nan"), 0.02], 3)
        self.assertTrue(math.isnan(out[1]))
        self.assertTrue(np.allclose(out[[0, 2]], [0.03, 0.04]))

    def test_substitution(self):
        # replacing test 1's p-value and adjusting again
        self.assertAlmostEqual(A.holm_with([0.01, 0.04, 0.03], 1, 0.001, 3), 0.003)


class TestBreakRatio(unittest.TestCase):
    def test_abs_log_ratio(self):
        # 2 of 10 broken against 4 of 40: ratio 2, |log 2|
        self.assertAlmostEqual(PS.abs_log_ratio(2, 10, 4, 40), math.log(2))
        self.assertEqual(PS.abs_log_ratio(0, 10, 4, 40), math.inf)
        self.assertEqual(PS.abs_log_ratio(2, 10, 0, 40), math.inf)
        self.assertIsNone(PS.abs_log_ratio(0, 10, 0, 40))

    def test_counts_match_tables(self):
        # two strains; strain a: 1 of 2 centromeres and 1 of 4 nulls broken,
        # strain b: 0 of 1 and 2 of 2. Pooled ratio (1/3) / (3/6) = 2/3.
        rows = [("a", "centromere", "split"), ("a", "centromere", "intact_called"),
                ("a", "null", "absent"), ("a", "null", "intact_uncalled"), ("a", "null", "intact_uncalled"),
                ("a", "null", "intact_uncalled"), ("b", "centromere", "intact_called"),
                ("b", "null", "n_run"), ("b", "null", "split")]
        el = pd.DataFrame(rows, columns=["strain", "kind", "status"])
        old_boot, old_perm = PS.N_BOOT, PS.N_PERM
        PS.N_BOOT, PS.N_PERM = 200, 200
        try:
            ratio, diff = PS.c4_breaks(el[el["kind"] == "centromere"], el[el["kind"] == "null"])
        finally:
            PS.N_BOOT, PS.N_PERM = old_boot, old_perm
        self.assertAlmostEqual(ratio["estimate"], (1 / 3) / (3 / 6))
        self.assertAlmostEqual(diff["estimate"], 1 / 3 - 3 / 6)
        self.assertEqual((ratio["k"], ratio["n_centromeres"], ratio["k_null"], ratio["n_null"]), (1, 3, 3, 6))

    def test_strain_counts_sorted(self):
        df = pd.DataFrame({"strain": ["b", "a", "b", "a", "a"]})
        self.assertEqual(PS.strain_counts(df, [1, 0, 1, 1, 0]), [(1, 3), (2, 2)])
        self.assertAlmostEqual(PS.pooled([(1, 3), (2, 2)]), 3 / 5)


if __name__ == "__main__":
    unittest.main()
