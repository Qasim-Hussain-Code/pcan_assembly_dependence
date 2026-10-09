"""Known-answer tests for the CDEII variant editor."""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "lib"))
import variants as V  # noqa: E402

# A toy chromosome: 10 bases of flank, CDEI (8), CDEII (20), CDEIII (26), flank.
FLANK_L = "GGGGGCCCCC"
CDEI = "TCACGTGA"
CDEII = "AATTAATTAACCTTAATTAA"
CDEIII = "TGTTTATGTTTTCCGAAATAAAAAAA"
FLANK_R = "CCCCCGGGGG"
SEQ = FLANK_L + CDEI + CDEII + CDEIII + FLANK_R
# Forward-strand call in 1-based coordinates
FWD = {"strand": "+", "start": 11, "end": 64, "cdei_start": 11, "cdei_end": 18,
       "cdeiii_start": 39, "cdeiii_end": 64}


class TestEditor(unittest.TestCase):
    def test_cdeii_interval(self):
        self.assertEqual(V.cdeii_interval(FWD), (19, 38))
        self.assertEqual(SEQ[19 - 1:38], CDEII)

    def test_insertion_is_a_tandem_copy(self):
        # duplicate the 5 bases ending at position 30 (CDEII bases 8 to 12,
        # "TTAAC"), inserting the copy after base 30
        new, rec = V.apply_edit(SEQ, FWD, 5, 30)
        self.assertEqual(rec["changed_bases"], SEQ[25:30])
        self.assertEqual(new, SEQ[:30] + SEQ[25:30] + SEQ[30:])
        self.assertEqual(len(new), len(SEQ) + 5)
        self.assertEqual(rec["expected_cdeii_len"], 25)
        self.assertEqual((rec["expected_start"], rec["expected_end"]), (11, 69))

    def test_insertion_leaves_cdei_and_cdeiii_untouched(self):
        new, rec = V.apply_edit(SEQ, FWD, 10, 37)
        self.assertEqual(new[10:18], CDEI)
        # the motif was at 1-based 39 to 64 and moves 10 bases right
        self.assertEqual(new[38 + 10:64 + 10], CDEIII)

    def test_deletion(self):
        new, rec = V.apply_edit(SEQ, FWD, -5, 21)
        self.assertEqual(rec["changed_bases"], SEQ[20:25])
        self.assertEqual(new, SEQ[:20] + SEQ[25:])
        self.assertEqual(rec["expected_cdeii_len"], 15)
        self.assertEqual(new[10:18], CDEI)
        self.assertEqual(new[33:59], CDEIII)

    def test_valid_positions(self):
        # insertion of 5: block [p-4, p] inside [19, 38], and p <= 37
        self.assertEqual(V.valid_positions(FWD, 5), (23, 37))
        # deletion of 5: [q, q+4] inside [19, 38]
        self.assertEqual(V.valid_positions(FWD, -5), (19, 34))
        # a deletion longer than CDEII is impossible
        self.assertIsNone(V.valid_positions(FWD, -21))

    def test_rejects_positions_that_touch_cdei_or_cdeiii(self):
        with self.assertRaises(ValueError):
            V.apply_edit(SEQ, FWD, 5, 22)      # block would start at 18, inside CDEI
        with self.assertRaises(ValueError):
            V.apply_edit(SEQ, FWD, 5, 38)      # copy would abut the CDEIII motif
        with self.assertRaises(ValueError):
            V.apply_edit(SEQ, FWD, -5, 35)     # deletion would reach base 39, the motif

    def test_reverse_strand_call(self):
        # the same element reverse-complemented: motif first in forward coordinates
        from padlib import revcomp
        seq = revcomp(SEQ)
        n = len(seq)
        rev = {"strand": "-", "start": n - 64 + 1, "end": n - 11 + 1, "cdeiii_start": n - 64 + 1,
               "cdeiii_end": n - 39 + 1, "cdei_start": n - 18 + 1, "cdei_end": n - 11 + 1}
        a, b = V.cdeii_interval(rev)
        self.assertEqual(seq[a - 1:b], revcomp(CDEII))
        new, rec = V.apply_edit(seq, rev, 10, b - 1)
        self.assertEqual(rec["expected_cdeii_len"], 30)
        self.assertEqual((rec["expected_start"], rec["expected_end"]), (rev["start"], rev["end"] + 10))
        self.assertEqual(revcomp(new)[10:18], CDEI)

    def test_shift_position(self):
        ins = {"size": 5, "position": 30}
        self.assertEqual(V.shift_position(30, ins), 30)
        self.assertEqual(V.shift_position(31, ins), 36)
        dele = {"size": -5, "position": 21}
        self.assertEqual(V.shift_position(20, dele), 20)
        self.assertIsNone(V.shift_position(23, dele))
        self.assertEqual(V.shift_position(26, dele), 21)

    def test_draw_position_is_seeded_and_valid(self):
        a = [V.draw_position(FWD, 5, np.random.default_rng(3)) for _ in range(5)]
        b = [V.draw_position(FWD, 5, np.random.default_rng(3)) for _ in range(5)]
        self.assertEqual(a, b)
        self.assertTrue(all(23 <= p <= 37 for p in a))


if __name__ == "__main__":
    unittest.main()
