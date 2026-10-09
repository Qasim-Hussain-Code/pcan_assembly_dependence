"""Known-answer tests for the five-way status classifier.

The target assembly is described by hand: contig A carries the element whole,
contig B carries it with an N run, contigs C and D each carry one flank.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "lib"))
import status as S  # noqa: E402

F = 2000          # flank length
E = 120           # element length on the long-read assembly
TARGETS = {
    "A": "C" * 5000 + "A" * E + "C" * 5000,
    "B": "C" * 5000 + "A" * 50 + "N" * 20 + "A" * 50 + "C" * 5000,
}


def aln(tname, tstart, tend, strand="+", identity=0.99, qcov=1.0, mapq=60, primary=True):
    qlen = F
    qend = int(qlen * qcov)
    alen = 1000
    return S.Aln("q", qlen, 0, qend, strand, tname, 20000, tstart, tend, int(alen * identity), alen, mapq, primary)


def seq(name):
    return TARGETS[name]


class TestPlacement(unittest.TestCase):
    def test_unique_good_alignment_is_placed(self):
        a = aln("A", 3000, 5000)
        self.assertIs(S.place_flank([a], 0.95, 0.9, 20), a)

    def test_thresholds_and_ambiguity(self):
        self.assertIsNone(S.place_flank([aln("A", 3000, 5000, identity=0.90)], 0.95, 0.9, 20))
        self.assertIsNone(S.place_flank([aln("A", 3000, 5000, qcov=0.5)], 0.95, 0.9, 20))
        self.assertIsNone(S.place_flank([aln("A", 3000, 5000, mapq=3)], 0.95, 0.9, 20))
        self.assertIsNone(S.place_flank([aln("A", 3000, 5000), aln("B", 3000, 5000)], 0.95, 0.9, 20))
        self.assertIsNone(S.place_flank([aln("A", 3000, 5000, primary=False)], 0.95, 0.9, 20))

    def test_parse_paf(self):
        line = "q\t2000\t0\t2000\t+\tA\t10120\t3000\t5000\t1990\t2000\t60\ttp:A:P\tNM:i:10\n"
        a = S.parse_paf([line])[0]
        self.assertEqual((a.tname, a.tstart, a.tend, a.mapq, a.primary), ("A", 3000, 5000, 60, True))
        self.assertAlmostEqual(a.identity, 0.995)
        self.assertEqual(S.parse_paf([line.replace("tp:A:P", "tp:A:S")])[0].primary, False)


class TestClassify(unittest.TestCase):
    def test_intact_called_forward(self):
        st, d = S.classify(aln("A", 3000, 5000), aln("A", 5120, 7120), E, seq, {"A": [(5001, 5120)]})
        self.assertEqual(st, "intact_called")
        self.assertEqual((d["region_start"], d["region_end"]), (5001, 5120))

    def test_intact_uncalled(self):
        st, _ = S.classify(aln("A", 3000, 5000), aln("A", 5120, 7120), E, seq, {"A": [(9000, 9100)]})
        self.assertEqual(st, "intact_uncalled")

    def test_intact_called_reverse(self):
        # reverse placement: the right flank comes first on the target
        st, d = S.classify(aln("A", 5120, 7120, "-"), aln("A", 3000, 5000, "-"), E, seq, {"A": [(5001, 5120)]})
        self.assertEqual(st, "intact_called")

    def test_n_run(self):
        st, _ = S.classify(aln("B", 3000, 5000), aln("B", 5120, 7120), E, seq, {"B": [(5001, 5120)]})
        self.assertEqual(st, "n_run")

    def test_split_on_two_contigs(self):
        st, d = S.classify(aln("C", 0, 2000), aln("D", 0, 2000), E, seq, {})
        self.assertEqual(st, "split")
        self.assertEqual(d["reason"], "flanks on different sequences")

    def test_split_one_flank(self):
        st, _ = S.classify(aln("C", 0, 2000), None, E, seq, {})
        self.assertEqual(st, "split")

    def test_split_opposite_orientation(self):
        st, _ = S.classify(aln("A", 3000, 5000, "+"), aln("A", 5120, 7120, "-"), E, seq, {})
        self.assertEqual(st, "split")

    def test_split_too_far_apart(self):
        st, _ = S.classify(aln("A", 0, 2000), aln("A", 2000 + E + 1001, 4000 + E + 1001), E, seq, {})
        self.assertEqual(st, "split")

    def test_absent(self):
        st, _ = S.classify(None, None, E, seq, {})
        self.assertEqual(st, "absent")

    def test_broken_set(self):
        self.assertEqual(S.BROKEN, {"split", "n_run", "absent"})


if __name__ == "__main__":
    unittest.main()
