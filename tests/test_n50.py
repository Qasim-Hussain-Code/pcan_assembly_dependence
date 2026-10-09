"""Known-answer tests for contiguity statistics.

Each case is small enough to count by hand. For lengths 40, 30, 20 and 10
(total 100), half the total is 50: 40 alone is short of it and 40 + 30 = 70
passes it, so N50 is 30 and L50 is 2.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "lib"))
import padlib as P  # noqa: E402


class TestN50(unittest.TestCase):
    def test_four_lengths(self):
        self.assertEqual(P.nx([10, 20, 30, 40]), 30)
        self.assertEqual(P.lx([10, 20, 30, 40]), 2)

    def test_order_does_not_matter(self):
        self.assertEqual(P.nx([40, 10, 30, 20]), 30)

    def test_exactly_half(self):
        # 50 + 50: the first sequence reaches exactly half, so N50 is 50 and L50 is 1
        self.assertEqual(P.nx([50, 50]), 50)
        self.assertEqual(P.lx([50, 50]), 1)

    def test_single_sequence(self):
        self.assertEqual(P.nx([12071326]), 12071326)
        self.assertEqual(P.lx([12071326]), 1)

    def test_n90(self):
        # 40 + 30 + 20 = 90 reaches 90 per cent of 100
        self.assertEqual(P.nx([10, 20, 30, 40], 90), 20)
        self.assertEqual(P.lx([10, 20, 30, 40], 90), 3)

    def test_empty_and_zero(self):
        self.assertEqual(P.nx([]), 0)
        self.assertEqual(P.nx([0, 0]), 0)

    def test_split_at_gaps(self):
        self.assertEqual(P.contig_lengths_split_at_gaps("ACGTNNNNACGT"), [4, 4])
        self.assertEqual(P.contig_lengths_split_at_gaps("NNACGTN"), [4])
        self.assertEqual(P.contig_lengths_split_at_gaps("ACGTnnAC"), [4, 2])
        self.assertEqual(P.contig_lengths_split_at_gaps("ACGT"), [4])
        self.assertEqual(P.contig_lengths_split_at_gaps("NNNN"), [])

    def test_contig_n50_of_a_gapped_scaffold(self):
        # one scaffold of 60 + 40 bases joined by 10 N: contigs 60 and 40,
        # contig N50 60 (60 alone is more than half of 100)
        seq = "A" * 60 + "N" * 10 + "C" * 40
        self.assertEqual(P.nx(P.contig_lengths_split_at_gaps(seq)), 60)


if __name__ == "__main__":
    unittest.main()
