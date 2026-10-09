"""End-to-end known-answer test of scripts/04_run_pcan.sh.

Three 20 kb pieces of the S288C genome that ships with PCAn, each centred on a
centromere, are written as a three-contig assembly and run through the
wrapper. The expected calls are fixed here from the S288C coordinates (which
agree with SGD apart from the 26th base of PCAn's CDEIII motif):

  CEN4   chromosome IV  449711-449822  forward  CDEII 78
  CEN9   chromosome IX  355629-355746  forward  CDEII 84
  CEN10  chromosome X   436306-436425  reverse  CDEII 86

Cut out from position 10,001 bases before each, a call at genome position x
sits at x - offset on its piece. The test needs PCAn installed by
scripts/01_install.sh and is skipped otherwise.
"""
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "lib"))
import padlib as P  # noqa: E402

PIECES = [  # contig, source, offset (bases before the piece), expected start, end, strand, CDEII
    ("t4", "NC_001136.10", 439000, 449711, 449822, "+", 78),
    ("t9", "NC_001141.2", 345000, 355629, 355746, "+", 84),
    ("t10", "NC_001142.9", 426000, 436306, 436425, "-", 86),
]


def installed():
    try:
        return os.path.exists(os.path.join(P.pcan_dir(), "AutomatedCENretrieval.py"))
    except SystemExit:
        return False


@unittest.skipUnless(os.path.exists(P.repo("project.conf")) and installed(), "PCAn not installed")
class TestWrapper(unittest.TestCase):
    def test_three_known_centromeres(self):
        genome = dict(P.read_fasta(os.path.join(P.pcan_dir(), "test_runs", "Saccharomyces_cerevisiae.fna")))
        work = tempfile.mkdtemp(dir=P.data_dir("work"))
        try:
            fa = os.path.join(work, "three.fna")
            P.write_fasta(fa, [(name, genome[src][off:off + 20000]) for name, src, off, *_ in PIECES])
            out = os.path.join(work, "calls.tsv")
            meta = P.run_pcan(fa, "Saccharomyces", out, label="wrapper_test")
            self.assertEqual(meta.get("status"), "ok", meta)
            self.assertEqual(meta.get("trace_matches_pcan"), "yes")
            import pandas as pd
            calls = pd.read_csv(out, sep="\t")
            self.assertEqual(len(calls), 3)
            got = {r.contig: (r.start, r.end, r.strand, r.cdeii_len) for r in calls.itertuples()}
            for name, src, off, s, e, strand, cdeii in PIECES:
                self.assertEqual(got[name], (s - off, e - off, strand, cdeii), name)
                seq = genome[src][s - 1:e]
                row = calls[calls.contig == name].iloc[0]
                expect = seq if strand == "+" else P.revcomp(seq)
                self.assertEqual(row.sequence.upper(), expect.upper(), name)
        finally:
            shutil.rmtree(work)


if __name__ == "__main__":
    unittest.main()
