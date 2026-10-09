"""End-to-end known-answer test of the arm 2 and 3 comparison.

S288C (GCA_000146045.2, cached by arm 0) is cut at chosen positions, which
makes a second "assembly" of the same strain whose breaks are known exactly.
Run through PCAn and the liftover against the intact S288C:

  - a centromere whose element or either 2 kb flank holds a cut must come out
    broken (split), and every other centromere intact and called;
  - every cut must be found again by the breakpoint locator, within 100 bp.

Needs PCAn installed and the S288C assembly downloaded; skipped otherwise.
"""
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "lib"))
import fragment as F  # noqa: E402
import liftover as L  # noqa: E402
import padlib as P  # noqa: E402

ACC = "GCA_000146045.2"


def available():
    try:
        return os.path.exists(P.data_dir("assemblies", ACC + ".fna.gz")) and \
            os.path.exists(os.path.join(P.pcan_dir(), "AutomatedCENretrieval.py"))
    except SystemExit:
        return False


@unittest.skipUnless(os.path.exists(P.repo("project.conf")) and available(), "S288C or PCAn not available")
class TestComparison(unittest.TestCase):
    def test_known_cuts(self):
        work = tempfile.mkdtemp(dir=P.data_dir("work"))
        try:
            long_fa = P.data_dir("assemblies", ACC + ".fna.gz")
            seqs = dict(P.read_fasta(long_fa))
            long_calls_path = os.path.join(work, "long.tsv")
            self.assertEqual(P.run_pcan(long_fa, "Saccharomyces", long_calls_path, label="long")["status"], "ok")
            long_calls = pd.read_csv(long_calls_path, sep="\t")
            self.assertEqual(len(long_calls), 16)
            c = long_calls.sort_values("contig").reset_index(drop=True)
            cuts = {}
            # cut inside the CDEII of the first centromere, 1 kb into the left
            # flank of the second, 1 kb beyond the right flank of the third
            # (so that one stays intact), and one cut far from any centromere
            cuts.setdefault(c.loc[0, "contig"], []).append(int(c.loc[0, "start"]) + 50)
            cuts.setdefault(c.loc[1, "contig"], []).append(int(c.loc[1, "start"]) - 1000)
            cuts.setdefault(c.loc[2, "contig"], []).append(int(c.loc[2, "end"]) + 3000)
            cuts.setdefault(c.loc[3, "contig"], []).append(int(c.loc[3, "end"]) + 60000)
            cuts = {k: np.array(sorted(v)) for k, v in cuts.items()}
            frags = F.apply_cuts(list(seqs.items()), cuts)
            # contigs from an assembler carry their own names, so the
            # fragments are renamed to hide their source coordinates
            q_fa = os.path.join(work, "query.fna")
            P.write_fasta(q_fa, [("ctg%04d" % i, f[1]) for i, f in enumerate(frags)])
            q_calls_path = os.path.join(work, "query.tsv")
            self.assertEqual(P.run_pcan(q_fa, "Saccharomyces", q_calls_path, label="query")["status"], "ok")
            q_calls = pd.read_csv(q_calls_path, sep="\t")

            cens = L.centromere_elements(long_calls, seqs, 2000)
            fq = os.path.join(work, "flanks.fa")
            L.flank_fasta(cens, seqs, fq, 2000)
            paf = os.path.join(work, "flanks.paf")
            L.run_minimap2(q_fa, fq, paf, "asm10", 4)
            by_seq = {}
            for r in q_calls.itertuples():
                by_seq.setdefault(r.contig, []).append((int(r.start), int(r.end)))
            st = L.classify_elements(cens, paf, dict(P.read_fasta(q_fa)), by_seq).set_index("contig")
            self.assertEqual(st.loc[c.loc[0, "contig"], "status"], "split")
            self.assertEqual(st.loc[c.loc[1, "contig"], "status"], "split")
            others = [k for k in st.index if k not in (c.loc[0, "contig"], c.loc[1, "contig"])]
            self.assertTrue((st.loc[others, "status"] == "intact_called").all(), st.loc[others, "status"].to_dict())

            bpaf = os.path.join(work, "contigs.paf")
            L.run_minimap2(long_fa, q_fa, bpaf, "asm5", 4)
            bps, _ = L.breakpoints(bpaf, {n: len(s) for n, s in seqs.items()})
            for contig, cp in cuts.items():
                for p in cp:
                    near = bps[(bps["target"] == contig) & ((bps["position"] - p).abs() <= 100)]
                    self.assertGreaterEqual(len(near), 1, "cut at %s:%d not found" % (contig, p))
        finally:
            shutil.rmtree(work)


if __name__ == "__main__":
    unittest.main()
