"""Five-way status of one long-read centromere (or one null window) in another
assembly of the same strain, from where its two flanks align.

The left flank is the F bases immediately before the element on the long-read
sequence and the right flank the F bases immediately after it. Each flank is
aligned to the other assembly with minimap2 and placed only if exactly one
primary alignment passes the identity, query-coverage and mapping-quality
thresholds. An unplaced flank is unplaced; nothing is forced.

Statuses
  intact_called    both flanks on one sequence, in order and orientation, at
                   most element length + max_extra apart, no N between them,
                   and a PCAn call overlaps the region between them
  intact_uncalled  as above, with no PCAn call there
  n_run            as above, but the region between the flanks holds an N
  split            one flank placed, or both placed but on different
                   sequences, in opposite orientations, or too far apart
  absent           neither flank placed

For a contig assembly n_run cannot occur, because contigs carry no N.
"""
from __future__ import annotations

from dataclasses import dataclass

STATUSES = ["intact_called", "intact_uncalled", "split", "n_run", "absent"]
BROKEN = {"split", "n_run", "absent"}


@dataclass
class Aln:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alen: int
    mapq: int
    primary: bool = True

    @property
    def identity(self):
        return self.nmatch / self.alen if self.alen else 0.0

    @property
    def qcov(self):
        return (self.qend - self.qstart) / self.qlen if self.qlen else 0.0


def parse_paf(lines):
    """minimap2 PAF records. Coordinates stay 0-based and half-open, as PAF
    writes them. An alignment is primary unless its tp tag says otherwise."""
    out = []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        tags = {t[:2]: t[5:] for t in f[12:]}
        out.append(Aln(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]),
                       int(f[9]), int(f[10]), int(f[11]), tags.get("tp", "P") == "P"))
    return out


def place_flank(alns, min_identity, min_qcov, min_mapq):
    good = [a for a in alns if a.primary and a.identity >= min_identity and a.qcov >= min_qcov
            and a.mapq >= min_mapq]
    return good[0] if len(good) == 1 else None


def classify(left, right, element_len, target_seq_fn, calls_by_seq, max_extra=1000, min_gap=-50, slack=10):
    """Status of one element.

    left, right      placements (Aln or None) of the left and right flank
    element_len      length of the element on the long-read assembly
    target_seq_fn    function(name) -> sequence of a target sequence
    calls_by_seq     {target name: [(start, end), ...]} PCAn calls, 1-based

    Returns (status, detail dict).
    """
    if left is None and right is None:
        return "absent", {"reason": "neither flank placed"}
    if left is None or right is None:
        return "split", {"reason": "one flank placed",
                         "placed_on": (left or right).tname}
    if left.tname != right.tname:
        return "split", {"reason": "flanks on different sequences"}
    if left.strand != right.strand:
        return "split", {"reason": "flanks in opposite orientations"}
    # On a forward placement the left flank comes first on the target; on a
    # reverse placement the right flank does. The element sits between them.
    if left.strand == "+":
        gap_start, gap_end = left.tend, right.tstart
    else:
        gap_start, gap_end = right.tend, left.tstart
    gap = gap_end - gap_start
    if gap < min_gap or gap > element_len + max_extra:
        return "split", {"reason": "flanks %d bp apart for a %d bp element" % (gap, element_len)}
    lo, hi = min(gap_start, gap_end), max(gap_start, gap_end)
    region = target_seq_fn(left.tname)[lo:hi]
    detail = {"target": left.tname, "region_start": lo + 1, "region_end": hi, "region_len": gap,
              "strand": left.strand}
    if "N" in region.upper():
        return "n_run", detail
    for s, e in calls_by_seq.get(left.tname, []):
        if s <= hi + slack and e >= lo + 1 - slack:
            detail["call_start"], detail["call_end"] = s, e
            return "intact_called", detail
    return "intact_uncalled", detail
