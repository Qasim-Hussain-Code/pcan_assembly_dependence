"""Planted CDEII variants for arm 1b, a controlled simulation on real genomes.

An insertion of s bp duplicates s bases of CDEII in tandem: the copy is the s
bases immediately before the insertion point, so the result reads ...XY + XY...
with both copies inside CDEII. A deletion removes s contiguous bases from
inside CDEII. Neither ever touches CDEI or the CDEIII motif. Duplicating
adjacent sequence keeps CDEII's AT content close to the original, whereas
random insertions would shift it, and PCAn scores AT content.

Coordinates are 1-based and inclusive, in the forward orientation of the
sequence that holds the call. CDEII is the stretch strictly between CDEI and
the CDEIII motif as PCAn reported them:
  forward-strand call: CDEI at [start, start + 7], motif ends at end,
                       CDEII = [cdei_end + 1, cdeiii_start - 1]
  reverse-strand call: motif at [start, start + 25], CDEI ends at end,
                       CDEII = [cdeiii_end + 1, cdei_start - 1]
Either way an edit of size s inside CDEII leaves the call's start where it was
and moves its end by s, and PCAn should report a CDEII length of the original
plus s.
"""
from __future__ import annotations


def cdeii_interval(call):
    if call["strand"] == "+":
        return call["cdei_end"] + 1, call["cdeiii_start"] - 1
    return call["cdeiii_end"] + 1, call["cdei_start"] - 1


def valid_positions(call, size):
    """Positions an edit may take.

    Insertion of size s > 0: the insertion point p is the last base of the
    duplicated block, which is [p - s + 1, p]; the copy goes in after base p.
    The block must lie in CDEII and p must leave at least one CDEII base after
    it, so the copy never abuts the CDEIII motif or CDEI.
    Deletion of size s < 0: the first deleted base q, removing [q, q + |s| - 1],
    which must lie in CDEII.
    """
    a, b = cdeii_interval(call)
    if size > 0:
        lo, hi = a + size - 1, b - 1
    else:
        lo, hi = a, b - (-size) + 1
    return (lo, hi) if hi >= lo else None


def apply_edit(seq, call, size, position):
    """Apply one edit to one sequence (a str, the whole chromosome or contig).
    Returns (new_sequence, record). The record carries the expected call after
    the edit."""
    rng = valid_positions(call, size)
    if rng is None or not (rng[0] <= position <= rng[1]):
        raise ValueError("position %d is not valid for a %+d bp edit of this CDEII" % (position, size))
    if size > 0:
        block = seq[position - size:position]          # bases position-size+1 .. position
        new = seq[:position] + block + seq[position:]
        changed = block
    else:
        k = -size
        changed = seq[position - 1:position - 1 + k]
        new = seq[:position - 1] + seq[position - 1 + k:]
    a, b = cdeii_interval(call)
    record = {"kind": "insertion" if size > 0 else "deletion", "size": size, "position": position,
              "changed_bases": changed, "cdeii_start": a, "cdeii_end": b,
              "expected_start": call["start"], "expected_end": call["end"] + size,
              "expected_cdeii_len": (b - a + 1) + size}
    return new, record


def shift_position(pos, edit):
    """Where a coordinate on the same sequence lands after an edit. Positions
    after the edit move by its size; positions before do not."""
    if edit["size"] > 0:
        return pos + edit["size"] if pos > edit["position"] else pos
    k = -edit["size"]
    if pos < edit["position"]:
        return pos
    if pos >= edit["position"] + k:
        return pos - k
    return None   # the base was deleted


def draw_position(call, size, rng):
    r = valid_positions(call, size)
    if r is None:
        return None
    return int(rng.integers(r[0], r[1] + 1))
