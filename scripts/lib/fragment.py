"""Controlled fragmentation of a real assembly (arm 1a, a simulation).

Breakpoints are clean cuts between two adjacent bases. No sequence is removed,
so every fragment maps back to its source by an offset alone. A cut at
position p (1-based) separates base p from base p + 1.

Two breakpoint models:

uniform      a homogeneous Poisson process along each sequence.
at_weighted  an inhomogeneous Poisson process whose intensity in each 500 bp
             window is proportional to exp(gamma * (a - a_genome)), where a is
             the window's AT fraction and a_genome the genome's. With the
             default gamma = ln(10) / 0.30, a window 30 percentage points more
             AT-rich than the genome is ten times as likely to break. The form
             and the constant are my assumptions, chosen to express "short-read
             libraries under-cover very AT-rich sequence" with one parameter;
             arms 2 and 3 measure where real assemblies break.

Under both models the genome-wide expected number of cuts for a given rate is
the same, because window weights are normalised to a length-weighted mean of 1.
"""
from __future__ import annotations

import math

import numpy as np

WINDOW = 500
GAMMA_DEFAULT = math.log(10) / 0.30


def window_weights(seq, genome_at, gamma, window=WINDOW):
    """Unnormalised breakpoint weight of each window of one sequence. Windows
    that are all N get the weight of a window at the genome's AT content."""
    n = len(seq)
    starts = np.arange(0, n, window)
    w = np.empty(len(starts))
    s = seq.upper()
    for i, a in enumerate(starts):
        chunk = s[a:a + window]
        at = chunk.count("A") + chunk.count("T")
        acgt = at + chunk.count("C") + chunk.count("G")
        frac = at / acgt if acgt else genome_at
        w[i] = math.exp(gamma * (frac - genome_at))
    return starts, w


def draw_cuts(records, rate, model, rng, genome_at=None, gamma=GAMMA_DEFAULT, weights=None):
    """Cut positions for every sequence: {name: sorted array of p}, with
    1 <= p <= len - 1.

    rate     expected cuts per base pair, averaged over the genome
    weights  precomputed {name: (starts, w)} for the AT-weighted model; the
             weights are normalised here so their length-weighted mean is 1
    """
    cuts = {}
    if model == "uniform":
        for name, seq in records:
            n = len(seq)
            if n < 2:
                cuts[name] = np.array([], dtype=np.int64)
                continue
            k = rng.poisson(rate * (n - 1))
            cuts[name] = np.unique(rng.integers(1, n, size=k))
        return cuts
    if model != "at_weighted":
        raise ValueError("unknown model %s" % model)
    if weights is None:
        weights = {name: window_weights(seq, genome_at, gamma) for name, seq in records}
    total_len = 0.0
    total_w = 0.0
    for name, seq in records:
        starts, w = weights[name]
        lens = np.minimum(starts + WINDOW, len(seq)) - starts
        total_len += lens.sum()
        total_w += (w * lens).sum()
    mean_w = total_w / total_len if total_len else 1.0
    for name, seq in records:
        n = len(seq)
        starts, w = weights[name]
        lens = np.minimum(starts + WINDOW, n) - starts
        lam = rate * lens * (w / mean_w)
        k = rng.poisson(lam)
        pos = []
        for a, ln, kk in zip(starts, lens, k):
            if kk:
                pos.append(a + rng.integers(0, ln, size=kk))
        p = np.concatenate(pos) if pos else np.array([], dtype=np.int64)
        # positions are 0-based bases inside the window; a cut after base
        # index i (0-based) is p = i + 1, and p = n would be no cut at all
        p = p + 1
        cuts[name] = np.unique(p[(p >= 1) & (p <= n - 1)])
    return cuts


def apply_cuts(records, cuts):
    """Fragments as [(fragment_name, sequence, source_name, source_start)]
    with source_start 1-based. Fragment names are <source>__<start>-<end> in
    source coordinates, so a call on a fragment maps back without a lookup."""
    out = []
    for name, seq in records:
        bounds = [0] + [int(p) for p in cuts.get(name, [])] + [len(seq)]
        for a, b in zip(bounds[:-1], bounds[1:]):
            if b > a:
                out.append(("%s__%d-%d" % (name, a + 1, b), seq[a:b], name, a + 1))
    return out


def parse_fragment_name(fragment):
    """Inverse of the naming in apply_cuts: (source, start, end), 1-based."""
    source, span = fragment.rsplit("__", 1)
    a, b = span.split("-")
    return source, int(a), int(b)


def fragment_lengths_split_at_gaps(fragments):
    from padlib import contig_lengths_split_at_gaps
    out = []
    for _, seq, _, _ in fragments:
        out.extend(contig_lengths_split_at_gaps(seq))
    return out


def calibrate_rate(records, target_n50, model, seed, genome_at=None, gamma=GAMMA_DEFAULT,
                   weights=None, draws=12, tol=0.02, max_iter=40):
    """Rate whose realised contig N50, averaged over `draws` calibration
    draws, matches the target. Bisection on the log of the rate. Calibration
    uses its own seed so that the replicate seeds are untouched by it.

    Returns (rate, mean realised N50 at that rate). A rate of 0 means the
    unfragmented assembly is already at or below the target."""
    from padlib import nx
    base = fragment_lengths_split_at_gaps(apply_cuts(records, {}))
    if nx(base) <= target_n50:
        return 0.0, float(nx(base))
    if weights is None and model == "at_weighted":
        weights = {name: window_weights(seq, genome_at, gamma) for name, seq in records}

    def mean_n50(rate):
        rng = np.random.default_rng(seed)
        vals = []
        for _ in range(draws):
            c = draw_cuts(records, rate, model, rng, genome_at, gamma, weights)
            vals.append(nx(fragment_lengths_split_at_gaps(apply_cuts(records, c))))
        return float(np.mean(vals))

    # For exponential fragment lengths with mean m, N50 is about 1.68 m, so
    # 1.68 / target is a good first guess for the rate.
    lo, hi = math.log(1.68 / target_n50) - 3, math.log(1.68 / target_n50) + 3
    best = (math.exp(lo), mean_n50(math.exp(lo)))
    for _ in range(max_iter):
        mid = (lo + hi) / 2
        val = mean_n50(math.exp(mid))
        best = (math.exp(mid), val)
        if abs(val - target_n50) / target_n50 <= tol:
            break
        if val > target_n50:
            lo = mid
        else:
            hi = mid
    return best


# ----------------------------------------------------------------------------
# geometry of one PCAn call

def extraction_window(call):
    """The window PCAn cuts around a CDEIII hit, in source coordinates.

    From AutomatedCENretrieval.py, lines 168 to 182, with CDEIIIup = 250:
      forward strand: genome[start - 250 : end] in Python slicing, which is
                      1-based start - 249 to end, 275 bp for a 26 bp motif;
      reverse strand: genome[start - 1 : end + 250], which is 1-based start to
                      end + 250, 276 bp, then reverse-complemented.
    The paper describes a 250 bp window (the motif plus 224 bp upstream)."""
    if call["strand"] == "+":
        return call["cdeiii_start"] - 249, call["cdeiii_end"]
    return call["cdeiii_start"], call["cdeiii_end"] + 250


def cut_inside(cut_positions, a, b):
    """True if any cut separates two bases that both lie in [a, b]."""
    cp = np.asarray(cut_positions)
    return bool(((cp >= a) & (cp < b)).any())


def geometric_loss(call, cut_positions):
    """Whether the cuts alone remove this call, before any filter sees it.

    Forward strand: a cut anywhere in the extraction window. A cut upstream of
    the CDEI makes PCAn's slice start negative, Python then returns an empty
    window, and FIMO has nothing to search for CDEI in.
    Reverse strand: a cut between the CDEI and the end of the CDEIII motif. A
    cut further upstream only shortens the window, because a slice running past
    the end of a sequence is truncated; the CDEI and the CDEII length survive.
    """
    if call["strand"] == "+":
        a, b = extraction_window(call)
    else:
        a, b = call["start"], call["end"]
    return cut_inside(cut_positions, a, b)


def window_cut(call, cut_positions):
    a, b = extraction_window(call)
    return cut_inside(cut_positions, a, b)


def synteny_checkable(call, fragment_start, fragment_end, flank=10000):
    """PCAn's synteny step extracts 10 kb either side of a call (Range = 10000).
    The evidence is complete only if the fragment holding the call extends that
    far on both sides."""
    return call["start"] - flank >= fragment_start and call["end"] + flank <= fragment_end
