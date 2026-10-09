"""Shared helpers for the Python stages.

Kept compatible with Python 3.8 because scripts/lib/pcan_runner.py imports it
from inside PCAn's own environment, which pins Python 3.8.17.
"""
from __future__ import annotations

import ast
import gzip
import math
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


# ----------------------------------------------------------------------------
# configuration, paths, logging

def read_conf():
    conf = {}
    path = os.path.join(ROOT, "project.conf")
    if not os.path.exists(path):
        sys.exit("project.conf not found; run: bash scripts/00_configure.sh --yes")
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                conf[k.strip()] = v.strip()
    conf.setdefault("DATA_DIR", os.path.join(ROOT, "data"))
    return conf


CONF = None


def conf():
    global CONF
    if CONF is None:
        CONF = read_conf()
    return CONF


def data_dir(*parts):
    return os.path.join(conf()["DATA_DIR"], *parts)


def repo(*parts):
    return os.path.join(ROOT, *parts)


def pcan_dir():
    return data_dir("pcan", "point-centromere-detection", "PCAn")


def rel(path):
    """Path as written in tracked files: relative to the repository root, with
    the data directory always shown as data/ wherever it really lives."""
    path = os.path.abspath(path)
    dd = os.path.abspath(conf()["DATA_DIR"])
    if path == dd or path.startswith(dd + os.sep):
        return "data" + path[len(dd):]
    if path.startswith(ROOT + os.sep):
        return path[len(ROOT) + 1:]
    return os.path.basename(path)


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(msg):
    sys.stderr.write("%s %s\n" % (utc_now(), msg))
    sys.stderr.flush()


def mark_done(name):
    os.makedirs(data_dir("state"), exist_ok=True)
    with open(data_dir("state", name + ".done"), "w") as fh:
        fh.write(utc_now() + "\n")


def is_done(name):
    return os.path.exists(data_dir("state", name + ".done"))


def atomic_write_text(path, text):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), prefix=".tmp_")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def atomic_write_tsv(df, path, **kw):
    """Write a pandas table so that an interrupted run never leaves half a file."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), prefix=".tmp_")
    os.close(fd)
    df.to_csv(tmp, sep="\t", index=False, **kw)
    os.replace(tmp, path)


# ----------------------------------------------------------------------------
# sequences

def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path):
    """Return [(name, sequence)] with the name cut at the first whitespace, the
    same rule FIMO and Biopython use for sequence identifiers."""
    records, name, chunks = [], None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def write_fasta(path, records, width=80):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt") as fh:
        for name, seq in records:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def at_fraction(seq):
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    return (s.count("A") + s.count("T")) / acgt if acgt else float("nan")


# ----------------------------------------------------------------------------
# contiguity

def nx(lengths, x=50):
    """Length L such that sequences of length >= L hold at least x per cent of
    the total. N50 is nx(lengths, 50)."""
    ls = sorted((int(v) for v in lengths if v > 0), reverse=True)
    total = sum(ls)
    if total == 0:
        return 0
    target = total * x / 100.0
    acc = 0
    for v in ls:
        acc += v
        if acc >= target:
            return v
    return ls[-1]


def lx(lengths, x=50):
    """Smallest number of sequences that together hold at least x per cent of
    the total. L50 is lx(lengths, 50)."""
    ls = sorted((int(v) for v in lengths if v > 0), reverse=True)
    total = sum(ls)
    if total == 0:
        return 0
    target = total * x / 100.0
    acc = 0
    for i, v in enumerate(ls, 1):
        acc += v
        if acc >= target:
            return i
    return len(ls)


def contig_lengths_split_at_gaps(seq, min_gap=1):
    """Lengths of the gap-free pieces of one sequence, splitting at every run
    of at least min_gap N. Every contig N50 in this repository, for published
    and fragmented assemblies alike, is computed this way with min_gap 1, so
    that all of them sit on one axis. NCBI's own contig N50, which follows the
    submitted AGP, is recorded beside it for comparison."""
    import re
    return [len(p) for p in re.split("[Nn]{%d,}" % min_gap, seq) if p]


# ----------------------------------------------------------------------------
# PCAn

def pcan_settings():
    """Read the genus list, the Kazachstania species list and the motif and
    threshold table straight from the checked-out PCAn script, so that this
    repository never restates PCAn's settings by hand."""
    path = os.path.join(pcan_dir(), "AutomatedCENretrieval.py")
    with open(path) as fh:
        tree = ast.parse(fh.read())
    found = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in ("available_genus_list", "available_kazachstania_species", "MotifsAndThresholds",
                        "CDEIIIup", "Nrows", "Range", "MinORFLength"):
                found[name] = ast.literal_eval(node.value)
    return found


def run_pcan(assembly, genus, out_tsv, species="", label="", keep_fimo="", timeout=None):
    """Run scripts/04_run_pcan.sh, the only place PCAn is invoked, and return
    the metadata it writes beside the call table."""
    cmd = ["bash", repo("scripts", "04_run_pcan.sh"), "--assembly", assembly, "--genus", genus,
           "--out", out_tsv]
    if species:
        cmd += ["--species", species]
    if label:
        cmd += ["--label", label]
    if keep_fimo:
        cmd += ["--keep-fimo", keep_fimo]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, timeout=timeout)
    meta_path = out_tsv + ".meta.tsv"
    meta = {"returncode": proc.returncode}
    if os.path.exists(meta_path):
        with open(meta_path) as fh:
            for line in fh:
                if "\t" in line:
                    k, v = line.rstrip("\n").split("\t", 1)
                    meta[k] = v
    if proc.returncode != 0:
        meta["stderr_tail"] = proc.stderr[-2000:]
    return meta


# ----------------------------------------------------------------------------
# intervals

def wilson(k, n, z=1.959963984540054):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def clopper_pearson(k, n, alpha=0.05):
    """Exact interval, used where counts are small."""
    from scipy.stats import beta
    if n == 0:
        return (float("nan"), float("nan"))
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k)
    return (float(lo), float(hi))


def cluster_bootstrap(units, statistic, n_boot=10000, seed=20261009, alpha=0.05):
    """Percentile bootstrap that resamples whole independent units (genomes or
    strains), never the centromeres inside them.

    units      list of per-unit data objects
    statistic  function(list of units) -> float
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    n = len(units)
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    est = statistic(units)
    reps = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        reps[b] = statistic([units[i] for i in idx])
    reps = reps[~np.isnan(reps)]
    if reps.size == 0:
        return (est, float("nan"), float("nan"))
    lo, hi = np.percentile(reps, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return (est, float(lo), float(hi))
