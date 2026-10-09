"""Run one stage, record its elapsed time, peak resident memory and peak disk
footprint, and stop it if the footprint crosses the ceiling in project.conf.

Usage (from a bash stage):
    python3 scripts/lib/measure.py --stage NAME -- command [args ...]

Standard library only, and compatible with Python 3.8, so that it runs under
the conda base interpreter before either project environment exists.

Peak memory is ru_maxrss of the waited-for children: the largest resident set
of any single process in the stage, which is the figure that decides whether a
stage fits in RAM. With --jobs above 1, concurrent processes add up and the
true peak is higher than this figure; logs/resources.tsv records the job count
beside every row for that reason.

Disk is sampled every --interval seconds with du on the data directory, and the
free space of the drive that holds it is sampled at the same time. Under WSL 2
that drive is the Windows drive holding the virtual disk file, because the
file grows there as data are written. A short-lived peak between two samples
can be missed; the interval is 15 s by default and the stages that write large
transient files (arm 3) log their own per-strain peaks as well.
"""
import argparse
import os
import resource
import shutil
import signal
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
HOST_FREE_FLOOR_GB = 8.0
COLUMNS = ["stage", "start_utc", "end_utc", "elapsed_s", "peak_rss_mb",
           "start_footprint_gb", "peak_footprint_gb", "min_free_drive_gb",
           "threads", "jobs", "exit_status", "note"]


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_conf():
    conf = {}
    path = os.path.join(ROOT, "project.conf")
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    k, v = line.split("=", 1)
                    conf[k.strip()] = v.strip()
    return conf


def drive_path(data_dir):
    try:
        with open("/proc/version") as fh:
            wsl = "microsoft" in fh.read().lower()
    except OSError:
        wsl = False
    if wsl and os.path.isdir("/mnt/c") and not data_dir.startswith("/mnt/"):
        return "/mnt/c"
    return data_dir


def footprint_gb(path):
    if not os.path.isdir(path):
        return 0.0
    try:
        out = subprocess.run(["du", "-sb", path], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, universal_newlines=True).stdout
        return int(out.split()[0]) / 1e9
    except (ValueError, IndexError):
        return float("nan")


class DiskWatch:
    """Samples the footprint in a background thread and enforces the ceiling."""

    def __init__(self, data_dir, ceiling_gb, interval, on_breach):
        self.data_dir = data_dir
        self.drive = drive_path(data_dir)
        self.ceiling = ceiling_gb
        self.interval = interval
        self.on_breach = on_breach
        self.start = footprint_gb(data_dir)
        self.peak = self.start
        self.min_free = shutil.disk_usage(self.drive).free / 1e9
        self.breach = ""
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def sample(self):
        fp = footprint_gb(self.data_dir)
        free = shutil.disk_usage(self.drive).free / 1e9
        self.peak = max(self.peak, fp)
        self.min_free = min(self.min_free, free)
        if self.ceiling and fp > self.ceiling and not self.breach:
            self.breach = "footprint %.1f GB above the %.1f GB ceiling" % (fp, self.ceiling)
            self.on_breach(self.breach)
        elif free < HOST_FREE_FLOOR_GB and not self.breach:
            self.breach = "drive free space %.1f GB below the %.0f GB floor" % (free, HOST_FREE_FLOOR_GB)
            self.on_breach(self.breach)

    def _run(self):
        while not self._stop.wait(self.interval):
            self.sample()

    def __enter__(self):
        self._thread.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._thread.join()
        self.sample()


def append_row(log_path, row):
    new = not os.path.exists(log_path)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "a") as fh:
        if new:
            fh.write("\t".join(COLUMNS) + "\n")
        fh.write("\t".join(str(row.get(c, "")) for c in COLUMNS) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--log", default=os.path.join(ROOT, "logs", "resources.tsv"))
    ap.add_argument("--interval", type=float, default=15.0)
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()
    cmd = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not cmd:
        ap.error("no command given")

    conf = read_conf()
    data_dir = conf.get("DATA_DIR", os.path.join(ROOT, "data"))
    ceiling = float(conf.get("DISK_CEILING_GB", "0") or 0)
    env = dict(os.environ, PAD_MEASURED="1")
    start_utc, t0 = utc_now(), time.time()
    proc = subprocess.Popen(cmd, env=env, start_new_session=True)

    def breach(reason):
        sys.stderr.write("%s refusing to continue %s: %s\n" % (utc_now(), args.stage, reason))
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass

    def forward(signum, _frame):
        try:
            os.killpg(proc.pid, signum)
        except ProcessLookupError:
            pass

    signal.signal(signal.SIGINT, forward)
    signal.signal(signal.SIGTERM, forward)
    with DiskWatch(data_dir, ceiling, args.interval, breach) as watch:
        status = proc.wait()
    rss_mb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0
    append_row(args.log, {
        "stage": args.stage, "start_utc": start_utc, "end_utc": utc_now(),
        "elapsed_s": "%.0f" % (time.time() - t0), "peak_rss_mb": "%.0f" % rss_mb,
        "start_footprint_gb": "%.2f" % watch.start, "peak_footprint_gb": "%.2f" % watch.peak,
        "min_free_drive_gb": "%.1f" % watch.min_free,
        "threads": conf.get("THREADS", ""), "jobs": os.environ.get("PAD_JOBS", conf.get("JOBS", "")),
        "exit_status": status, "note": watch.breach,
    })
    sys.exit(3 if watch.breach else status)


if __name__ == "__main__":
    main()
