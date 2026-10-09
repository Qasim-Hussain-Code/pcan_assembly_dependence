"""Run PCAn v1.0 once, unchanged, on one assembly and turn its output into a
call table. Called only by scripts/04_run_pcan.sh, inside PCAn's environment
(Python 3.8.17, MEME 4.11.2).

PCAn asks three questions on standard input: the path to the assembly, the
index of the genus in its list, and whether to run the BLAST synteny check
(plus a fourth, the species, for Kazachstania). The answers are written to its
standard input here. The synteny check is always declined: PCAn writes its call
table, CENsequences.txt, before the BLAST step begins (lines 300 to 316 of the
script), so the calls do not depend on the answer, and the check produces only
a plot.

FIMO is reached through a pass-through wrapper placed first on PATH. PCAn runs
FIMO twice with the same --oc directory, so the second run (CDEI) overwrites
the first (CDEIII). The wrapper copies each run's fimo.txt and standard error
aside. Nothing PCAn sees changes. The first-pass output is what links each call
to its CDEIII hit, and with it the strand, which PCAn does not report.
"""
import argparse
import gzip
import os
import resource
import shutil
import stat
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import padlib as P  # noqa: E402

SHIM = """#!/usr/bin/env bash
# Pass-through FIMO wrapper written by scripts/lib/pcan_runner.py.
n=$(( $(cat "$PAD_SHIM_DIR/count" 2>/dev/null || echo 0) + 1 ))
echo "$n" > "$PAD_SHIM_DIR/count"
printf '%s\\n' "$*" > "$PAD_SHIM_DIR/pass$n.args"
"$PAD_REAL_FIMO" "$@" 2> "$PAD_SHIM_DIR/pass$n.stderr"
status=$?
if [[ -f FIMOresults/fimo.txt ]]; then cp FIMOresults/fimo.txt "$PAD_SHIM_DIR/pass$n.fimo.txt"; fi
exit $status
"""

CALL_COLUMNS = ["label", "contig", "contig_hit", "start", "end", "strand", "cdei_start", "cdei_end",
                "cdeiii_start", "cdeiii_end", "cdeiii_pvalue", "cdeii_len", "cdeii_at", "fimo_i_score",
                "fimo_ii_score", "overall_score", "sequence"]
PCAN_COLUMNS = ["Contig", "ContigHit", "Start", "End", "Sequence", "CDEIIlen", "CDEIIAT",
                "FIMOIscore", "FIMOIIscore", "OverallScore"]


def read_fimo(path):
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        header = fh.readline()
        if not header.startswith("#pattern name"):
            raise SystemExit("unexpected FIMO output format in %s" % path)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                rows.append({"seq": f[1], "start": int(f[2]), "stop": int(f[3]), "strand": f[4],
                             "score": float(f[5]), "p": f[6]})
    return rows


def parse_censequences(path):
    """The table above the line of asterisks, and the settings PCAn reports below it."""
    with open(path) as fh:
        text = fh.read()
    table, _, footer = text.partition("\n****************")
    lines = [l for l in table.split("\n") if l.strip()]
    if not lines or lines[0].split("\t") != PCAN_COLUMNS:
        raise SystemExit("unexpected CENsequences.txt header: %r" % (lines[:1],))
    rows = [dict(zip(PCAN_COLUMNS, l.split("\t"))) for l in lines[1:]]
    settings = {}
    for l in footer.split("\n"):
        if l.startswith("CDEI motif used:"):
            settings["cdei_motif"] = os.path.basename(l.split(":", 1)[1].strip())
        elif l.startswith("CDEIII motif used:"):
            settings["cdeiii_motif"] = os.path.basename(l.split(":", 1)[1].strip())
        elif l.startswith("FIMO significance threshold for the CDEI motif"):
            settings["cdei_threshold"] = l.rsplit(":", 1)[1].strip()
        elif l.startswith("FIMO significance threshold for the CDEIII motif"):
            settings["cdeiii_threshold"] = l.rsplit(":", 1)[1].strip()
        elif l.startswith("Genus/species selected:"):
            settings["selected"] = l.split(":", 1)[1].strip()
    return rows, settings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--genus", required=True)
    ap.add_argument("--species", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--keep-fimo", default="")
    args = ap.parse_args()

    st = P.pcan_settings()
    genera = st["available_genus_list"]
    if args.genus not in genera:
        raise SystemExit("genus %s is not in PCAn's list" % args.genus)
    answers = [None, str(genera.index(args.genus))]
    if args.genus == "Kazachstania":
        kaz = st["available_kazachstania_species"]
        if args.species not in kaz:
            raise SystemExit("Kazachstania %s is not in PCAn's list %s" % (args.species, kaz))
        answers.append(str(kaz.index(args.species)))
    answers.append("No")

    work = os.path.abspath(args.workdir)
    run_dir = os.path.join(work, "run")
    shim_dir = os.path.join(work, "shim")
    os.makedirs(run_dir)
    os.makedirs(shim_dir)
    genome = os.path.join(work, "genome.fna")
    if args.assembly.endswith(".gz"):
        with gzip.open(args.assembly, "rb") as src, open(genome, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
    else:
        shutil.copyfile(args.assembly, genome)
    answers[0] = genome

    real_fimo = shutil.which("fimo")
    if not real_fimo:
        raise SystemExit("fimo not found in the PCAn environment")
    shim = os.path.join(shim_dir, "fimo")
    with open(shim, "w") as fh:
        fh.write(SHIM)
    os.chmod(shim, os.stat(shim).st_mode | stat.S_IXUSR)
    env = dict(os.environ, PATH=shim_dir + os.pathsep + os.environ["PATH"], PAD_SHIM_DIR=shim_dir,
               PAD_REAL_FIMO=real_fimo, MPLBACKEND="Agg")

    script = os.path.join(P.pcan_dir(), "AutomatedCENretrieval.py")
    t0 = time.time()
    with open(os.path.join(work, "pcan.stdout"), "w") as out, open(os.path.join(work, "pcan.stderr"), "w") as err:
        proc = subprocess.run([sys.executable, script], input="\n".join(answers) + "\n", cwd=run_dir,
                              env=env, stdout=out, stderr=err, universal_newlines=True)
    elapsed = time.time() - t0
    rss_mb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0

    pass1 = read_fimo(os.path.join(shim_dir, "pass1.fimo.txt"))
    pass2 = read_fimo(os.path.join(shim_dir, "pass2.fimo.txt"))
    cap_warning = "no"
    for n in (1, 2):
        e = os.path.join(shim_dir, "pass%d.stderr" % n)
        if os.path.exists(e):
            with open(e) as fh:
                if "max stored scores" in fh.read().lower():
                    cap_warning = "yes"
    meta = {"label": args.label, "assembly": P.rel(args.assembly), "genus": args.genus,
            "species": args.species, "pcan_exit_status": proc.returncode, "elapsed_s": "%.1f" % elapsed,
            "max_rss_mb": "%.0f" % rss_mb, "fimo_pass1_hits": len(pass1), "fimo_pass2_hits": len(pass2),
            "fimo_max_stored_scores_warning": cap_warning, "finished_utc": P.utc_now()}

    cen = os.path.join(run_dir, "CENsequences.txt")
    status = "ok"
    if proc.returncode != 0 or not os.path.exists(cen):
        status = "failed"
        with open(os.path.join(work, "pcan.stderr")) as fh:
            tail = fh.read()[-600:].replace("\n", " | ").replace("\t", " ")
        meta["error"] = tail.replace(work, "<workdir>")
    else:
        rows, settings = parse_censequences(cen)
        meta.update(settings)
        calls = []
        for r in rows:
            k = int(r["ContigHit"].rsplit("_", 1)[1])
            hit = pass1[k - 1]
            if hit["seq"] != r["Contig"]:
                raise SystemExit("call %s does not match FIMO row %d (%s)" % (r["ContigHit"], k, hit["seq"]))
            start, end = int(r["Start"]), int(r["End"])
            if hit["strand"] == "+":
                cdei = (start, start + 7)
                if end != hit["stop"]:
                    raise SystemExit("call %s ends at %d, CDEIII hit at %d" % (r["ContigHit"], end, hit["stop"]))
            else:
                cdei = (end - 7, end)
                if start != hit["start"]:
                    raise SystemExit("call %s starts at %d, CDEIII hit at %d" % (r["ContigHit"], start, hit["start"]))
            calls.append([args.label, r["Contig"], r["ContigHit"], start, end, hit["strand"], cdei[0], cdei[1],
                          hit["start"], hit["stop"], hit["p"], r["CDEIIlen"], r["CDEIIAT"], r["FIMOIscore"],
                          r["FIMOIIscore"], r["OverallScore"], r["Sequence"]])
        # Replay PCAn's candidate table and filters from the two FIMO passes, to
        # say later which filter removed a candidate. The replay must leave
        # exactly PCAn's calls standing, or it is not used.
        import pcan_trace
        genome_seqs = dict(P.read_fasta(genome))
        labelled, anchors, matches = pcan_trace.trace(
            genome_seqs, os.path.join(shim_dir, "pass1.fimo.txt"), os.path.join(shim_dir, "pass2.fimo.txt"), rows)
        labelled.to_csv(args.out + ".candidates.tsv.gz", sep="\t", index=False, compression="gzip")
        meta["trace_matches_pcan"] = "yes" if matches else "no"
        meta["candidates"] = len(labelled)
        for k, v in anchors.items():
            meta[k] = v
        tmp = args.out + ".part"
        with open(tmp, "w") as fh:
            fh.write("\t".join(CALL_COLUMNS) + "\n")
            for c in calls:
                fh.write("\t".join(str(v) for v in c) + "\n")
        os.replace(tmp, args.out)
        meta["n_calls"] = len(calls)
    meta["status"] = status

    if args.keep_fimo and status == "ok":
        os.makedirs(args.keep_fimo, exist_ok=True)
        for n in (1, 2):
            src = os.path.join(shim_dir, "pass%d.fimo.txt" % n)
            if os.path.exists(src):
                with open(src, "rb") as s, gzip.open(os.path.join(args.keep_fimo, "pass%d.fimo.txt.gz" % n), "wb") as d:
                    shutil.copyfileobj(s, d)
            a = os.path.join(shim_dir, "pass%d.args" % n)
            if os.path.exists(a):
                with open(a) as s, open(os.path.join(args.keep_fimo, "pass%d.args" % n), "w") as d:
                    d.write(s.read().replace(work, "<workdir>").replace(P.pcan_dir(), "<PCAn>"))

    with open(args.out + ".meta.tsv", "w") as fh:
        for k, v in meta.items():
            fh.write("%s\t%s\n" % (k, v))
    sys.exit(0 if status == "ok" else 2)


if __name__ == "__main__":
    main()
