#!/usr/bin/env bash
# Arm 3: assemble real Illumina reads at known depths with two assemblers, one
# strain at a time, by the design registered in config/analysis_plan.md.
# For each strain: download both read files from ENA, resuming any transfer
# that breaks, check them against ENA's MD5, draw every depth and seed in one
# pass over each file, delete the files, check that mates still pair, trim
# with fastp, assemble every subsample with SPAdes and MEGAHIT, keep the
# contigs and reports, delete the reads. The next strain's download runs
# while the current strain assembles; assemblies never overlap.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash scripts/11_assemble_reads.sh [--max-strains N] [--force]

Strains come from config/arm3_runs.tsv (arm3_eligible = yes), in ascending
order_key, the order fixed in the analysis plan. Per strain, kept in
data/arm3/<code>/: assemblies/*.fna.gz (contigs, gzip), fastp/*.json,
manifest.tsv. Reads and assembler working directories are deleted when the
strain is done, and by the cleanup trap if the run is interrupted.

Strains are assembled one at a time. While one strain assembles, the reads of
the next one download, so that the link and the processors are both in use.
A strain whose reads could not be downloaded intact is logged as such and
left without a done marker, so that a later run tries it again.

Before a strain's download starts, the arm's wall time is projected to the end
of that strain from the longest download and the longest assembly seen so far
(priors of 1.5 h each before the first strain), and its disk from the largest
strain so far (prior 12 GB). A strain that would end past HOURS_ARM3, or take
the data directory past the disk ceiling in project.conf, is refused with the
shortfall. Strains are cut from the end of the order; depths, seeds and
controls never are.

  --max-strains N   stop after N strains (default: as many as the budget allows)
  --force           redo strains that are already complete
  --runs FILE       a runs table other than config/arm3_runs.tsv, for testing; with it
                    the stage is never marked complete
  --outdir DIR      where strains go instead of data/arm3 (testing)
  --strain-log FILE per-strain log instead of logs/arm3_strains.tsv (testing)
EOF
}

# shellcheck source=lib/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
# The arguments as given, kept for the measured re-run: the parsing loop below
# consumes them with shift.
args=("$@")
max_strains="" force=0 runs_arg="" outdir_arg="" log_arg=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-strains) max_strains="$2"; shift 2 ;;
        --force) force=1; shift ;;
        --runs) runs_arg="$2"; shift 2 ;;
        --outdir) outdir_arg="$2"; shift 2 ;;
        --strain-log) log_arg="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; pad_die "unknown option $1" ;;
    esac
done
pad_load_conf
stage="assemble_reads"
[[ -z "$runs_arg" ]] || stage="assemble_reads_test"
pad_skip_if_done "$stage" "$force"
pad_measure_self "$stage" "${args[@]}"
pad_activate tools

: "${SPADES_MEM_GB:=4}" "${HOURS_ARM3:=48}" "${DISK_CEILING_GB:=55}"
# SPAdes' -m caps its address space, not its resident memory, and each thread
# reserves address space of its own. With 16 threads a 4 GB cap ran out at 40x
# (mmap failed while SPAdes held 1.5 GB resident) for the first three strains;
# with 4 threads the same 40x reads of BAM assembled within the cap, at 1.5 GB
# resident, into contigs identical to a 16-thread run without any cap. One
# SPAdes thread per GB of the cap is the rule from the fourth strain on
# (config/analysis_plan.md, 9 October 2026). MEGAHIT keeps every thread.
SPADES_THREADS=$(( THREADS < ${SPADES_MEM_GB%.*} ? THREADS : ${SPADES_MEM_GB%.*} ))
(( SPADES_THREADS >= 1 )) || SPADES_THREADS=1
RUNS="${runs_arg:-$PAD_ROOT/config/arm3_runs.tsv}"
[[ -f "$RUNS" ]] || pad_die "$(pad_rel "$RUNS") not found; run scripts/09_build_pairs.py first"
STRAIN_LOG="${log_arg:-$PAD_ROOT/logs/arm3_strains.tsv}"
ARM3_DIR="${outdir_arg:-$DATA_DIR/arm3}"
LOG_HEADER=$'strain\trun\tstart_utc\tend_utc\telapsed_s\tdownload_s\tassembly_s\tpeak_strain_disk_gb\tmin_free_drive_gb\tdownload_bytes\tmd5_ok\tassemblies_ok\tassemblies_failed\tmax_assembler_rss_mb\tstatus'
[[ -f "$STRAIN_LOG" ]] || printf '%s\n' "$LOG_HEADER" > "$STRAIN_LOG"
[[ "$(head -1 "$STRAIN_LOG")" == "$LOG_HEADER" ]] || pad_die "$(pad_rel "$STRAIN_LOG") has other columns than this script writes"
PRIOR_DOWNLOAD_S=5400
PRIOR_ASSEMBLY_S=5400
PRIOR_DISK_GB=12
SEED_MAIN=11
SEEDS_10X=(22 33)
DEPTHS=(5 10 20 40 full)
FULL_CAP=80
mkdir -p "$ARM3_DIR"

# Column lookup by header name, so the order of columns in the TSV can change.
col() { awk -F'\t' -v k="$1" 'NR==1 {for (i=1;i<=NF;i++) if ($i==k) c=i; next} {print $c}' "$RUNS"; }
mapfile -t codes < <(paste <(col code) <(col arm3_eligible) <(col order_key) | awk -F'\t' '$2=="yes"' | sort -t$'\t' -k3,3 | cut -f1)
pad_log "arm 3: ${#codes[@]} eligible strains in the registered order; SPAdes with $SPADES_THREADS threads and a ${SPADES_MEM_GB} GB cap, MEGAHIT with $THREADS threads"

field() {   # field <code> <column>
    awk -F'\t' -v c="$1" -v k="$2" 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next} $h["code"]==c {print $h[k]}' "$RUNS"
}

# du exits non-zero when a file vanishes while it walks the tree, which other
# stages writing to the data directory make likely; the total is still right
# to within that file, so the status is ignored.
footprint_gb() { { du -sb "$DATA_DIR" 2>/dev/null || true; } | awk '{printf "%.2f", $1 / 1e9}'; }
drive_free_gb() {
    if grep -qi microsoft /proc/version && [[ -d /mnt/c ]] && [[ "$DATA_DIR" != /mnt/* ]]; then
        df -B1 --output=avail /mnt/c | tail -1 | awk '{printf "%.1f", $1 / 1e9}'
    else
        df -B1 --output=avail "$DATA_DIR" | tail -1 | awk '{printf "%.1f", $1 / 1e9}'
    fi
}

# Contiguity of one contig file: N50, L50, contig count, total length.
contig_stats() {
    python - "$1" <<'EOF'
import os
import sys
sys.path.insert(0, os.path.join(os.environ["PAD_ROOT"], "scripts", "lib"))
import padlib as P
ls = [len(s) for _, s in P.read_fasta(sys.argv[1])]
print("%d\t%d\t%d\t%d" % (P.nx(ls), P.lx(ls), len(ls), sum(ls)))
EOF
}

bases_of() {   # total bases in one or more FASTQ files
    seqkit stats -T -j 2 "$@" 2>/dev/null | awk -F'\t' 'NR>1 {s += $5} END {printf "%d", s}'
}

names_match() {   # prints "<pairs> <mismatches>" for two FASTQ files
    paste <(zcat "$1" | awk 'NR % 4 == 1 {n = $1; sub(/\/[12]$/, "", n); print n}') \
          <(zcat "$2" | awk 'NR % 4 == 1 {n = $1; sub(/\/[12]$/, "", n); print n}') |
        awk -F'\t' '$1 != $2 {bad++} END {printf "%d %d\n", NR, bad + 0}'
}

specs_of() {   # specs_of <G> <B>: one "target:seed:fraction" line per subsample
    local d f s
    for d in "${DEPTHS[@]}"; do
        if [[ "$d" == "full" ]]; then
            f=$(awk -v g="$1" -v b="$2" -v c="$FULL_CAP" 'BEGIN {x = c * g / b; if (x > 1) x = 1; printf "%.6f", x}')
        else
            f=$(awk -v g="$1" -v b="$2" -v d="$d" 'BEGIN {x = d * g / b; if (x > 1) x = 1; printf "%.6f", x}')
        fi
        echo "$d:$SEED_MAIN:$f"
        if [[ "$d" == "10" ]]; then
            for s in "${SEEDS_10X[@]}"; do echo "$d:$s:$f"; done
        fi
    done
}

# download_mate <m> <url> <md5> <bytes>: one read file into the strain's reads
# directory, checked against ENA's MD5. On the link these results came from,
# transfers of a few GB broke after 7 to 40 minutes, once with a TLS record
# that failed its integrity check ("bad record mac"), so a stream that had to
# start again from its first byte never finished. A file can be resumed from
# its last good byte, which TLS has verified; curl -C - does that, as often as
# needed while each try still adds bytes. The MD5 then checks the whole file,
# and a mismatch starts the file again, at most twice. Called by fetch, whose
# code and sdir it uses; writes reads/result_<m>: "ok <bytes>" or "failed 0".
download_mate() {
    local m="$1" url="$2" want="$3" size="$4" f="$sdir/reads/R$1.fastq.gz" whole tries have before got
    if ! [[ "$size" =~ ^[0-9]+$ ]]; then
        # no size in the runs table: ask the server for it
        size=$(curl -fsSIL "$url" | tr -d '\r' | awk 'tolower($1) == "content-length:" {n = $2} END {print n + 0}')
    fi
    for whole in 1 2 3; do
        rm -f "$f"
        tries=0
        have=0
        while (( have < size && tries < 200 )); do
            tries=$((tries + 1))
            before=$have
            curl -fsSL -C - -o "$f" "$url" 2>> "$sdir/reads/curl_$m.log" || true
            have=$(stat -c %s "$f" 2>/dev/null || echo 0)
            if (( have < size )); then
                pad_log "$code: mate $m stopped at $have of $size bytes; resuming"
                # a try that added nothing is counted ten times over, so that
                # a link that has gone for good gives up within 20 tries
                (( have > before )) || tries=$((tries + 9))
                sleep 10
            fi
        done
        got=$(md5sum "$f" 2>/dev/null | awk '{print $1}')
        if [[ "$have" == "$size" && "$got" == "$want" ]]; then
            echo "ok $have" > "$sdir/reads/result_$m"
            return 0
        fi
        pad_log "$code: mate $m, download $whole of 3: $have of $size bytes, MD5 $got, expected $want"
    done
    echo "failed 0" > "$sdir/reads/result_$m"
}

# subsample_mate <m>: one pass over the downloaded file of mate m, fanned out
# through named pipes to a seqtk for every depth and seed. The exit status of
# every process in the fan-out is checked, because a reader that stops early
# leaves a short subsample without any error message. I suspect that is how
# the streamed downloads broke here, with seqtk stopping at damaged data; I
# did not confirm it. A failed pass is run again from the file, which is
# deleted once a pass has succeeded. Writes reads/subsampled_<m>.
subsample_mate() {
    local m="$1" f="$sdir/reads/R$1.fastq.gz" pass spec t s fr fifo out ok
    local -a pids fifos
    for pass in 1 2 3; do
        rm -f "$sdir"/reads/fifo_*_"$m" "$sdir"/reads/d*_R"$m".fq.gz
        pids=()
        fifos=()
        for spec in "${specs[@]}"; do
            IFS=: read -r t s fr <<< "$spec"
            fifo="$sdir/reads/fifo_${t}_${s}_$m"
            mkfifo "$fifo"
            fifos+=("$fifo")
            out="$sdir/reads/d${t}_s${s}_R$m.fq.gz"
            if awk -v f="$fr" 'BEGIN {exit !(f >= 1)}'; then
                # seqtk reads a value of 1 or more as a read count, not a
                # fraction, so the whole run is copied as it comes
                cat "$fifo" > "$out" &
            else
                ( set -o pipefail; seqtk sample -s "$s" "$fifo" "$fr" | pigz -p 2 > "$out" ) &
            fi
            pids+=($!)
        done
        ok=1
        tee "${fifos[@]}" < "$f" > /dev/null || ok=0
        for p in "${pids[@]}"; do wait "$p" || ok=0; done
        if [[ "$ok" == "1" ]]; then
            rm -f "$f" "$sdir"/reads/fifo_*_"$m"
            touch "$sdir/reads/subsampled_$m"
            return 0
        fi
        pad_log "$code: subsampling pass $pass over mate $m failed; repeating it from the downloaded file"
    done
    return 1
}

# fetch <code>: both mates of one strain, downloaded at the same time (on this
# link two transfers together moved about 30 per cent more than one), then
# each subsampled in one pass. Ends by writing reads/fetched:
# "<reads intact: yes or no> <bytes> <seconds>".
fetch() {
    local code="$1" sdir="$ARM3_DIR/$1" t0 m res nbytes md5_ok="yes" dl_bytes=0 pid1 pid2
    local -a specs sizes
    mapfile -t specs < <(specs_of "$(field "$code" long_read_length)" "$(field "$code" base_count)")
    IFS=';' read -r -a sizes <<< "$(field "$code" fastq_bytes)"
    t0=$(date +%s)
    download_mate 1 "https://$(field "$code" fastq_1)" "$(field "$code" md5_1)" "${sizes[0]:-}" &
    pid1=$!
    download_mate 2 "https://$(field "$code" fastq_2)" "$(field "$code" md5_2)" "${sizes[1]:-}" &
    pid2=$!
    wait "$pid1" || true
    wait "$pid2" || true
    for m in 1 2; do
        read -r res nbytes < "$sdir/reads/result_$m" || { res="failed"; nbytes=0; }
        [[ "$res" == "ok" ]] || md5_ok="no"
        dl_bytes=$((dl_bytes + nbytes))
    done
    if [[ "$md5_ok" == "yes" ]]; then
        for m in 1 2; do
            subsample_mate "$m" || md5_ok="no"
        done
    fi
    rm -f "$sdir"/reads/result_* "$sdir"/reads/R?.fastq.gz
    printf '%s %s %s\n' "$md5_ok" "$dl_bytes" "$(( $(date +%s) - t0 ))" > "$sdir/reads/fetched"
}

fetch_pid="" fetch_code="" assembling_dir="" sampler_pid=""

# start_strain <code>: fresh directories, then the download in the background.
start_strain() {
    local sdir="$ARM3_DIR/$1"
    rm -rf "$sdir"
    mkdir -p "$sdir/reads" "$sdir/work" "$sdir/assemblies" "$sdir/fastp"
    pad_utc > "$sdir/work/start_utc"
    pad_log "$1: download started, run $(field "$1" run_accession), $(field "$1" base_count) bases, genome $(field "$1" long_read_length) bp"
    fetch "$1" &
    fetch_pid=$!
    fetch_code="$1"
}

# A process and everything it started. Each process is stopped before its
# children are listed, so none can start another in between.
kill_tree() {
    local c
    kill -STOP "$1" 2>/dev/null || return 0
    for c in $(pgrep -P "$1" 2>/dev/null); do kill_tree "$c"; done
    kill -TERM "$1" 2>/dev/null || true
    kill -CONT "$1" 2>/dev/null || true
}

cleanup() {
    local status=$? d
    [[ -n "$sampler_pid" ]] && kill "$sampler_pid" 2>/dev/null || true
    [[ -n "$fetch_pid" ]] && kill_tree "$fetch_pid"
    for d in "$assembling_dir" "${fetch_code:+$ARM3_DIR/$fetch_code}"; do
        [[ -n "$d" ]] || continue
        # the same strain can be both the one assembling and the one being
        # fetched; it is cleaned, and logged, once
        [[ -d "$d/reads" || -d "$d/work" ]] || continue
        rm -rf "$d/reads" "$d/work"
        if [[ $status -ne 0 ]]; then
            pad_log "interrupted; removed the reads and working directories of $(pad_rel "$d")"
        fi
    done
}
trap cleanup EXIT
trap 'exit 130' INT TERM

# The wall time arm 3 has used: the union of the logged strains' intervals and
# of this run so far, so that overlapping strains are not counted twice and
# earlier runs of this script are.
wall_spent() {
    {
        awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} {print $h["start_utc"], $h["end_utc"]}' "$STRAIN_LOG" |
            while read -r a b; do echo "$(date -d "$a" +%s) $(date -d "$b" +%s)"; done
        echo "$run_start $(date +%s)"
    } | sort -n -k1,1 | awk 'NR == 1 {s = $1; e = $2; next}
        $1 > e {t += e - s; s = $1; e = $2; next}
        $2 > e {e = $2}
        END {printf "%d\n", t + e - s}'
}

# The longest download and assembly, and the largest strain on disk, among the
# strains complete so far, or the priors before the first.
estimates() {
    awk -F'\t' -v pf="$PRIOR_DOWNLOAD_S" -v pa="$PRIOR_ASSEMBLY_S" -v pd="$PRIOR_DISK_GB" '
        NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
        $h["status"] == "complete" {
            n++
            if ($h["download_s"] + 0 > f) f = $h["download_s"] + 0
            if ($h["assembly_s"] + 0 > a) a = $h["assembly_s"] + 0
            if ($h["peak_strain_disk_gb"] + 0 > d) d = $h["peak_strain_disk_gb"] + 0
        }
        END {printf "%d %d %.2f\n", (n ? f : pf), (n ? a : pa), (n ? d : pd)}' "$STRAIN_LOG"
}

# fits <code> <overlapped: 0 or 1>: whether the strain's download may start
# now. Overlapped means another strain is about to assemble while it
# downloads: that strain finishes after max(download, assembly) and this one
# an assembly later, and both need their disk at once.
fits() {
    local code="$1" overlapped="$2" f a pk spent need fp short
    read -r f a pk < <(estimates)
    spent=$(wall_spent)
    if [[ "$overlapped" == "1" ]]; then
        need=$(( (f > a ? f : a) + a ))
    else
        need=$(( f + a ))
    fi
    if (( spent + need > budget )); then
        short=$(awk -v s="$spent" -v n="$need" -v b="$budget" 'BEGIN {printf "%.1f", (s + n - b) / 3600}')
        pad_log "refusing to start $code: it would take arm 3 $short h past its ${HOURS_ARM3} h budget; strains from here on are cut"
        return 1
    fi
    fp=$(footprint_gb)
    if awk -v f="$fp" -v p="$pk" -v n="$overlapped" -v c="$DISK_CEILING_GB" 'BEGIN {exit !(f + (n + 1) * p > c)}'; then
        short=$(awk -v f="$fp" -v p="$pk" -v n="$overlapped" -v c="$DISK_CEILING_GB" 'BEGIN {printf "%.1f", f + (n + 1) * p - c}')
        pad_log "refusing to start $code: projected footprint is $short GB over the ${DISK_CEILING_GB} GB ceiling; strains from here on are cut"
        return 1
    fi
}

done_strains=0
todo=()
for code in "${codes[@]}"; do
    if [[ -f "$ARM3_DIR/$code/done" && "$force" != "1" ]]; then
        pad_log "$code: already assembled, skipping"
        done_strains=$((done_strains + 1))
    else
        todo+=("$code")
    fi
done
if [[ -n "$max_strains" ]]; then
    allowed=$(( max_strains > done_strains ? max_strains - done_strains : 0 ))
    (( allowed >= ${#todo[@]} )) || pad_log "stopping after $max_strains strains as asked"
    todo=("${todo[@]:0:allowed}")
fi
budget=$(awk -v h="$HOURS_ARM3" 'BEGIN {printf "%.0f", h * 3600}')
run_start=$(date +%s)

# Disk sampler: every 10 s, each strain that has a working directory gets the
# size of its whole directory, and its peak is kept in work/.peak_bytes.
( while true; do
      for s in "$ARM3_DIR"/*/; do
          [[ -d "$s/work" ]] || continue
          b=$({ du -sb "$s" 2>/dev/null || true; } | awk '{print $1}')
          p=$(cat "$s/work/.peak_bytes" 2>/dev/null || echo 0)
          if [[ -n "$b" && "$b" -gt "${p:-0}" ]]; then
              { echo "$b" > "$s/work/.peak_bytes"; } 2>/dev/null || true
          fi
      done
      sleep 10
  done ) &
sampler_pid=$!

next=0
if (( ${#todo[@]} > 0 )) && fits "${todo[0]}" 0; then
    start_strain "${todo[0]}"
    next=1
fi
done_this_run=0
while [[ -n "$fetch_code" ]]; do
    code="$fetch_code"
    sdir="$ARM3_DIR/$code"
    assembling_dir="$sdir"
    wait "$fetch_pid" || true
    fetch_pid="" fetch_code=""
    read -r md5_ok dl_bytes dl_s < "$sdir/reads/fetched" || { md5_ok="no"; dl_bytes=0; dl_s=0; }
    start_utc=$(cat "$sdir/work/start_utc")
    pad_log "$code: download finished in ${dl_s} s, md5 ${md5_ok}"
    if (( next < ${#todo[@]} )) && fits "${todo[next]}" 1; then
        start_strain "${todo[next]}"
        next=$((next + 1))
    else
        next=${#todo[@]}
    fi

    run=$(field "$code" run_accession)
    G=$(field "$code" long_read_length)
    mapfile -t specs < <(specs_of "$G" "$(field "$code" base_count)")
    t_asm=$(date +%s)
    printf 'assembler\ttarget_depth\tseed\tfraction\tread_pairs\tbases_raw\trealised_depth_raw\tbases_trimmed\trealised_depth_trimmed\tnames_checked\tname_mismatches\tstatus\terror\telapsed_s\tpeak_rss_mb\tcontig_n50\tcontig_l50\tn_contigs\ttotal_length\tthreads\n' > "$sdir/manifest.tsv"
    n_ok=0 n_fail=0 max_rss=0
    for spec in "${specs[@]}"; do
        IFS=: read -r t s f <<< "$spec"
        r1="$sdir/reads/d${t}_s${s}_R1.fq.gz"
        r2="$sdir/reads/d${t}_s${s}_R2.fq.gz"
        if [[ "$md5_ok" != "yes" ]]; then
            for a in spades megahit; do
                printf '%s\t%s\t%s\t%s\t\t\t\t\t\t\t\tfailed\treads not downloaded intact\t\t\t\t\t\t\t\n' "$a" "$t" "$s" "$f" >> "$sdir/manifest.tsv"
                n_fail=$((n_fail + 1))
            done
            continue
        fi
        read -r pairs mism < <(names_match "$r1" "$r2")
        raw=$(bases_of "$r1" "$r2")
        tag="d${t}_s${s}"
        t1="$sdir/work/${tag}_t1.fq.gz"
        t2="$sdir/work/${tag}_t2.fq.gz"
        fastp -i "$r1" -I "$r2" -o "$t1" -O "$t2" --detect_adapter_for_pe -w 8 \
            -j "$sdir/fastp/$tag.json" -h /dev/null 2> "$sdir/fastp/$tag.log"
        trimmed=$(bases_of "$t1" "$t2")
        for a in spades megahit; do
            wd="$sdir/work/${a}_$tag"
            tm="$sdir/work/${a}_$tag.time"
            err=""
            rm -rf "$wd"
            if [[ "$a" == "spades" ]]; then
                /usr/bin/time -f "%e %M" -o "$tm" spades.py --isolate -1 "$t1" -2 "$t2" -o "$wd" -t "$SPADES_THREADS" \
                    -m "$SPADES_MEM_GB" --tmp-dir "$sdir/work/tmp_spades" > "$sdir/work/$a.$tag.log" 2>&1 || err="SPAdes exited with an error"
                contigs="$wd/contigs.fasta"
            else
                /usr/bin/time -f "%e %M" -o "$tm" megahit -1 "$t1" -2 "$t2" -o "$wd" -t "$THREADS" \
                    -m "$((SPADES_MEM_GB * 1000000000))" > "$sdir/work/$a.$tag.log" 2>&1 || err="MEGAHIT exited with an error"
                contigs="$wd/final.contigs.fa"
            fi
            read -r el rss < <(tail -1 "$tm" 2>/dev/null || echo "0 0")
            rss_mb=$((rss / 1024))
            (( rss_mb > max_rss )) && max_rss=$rss_mb
            if [[ -z "$err" && -s "$contigs" ]]; then
                pigz -c "$contigs" > "$sdir/assemblies/${a}_$tag.fna.gz"
                read -r n50 l50 nc tot < <(contig_stats "$sdir/assemblies/${a}_$tag.fna.gz")
                st="ok"
                n_ok=$((n_ok + 1))
            else
                [[ -n "$err" ]] || err="no contigs written"
                # SPAdes ends with a generic line; the first ERROR line it
                # logged says what went wrong
                why=$(grep -m1 ' ERROR ' "$sdir/work/$a.$tag.log" | sed 's/.*) *//')
                [[ -n "$why" ]] || why=$(grep -iE 'error|memory|killed' "$sdir/work/$a.$tag.log" | tail -1)
                err="$err: $(printf '%s' "$why" | pad_scrub | tr '\t' ' ' | cut -c1-200)"
                n50="" l50="" nc="" tot="" st="failed"
                n_fail=$((n_fail + 1))
            fi
            threads_used=$THREADS
            [[ "$a" == "spades" ]] && threads_used=$SPADES_THREADS
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$a" "$t" "$s" "$f" \
                "$pairs" "$raw" "$(awk -v x="$raw" -v g="$G" 'BEGIN {printf "%.2f", x / g}')" "$trimmed" \
                "$(awk -v x="$trimmed" -v g="$G" 'BEGIN {printf "%.2f", x / g}')" "$pairs" "$mism" "$st" "$err" \
                "$el" "$rss_mb" "$n50" "$l50" "$nc" "$tot" "$threads_used" >> "$sdir/manifest.tsv"
            rm -rf "$wd" "$sdir/work/tmp_spades"
            pad_log "$code: $a $tag $st (${el}s, ${rss_mb} MB)"
        done
        rm -f "$t1" "$t2" "$r1" "$r2"
    done

    asm_s=$(( $(date +%s) - t_asm ))
    peak_gb=$(awk '{printf "%.2f", $1 / 1e9}' "$sdir/work/.peak_bytes" 2>/dev/null || echo "")
    rm -rf "$sdir/reads" "$sdir/work"
    assembling_dir=""
    status="complete"
    [[ "$md5_ok" == "yes" ]] || status="download failed"
    end_utc=$(pad_utc)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$code" "$run" "$start_utc" "$end_utc" \
        "$(( $(date -d "$end_utc" +%s) - $(date -d "$start_utc" +%s) ))" "$dl_s" "$asm_s" "$peak_gb" "$(drive_free_gb)" \
        "$dl_bytes" "$md5_ok" "$n_ok" "$n_fail" "$max_rss" "$status" >> "$STRAIN_LOG"
    pad_log "$code: $n_ok assemblies, $n_fail failed; download ${dl_s} s, assembly ${asm_s} s, peak ${peak_gb} GB"
    if [[ "$status" != "complete" ]]; then
        # no done marker: the reads never arrived intact, which says nothing
        # about the strain, so a later run tries it again
        pad_log "$code: reads not downloaded intact; left for a later run"
        continue
    fi
    pad_utc > "$sdir/done"
    done_strains=$((done_strains + 1))
    done_this_run=$((done_this_run + 1))
    if [[ "$done_this_run" == "1" ]]; then
        fit=$(awk -v b="$budget" -v f="$dl_s" -v a="$asm_s" -v n="${#codes[@]}" 'BEGIN {
            m = (f > a ? f : a); k = (f + a > b) ? 0 : 1 + int((b - f - a) / (m > 0 ? m : 1)); if (k > n) k = n; printf "%d", k}')
        pad_log "at the pace of the first strain, with downloads overlapping assembly, the ${HOURS_ARM3} h budget fits about $fit of ${#codes[@]} strains"
    fi
done
kill "$sampler_pid" 2>/dev/null || true
wait "$sampler_pid" 2>/dev/null || true
sampler_pid=""
pad_log "arm 3 assembly finished: $done_strains strains"
if [[ -z "$runs_arg" ]] && (( done_strains == ${#codes[@]} )); then pad_mark_done assemble_reads; fi
