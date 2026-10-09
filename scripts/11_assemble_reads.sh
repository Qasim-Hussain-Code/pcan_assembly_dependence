#!/usr/bin/env bash
# Arm 3: assemble real Illumina reads at known depths with two assemblers, one
# strain at a time, by the design registered in config/analysis_plan.md.
# For each strain: stream both read files from ENA once, subsample during the
# download (all depths and seeds in the same pass), check that mates still
# pair, trim with fastp, assemble every subsample with SPAdes and MEGAHIT, keep
# the contigs and reports, delete the reads.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash scripts/11_assemble_reads.sh [--max-strains N] [--force]

Strains come from config/arm3_runs.tsv (arm3_eligible = yes), in ascending
order_key, the order fixed in the analysis plan. Per strain, kept in
data/arm3/<code>/: assemblies/*.fna.gz (contigs, gzip), fastp/*.json,
manifest.tsv. Reads and assembler working directories are deleted before the
next strain starts, and by the cleanup trap if a strain is interrupted.

Before each strain the time and disk it will need are projected from the
strains already done (from priors of 3 h and 12 GB before the first) and the
strain is refused, with the shortfall, if it would exceed HOURS_ARM3 or the
disk ceiling in project.conf. Strains are cut from the end of the order;
depths, seeds and controls never are.

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
RUNS="${runs_arg:-$PAD_ROOT/config/arm3_runs.tsv}"
[[ -f "$RUNS" ]] || pad_die "$(pad_rel "$RUNS") not found; run scripts/09_build_pairs.py first"
STRAIN_LOG="${log_arg:-$PAD_ROOT/logs/arm3_strains.tsv}"
ARM3_DIR="${outdir_arg:-$DATA_DIR/arm3}"
[[ -f "$STRAIN_LOG" ]] || printf 'strain\trun\tstart_utc\tend_utc\telapsed_s\tpeak_strain_disk_gb\tmin_free_drive_gb\tdownload_bytes\tmd5_ok\tassemblies_ok\tassemblies_failed\tmax_assembler_rss_mb\tstatus\n' > "$STRAIN_LOG"
PRIOR_HOURS=3
PRIOR_DISK_GB=12
SEED_MAIN=11
SEEDS_10X=(22 33)
DEPTHS=(5 10 20 40 full)
FULL_CAP=80

# Column lookup by header name, so the order of columns in the TSV can change.
col() { awk -F'\t' -v k="$1" 'NR==1 {for (i=1;i<=NF;i++) if ($i==k) c=i; next} {print $c}' "$RUNS"; }
mapfile -t codes < <(paste <(col code) <(col arm3_eligible) <(col order_key) | awk -F'\t' '$2=="yes"' | sort -t$'\t' -k3,3 | cut -f1)
pad_log "arm 3: ${#codes[@]} eligible strains in the registered order"

field() {   # field <code> <column>
    awk -F'\t' -v c="$1" -v k="$2" 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next} $h["code"]==c {print $h[k]}' "$RUNS"
}

current_dir=""
sampler_pid=""
cleanup() {
    local status=$?
    [[ -n "$sampler_pid" ]] && kill "$sampler_pid" 2>/dev/null || true
    if [[ -n "$current_dir" ]]; then
        rm -rf "$current_dir/reads" "$current_dir/work"
        if [[ $status -ne 0 ]]; then
            pad_log "interrupted; removed the reads and working directories of $(pad_rel "$current_dir")"
        fi
    fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM

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

done_strains=0
for code in "${codes[@]}"; do
    if [[ -n "$max_strains" && "$done_strains" -ge "$max_strains" ]]; then
        pad_log "stopping after $max_strains strains as asked"
        break
    fi
    sdir="$ARM3_DIR/$code"
    if [[ -f "$sdir/done" && "$force" != "1" ]]; then
        pad_log "$code: already assembled, skipping"
        done_strains=$((done_strains + 1))
        continue
    fi

    # Budget: time from the strains logged so far, disk from their peaks.
    read -r spent per_strain peak_disk < <(awk -F'\t' -v ph="$PRIOR_HOURS" -v pd="$PRIOR_DISK_GB" '
        NR > 1 && $13 == "complete" {s += $5; if ($5 > m) m = $5; if ($6 > d) d = $6; n++}
        END {printf "%.0f %.0f %.2f\n", s, (n ? m : ph * 3600), (n ? d : pd)}' "$STRAIN_LOG")
    budget=$(awk -v h="$HOURS_ARM3" 'BEGIN {printf "%.0f", h * 3600}')
    if (( spent + per_strain > budget )); then
        short=$(awk -v s="$spent" -v p="$per_strain" -v b="$budget" 'BEGIN {printf "%.1f", (s + p - b) / 3600}')
        pad_log "refusing to start $code: it would take arm 3 $short h past its ${HOURS_ARM3} h budget; strains from here on are cut"
        break
    fi
    fp=$(footprint_gb)
    if awk -v f="$fp" -v p="$peak_disk" -v c="$DISK_CEILING_GB" 'BEGIN {exit !(f + p > c)}'; then
        short=$(awk -v f="$fp" -v p="$peak_disk" -v c="$DISK_CEILING_GB" 'BEGIN {printf "%.1f", f + p - c}')
        pad_log "refusing to start $code: projected footprint is $short GB over the ${DISK_CEILING_GB} GB ceiling"
        break
    fi

    run=$(field "$code" run_accession)
    G=$(field "$code" long_read_length)
    B=$(field "$code" base_count)
    url1="https://$(field "$code" fastq_1)"
    url2="https://$(field "$code" fastq_2)"
    md5_1=$(field "$code" md5_1)
    md5_2=$(field "$code" md5_2)
    pad_log "$code: run $run, $B bases, genome $G bp ($(awk -v b="$B" -v g="$G" 'BEGIN {printf "%.0f", b / g}')x available)"

    rm -rf "$sdir"
    mkdir -p "$sdir/reads" "$sdir/work" "$sdir/assemblies" "$sdir/fastp"
    current_dir="$sdir"
    start_utc=$(pad_utc)
    t0=$(date +%s)
    # Disk sampler for this strain's peak, every 10 s.
    ( peak=0; while true; do
          b=$({ du -sb "$sdir" 2>/dev/null || true; } | awk '{print $1}'); [[ -n "$b" && "$b" -gt "$peak" ]] && peak=$b
          echo "$peak" > "$sdir/peak_bytes"; sleep 10
      done ) &
    sampler_pid=$!

    # Subsample specifications: target, seed, fraction of the run.
    specs=()
    for d in "${DEPTHS[@]}"; do
        if [[ "$d" == "full" ]]; then
            f=$(awk -v g="$G" -v b="$B" -v c="$FULL_CAP" 'BEGIN {x = c * g / b; if (x > 1) x = 1; printf "%.6f", x}')
        else
            f=$(awk -v g="$G" -v b="$B" -v d="$d" 'BEGIN {x = d * g / b; if (x > 1) x = 1; printf "%.6f", x}')
        fi
        specs+=("$d:$SEED_MAIN:$f")
        if [[ "$d" == "10" ]]; then
            for s in "${SEEDS_10X[@]}"; do specs+=("$d:$s:$f"); done
        fi
    done

    # One download per mate, fanned out to every subsampler through named
    # pipes. curl does not retry inside a stream (a retried transfer would
    # repeat data), so a failed download is repeated whole, up to three times.
    # The two mates stream at the same time: on the link these results came
    # from, two streams together moved about 30 per cent more than one.
    stream_mate() {
        local m="$1" url="$2" want="$3" attempt spec t s f fifo out ok got
        local -a pids fifos
        for attempt in 1 2 3; do
            rm -f "$sdir"/reads/fifo_*_"$m" "$sdir"/reads/*_R"$m".fq.gz
            pids=()
            fifos=()
            for spec in "${specs[@]}"; do
                IFS=: read -r t s f <<< "$spec"
                fifo="$sdir/reads/fifo_${t}_${s}_$m"
                mkfifo "$fifo"
                fifos+=("$fifo")
                out="$sdir/reads/d${t}_s${s}_R$m.fq.gz"
                if awk -v f="$f" 'BEGIN {exit !(f >= 1)}'; then
                    # seqtk reads a value of 1 or more as a read count, not a
                    # fraction, so the whole run is copied as it comes
                    cat "$fifo" > "$out" &
                else
                    ( seqtk sample -s "$s" "$fifo" "$f" | pigz -p 2 > "$out" ) &
                fi
                pids+=($!)
            done
            mkfifo "$sdir/reads/fifo_md5_$m"
            ( md5sum < "$sdir/reads/fifo_md5_$m" | awk '{print $1}' > "$sdir/reads/md5_$m" ) &
            pids+=($!)
            mkfifo "$sdir/reads/fifo_count_$m"
            ( wc -c < "$sdir/reads/fifo_count_$m" > "$sdir/reads/bytes_$m" ) &
            pids+=($!)
            ok=1
            curl -fsSL "$url" | tee "${fifos[@]}" "$sdir/reads/fifo_md5_$m" > "$sdir/reads/fifo_count_$m" || ok=0
            for p in "${pids[@]}"; do wait "$p" || ok=0; done
            got=$(cat "$sdir/reads/md5_$m" 2>/dev/null || echo none)
            if [[ "$ok" == "1" && "$got" == "$want" ]]; then
                echo "ok $(cat "$sdir/reads/bytes_$m")" > "$sdir/reads/result_$m"
                return 0
            fi
            pad_log "$code: mate $m attempt $attempt failed (md5 $got, expected $want)"
        done
        echo "failed 0" > "$sdir/reads/result_$m"
    }
    stream_mate 1 "$url1" "$md5_1" &
    pid1=$!
    stream_mate 2 "$url2" "$md5_2" &
    pid2=$!
    wait "$pid1" || true
    wait "$pid2" || true
    md5_ok="yes"
    dl_bytes=0
    for m in 1 2; do
        read -r res nbytes < "$sdir/reads/result_$m" || { res="failed"; nbytes=0; }
        [[ "$res" == "ok" ]] || md5_ok="no"
        dl_bytes=$((dl_bytes + nbytes))
    done
    rm -f "$sdir"/reads/fifo_* "$sdir"/reads/result_*

    printf 'assembler\ttarget_depth\tseed\tfraction\tread_pairs\tbases_raw\trealised_depth_raw\tbases_trimmed\trealised_depth_trimmed\tnames_checked\tname_mismatches\tstatus\terror\telapsed_s\tpeak_rss_mb\tcontig_n50\tcontig_l50\tn_contigs\ttotal_length\n' > "$sdir/manifest.tsv"
    n_ok=0 n_fail=0 max_rss=0
    for spec in "${specs[@]}"; do
        IFS=: read -r t s f <<< "$spec"
        r1="$sdir/reads/d${t}_s${s}_R1.fq.gz"
        r2="$sdir/reads/d${t}_s${s}_R2.fq.gz"
        if [[ "$md5_ok" != "yes" ]]; then
            for a in spades megahit; do
                printf '%s\t%s\t%s\t%s\t\t\t\t\t\t\t\tfailed\tdownload failed its md5 check\t\t\t\t\t\t\n' "$a" "$t" "$s" "$f" >> "$sdir/manifest.tsv"
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
                /usr/bin/time -f "%e %M" -o "$tm" spades.py --isolate -1 "$t1" -2 "$t2" -o "$wd" -t "$THREADS" \
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
                err="$err: $(grep -iE 'error|memory|killed' "$sdir/work/$a.$tag.log" | tail -1 | pad_scrub | tr '\t' ' ' | cut -c1-200)"
                n50="" l50="" nc="" tot="" st="failed"
                n_fail=$((n_fail + 1))
            fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$a" "$t" "$s" "$f" \
                "$pairs" "$raw" "$(awk -v x="$raw" -v g="$G" 'BEGIN {printf "%.2f", x / g}')" "$trimmed" \
                "$(awk -v x="$trimmed" -v g="$G" 'BEGIN {printf "%.2f", x / g}')" "$pairs" "$mism" "$st" "$err" \
                "$el" "$rss_mb" "$n50" "$l50" "$nc" "$tot" >> "$sdir/manifest.tsv"
            rm -rf "$wd" "$sdir/work/tmp_spades"
            pad_log "$code: $a $tag $st (${el}s, ${rss_mb} MB)"
        done
        rm -f "$t1" "$t2" "$r1" "$r2"
    done

    kill "$sampler_pid" 2>/dev/null || true
    wait "$sampler_pid" 2>/dev/null || true
    sampler_pid=""
    peak_gb=$(awk '{printf "%.2f", $1 / 1e9}' "$sdir/peak_bytes" 2>/dev/null || echo "")
    rm -rf "$sdir/reads" "$sdir/work" "$sdir/peak_bytes"
    current_dir=""
    status="complete"
    [[ "$md5_ok" == "yes" ]] || status="download failed"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$code" "$run" "$start_utc" "$(pad_utc)" \
        "$(( $(date +%s) - t0 ))" "$peak_gb" "$(drive_free_gb)" "$dl_bytes" "$md5_ok" "$n_ok" "$n_fail" "$max_rss" "$status" >> "$STRAIN_LOG"
    date -u +%Y-%m-%dT%H:%M:%SZ > "$sdir/done"
    done_strains=$((done_strains + 1))
    pad_log "$code: $n_ok assemblies, $n_fail failed, $(( $(date +%s) - t0 )) s, peak ${peak_gb} GB"
    if [[ "$done_strains" == "1" ]]; then
        el=$(( $(date +%s) - t0 ))
        fit=$(awk -v b="$budget" -v e="$el" 'BEGIN {printf "%d", b / e}')
        pad_log "first strain took ${el} s; at that pace the ${HOURS_ARM3} h budget fits about $fit of ${#codes[@]} strains"
    fi
done
pad_log "arm 3 assembly finished: $done_strains strains"
if [[ -z "$runs_arg" ]] && (( done_strains == ${#codes[@]} )); then pad_mark_done assemble_reads; fi
