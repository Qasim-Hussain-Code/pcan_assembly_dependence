#!/usr/bin/env bash
# Download assemblies by exact accession and version, cache them compressed with
# their NCBI sequence reports, and record contiguity statistics for each.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash scripts/03_fetch_assemblies.sh [--set arm0|arm2] [--force]

  --set arm0   every assembly in config/species_arm0.tsv whose arm 0 role is not
               "excluded" (default)
  --set arm2   every long-read assembly in config/arm2_pairs.tsv

For each accession the NCBI Datasets report is read first. A version NCBI
marks "previous" (replaced by a newer version) or "suppressed" is not
downloaded and nothing is used in its place; it is listed in
results/<set>/excluded_assemblies.tsv with its status.

Writes
  data/assemblies/<accession>.fna.gz                 genome, gzip
  data/assemblies/<accession>.sequence_report.jsonl  NCBI sequence roles and names
  data/assemblies/<accession>.report.json            NCBI dataset report
  logs/assemblies_<set>.tsv                          checksum and date per assembly
  results/<set>/assembly_stats.tsv                   contig N50, L50, counts, NCBI level
EOF
}

# shellcheck source=lib/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
# The arguments as given, kept for the measured re-run: the parsing loop below
# consumes them with shift.
args=("$@")
set_name="arm0"
force=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --set) set_name="$2"; shift 2 ;;
        --force) force=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; pad_die "unknown option $1" ;;
    esac
done
case "$set_name" in arm0|arm2) ;; *) pad_die "--set must be arm0 or arm2" ;; esac
pad_load_conf
pad_skip_if_done "fetch_assemblies_$set_name" "$force"
pad_measure_self "fetch_assemblies_$set_name" "${args[@]}"
pad_activate tools

asm="$DATA_DIR/assemblies"
mkdir -p "$DATA_DIR/work"
work="$(mktemp -d "$DATA_DIR/work/fetch.XXXXXX")"
trap 'rm -rf "$work"' EXIT
mkdir -p "$asm" "$PAD_ROOT/results/$set_name"

case "$set_name" in
    arm0) list="$PAD_ROOT/config/species_arm0.tsv"
          mapfile -t accs < <(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next} $h["arm0_role"] != "excluded" {print $h["accession"]}' "$list") ;;
    # The short-read side of every pair comes from the Peter et al. 2018
    # archive, which 09_build_pairs.py fetches and unpacks; only the long-read
    # assemblies have accessions to fetch here.
    arm2) list="$PAD_ROOT/config/arm2_pairs.tsv"
          [[ -f "$list" ]] || pad_die "config/arm2_pairs.tsv not found; run 09_build_pairs.py first"
          mapfile -t accs < <(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next} {print $h["long_read_accession"]}' "$list" | sort -u) ;;
esac
pad_log "$set_name: ${#accs[@]} accessions to check"

manifest="$PAD_ROOT/logs/assemblies_${set_name}.tsv"
excluded="$PAD_ROOT/results/$set_name/excluded_assemblies.tsv"
printf 'accession\tstatus\tcurrent_accession\tassembly_level\tbytes_gz\tsha256_gz\tdownloaded_utc\n' > "$manifest.tmp"
printf 'accession\tstatus\tcurrent_accession\treason\n' > "$excluded.tmp"

API="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession"
for acc in "${accs[@]}"; do
    report="$asm/$acc.report.json"
    if [[ ! -s "$report" || "$force" == "1" ]]; then
        # filters.assembly_version=all_assemblies makes the API answer for this
        # exact version even when a newer one exists, instead of the latest.
        curl -fsS --retry 5 --retry-delay 10 "$API/$acc/dataset_report?filters.assembly_version=all_assemblies" > "$work/report.json"
        mv "$work/report.json" "$report"
        sleep 0.4   # stay under NCBI's limit of three requests per second without a key
    fi
    read -r status current level < <(python - "$report" "$acc" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
rep = [r for r in d.get("reports", []) if r.get("accession") == sys.argv[2]]
if not rep:
    print("not_found", "-", "-")
else:
    r = rep[0]
    ai = r.get("assembly_info", {})
    print(ai.get("assembly_status", "unknown"), r.get("current_accession", "-"), ai.get("assembly_level", "-").replace(" ", "_"))
EOF
)
    if [[ "$status" != "current" ]]; then
        case "$status" in
            previous)  why="replaced by $current; nothing used in its place" ;;
            suppressed) why="suppressed by NCBI; nothing used in its place" ;;
            *)         why="NCBI status $status; nothing used in its place" ;;
        esac
        printf '%s\t%s\t%s\t%s\n' "$acc" "$status" "$current" "$why" >> "$excluded.tmp"
        pad_log "$acc: $why"
        continue
    fi
    fna="$asm/$acc.fna.gz"
    if [[ ! -s "$fna" || "$force" == "1" ]]; then
        rm -rf "$work/dl" && mkdir -p "$work/dl"
        datasets download genome accession "$acc" --include genome,seq-report \
            --filename "$work/dl/pkg.zip" --no-progressbar > /dev/null 2>&1 \
            || { sleep 30; datasets download genome accession "$acc" --include genome,seq-report \
                    --filename "$work/dl/pkg.zip" --no-progressbar > /dev/null; }
        unzip -q -o "$work/dl/pkg.zip" -d "$work/dl"
        mapfile -t fnas < <(find "$work/dl/ncbi_dataset/data/$acc" -name '*_genomic.fna')
        [[ ${#fnas[@]} -eq 1 ]] || pad_die "$acc: expected one genomic FASTA in the package, found ${#fnas[@]}"
        cp "$work/dl/ncbi_dataset/data/$acc/sequence_report.jsonl" "$asm/$acc.sequence_report.jsonl.part"
        pigz -p "$THREADS" -c "${fnas[0]}" > "$fna.part"
        mv "$asm/$acc.sequence_report.jsonl.part" "$asm/$acc.sequence_report.jsonl"
        mv "$fna.part" "$fna"
        rm -rf "$work/dl"
        stamp="$(pad_utc)"
        echo "$stamp" > "$asm/$acc.downloaded_utc"
        pad_log "$acc: downloaded ($(stat -c %s "$fna") bytes gzip)"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$acc" "$status" "$current" "$level" "$(stat -c %s "$fna")" \
        "$(sha256sum "$fna" | awk '{print $1}')" "$(cat "$asm/$acc.downloaded_utc" 2>/dev/null || echo unknown)" >> "$manifest.tmp"
done
mv "$manifest.tmp" "$manifest"
mv "$excluded.tmp" "$excluded"

python "$PAD_ROOT/scripts/lib/assembly_stats.py" --manifest "$manifest" --out "$PAD_ROOT/results/$set_name/assembly_stats.tsv"
pad_mark_done "fetch_assemblies_$set_name"
pad_log "$set_name: $(($(wc -l < "$manifest") - 1)) assemblies cached, $(($(wc -l < "$excluded") - 1)) excluded"
