#!/usr/bin/env bash
# One assembly in, one call table out. This is the only place PCAn is invoked.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash scripts/04_run_pcan.sh --assembly FILE --genus GENUS --out CALLS.tsv
                                   [--species EPITHET] [--label ID] [--keep-fimo DIR] [--force]

  --assembly   FASTA, plain or gzip
  --genus      one of the genera in PCAn v1.0's list, e.g. Saccharomyces
  --species    the species epithet, required when --genus is Kazachstania
  --out        call table to write (TSV); CALLS.tsv.meta.tsv records the run
  --label      identifier copied into every row of the call table
  --keep-fimo  directory to keep both FIMO passes in, gzip
  --force      run again even if --out exists

PCAn runs exactly as released at commit a5fa46f with patches/ applied, with
the motifs and thresholds its own table assigns to the genus. The BLAST
synteny step is declined; it does not change the calls (see
scripts/lib/pcan_runner.py). Exit status 2 means PCAn ran and failed; the
meta file then holds the end of its error output.
EOF
}

# shellcheck source=lib/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
assembly="" genus="" species="" out="" label="" keep="" force=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)  assembly="$2"; shift 2 ;;
        --genus)     genus="$2"; shift 2 ;;
        --species)   species="$2"; shift 2 ;;
        --out)       out="$2"; shift 2 ;;
        --label)     label="$2"; shift 2 ;;
        --keep-fimo) keep="$2"; shift 2 ;;
        --force)     force=1; shift ;;
        -h|--help)   usage; exit 0 ;;
        *)           usage >&2; pad_die "unknown option $1" ;;
    esac
done
[[ -n "$assembly" && -n "$genus" && -n "$out" ]] || { usage >&2; pad_die "--assembly, --genus and --out are required"; }
[[ -s "$assembly" ]] || pad_die "assembly not found: $assembly"
pad_load_conf

if [[ -s "$out" && -s "$out.meta.tsv" && "$force" != "1" ]]; then
    pad_log "$(pad_rel "$out") exists, skipping"
    exit 0
fi

pad_activate pcan
mkdir -p "$DATA_DIR/work" "$(dirname "$out")"
work="$(mktemp -d "$DATA_DIR/work/pcan.XXXXXX")"
# The working directory holds the decompressed assembly and PCAn's FIMO output.
# It goes whatever happens, and a half-written call table goes with it.
trap 'rm -rf "$work" "$out.part"' EXIT
trap 'exit 130' INT TERM

python "$PAD_ROOT/scripts/lib/pcan_runner.py" --assembly "$assembly" --genus "$genus" \
    --species "$species" --out "$out" --workdir "$work" --label "$label" --keep-fimo "$keep"
