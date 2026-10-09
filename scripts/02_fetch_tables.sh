#!/usr/bin/env bash
# Fetch the published tables this analysis compares against or draws from, each
# with a checksum and download date: Supplementary Data 2, 4, 5 and 6 of the
# PCAn paper from the publisher, the matching Figshare files, and the dated SGD
# S288C genome release. Then convert them to TSV and build config/species_arm0.tsv.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash scripts/02_fetch_tables.sh [--force]

Writes
  data/tables/                raw downloads and their TSV conversions (untracked)
  logs/downloads.tsv          URL, bytes, SHA-256, source checksum and UTC date per file
  config/species_arm0.tsv     every Supplementary Data 5 row with its PCAn genus
                              selection and arm 0 role, with reasons
  results/arm0/source_concordance.tsv
                              publisher versus Figshare copies of the same tables
EOF
}

# shellcheck source=lib/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
force=0
for a in "$@"; do
    case "$a" in
        --force) force=1 ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; pad_die "unknown option $a" ;;
    esac
done
pad_load_conf
pad_skip_if_done fetch_tables "$force"
pad_measure_self fetch_tables "$@"
pad_activate tools

tables="$DATA_DIR/tables"
mkdir -p "$tables"
manifest="$PAD_ROOT/logs/downloads.tsv"
partial=""
trap '[[ -n "$partial" ]] && rm -f "$partial"' EXIT

ESM="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09779-1/MediaObjects"
# SGD's archive host resets HTTPS connections; the same bucket answers over
# HTTPS at its path-style S3 address.
SGD="https://s3-us-west-2.amazonaws.com/sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases"

# name|url|checksum published by the source (figshare md5, or none)
SOURCES=(
    "supp_data2_publisher.xlsx|$ESM/41586_2025_9779_MOESM3_ESM.xlsx|"
    "supp_data4_publisher.xlsx|$ESM/41586_2025_9779_MOESM5_ESM.xlsx|"
    "supp_data5_publisher.xlsx|$ESM/41586_2025_9779_MOESM6_ESM.xlsx|"
    "supp_data6_publisher.xlsx|$ESM/41586_2025_9779_MOESM7_ESM.xlsx|"
    "supp_data2_figshare_28225061v2.xlsx|https://ndownloader.figshare.com/files/51725942|md5:3c2d50afc29fe45c317d4fc19452a107"
    "supp_data5_figshare_28225067v3.xlsx|https://ndownloader.figshare.com/files/56861582|md5:d4dc4fc776d69ed9ae53857fa6563dd6"
    "supp_data6_figshare_28225076v1.xlsx|https://ndownloader.figshare.com/files/51725960|md5:738bfa0a1eb329121a3753803dacb39f"
    "figshare_collection_7630151.json|https://api.figshare.com/v2/collections/7630151/articles?page_size=100|"
    "S288C_reference_genome_R64-5-1_20240529.tgz|$SGD/S288C_reference_genome_R64-5-1_20240529.tgz|"
)

[[ -f "$manifest" ]] || printf 'file\turl\tbytes\tsha256\tsource_checksum\tsource_checksum_ok\tlast_modified\tdownloaded_utc\n' > "$manifest"

for entry in "${SOURCES[@]}"; do
    IFS='|' read -r name url expect <<< "$entry"
    dest="$tables/$name"
    if [[ -s "$dest" && "$force" != "1" ]]; then
        pad_log "have $name"
        continue
    fi
    partial="$dest.part"
    # Figshare refuses HEAD requests, so a missing Last-Modified is expected there.
    headers="$(curl -fsSIL --retry 2 "$url" 2>/dev/null | tr -d '\r' | awk 'tolower($1)=="last-modified:" {sub(/^[^:]*: /,""); lm=$0} END {print lm}')" || headers=""
    curl -fsSL --retry 5 --retry-delay 10 -o "$partial" "$url"
    mv "$partial" "$dest"
    partial=""
    sha="$(sha256sum "$dest" | awk '{print $1}')"
    ok="not published"
    if [[ "$expect" == md5:* ]]; then
        got="$(md5sum "$dest" | awk '{print $1}')"
        if [[ "$got" == "${expect#md5:}" ]]; then ok="yes"; else ok="no"; fi
        [[ "$ok" == "yes" ]] || pad_die "$name: md5 $got does not match the published $expect"
    fi
    # The Figshare listing changes as the collection is updated, so it is kept
    # as a dated snapshot and has no fixed checksum to compare against.
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "data/tables/$name" "$url" "$(stat -c %s "$dest")" "$sha" \
        "${expect:-none}" "$ok" "${headers:-not given}" "$(pad_utc)" >> "$manifest"
    pad_log "fetched $name ($(stat -c %s "$dest") bytes)"
done

# The SGD release is one archive. It is stored once, compressed. Two members are
# extracted: the GFF, which holds the centromere features, and the reference
# sequence, which arm 0 compares with NCBI's copy so that SGD's coordinates
# can be applied to the assembly PCAn is run on. The archive was made on a Mac
# and carries AppleDouble companions (._name) beside each file; they are skipped.
sgd_tgz="$tables/S288C_reference_genome_R64-5-1_20240529.tgz"
extract_member() {
    local pattern="$1" dest="$2" member
    member="$(tar --warning=no-unknown-keyword -tzf "$sgd_tgz" | grep -E "$pattern" | grep -v '/\._' | head -1)"
    [[ -n "$member" ]] || pad_die "no member matching $pattern in the SGD release archive"
    tar --warning=no-unknown-keyword -xzOf "$sgd_tgz" "$member" | gzip -dc > "$dest.part"
    mv "$dest.part" "$dest"
    pad_log "extracted $(basename "$member")"
}
if [[ ! -s "$tables/sgd_R64-5-1.gff" || "$force" == "1" ]]; then
    extract_member '\.gff\.gz$' "$tables/sgd_R64-5-1.gff"
fi
if [[ ! -s "$tables/sgd_R64-5-1_reference.fsa" || "$force" == "1" ]]; then
    extract_member 'S288C_reference_sequence_.*\.fsa\.gz$' "$tables/sgd_R64-5-1_reference.fsa"
fi

python "$PAD_ROOT/scripts/lib/tables.py"

pad_mark_done fetch_tables
pad_log "tables complete"
