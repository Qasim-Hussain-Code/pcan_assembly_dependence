# shellcheck shell=bash
# Sourced by every bash stage. Finds the repository root from this file's own
# position, loads project.conf, and provides logging, path scrubbing, conda
# activation and the skip check that makes each stage safe to re-run.

PAD_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PAD_CONF="$PAD_ROOT/project.conf"
export PAD_ROOT

pad_log() {
    printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2
}

pad_die() {
    pad_log "error: $*"
    exit 1
}

# project.conf is machine-specific and untracked. Every stage after
# 00_configure.sh refuses to run without it, so that thread counts, memory caps
# and disk budgets always come from one place.
pad_load_conf() {
    [[ -f "$PAD_CONF" ]] || pad_die "project.conf not found; run: bash scripts/00_configure.sh --yes"
    # shellcheck source=/dev/null
    source "$PAD_CONF"
    : "${DATA_DIR:=$PAD_ROOT/data}"
    : "${THREADS:=1}"
    : "${JOBS:=1}"
    : "${CONDA_BASE:?project.conf has no CONDA_BASE}"
    export DATA_DIR THREADS JOBS CONDA_BASE
    export PAD_ENV_PCAN="$DATA_DIR/env/pcan"
    export PAD_ENV_TOOLS="$DATA_DIR/env/tools"
    export PAD_PCAN_DIR="$DATA_DIR/pcan/point-centromere-detection/PCAn"
    mkdir -p "$DATA_DIR/state" "$PAD_ROOT/logs"
}

# Activate one of the two environments built by 01_install.sh. PCAn keeps its
# own environment, exactly as its pcan_specs.yml pins it (Python 3.8.17, MEME
# 4.11.2); everything else runs in the tools environment.
pad_activate() {
    local which="$1"
    # shellcheck source=/dev/null
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    set +u
    case "$which" in
        pcan)  conda activate "$PAD_ENV_PCAN" ;;
        tools) conda activate "$PAD_ENV_TOOLS" ;;
        *)     set -u; pad_die "unknown environment $which" ;;
    esac
    set -u
}

# Paths written to logs/ and results/ must be relative to the repository root.
# The data directory is written as data/ whatever its real location.
pad_rel() {
    local p="$1"
    p="${p/#$PAD_ROOT\//}"
    p="${p/#$DATA_DIR/data}"
    printf '%s' "$p"
}

# Remove machine-identifying strings from tool output before it is tracked.
pad_scrub() {
    local user host
    user="$(id -un)"
    host="$(hostname)"
    sed -e "s#${PAD_ROOT}#.#g" \
        -e "s#${DATA_DIR}#data#g" \
        -e "s#${CONDA_BASE}#<conda>#g" \
        -e "s#${HOME}#<home>#g" \
        -e "s#/mnt/[a-z]/Users/[^/[:space:]]*#<home>#g" \
        -e "s#/home/[^/[:space:]]*#<home>#g" \
        -e "s#${host}#<host>#g" \
        -e "s#\\b${user}\\b#<user>#g"
}

# A stage writes data/state/<name>.done when it finishes. Re-running it then
# skips and says so, unless --force is given.
pad_is_done() {
    [[ -f "$DATA_DIR/state/$1.done" ]]
}

pad_mark_done() {
    date -u +%Y-%m-%dT%H:%M:%SZ > "$DATA_DIR/state/$1.done"
}

pad_skip_if_done() {
    local name="$1" force="${2:-0}"
    if [[ "$force" != "1" ]] && pad_is_done "$name"; then
        pad_log "$name: already complete on $(cat "$DATA_DIR/state/$name.done"), skipping (use --force to redo)"
        exit 0
    fi
}

pad_utc() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

# Re-run the calling stage under scripts/lib/measure.py unless it is already
# being measured, so that every stage appends its elapsed time, peak memory and
# peak disk to logs/resources.tsv whether run_all.sh started it or a person did.
pad_measure_self() {
    local stage="$1"
    shift
    if [[ -z "${PAD_MEASURED:-}" ]]; then
        local py
        py="$(command -v python3 || true)"
        [[ -n "$py" ]] || py="$CONDA_BASE/bin/python"
        exec "$py" "$PAD_ROOT/scripts/lib/measure.py" --stage "$stage" -- bash "$0" "$@"
    fi
}
