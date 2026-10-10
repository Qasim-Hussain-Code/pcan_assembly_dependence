#!/usr/bin/env bash
# Check out PCAn at the pinned commit and apply patches/, build the PCAn and
# tools environments, record every version, export scrubbed lock files to env/,
# then clean the conda package cache and log the space it freed.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash scripts/01_install.sh [--force]

Needs only a conda base installation (CONDA_BASE in project.conf) and git.

  data/pcan/       PCAn at commit a5fa46f (the v1.0 release on Zenodo), with patches/*.patch applied
  data/env/pcan    PCAn's own environment, pinned by its pcan_specs.yml
  data/env/tools   minimap2, seqkit, seqtk, fastp, SPAdes, MEGAHIT, NCBI datasets,
                   shellcheck, pigz and the Python packages for the analysis
  env/             explicit lock files for both environments (tracked)
  logs/versions.tsv, logs/conda_clean.tsv

If env/*.explicit.txt exist, the environments are rebuilt from them exactly.
Otherwise they are solved from the pinned specifications below and the lock
files are written. --force rebuilds everything.
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
pad_skip_if_done install "$force"
pad_measure_self install "$@"

PCAN_URL="https://github.com/JHelsen/point-centromere-detection.git"
# The Zenodo record for PCAn v1.0 (10.5281/zenodo.17293587; all versions:
# 10.5281/zenodo.17293586) archives exactly this commit. The v1.0 tag on GitHub points one commit later
# (86184e5), and HEAD on 9 October 2026 (1894bcc) later still; both differ
# from a5fa46f only in README files.
PCAN_COMMIT="a5fa46f0cb971e7d499fd52b38b0c29166b33ef8"

# Tools environment. Versions are pinned to what bioconda and conda-forge served
# on 9 October 2026; the explicit lock files in env/ pin every dependency.
# xlrd reads the .xls supplementary tables of Peter et al. 2018 (arm 2); it was
# added on 9 October 2026, after the first install, with --freeze-installed so
# that nothing else in the environment changed.
TOOLS_SPEC=(
    "python=3.12" "minimap2=2.31" "seqkit=2.14.0" "seqtk=1.5" "fastp=1.4.0"
    "spades=4.3.0" "megahit=1.2.9" "ncbi-datasets-cli=18.38.0"
    "pandas>=2.2,<3" "numpy>=2" "scipy" "matplotlib" "openpyxl" "xlrd=2.0.1" "shellcheck" "pigz"
)

building=""
cleanup() {
    local status=$?
    if [[ $status -ne 0 && -n "$building" && -d "$building" ]]; then
        pad_log "install interrupted; removing the partial environment $(pad_rel "$building")"
        rm -rf "$building"
    fi
}
trap cleanup EXIT

# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
mkdir -p "$DATA_DIR/pcan" "$DATA_DIR/env" "$DATA_DIR/logs" "$PAD_ROOT/env"

# 1. PCAn at the pinned commit
checkout="$DATA_DIR/pcan/point-centromere-detection"
if [[ ! -d "$checkout/.git" ]]; then
    pad_log "cloning PCAn"
    git clone -q "$PCAN_URL" "$checkout"
fi
git -C "$checkout" fetch -q --tags origin
if [[ "$(git -C "$checkout" rev-parse HEAD)" != "$PCAN_COMMIT" ]]; then
    git -C "$checkout" checkout -q --force "$PCAN_COMMIT"
fi
[[ "$(git -C "$checkout" rev-parse HEAD)" == "$PCAN_COMMIT" ]] || pad_die "PCAn checkout is not at $PCAN_COMMIT"
pad_log "PCAn at $(git -C "$checkout" rev-parse --short HEAD), the commit archived as v1.0 on Zenodo"

for p in "$PAD_ROOT"/patches/*.patch; do
    [[ -e "$p" ]] || continue
    if git -C "$checkout" apply --reverse --check "$p" 2>/dev/null; then
        pad_log "patch $(basename "$p") already applied"
    else
        git -C "$checkout" apply --check "$p" || pad_die "patch $(basename "$p") does not apply to $PCAN_COMMIT"
        git -C "$checkout" apply "$p"
        pad_log "applied patch $(basename "$p")"
    fi
done

# 2. PCAn's environment
# The pins are read from pcan_specs.yml at the pinned commit and passed to
# conda create with only bioconda and conda-forge, in the order the file lists
# them. Two departures from running `conda env create -f pcan_specs.yml`:
#   - --override-channels, so that a machine whose conda is configured with
#     the defaults channel does not mix it in;
#   - flexible channel priority. Under strict priority, which recent conda
#     uses by default, the specification does not solve: biopython 1.83 and
#     tqdm 4.67.1 for Python 3.8 exist only on conda-forge, and the file lists
#     bioconda first. Flexible priority is the setting under which the file
#     solves.
spec="$checkout/PCAn/pcan_specs.yml"
mapfile -t pcan_conda < <(awk '/^dependencies:/{d=1;next} d && /^  - [A-Za-z]/ && !/pip:/ {sub(/^  - /,""); print}' "$spec")
mapfile -t pcan_pip < <(awk '/^  - pip:/{p=1;next} p && /^      - / {sub(/^      - /,""); print}' "$spec")
[[ ${#pcan_conda[@]} -ge 5 ]] || pad_die "could not read the conda pins from pcan_specs.yml"
if [[ ! -x "$PAD_ENV_PCAN/bin/python" || "$force" == "1" ]]; then
    rm -rf "$PAD_ENV_PCAN"
    building="$PAD_ENV_PCAN"
    if [[ -f "$PAD_ROOT/env/pcan.explicit.txt" ]]; then
        pad_log "building the PCAn environment from env/pcan.explicit.txt"
        conda create -q -y -p "$PAD_ENV_PCAN" --file "$PAD_ROOT/env/pcan.explicit.txt" > "$DATA_DIR/logs/install_pcan_env.log" 2>&1
    else
        pad_log "solving the PCAn environment: ${pcan_conda[*]}"
        CONDA_CHANNEL_PRIORITY=flexible conda create -q -y -p "$PAD_ENV_PCAN" --override-channels \
            -c bioconda -c conda-forge "${pcan_conda[@]}" > "$DATA_DIR/logs/install_pcan_env.log" 2>&1
    fi
    "$PAD_ENV_PCAN/bin/python" -m pip install -q --no-deps "${pcan_pip[@]}" >> "$DATA_DIR/logs/install_pcan_env.log" 2>&1
    building=""
fi

# 3. Tools environment
if [[ ! -x "$PAD_ENV_TOOLS/bin/python" || "$force" == "1" ]]; then
    rm -rf "$PAD_ENV_TOOLS"
    building="$PAD_ENV_TOOLS"
    if [[ -f "$PAD_ROOT/env/tools.explicit.txt" ]]; then
        pad_log "building the tools environment from env/tools.explicit.txt"
        conda create -q -y -p "$PAD_ENV_TOOLS" --file "$PAD_ROOT/env/tools.explicit.txt" > "$DATA_DIR/logs/install_tools_env.log" 2>&1
    else
        pad_log "solving the tools environment"
        conda create -q -y -p "$PAD_ENV_TOOLS" --override-channels --strict-channel-priority \
            -c conda-forge -c bioconda "${TOOLS_SPEC[@]}" > "$DATA_DIR/logs/install_tools_env.log" 2>&1
    fi
    building=""
fi

# 4. Versions, as the tools themselves report them
P="$PAD_ENV_PCAN/bin"
T="$PAD_ENV_TOOLS/bin"
pyver() { "$1" -c "import importlib.metadata as m; print(m.version('$2'))"; }
{
    printf 'tool\tversion\tenvironment\n'
    printf 'PCAn commit\t%s\tdata/pcan\n' "$PCAN_COMMIT"
    printf 'PCAn patches\t%s\tdata/pcan\n' "$(find "$PAD_ROOT/patches" -maxdepth 1 -name '*.patch' -printf '%f\n' | sort | paste -sd, -)"
    printf 'fimo (MEME suite)\t%s\tpcan\n' "$("$P/fimo" --version 2>&1 | head -1)"
    printf 'blastp\t%s\tpcan\n' "$("$P/blastp" -version 2>&1 | head -1 | awk '{print $2}')"
    printf 'python\t%s\tpcan\n' "$("$P/python" -c 'import platform; print(platform.python_version())')"
    for m in biopython pandas matplotlib tqdm orffinder numpy; do
        printf '%s\t%s\tpcan\n' "$m" "$(pyver "$P/python" "$m")"
    done
    printf 'python\t%s\ttools\n' "$("$T/python" -c 'import platform; print(platform.python_version())')"
    for m in pandas numpy scipy matplotlib openpyxl xlrd; do
        printf '%s\t%s\ttools\n' "$m" "$(pyver "$T/python" "$m")"
    done
    printf 'minimap2\t%s\ttools\n' "$("$T/minimap2" --version 2>&1)"
    printf 'seqkit\t%s\ttools\n' "$("$T/seqkit" version 2>&1 | awk '{print $2}')"
    printf 'seqtk\t%s\ttools\n' "$(conda list -p "$PAD_ENV_TOOLS" '^seqtk$' | awk '!/^#/ {print $2}')"
    printf 'fastp\t%s\ttools\n' "$("$T/fastp" --version 2>&1 | awk '{print $2}')"
    printf 'spades\t%s\ttools\n' "$("$T/spades.py" --version 2>&1 | awk '{print $NF}')"
    printf 'megahit\t%s\ttools\n' "$("$T/megahit" --version 2>&1 | awk '{print $NF}')"
    printf 'ncbi datasets\t%s\ttools\n' "$("$T/datasets" --version 2>&1 | awk '{print $NF}')"
    printf 'shellcheck\t%s\ttools\n' "$("$T/shellcheck" --version | awk '/^version:/ {print $2}')"
    printf 'pigz\t%s\ttools\n' "$("$T/pigz" --version 2>&1 | awk '{print $2}')"
    printf 'conda\t%s\tbase\n' "$(conda --version | awk '{print $2}')"
} > "$PAD_ROOT/logs/versions.tsv.tmp"
mv "$PAD_ROOT/logs/versions.tsv.tmp" "$PAD_ROOT/logs/versions.tsv"
pad_log "versions written to logs/versions.tsv"

# 5. Lock files. conda list --explicit records package URLs and checksums only.
# pip-installed packages are listed separately, by name and version, because
# pip freeze would record build paths for packages that conda installed.
conda list -p "$PAD_ENV_PCAN" --explicit --md5 > "$PAD_ROOT/env/pcan.explicit.txt"
conda list -p "$PAD_ENV_PCAN" | awk '$4 == "pypi" {print $1 "==" $2}' > "$PAD_ROOT/env/pcan.pip.txt"
conda list -p "$PAD_ENV_TOOLS" --explicit --md5 > "$PAD_ROOT/env/tools.explicit.txt"
if grep -nE 'file://|/home/|/mnt/|/Users/|^prefix:' "$PAD_ROOT"/env/*.txt; then
    pad_die "a lock file in env/ contains a local path"
fi
pad_log "lock files written to env/"

# 6. Package cache. Removes tarballs and unpacked packages that no environment
# links to; environments, including ones from other projects, are untouched.
before="$(du -sb "$CONDA_BASE/pkgs" 2>/dev/null | awk '{print $1}')"
conda clean --all --yes > "$DATA_DIR/logs/conda_clean.log" 2>&1
after="$(du -sb "$CONDA_BASE/pkgs" 2>/dev/null | awk '{print $1}')"
{
    printf 'date_utc\tpackage_cache_before_gb\tpackage_cache_after_gb\tfreed_gb\n'
    awk -v d="$(pad_utc)" -v b="$before" -v a="$after" 'BEGIN {printf "%s\t%.2f\t%.2f\t%.2f\n", d, b/1e9, a/1e9, (b-a)/1e9}'
} > "$PAD_ROOT/logs/conda_clean.tsv"
pad_log "conda clean freed $(awk -v b="$before" -v a="$after" 'BEGIN {printf "%.2f", (b-a)/1e9}') GB"

pad_mark_done install
pad_log "install complete"
