#!/usr/bin/env bash
# Detect or accept the machine's resources and budgets, project the footprint
# of a full run, and write project.conf, which every later stage sources.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash scripts/00_configure.sh [options] --yes

Writes project.conf (untracked) at the repository root.

Options
  --threads N        threads given to FIMO, minimap2, fastp and the assemblers (default: nproc)
  --ram GB           RAM available to this run (default: MemTotal of this system)
  --disk GB          working disk budget for the data directory (default: 50)
  --hours-arm1 H     time budget for arm 1 (default: 24)
  --hours-arm3 H     time budget for arm 3 (default: 48)
  --spades-ram GB    SPAdes memory cap (default: 75 per cent of --ram, rounded down)
  --jobs N           PCAn runs in parallel in arms 0 and 1 (default: 1, see below)
  --data-dir PATH    where assemblies, reads, environments and the PCAn checkout live (default: data/)
  --conda-base PATH  conda installation to build the environments from (default: detected)
  --yes              write project.conf, replacing an existing one
  -h, --help         show this help

The hard disk ceiling is set 10 per cent above the working budget. Stages that
measure a projected footprint above the budget, or above the free space of the
drive, refuse to start and print the shortfall.

--jobs defaults to 1 so that the per-run timings in logs/ stay comparable
between runs on different machines. Memory is not the reason: one PCAn run on
a 12 Mb yeast genome peaks well under 1 GB.
EOF
}

# shellcheck source=lib/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

threads="" ram="" disk="50" h1="24" h3="48" spades="" jobs="1" data_dir="" conda_base="" yes=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)     threads="$2"; shift 2 ;;
        --ram)         ram="$2"; shift 2 ;;
        --disk)        disk="$2"; shift 2 ;;
        --hours-arm1)  h1="$2"; shift 2 ;;
        --hours-arm3)  h3="$2"; shift 2 ;;
        --spades-ram)  spades="$2"; shift 2 ;;
        --jobs)        jobs="$2"; shift 2 ;;
        --data-dir)    data_dir="$2"; shift 2 ;;
        --conda-base)  conda_base="$2"; shift 2 ;;
        --yes)         yes=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *)             usage >&2; pad_die "unknown option $1" ;;
    esac
done

if [[ -f "$PAD_CONF" && "$yes" != "1" ]]; then
    pad_log "project.conf already exists; pass --yes to replace it"
    exit 0
fi
[[ "$yes" == "1" ]] || { usage >&2; pad_die "pass --yes to write project.conf"; }

[[ -n "$threads" ]] || threads="$(nproc)"
if [[ -z "$ram" ]]; then
    # MemTotal is what this kernel can use. Under WSL 2 that is the VM limit
    # set in .wslconfig, not the physical RAM of the host.
    ram="$(awk '/^MemTotal:/ {printf "%.1f", $2 / 1048576}' /proc/meminfo)"
fi
[[ -n "$spades" ]] || spades="$(awk -v r="$ram" 'BEGIN {printf "%d", r * 0.75}')"
[[ -n "$data_dir" ]] || data_dir="$PAD_ROOT/data"
mkdir -p "$data_dir"
data_dir="$(cd "$data_dir" && pwd)"

if [[ -z "$conda_base" ]]; then
    if command -v conda >/dev/null 2>&1; then
        conda_base="$(conda info --base 2>/dev/null)"
    else
        for c in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/mambaforge" "$HOME/anaconda3" /opt/conda; do
            if [[ -f "$c/etc/profile.d/conda.sh" ]]; then conda_base="$c"; break; fi
        done
    fi
fi
[[ -n "$conda_base" && -f "$conda_base/etc/profile.d/conda.sh" ]] \
    || pad_die "no conda installation found; install Miniforge or pass --conda-base"

# Free space. Under WSL 2 the Linux filesystem is a virtual disk file on the
# Windows drive. df on that filesystem reports the virtual disk's nominal size
# (about 1 TB), which is meaningless here: the file grows on C: as data are
# written and does not shrink when files are deleted. The real limit is the
# free space of the Windows drive, so on WSL the smaller of the two is used.
free_gb="$(df -B1 --output=avail "$data_dir" | tail -1 | awk '{printf "%.1f", $1 / 1e9}')"
host_note=""
if grep -qi microsoft /proc/version 2>/dev/null && [[ -d /mnt/c ]]; then
    host_gb="$(df -B1 --output=avail /mnt/c | tail -1 | awk '{printf "%.1f", $1 / 1e9}')"
    if awk -v a="$host_gb" -v b="$free_gb" 'BEGIN {exit !(a < b)}'; then
        free_gb="$host_gb"
        host_note=" (free space of the Windows drive that holds the WSL virtual disk)"
    fi
    case "$data_dir" in
        /mnt/*) pad_log "warning: the data directory is on a Windows drive mounted in WSL. File creation there was 300 times slower than on the WSL disk in my tests, and conda environments do not install reliably on it. Consider --data-dir on the Linux filesystem." ;;
    esac
fi
ceiling="$(awk -v d="$disk" 'BEGIN {printf "%.1f", d * 1.1}')"

# Footprint projected before anything has run. These figures are priors; the
# stages replace them with measurements after the first ten assemblies and
# after the first arm 3 strain, and refuse to continue if the measured
# projection does not fit.
#   environments          4.0 GB  (PCAn environment plus the tools environment)
#   arm 0 assemblies      1.0 GB  (138 assemblies, gzip, with sequence reports)
#   arm 1 working space   0.1 GB  per parallel job (one perturbed FASTA at a time)
#   arm 3, one strain    12.0 GB  (subsampled reads, fastp output, SPAdes working directory)
#   arm 3, kept outputs   0.2 GB  per strain
proj="$(awk -v j="$jobs" 'BEGIN {printf "%.1f", 4.0 + 1.0 + 0.1 * j + 12.0 + 0.2 * 10}')"

if awk -v c="$ceiling" -v f="$free_gb" 'BEGIN {exit !(c > f - 5)}'; then
    pad_die "the disk ceiling of ${ceiling} GB leaves less than 5 GB free on a drive with ${free_gb} GB free${host_note}; pass a smaller --disk"
fi
if awk -v p="$proj" -v d="$disk" 'BEGIN {exit !(p > d)}'; then
    short="$(awk -v p="$proj" -v d="$disk" 'BEGIN {printf "%.1f", p - d}')"
    pad_die "projected peak footprint ${proj} GB exceeds the ${disk} GB budget by ${short} GB"
fi

cat > "$PAD_CONF" <<EOF
# Written by scripts/00_configure.sh on $(pad_utc). Machine-specific; not tracked.
THREADS=$threads
JOBS=$jobs
RAM_GB=$ram
SPADES_MEM_GB=$spades
DISK_BUDGET_GB=$disk
DISK_CEILING_GB=$ceiling
FREE_DISK_GB_AT_CONFIGURE=$free_gb
HOURS_ARM1=$h1
HOURS_ARM3=$h3
DATA_DIR=$data_dir
CONDA_BASE=$conda_base
EOF

pad_log "wrote project.conf: threads $threads, jobs $jobs, RAM ${ram} GB, SPAdes cap ${spades} GB"
pad_log "disk: budget ${disk} GB, ceiling ${ceiling} GB, free ${free_gb} GB${host_note}"
pad_log "prior projected peak footprint ${proj} GB; replaced by measurements once stages run"
pad_log "time budgets: arm 1 ${h1} h, arm 3 ${h3} h"
