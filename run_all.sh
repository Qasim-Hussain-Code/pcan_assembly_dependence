#!/usr/bin/env bash
# Run the whole analysis, or part of it, in order. Every stage skips work it has
# already done and says so, so re-running this is safe.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash run_all.sh [--arm LIST] [--from STAGE] [--jobs N] [--test] [--help]

  --arm LIST    arms to run, comma separated from 0,1,2,3 (default: 0,1,2,3).
                The shared stages (configure to assemblies) always run first.
  --from STAGE  start at this stage, skipping the ones before it
  --jobs N      PCAn runs in parallel in arms 0 and 1 (default: JOBS in
                project.conf, which defaults to 1 so that timings stay
                comparable between machines; memory is not the reason)
  --test        run the known-answer tests in tests/ and stop
  --help        show this help

Stages, in order:
  configure    scripts/00_configure.sh (only when project.conf is missing)
  install      scripts/01_install.sh
  tables       scripts/02_fetch_tables.sh
  assemblies   scripts/03_fetch_assemblies.sh --set arm0
  reproduce    scripts/05_reproduce.py                arm 0
  fragment     scripts/06_fragment.py                 arm 1 (simulation)
  variants     scripts/07_plant_variants.py           arm 1 (simulation)
  score        scripts/08_score_perturbations.py      arm 1
  pairs        scripts/09_build_pairs.py              arms 2 and 3
  assemblies2  scripts/03_fetch_assemblies.sh --set arm2
  compare      scripts/10_compare_pairs.py --arm 2    arm 2
  reads        scripts/11_assemble_reads.sh           arm 3
  compare3     scripts/10_compare_pairs.py --arm 3    arm 3
  analyse      scripts/12_analyse.py
  figures      scripts/13_figures.py

Arm 0 is a gate. If fewer than 95 per cent of published calls are reproduced
exactly (config/arm0_reproduction.md), arms 1 to 3 do not start.
EOF
}

# shellcheck source=scripts/lib/common.sh
source "$(dirname "${BASH_SOURCE[0]}")/scripts/lib/common.sh"
arms="0,1,2,3" from="" jobs="" test=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --arm)  arms="$2"; shift 2 ;;
        --from) from="$2"; shift 2 ;;
        --jobs) jobs="$2"; shift 2 ;;
        --test) test=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; pad_die "unknown option $1" ;;
    esac
done

STAGES=(configure install tables assemblies reproduce fragment variants score pairs assemblies2 compare reads compare3 analyse figures)
declare -A ARM_OF=([reproduce]="0" [fragment]="1" [variants]="1" [score]="1" [pairs]="2 3" [assemblies2]="2"
                   [compare]="2" [reads]="3" [compare3]="3")

if [[ ! -f "$PAD_CONF" ]]; then
    pad_log "no project.conf; configuring with detected defaults"
    bash "$PAD_ROOT/scripts/00_configure.sh" --yes
fi
pad_load_conf
[[ -n "$jobs" ]] || jobs="$JOBS"
export PAD_JOBS="$jobs"

if [[ "$test" == "1" ]]; then
    pad_activate tools
    cd "$PAD_ROOT"
    python -m unittest discover -s tests -p 'test_*.py' -v
    exit $?
fi

wanted() {
    local stage="$1" a
    [[ -z "${ARM_OF[$stage]:-}" ]] && return 0
    read -ra stage_arms <<< "${ARM_OF[$stage]}"
    for a in "${stage_arms[@]}"; do
        [[ ",$arms," == *",$a,"* ]] && return 0
    done
    return 1
}

python_stage() {
    local name="$1"
    shift
    pad_activate tools
    python3 "$PAD_ROOT/scripts/lib/measure.py" --stage "$name" -- python "$@"
}

started=0
[[ -z "$from" ]] && started=1
for stage in "${STAGES[@]}"; do
    if [[ "$started" == "0" ]]; then
        if [[ "$stage" == "$from" ]]; then started=1; else continue; fi
    fi
    wanted "$stage" || continue
    pad_log "stage $stage"
    case "$stage" in
        configure)   : ;;
        install)     bash "$PAD_ROOT/scripts/01_install.sh" ;;
        tables)      bash "$PAD_ROOT/scripts/02_fetch_tables.sh" ;;
        assemblies)  bash "$PAD_ROOT/scripts/03_fetch_assemblies.sh" --set arm0 ;;
        reproduce)   python_stage reproduce "$PAD_ROOT/scripts/05_reproduce.py" --jobs "$jobs" ;;
        fragment|variants|score|pairs|compare|compare3|analyse|figures)
            if [[ ! -f "$PAD_ROOT/results/arm0/gate.tsv" ]] \
               || [[ "$(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next} {print $h["decision"]}' "$PAD_ROOT/results/arm0/gate.tsv")" != "proceed" ]]; then
                pad_die "arm 0 has not passed its gate (results/arm0/gate.tsv); later arms do not start"
            fi
            case "$stage" in
                fragment) python_stage fragment "$PAD_ROOT/scripts/06_fragment.py" --jobs "$jobs" ;;
                variants) python_stage variants "$PAD_ROOT/scripts/07_plant_variants.py" --jobs "$jobs" ;;
                score)    python_stage score "$PAD_ROOT/scripts/08_score_perturbations.py" ;;
                pairs)    python_stage pairs "$PAD_ROOT/scripts/09_build_pairs.py" ;;
                compare)  python_stage compare "$PAD_ROOT/scripts/10_compare_pairs.py" --arm 2 --jobs "$jobs" ;;
                compare3) python_stage compare3 "$PAD_ROOT/scripts/10_compare_pairs.py" --arm 3 --jobs "$jobs" ;;
                analyse)  python_stage analyse "$PAD_ROOT/scripts/12_analyse.py" ;;
                figures)  python_stage figures "$PAD_ROOT/scripts/13_figures.py" ;;
            esac ;;
        assemblies2) bash "$PAD_ROOT/scripts/03_fetch_assemblies.sh" --set arm2 ;;
        reads)       bash "$PAD_ROOT/scripts/11_assemble_reads.sh" ;;
    esac
done
pad_log "run_all finished"
