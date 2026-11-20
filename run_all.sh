#!/usr/bin/env bash
# =============================================================================
# run_all.sh  —  all standard steps for every simulation listed in config.ini
#
# Usage:
#   bash run_all.sh                      # every step
#   bash run_all.sh --only analysis,run_comparison
#   bash run_all.sh --skip bns_all
#
# Steps, in this order (each script keeps its own settings, see its header):
#   bns_all         BNS_all.sh for each simulation: history plots, frames, movies
#   analysis        sph_data_analysis/run_analysis.sh for each simulation
#   bns_compare     BNS_compare.sh once for all simulations (needs only .hst files)
#   run_comparison  sph_data_analysis/run_comparison.sh once for all simulations
#                   with an analysis directory (written by the analysis step)
# The simulations are the `simpath` lines of config.ini; the two comparisons
# need at least two of them.  A failing simulation does not stop the others;
# the summary at the end lists what ran, was skipped or failed.
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALL_STEPS=(bns_all analysis bns_compare run_comparison)

usage() { echo "Usage: $0 [--only STEP,...] [--skip STEP,...]  (steps: ${ALL_STEPS[*]})"; exit 1; }

STEPS=("${ALL_STEPS[@]}")
while [[ $# -gt 0 ]]; do
    case "$1" in
        --only) IFS=, read -ra STEPS <<< "${2:?}"; shift 2 ;;
        --skip) IFS=, read -ra skip <<< "${2:?}"; shift 2
                keep=()
                for s in "${STEPS[@]}"; do
                    [[ " ${skip[*]} " == *" ${s} "* ]] || keep+=("${s}")
                done
                STEPS=("${keep[@]}") ;;
        *) usage ;;
    esac
done
for s in "${STEPS[@]}"; do
    [[ " ${ALL_STEPS[*]} " == *" ${s} "* ]] || { echo "Unknown step: ${s}"; usage; }
done
wanted() { [[ " ${STEPS[*]} " == *" $1 "* ]]; }

# config.ini: the simulations (_cfg_sims) and the environment, set up once here.
source "${SCRIPT_DIR}/load_config.sh"
mapfile -t SIMS < <(_cfg_sims)
if [[ ${#SIMS[@]} -eq 0 ]]; then
    echo "ERROR: no simpath in ${CONFIG_FILE}."; exit 1
fi
if [[ ${#SIMS[@]} -gt 1 ]]; then
    _cfg_check_multi_sim || exit 1
fi

SUMMARY=()
FAILED=0
note() {   # STEP STATUS WHAT
    SUMMARY+=("$(printf '%-15s %-8s %s' "$1" "$2" "$3")")
    [[ "$2" == FAILED ]] && FAILED=1
    return 0
}
banner() { echo; echo "################ $*"; }

# One script, once per simulation.
per_sim() {   # STEP SCRIPT
    local sim
    for sim in "${SIMS[@]}"; do
        banner "$1: ${sim#*|}  (${sim%%|*})"
        if _BNS_ENV_READY=1 _SIM_CURRENT="${sim}" bash "$2"; then
            note "$1" ok "${sim#*|}"
        else
            note "$1" FAILED "${sim#*|}"
        fi
    done
}

# One script, once for all the runs given as DIR LABEL pairs.
together() {   # STEP SCRIPT DIR LABEL [DIR LABEL ...]
    local step="$1" script="$2"; shift 2
    local n=$(( $# / 2 ))
    if (( n < 2 )); then
        note "${step}" skipped "needs at least two simulations (have ${n})"
        return
    fi
    banner "${step}: ${n} simulations"
    if _BNS_ENV_READY=1 bash "${script}" "$@"; then
        note "${step}" ok "${n} simulations"
    else
        note "${step}" FAILED "${n} simulations"
    fi
}

wanted bns_all  && per_sim bns_all  "${SCRIPT_DIR}/BNS_all.sh"
wanted analysis && per_sim analysis "${SCRIPT_DIR}/sph_data_analysis/run_analysis.sh"

if wanted bns_compare; then
    runs=()
    for sim in "${SIMS[@]}"; do runs+=("${sim%%|*}" "${sim#*|}"); done
    together bns_compare "${SCRIPT_DIR}/BNS_compare.sh" "${runs[@]}"
fi

if wanted run_comparison; then
    models=()
    for sim in "${SIMS[@]}"; do
        if [[ -f "${sim%%|*}/analysis/merger_time.txt" ]]; then
            models+=("${sim%%|*}/analysis" "${sim#*|}")
        else
            note run_comparison "left out" "${sim#*|}: no analysis yet"
        fi
    done
    together run_comparison "${SCRIPT_DIR}/sph_data_analysis/run_comparison.sh" "${models[@]}"
fi

banner "summary"
printf '%s\n' "${SUMMARY[@]}"
exit "${FAILED}"
