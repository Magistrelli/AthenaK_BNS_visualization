#!/usr/bin/env bash
# =============================================================================
# run_comparison.sh  —  overlay neutrino and ejecta diagnostics from 2+ models
#
# Wraps comparison_plot.py (neutrinos, and the ejecta figures as fig_ejecta.pdf
# and fig_ejecta_time.pdf of a single run).  Reuses the .txt outputs already
# produced by run_analysis.sh in each model's analysis directory (no data
# regeneration).
#
# Usage:
#   bash run_comparison.sh                 # models from config.ini (or MODELS)
#   bash run_comparison.sh DIR1 LABEL1 DIR2 LABEL2 [...]   # override on CLI
#
# By default the models are the simulations listed in config.ini (`simpath`
# lines), each with its <simpath>/analysis directory; those without one yet
# are left out.  Or fill the MODELS list below (one "analysis_dir|label" each).
# First model is drawn solid, second dashed, third dotted, fourth dash-dot.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS  ← edit per run
# ─────────────────────────────────────────────────────────────────────────────

# Models to overlay: one "analysis_dir|legend label" per line (2 to 4 models).
# Each analysis_dir must contain the merger_time.txt / Lnu_E_*.txt / Eav_*.txt
# produced by run_analysis.sh.  Empty = the simulations listed in config.ini, e.g.
#   "/lus/flare/projects/CompactBinaryMerger/yiqiu/Large_sim/bundle_run/VVLR_eq/analysis|baseline"
MODELS=()

# Alternate Perlmutter comparison (or supply these DIR/LABEL pairs on the CLI):
# "/pscratch/sd/y/yiqiu/TOV_test_athenaK/lorene_bns/DD2_m1_nn_opacity/VVLR/analysis|Baseline"
# "/pscratch/sd/y/yiqiu/TOV_test_athenaK/lorene_bns/DD2_m1_nn_opacity/VVLR_nn_new/analysis|NN"

# Output figure filename.
OUT="fig_neutrino_comparison.png"

# Where to write the figure (empty = <first analysis dir>/comparison_plot/).
OUTDIR="${SCRATCH}/athenak/comparisons"

# Override x-range [ms] about merger, e.g. XLIM=(-2 20).  Empty = all models.
XLIM=()

# Output figure filenames of the ejecta comparisons (as fig_ejecta.pdf and
# fig_ejecta_time.pdf).  The distributions cover each model's post-merger
# window (T_POST_MS in run_analysis.sh), which must be the same for all models.
OUT_EJ="fig_ejecta_comparison.png"
OUT_EJ_TIME="fig_ejecta_time_comparison.png"

# ─────────────────────────────────────────────────────────────────────────────
# Allow overriding the MODELS list from the command line (DIR LABEL pairs).
# ─────────────────────────────────────────────────────────────────────────────
if [[ $# -gt 0 ]]; then
    if [[ $(( $# % 2 )) -ne 0 ]]; then
        echo "ERROR: command-line args must be DIR LABEL pairs (got $# args)."
        echo "Usage: $0 [DIR1 LABEL1 DIR2 LABEL2 ...]"
        exit 1
    fi
    MODELS=()
    while [[ $# -gt 0 ]]; do
        MODELS+=("$1|$2")
        shift 2
    done
fi

# ─────────────────────────────────────────────────────────────────────────────
# SETUP  — machine environment from config.ini (env_setup: module/conda/venv)
# ─────────────────────────────────────────────────────────────────────────────
source "${SCRIPT_DIR}/../load_config.sh"

if [[ ${#MODELS[@]} -eq 0 ]]; then
    while IFS= read -r sim; do
        if [[ -d "${sim%%|*}/analysis" ]]; then
            MODELS+=("${sim%%|*}/analysis|${sim#*|}")
        else
            echo "  [skip] ${sim#*|}: no analysis directory yet (${sim%%|*}/analysis)"
        fi
    done < <(_cfg_sims)
fi
if [[ ${#MODELS[@]} -lt 2 ]]; then
    echo "ERROR: need at least two models with an analysis directory (set MODELS," \
         "pass them on the command line, or list simulations in ${CONFIG_FILE})."
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Build comparison_plot.py args from MODELS
# ─────────────────────────────────────────────────────────────────────────────
ARGS=()
for entry in "${MODELS[@]}"; do
    dir="${entry%%|*}"
    label="${entry#*|}"
    if [[ ! -d "${dir}" ]]; then
        echo "ERROR: analysis directory not found: ${dir}"; exit 1
    fi
    ARGS+=(--dir "${dir}" --label "${label}")
done

ARGS+=(--out "${OUT}" --ejecta-out "${OUT_EJ}" --ejecta-time-out "${OUT_EJ_TIME}")
[[ -n "${OUTDIR}" ]] && ARGS+=(--outdir "${OUTDIR}")
[[ ${#XLIM[@]} -eq 2 ]] && ARGS+=(--xlim-ms "${XLIM[0]}" "${XLIM[1]}")

echo "============================================================"
echo "  BNS neutrino + ejecta comparison plots"
echo "  Models : ${#MODELS[@]}"
for entry in "${MODELS[@]}"; do echo "    - ${entry%%|*}  (${entry#*|})"; done
echo "============================================================"

python3 "${SCRIPT_DIR}/comparison_plot.py" "${ARGS[@]}"
