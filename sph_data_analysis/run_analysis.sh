#!/usr/bin/env bash
# =============================================================================
# run_analysis.sh  —  BNS SPH ejecta + neutrino analysis wrapper
#
# Usage:
#   bash run_analysis.sh [--ejecta] [--neutrino] [--plot] [--accretion] [--all]
#
# With no flags all steps (merger time, ejecta, neutrino, plot, accretion) run.
#
# Machine-specific PATHS (simpath, eos_table, optional pythonpath_extra) are read
# from ../config.ini.  The run settings below stay in this script — edit per run.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
# Honor an exported CONFIG_FILE (same rule as load_config.sh) so a driver can
# point this script at a different sim without touching the repo's config.ini.
CONFIG_FILE="${CONFIG_FILE:-${REPO_DIR}/config.ini}"
PLOTTOOLS="${REPO_DIR}/plot-tools"

# All machine-specific settings (paths + environment) come from config.ini.
# load_config.sh reads it, defines _cfg_get, and activates the environment.
source "${REPO_DIR}/load_config.sh"
# Several simulations in config.ini: run this script once for each of them.
_for_each_sim "$@"

SIMPATH="$(_cfg_get simpath)"
EOS_TABLE="$(_cfg_get eos_table)"
PYTHONPATH_EXTRA="$(_cfg_get pythonpath_extra)"   # optional (e.g. a vtk install)

if [[ -z "${SIMPATH}" ]]; then
    echo "ERROR: 'simpath' is not set in ${CONFIG_FILE}."; exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# RUN SETTINGS  ← edit per run
# ─────────────────────────────────────────────────────────────────────────────

# Directory where all output txt/npy/pdf files will be written.
OUTPUT_DIR="${SIMPATH}/analysis"

# AthenaK job name prefixing the .hst files (used by the accretion step).
JOBNAME=bns

# SPH extraction radii [M_sun]. Ejecta (and the accretion bookkeeping, which
# subtracts the ejecta flux) use EJ_RADIUS = 200 M_sun ~ 295 km; neutrinos keep
# RADIUS = 300. On multi-radius sph output the nearest shell is used.
# Must have mhd/adm/z4c sph output at EJ_RADIUS and rad_m1_* at RADIUS.
RADIUS=${RADIUS:-300.0}
EJ_RADIUS=${EJ_RADIUS:-200.0}

# Ejecta criterion the accretion step uses ('geo' or 'Ber').
ACCR_CRIT=Ber

# GW extraction radius for the merger-time estimate (rpsi4_*_RRRR.txt).
WF_RADIUS=300

# Neutrino plot x-axis: 0 = absolute time [ms]; 1 = t − t_merger [ms] (xlim left=0).
PLOT_FROM_MERGER=0

# How many ms post-merger to include in the 2D ejecta histogram.
T_POST_MS=25.0

# Density floor used in the simulation [code units].
DFLOOR=3e-15

# Number of parallel workers for the snapshot loop.
# Login node: keep ≤ 4.  Compute node (salloc/srun): set to the core count.
N_WORKERS=${N_WORKERS:-8}

# Process every STRIDE-th SPH snapshot (1 = all).  Long runs (many output-*
# segments) can have thousands of snapshots — the neutrino step reads them ALL,
# which stalls/gets killed on a login node.  Raise this (e.g. 5-10) to keep the
# ejecta + neutrino loops tractable; luminosity/rate curves stay smooth.
STRIDE=${STRIDE:-2}

# Reuse per-snapshot results cached under ${OUTPUT_DIR}/cache/ by earlier runs
# and process only snapshots never seen before (the SPH analogue of BNS_all.sh's
# --skip-existing). The reduction into the .txt/.npy outputs always re-runs over
# cached + new results together, so the outputs are the same either way -- only
# the runtime changes. Set 0 to recompute every snapshot from scratch.
SPH_CACHE=${SPH_CACHE:-1}
if [[ "${SPH_CACHE}" == "0" ]]; then _USE_CACHE=False; else _USE_CACHE=True; fi

# ── auto-discover all output-*/sph directories (sorted) ──────────────────────
SPH_DIRS=()
for d in "${SIMPATH}"/output-*/sph; do
    [[ -d "$d" ]] && SPH_DIRS+=("$d")
done
if [[ ${#SPH_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: no output-*/sph directories found under ${SIMPATH}"; exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# ARG PARSING
# ─────────────────────────────────────────────────────────────────────────────
RUN_EJECTA=0; RUN_NEUTRINO=0; RUN_PLOT=0; RUN_ACCRETION=0
if [[ $# -eq 0 ]]; then
    RUN_EJECTA=1; RUN_NEUTRINO=1; RUN_PLOT=1; RUN_ACCRETION=1
fi
for arg in "$@"; do
    case "$arg" in
        --ejecta)    RUN_EJECTA=1    ;;
        --neutrino)  RUN_NEUTRINO=1  ;;
        --plot)      RUN_PLOT=1      ;;
        --accretion) RUN_ACCRETION=1 ;;
        --all)       RUN_EJECTA=1; RUN_NEUTRINO=1; RUN_PLOT=1; RUN_ACCRETION=1 ;;
        *) echo "Unknown flag: $arg"; echo "Usage: $0 [--ejecta] [--neutrino] [--plot] [--accretion] [--all]"; exit 1 ;;
    esac
done

# ─────────────────────────────────────────────────────────────────────────────
# SETUP  — environment already activated by load_config.sh (config.ini env_setup)
# ─────────────────────────────────────────────────────────────────────────────
# The module Python may have no pip/user-site; provide extra packages via
# pythonpath_extra in config.ini.  (The sph VTK reader, sph_vtk.py, is pure
# numpy and needs no vtk package.)
if [[ -n "${PYTHONPATH_EXTRA}" ]]; then
    export PYTHONPATH="${PYTHONPATH_EXTRA}:${PYTHONPATH:-}"
fi

mkdir -p "${OUTPUT_DIR}"

# Build a Python list literal from SPH_DIRS, e.g. ['/a', '/b']
SPH_LIST_PY="["
for d in "${SPH_DIRS[@]}"; do SPH_LIST_PY+="'$d',"; done
SPH_LIST_PY="${SPH_LIST_PY%,}]"

echo "============================================================"
echo "  BNS SPH analysis"
echo "  SPH dirs   : ${#SPH_DIRS[@]} segment(s)"
echo "  Radius     : ejecta ${EJ_RADIUS}, neutrino ${RADIUS} M_sun"
echo "  Cache      : $([[ "${SPH_CACHE}" == "0" ]] && echo 'disabled (recompute all)' || echo "${OUTPUT_DIR}/cache")"
echo "  Output     : ${OUTPUT_DIR}"
echo "============================================================"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: MERGER TIME  (peak |h22| of the (2,2) GW strain)
# ─────────────────────────────────────────────────────────────────────────────
MERGER_FILE="${OUTPUT_DIR}/merger_time.txt"
echo ""
echo "── Merger time (GW (2,2) mode) ──────────────────────────────"
# Peak of |h22| (strain from psi4 by FFI); a run whose own waveform stops
# before merger falls back to its restart parent's. t_merger in config.ini (or
# T_MERGER in the environment) overrides the computation.
# The peak time is shifted to retarded time with the ADM mass (m_adm in
# config.ini, default 2.7 as in THC).
T_MERGER_OVERRIDE="${T_MERGER:-$(_cfg_get t_merger)}"
M_ADM="${M_ADM:-$(_cfg_get m_adm)}"
python3 "${SCRIPT_DIR}/merger_time.py" \
    --simpath "${SIMPATH}" --radius "${WF_RADIUS}" --out "${MERGER_FILE}" \
    ${T_MERGER_OVERRIDE:+--t-merger "${T_MERGER_OVERRIDE}"} \
    ${M_ADM:+--m-adm "${M_ADM}"} || exit 1
# Read the merger time (first non-comment numeric line) from the file.
T_MERGER_MSUN="$(awk '!/^#/ && NF {print $1; exit}' "${MERGER_FILE}")"
if [[ -z "${T_MERGER_MSUN}" ]]; then
    echo "ERROR: could not read merger time from ${MERGER_FILE}"; exit 1
fi
echo "  T_merger = ${T_MERGER_MSUN} M_sun"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 + 2: EJECTA and NEUTRINO  (run concurrently — they share no state)
# ─────────────────────────────────────────────────────────────────────────────
EJECTA_PID=""; NEUTRINO_PID=""

if [[ $RUN_EJECTA -eq 1 ]]; then
    echo ""
    echo "── Ejecta analysis (background) ─────────────────────────────"
    # Which EOS? Don't make the user keep eos_table in step with simpath by
    # hand: every AthenaK .bin output embeds the parfile it ran with, including
    # <mhd> `table`, so the run tells us whether it is SFHo, DD2, ... A non-empty
    # eos_table in config.ini still wins, as an explicit override.
    if [[ -z "${EOS_TABLE}" ]]; then
        EOS_TABLE="$(python3 "${SCRIPT_DIR}/eos_from_run.py" \
                         --simpath "${SIMPATH}" --jobname "${JOBNAME}")" || exit 1
        echo "  EOS table : ${EOS_TABLE}  (auto-detected from the run's own output)"
    else
        echo "  EOS table : ${EOS_TABLE}  (pinned in ${CONFIG_FILE})"
    fi
    if [[ ! -f "${EOS_TABLE}" ]]; then
        echo "ERROR: eos_table '${EOS_TABLE}' (from ${CONFIG_FILE}) does not exist."
        echo "       It must be the PyCompOSE HDF5 twin of the .athtab the run used."
        echo "       Build it from the .athtab with:"
        echo "         python3 ${SCRIPT_DIR}/athtab_to_h5.py <table>.athtab <table>.h5"
        exit 1
    fi
    python3 - << PYEOF &
import sys
sys.path.insert(0, '${PLOTTOOLS}')
sys.path.insert(0, '${SCRIPT_DIR}')
import analyze_ejecta as ae

ae.SPH_DIR_LIST   = ${SPH_LIST_PY}
ae.EOS_TABLE_PATH = '${EOS_TABLE}'
ae.RADIUS         = ${EJ_RADIUS}
ae.STRIDE         = ${STRIDE}
ae.DFLOOR         = ${DFLOOR}
ae.T_MERGER_MSUN  = ${T_MERGER_MSUN}
ae.T_POST_MS      = ${T_POST_MS}
ae.N_WORKERS      = ${N_WORKERS}
ae.USE_CACHE      = ${_USE_CACHE}
MSUN_TO_S = 6.67430e-11 * 1.98892e30 / 2.99792458e8**3
ae.T_STOP  = ae.T_MERGER_MSUN + ae.T_POST_MS * 1e-3 / MSUN_TO_S
ae.OUTPUT_DIR = '${OUTPUT_DIR}'
ae.main()
PYEOF
    EJECTA_PID=$!
fi

if [[ $RUN_NEUTRINO -eq 1 ]]; then
    echo "── Neutrino analysis (background) ───────────────────────────"
    python3 - << PYEOF &
import sys
sys.path.insert(0, '${PLOTTOOLS}')
sys.path.insert(0, '${SCRIPT_DIR}')
import analyze_neutrinos as an

an.SPH_DIR_LIST = ${SPH_LIST_PY}
an.RADIUS       = ${RADIUS}
an.STRIDE       = ${STRIDE}
an.N_WORKERS    = ${N_WORKERS}
an.USE_CACHE    = ${_USE_CACHE}
an.OUTPUT_DIR   = '${OUTPUT_DIR}'
an.main()
PYEOF
    NEUTRINO_PID=$!
fi

if [[ -n "${EJECTA_PID}" ]]; then
    wait "${EJECTA_PID}" && echo "── Ejecta done ──" || { echo "ERROR: ejecta failed"; exit 1; }
fi
if [[ -n "${NEUTRINO_PID}" ]]; then
    wait "${NEUTRINO_PID}" && echo "── Neutrino done ──" || { echo "ERROR: neutrino failed"; exit 1; }
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: PLOT
# ─────────────────────────────────────────────────────────────────────────────
if [[ $RUN_PLOT -eq 1 ]]; then
    echo ""
    echo "── Plotting ─────────────────────────────────────────────────"
    PLOT_ARGS=(--output-dir "${OUTPUT_DIR}" --t-merger "${T_MERGER_MSUN}" --radius "${RADIUS}" --ej-radius "${EJ_RADIUS}")
    [[ "${PLOT_FROM_MERGER}" == "1" ]] && PLOT_ARGS+=(--from-merger)
    python3 "${SCRIPT_DIR}/plot_analysis.py" "${PLOT_ARGS[@]}"
    echo "── Plotting done ────────────────────────────────────────────"
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: ACCRETION  (post-merger accretion rate from baryon bookkeeping)
# Needs the ejecta outputs (Mej_rate_*/vinf_avg_*) — run ejecta first, or in the
# same invocation (they are produced above).
# ─────────────────────────────────────────────────────────────────────────────
if [[ $RUN_ACCRETION -eq 1 ]]; then
    echo ""
    echo "── Accretion rate (baryon bookkeeping) ──────────────────────"
    if [[ -f "${OUTPUT_DIR}/Mej_rate_${ACCR_CRIT}.txt" ]]; then
        python3 "${SCRIPT_DIR}/analyze_accretion.py" \
            --simpath "${SIMPATH}" --analysis-dir "${OUTPUT_DIR}" \
            --jobname "${JOBNAME}" --r-ej "${EJ_RADIUS}" --crit "${ACCR_CRIT}"
        echo "── Accretion done ───────────────────────────────────────────"
    else
        echo "  Skipping: ${OUTPUT_DIR}/Mej_rate_${ACCR_CRIT}.txt not found."
        echo "  Run the ejecta step first (bash run_analysis.sh --ejecta)."
    fi
fi

echo ""
echo "All done. Results in: ${OUTPUT_DIR}"
