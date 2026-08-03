#!/usr/bin/env bash
# BNS_all visualization driver
# -----------------------------
# Machine-specific PATHS live in config.ini (copy it from config.example.ini):
#   cp config.example.ini config.ini   # then edit config.ini for your machine
#
# config.ini is a plain INI file (read, not executed) holding simpath and an
# optional figpath. It is gitignored, so those paths never get committed. All
# the run settings below (which diagnostics to plot, units, movie options) stay
# here in the script — edit them per run.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Paths — read from config.ini
# simpath should be the parent directory that contains the output-XXXX
# restart-segment subdirectories, e.g.:
#   DD2_m1_VVLR/
#     output-0000/
#     output-0001/
#     ...
# ---------------------------------------------------------------------------
# All machine-specific settings (paths + environment) come from config.ini.
# load_config.sh reads it, defines _cfg_get, and activates the environment.
source "${SCRIPT_DIR}/load_config.sh"
# Several simulations in config.ini: run this script once for each of them.
_for_each_sim "$@"

SIMPATH="$(_cfg_get simpath)"
FIGPATH="$(_cfg_get figpath)"
if [[ -z "${SIMPATH}" ]]; then
    echo "ERROR: 'simpath' is not set in ${CONFIG_FILE}."
    exit 1
fi
: "${FIGPATH:=${SIMPATH}/Figs}"      # default: <simpath>/Figs

# ---------------------------------------------------------------------------
# Run settings — edit per run
# ---------------------------------------------------------------------------
# Set to "true" to show NS tracker markers on 2-D plots
SHOW_TRACKERS=false

# Set to "true" to skip frames whose PNG already exists in the Figs directory
# (the tests below lowercase these, so True/TRUE/true all work)
SKIP_EXISTING=true

# Set to "true" to show the full simulation domain instead of the ±60 Msun window
FULL_DOMAIN=false

# Unit system for physical quantities: 'code', 'cgs', or 'ngs'  (default: cgs)
UNITS=cgs

# AthenaK job name prefixing all input/output files (bns.*, bns_frame_*, ...).
# Set to bhns/bbh/... to visualize other systems with the same driver.
JOBNAME=bns

# Movie encoding (passed to make_movies.sh)
MOVIE_FPS=10        # frames per second
MOVIE_STRIDE=10     # use every Nth frame (1 = all)
PACK_FRAMES=1       # 1 = after its movie is made, pack a section's frames into frames.tar

# ---------------------------------------------------------------------------
# SECTIONS — comment out any line to disable that diagnostic
# ---------------------------------------------------------------------------
SECTIONS=(
    history       # baryonic mass & Hamiltonian constraint norm vs time
    rho_max       # central density & min lapse vs time
    # trackers      # NS trajectory plot (both stars, smoothed)
    run_speed     # AthenaK cycles/sec from log file

    # density       # density (dens) slice frames  [+ trackers + AMR grid]
    # temperature   # temperature slice frames      [+ trackers + AMR grid]
    # s00           # electron fraction (s_00) frames
    # s00_avg       # area-weighted average Y_e vs time (1-D timeseries)
    ynue          # neutrino fraction Y_nue  = n_nue/n_b  [bns.prim_xy + bns.m1ff_xy n:0]
    # ynua          # neutrino fraction Y_anue = n_anue/n_b [bns.prim_xy + bns.m1ff_xy n:1]
    # ynux          # neutrino fraction Y_nux  (species 2)  [bns.prim_xy + bns.m1ff_xy n:2]
    # ynuax         # neutrino fraction Y_anux (species 3, 4-species runs only)
    nu_y          # bns_nurates degeneracy y = n(n/J)^3/(4pi/(hc)^3), nu_e; eq fallback at y >= 64/81  [bns.m1ff_xy]
    # nu_y1         # same, anti-nu_e
    # nu_y2         # same, nu_x  (bns_nurates gets n/2, J/2)
    # nu_y3         # same, anti-nu_x (4-species runs only)
    # bnorm         # magnetic field magnitude |B| slice frames   [bns.prim_xy bcc1-3]
    # bmax          # global max |B| & location vs time (1-D timeseries) [bns.prim3D]
    # z4c
    # rad_E         # radiation energy E:0 frames          [bns.rad_m1_E]
    # rad_E2         # radiation energy E:0 frames          [bns.rad_m1_E]
    # rad_N         # radiation number density N:0 frames  [bns.rad_m1_N]
    # rad_N2        # radiation number density N:2 (heavy lepton) frames
    # rad_abs       # absorption opacity abs_1:0 frames    [bns.rad_m1_abs_1]
     rad_abs0    # absorption opacity abs_0:0 frames    [bns.rad_m1_abs_0]
     rad_eta0    # emissivity eta_0:0 frames            [bns.rad_m1_eta_0]
    # rad_eta1    # emissivity eta_1:0 frames            [bns.rad_m1_eta_1]
    # avg_energy  # derived average neutrino energy frames
    # nu_energy   # neutrino energy E/N species 0 (ν_e)    [MeV in cgs]
    # nu_energy1  # neutrino energy E/N species 1 (ν̄_e)    [MeV in cgs]
    # nu_energy2  # neutrino energy E/N species 2 (ν_x)    [MeV in cgs]
    # --- fluid-frame moments & bns_nurates reconstruction [need bns.rad_J / bns.rad_n] ---
    # rad_J       # fluid-frame energy density J:0 frames  [bns.rad_J]  [MeV/nm^3 in ngs]
    # rad_n       # fluid-frame number density n:0 frames  [bns.rad_n]  [nm^-3 in ngs]
    # nu_avge     # fluid-frame mean energy <e>=J/n (ν_e)  [MeV]
    # nu_eta_t    # reconstructed degeneracy η_t = μ/T (clamped at 20)
    # nu_temp_t   # reconstructed spectral temperature T_ν (ν_e)  [MeV]
    # nu_Ynu      # neutrino fraction per baryon Y_ν = n_ν/n_b (ν_e)  [needs bns.prim_xy]
    # --- heavy-lepton (ν_x, species 2) reconstruction: 0.5 μ+τ weighting applied ---
    # nux_avge    # fluid-frame mean energy <e>=J/n (ν_x)  [MeV]
    # nux_y       # reconstructed degeneracy invariant y (ν_x, ceiling 0.5·64/81)
    # nux_eta_t   # reconstructed degeneracy η_t = μ/T (ν_x, uncapped)
    # nux_temp_t  # reconstructed spectral temperature T_ν (ν_x)  [MeV]
)
# ---------------------------------------------------------------------------

# (Environment already activated by load_config.sh, above.)

ARGS=(--simpath "${SIMPATH}")
[[ -n "${FIGPATH:-}"  ]] && ARGS+=(--figpath  "${FIGPATH}")
[[ -n "${JOBNAME:-}"  ]] && ARGS+=(--jobname  "${JOBNAME}")
[[ "${SHOW_TRACKERS,,}" == "false" ]] && ARGS+=(--no-trackers)
[[ "${SKIP_EXISTING,,}" == "true"  ]] && ARGS+=(--skip-existing)
[[ "${FULL_DOMAIN,,}"   == "true"  ]] && ARGS+=(--full-domain)
ARGS+=(--units "${UNITS}")
ARGS+=(--sections "${SECTIONS[@]}")

python3 "${SCRIPT_DIR}/BNS_all.py" "${ARGS[@]}"

# ---------------------------------------------------------------------------
# Movie generation — automatically create MP4s from rendered frames
# ---------------------------------------------------------------------------
echo ""
echo "=== Making movies from ${FIGPATH} ==="
PACK_FRAMES="${PACK_FRAMES}" \
    bash "${SCRIPT_DIR}/make_movies.sh" "${MOVIE_FPS}" "${MOVIE_STRIDE}" "${FIGPATH}"
