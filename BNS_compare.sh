#!/usr/bin/env bash
# =============================================================================
# BNS_compare.sh  —  compare opacity / radiation slices between TWO runs
#
# Wraps BNS_compare.py.  For each variable and matched output time it renders a
# 3-panel figure  [run A] [run B] [ratio A/B]  on the xz slice, reusing the same
# units / colormaps as BNS_all.py.  Mirrors the two-folder convention of
# sph_data_analysis/run_comparison.sh.
#
# Usage:
#   bash BNS_compare.sh                          # uses RUNS / VARS below
#   bash BNS_compare.sh DIRA LABELA DIRB LABELB  # override the two runs on CLI
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS  ← edit per comparison
# ─────────────────────────────────────────────────────────────────────────────
# Exactly two runs: "run_dir|label".  Each run_dir contains output-*/bin/.
RUNS=(
    "/lus/flare/projects/CompactBinaryMerger/yiqiu/Large_sim/bundle_run/VVLR_test_2|eq"
    "/lus/flare/projects/CompactBinaryMerger/yiqiu/Large_sim/bundle_run/VVLR_test_3|non-eq"
)

# Variables to compare, as base:species  (nue=0, anue=1, nux=2, anux=3).
# bases: eta_0 eta_1 abs_0 abs_1 E N
VARS=(eta_1:2 abs_1:2 eta_0:2 eta_1:0 abs_1:0)

UNITS=cgs              # code | cgs | ngs   (cgs/ngs -> axes in km, opacities scaled)
XLIM=(-60 60)          # x-range [M_sun]
ZLIM=(-60 60)          # z-range [M_sun]
PLANE=xz               # xz | xy  (xy needs genuine *_xy dumps of each variable)
NX=400; NZ=400         # rasterization resolution of the common grid
RMAX=3.0               # ratio panel spans 1/RMAX .. RMAX  (white = equal)
EVERY=1                # use every Nth matched time
OUTDIR=""              # empty = <run A>/compare_<PLANE>/

MAKE_MOVIES=1          # 1 = also encode one MP4 per variable from the frames
FRAMERATE=10           # movie frames per second
VCODEC=""              # ffmpeg video encoder; empty = auto (libx264, libopenh264)

# ─────────────────────────────────────────────────────────────────────────────
# CLI override:  DIRA LABELA DIRB LABELB
# ─────────────────────────────────────────────────────────────────────────────
if [[ $# -eq 4 ]]; then
    RUNS=("$1|$2" "$3|$4")
elif [[ $# -ne 0 ]]; then
    echo "Usage: $0 [DIRA LABELA DIRB LABELB]"; exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# SETUP  — machine environment from config.ini (env_setup: module/conda/venv)
# ─────────────────────────────────────────────────────────────────────────────
source "${SCRIPT_DIR}/load_config.sh"

# ─────────────────────────────────────────────────────────────────────────────
ARGS=()
for entry in "${RUNS[@]}"; do
    dir="${entry%%|*}"; label="${entry#*|}"
    [[ -d "${dir}" ]] || { echo "ERROR: run directory not found: ${dir}"; exit 1; }
    ARGS+=(--run "${dir}" --label "${label}")
done
ARGS+=(--vars "${VARS[@]}")
ARGS+=(--units "${UNITS}")
ARGS+=(--xlim "${XLIM[0]}" "${XLIM[1]}")
ARGS+=(--zlim "${ZLIM[0]}" "${ZLIM[1]}")
ARGS+=(--nx "${NX}" --nz "${NZ}" --rmax "${RMAX}" --every "${EVERY}")
ARGS+=(--plane "${PLANE}")
[[ -n "${OUTDIR}" ]] && ARGS+=(--outdir "${OUTDIR}")

echo "============================================================"
echo "  BNS opacity slice comparison (A vs B, + ratio)"
for entry in "${RUNS[@]}"; do echo "    - ${entry%%|*}  (${entry#*|})"; done
echo "  vars : ${VARS[*]}"
echo "============================================================"

python3 "${SCRIPT_DIR}/BNS_compare.py" "${ARGS[@]}"

# ─────────────────────────────────────────────────────────────────────────────
# Movies: one MP4 per variable from the cmp_<var>_<time>.png frames
# ─────────────────────────────────────────────────────────────────────────────
if [[ "${MAKE_MOVIES}" == 1 ]]; then
    FRAMEDIR="${OUTDIR:-${RUNS[0]%%|*}/compare_${PLANE}}"
    # H.264 MP4 plays in browsers and chats; mpeg4 only as a last resort.
    if [[ -z "${VCODEC}" ]]; then
        _encoders="$(ffmpeg -hide_banner -encoders 2>/dev/null || true)"
        for _c in libx264 libopenh264 mpeg4; do
            if grep -qw "${_c}" <<< "${_encoders}"; then
                VCODEC="${_c}"; break
            fi
        done
        if [[ "${VCODEC}" == mpeg4 ]]; then
            echo "WARNING: no H.264 encoder in ffmpeg; mpeg4 movies may not play in browsers"
        fi
    fi
    for var in "${VARS[@]}"; do
        name="${var//:/_}"
        ffmpeg -y -loglevel error -framerate "${FRAMERATE}" \
               -pattern_type glob -i "${FRAMEDIR}/cmp_${name}_*.png" \
               -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v "${VCODEC}" \
               -pix_fmt yuv420p -movflags +faststart \
               "${FRAMEDIR}/${name}.mp4" \
            && echo "  -> ${FRAMEDIR}/${name}.mp4" \
            || echo "WARNING: no movie for ${var} (no frames or ffmpeg error)"
    done
fi
