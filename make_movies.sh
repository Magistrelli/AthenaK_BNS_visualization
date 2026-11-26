#!/usr/bin/env bash
# make_movies.sh
#
# Creates one MP4 per frame sequence found in FIGPATH.  Handles both layouts:
#   * per-diagnostic subfolders  <name>_xy/<...>_NNNNN.png   (2-D BNS_all driver)
#   * flat sequences             <prefix>_NNNNN.png          (3-D BNS_3d: iso_*,
#                                                             ortho_*; legacy 2-D)
# STRIDE controls sub-sampling: 1 = every frame, 10 = every 10th frame, etc.
# In the subfolder layout, a movie is only re-encoded when its frames or
# settings changed, and afterwards its frames are packed into
# <subfolder>/frames.tar (PACK_FRAMES; see movie_tools.sh).
#
# Usage:
#   bash make_movies.sh               # use settings below
#   bash make_movies.sh 25 1          # FPS=25, every frame
#   bash make_movies.sh 25 10         # FPS=25, every 10th frame
#   bash make_movies.sh 25 10 /path   # also override FIGPATH

# Log everything to a file so errors are readable even if the terminal closes.
_LOGFILE="${HOME}/make_movies.log"
exec > >(tee -a "${_LOGFILE}") 2>&1
echo "===== $(date) ====="
echo "Log: ${_LOGFILE}"

# ---------------------------------------------------------------------------
# Configuration — edit these
# ---------------------------------------------------------------------------
# Directory containing the rendered frames (subfolders or flat *_NNNNN.png).
# Usually passed as the 3rd argument by BNS_all.sh / BNS_3d.sh; edit here only
# if you run make_movies.sh standalone.
FIGPATH="/path/to/your/simulation/run/Figs"

# Frames per second in the output movie.
# 10–15  → slow, good for inspecting individual frames
# 24–30  → cinematic / standard
FPS=10

# Use every STRIDEth frame from the available PNGs.
# 1  → all frames (00000, 00001, 00002, ...)
# 10 → every 10th (00000, 00010, 00020, ...)
STRIDE=1

# H.264 quality: 0 (lossless) – 51 (worst).  18 ≈ visually lossless.
# Only used by encoders that support -crf (libx264, libx265, h264_nvenc).
CRF=18

# Video encoder.  Leave empty for auto-detection (recommended).
# Override if needed, e.g.:  VCODEC="mpeg4"
VCODEC=""

# 1 = once a subfolder's movie is made, pack its frames into frames.tar (exact
# copies, `tar -xf frames.tar` restores them) instead of keeping thousands of
# PNGs; reruns unpack them when new frames need a new movie.  Env overrides.
PACK_FRAMES="${PACK_FRAMES:-1}"
# ---------------------------------------------------------------------------

# Command-line overrides
FPS="${1:-${FPS}}"
STRIDE="${2:-${STRIDE}}"
FIGPATH="${3:-${FIGPATH}}"

# ---------------------------------------------------------------------------
# ffmpeg, encoder and the per-folder movie logic are shared with BNS_compare.sh.
# ---------------------------------------------------------------------------
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/movie_tools.sh"
_movie_setup || { return 1 2>/dev/null || exit 1; }

if [[ ! -d "${FIGPATH}" ]]; then
    echo "ERROR: FIGPATH not found: ${FIGPATH}"
    echo "  Edit FIGPATH at the top of make_movies.sh to point to your Figs/ directory."
    return 1 2>/dev/null || exit 1
fi

# Absolutize so the ffmpeg concat lists are unambiguous regardless of CWD.
FIGPATH="$(cd "${FIGPATH}" && pwd)"

MOVIEDIR="${FIGPATH}/all_movies"
mkdir -p "${MOVIEDIR}"

echo "FIGPATH : ${FIGPATH}"
echo "FPS     : ${FPS}"
echo "STRIDE  : ${STRIDE}"

_n_made=0

# ---------------------------------------------------------------------------
# 1) Per-diagnostic subfolders (2-D driver layout):  <name>_xy/*_NNNNN.png
#    One movie per subfolder, written straight into all_movies/.
# ---------------------------------------------------------------------------
shopt -s nullglob
SUBDIRS=( "${FIGPATH}"/*/ )
shopt -u nullglob

for d in "${SUBDIRS[@]}"; do
    d="${d%/}"
    [[ "$(basename "${d}")" == "all_movies" ]] && continue
    shopt -s nullglob
    _frames=( "${d}"/*.png "${d}"/frames.tar )
    shopt -u nullglob
    [[ ${#_frames[@]} -eq 0 ]] && continue
    _movie_folder "${d}" "${MOVIEDIR}/$(basename "${d}").${_EXT}"
    _n_made=$((_n_made + 1))
    echo ""
done

# ---------------------------------------------------------------------------
# 2) Flat frame sequences (3-D BNS_3d Figs_3d, or legacy flat 2-D output):
#    <prefix>_NNNNN.png directly in FIGPATH.  Only used when no subfolder
#    movies were made, so the two layouts never double-count.
# ---------------------------------------------------------------------------
if [[ ${_n_made} -eq 0 ]]; then
    shopt -s nullglob
    PNGS=( "${FIGPATH}"/bns_frame_*_?????.png \
           "${FIGPATH}"/iso_*_?????.png \
           "${FIGPATH}"/ortho_*_?????.png )
    shopt -u nullglob

    unset _seen; declare -A _seen
    PREFIXES=()
    for f in "${PNGS[@]}"; do
        b="$(basename "${f}")"; p="${b%_?????.png}"
        if [[ -z "${_seen[$p]+x}" ]]; then
            _seen[$p]=1
            PREFIXES+=("$p")
        fi
    done

    for p in "${PREFIXES[@]}"; do
        shopt -s nullglob
        ALL_FRAMES=( "${FIGPATH}/${p}"_?????.png )
        shopt -u nullglob
        [[ ${#ALL_FRAMES[@]} -eq 0 ]] && continue

        IFS=$'\n' ALL_FRAMES=( $(printf '%s\n' "${ALL_FRAMES[@]}" | sort) ); unset IFS
        IFS=$'\n' FRAMES=( $(_subsample "${ALL_FRAMES[@]}") ); unset IFS

        echo "[${p}]  ${#ALL_FRAMES[@]} available, ${#FRAMES[@]} used  ->  ${p}.${_EXT}"
        _encode_movie "${MOVIEDIR}/${p}.${_EXT}" "${FRAMES[@]}"
        _n_made=$((_n_made + 1))
        echo ""
    done
fi

if [[ ${_n_made} -eq 0 ]]; then
    echo "No frame sequences found under ${FIGPATH}"
    echo "  (looked for <name>/*_NNNNN.png subfolders and flat"
    echo "   bns_frame_*/iso_*/ortho_*_NNNNN.png)"
else
    echo "Made ${_n_made} movie(s) in ${MOVIEDIR}/"
fi
echo "Done."
