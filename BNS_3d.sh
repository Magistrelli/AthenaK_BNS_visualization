#!/usr/bin/env bash
# BNS_3d — 3D isosurface and orthoslice driver  (Aurora edition)
# --------------------------------------------------------------
# Renders the 3D prim snapshots (bns.prim3D.NNNNN.bin) found in every
# output-XXXX/ segment as Matplotlib orthoslices and/or PyVista isosurfaces,
# then stitches them into movies.
#
# Machine-specific PATHS live in config.ini (same file as BNS_all.sh):
#   cp config.example.ini config.ini   # then edit simpath / pythonpath_extra
#
# config.ini is a plain INI file (read, not executed). It is gitignored, so
# those paths never get committed. All run settings below stay in this script.
#
# Timing (256^3, crop 80, per snapshot, single core):
#   orthoslice  ~7 s       isosurface  ~17 s

set -uo pipefail

# ======================== internals (paths/env) ========================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOT3D_DIR="${SCRIPT_DIR}/plot_3d_data"
# All machine-specific settings (paths + environment) come from config.ini.
# load_config.sh reads it, defines _cfg_get, and activates the environment.
source "${SCRIPT_DIR}/load_config.sh"
# Several simulations in config.ini: run this script once for each of them.
_for_each_sim "$@"

SIMPATH="$(_cfg_get simpath)"
PYTHONPATH_EXTRA="$(_cfg_get pythonpath_extra)"
if [[ -z "${SIMPATH}" ]]; then
    echo "ERROR: 'simpath' is not set in ${CONFIG_FILE}."
    exit 1
fi

# 3D frames go to a separate directory so they don't collide with the 2D Figs.
FIGPATH="${SIMPATH}/Figs_3d"

# Python environment was activated by load_config.sh (config.ini's env_setup).
# pyvista / h5py / vtk / imageio-ffmpeg are not module-provided — install them
# and point pythonpath_extra at them in config.ini (see below).
PYTHON="${PYTHON:-python3.12}"
command -v "${PYTHON}" &>/dev/null || PYTHON=python3

# pythonpath_extra (vtk/pyvista/h5py/imageio-ffmpeg) + repo plot-tools (athplot),
# prepended so they take precedence over any module-provided copies.
export PYTHONPATH="${PYTHONPATH_EXTRA:+${PYTHONPATH_EXTRA}:}${SCRIPT_DIR}/plot-tools:${PYTHONPATH:-}"

# ======================== what to plot =================
PLOTS=(
    # orthoslice          # xy / xz / yz Matplotlib slices
    isosurface          # PyVista 3D nested isosurfaces
    # isosurface_colorT   # density isosurface colored by temperature
    # isosurface_colorYe  # density isosurface colored by Y_e  (composition)
    # volume              # density volume rendering (see VOLUME settings)
)

SCALARS=(
    log_dens            # log10(density)
    # temperature
    # log_bmag            # log10(|B|)  — Gauss in cgs, code-log10 in code units
    # current           # |curl B| (B-field "current" proxy)
    # vorticity         # |curl v|
    # qcriterion        # Q-criterion
)

# Overlay magnetic field lines on every isosurface: set to "B" (or "velocity"),
# empty to disable.  Adds noticeable render time.
STREAMLINES=""

# ======================== frame selection ==============
STRIDE=1              # render every Nth snapshot (1 = all)
START=0              # uncomment to override first snapshot index
# END=80               # uncomment to override last  snapshot index

# ======================== parallelism ==================
JOBS=4                 # 320^3 rendering is memory-heavy; raise only on a compute node

# ======================== units ========================
UNITS=cgs              # 'code' or 'cgs'  (cgs: density g/cm^3, |B| Gauss)
# UNITS=code

# ======================== grid / crop ==================
RESOLUTION=320         # uniform reconstruction cells per longest axis
INTERPOLATION=linear   # linear for smooth rho/T/Ye; nearest for native cells
IMAGE_INTERPOLATION=bilinear  # orthoslice display: nearest | bilinear | bicubic
CROP=80                # half-width crop: domain -> [-CROP, CROP]^3

# ======================== isosurface settings ==========
WINDOW_SIZE="1600,1600"
ANTI_ALIASING=ssaa     # none | ssaa | msaa | fxaa
CAMERA_ELEVATION=30
CAMERA_AZIMUTH=25
VIEW_ANGLE=35
BACKGROUND=black

# Per-scalar isosurface levels and color limits.
#   cgs log_dens : log10(6.178e17 * rho_code), so ~8-14 in g/cm^3.
#   cgs log_bmag : log10(|B| in Gauss), so ~9-15 (peak ~1e15 G).
#   temperature  : MeV in both code and cgs.
# NOTE: these levels assume UNITS=cgs.  For UNITS=code use the code-unit ranges
# noted in comments (e.g. log_bmag -13..-7, log_dens via 6e8..2e15 / 6.178e17).
declare -A ISO_LEVELS ISO_OPACITY ISO_CLIM ISO_CMAP
ISO_LEVELS[log_dens]="8,10,12,14"
ISO_OPACITY[log_dens]="0.08,0.15,0.30,0.70"
ISO_CLIM[log_dens]="6,15"
ISO_CMAP[log_dens]="inferno"

ISO_LEVELS[temperature]="5,15,30,45"
ISO_OPACITY[temperature]="0.08,0.18,0.35,0.70"
ISO_CLIM[temperature]="0,55"
ISO_CMAP[temperature]="viridis"

# cgs (Gauss), 1e12-1e15 G to match athenak_units.plot_params('bmag').
# For UNITS=code use ISO_LEVELS="-12,-10.5,-9,-7.5" / ISO_CLIM="-13,-7".
ISO_LEVELS[log_bmag]="12,13,14"
ISO_OPACITY[log_bmag]="0.2,0.35,0.8"
ISO_CLIM[log_bmag]="12,15"
ISO_CMAP[log_bmag]="plasma"

# current (|curl B|): levels/clim left empty -> auto percentiles per frame.
ISO_CMAP[current]="inferno"

ISO_CMAP[vorticity]="magma"
ISO_CMAP[qcriterion]="magma"

# ======================== timing / annotation ==========
# Merger time in CODE units.  When set, every frame is labelled
# "t - t_merger = +X.XXX ms" instead of absolute simulation time.
# Find it with sph_data_analysis/merger_time.py.
T0=""                  # e.g. T0=5310

# Orthographic projection makes the on-screen scale exact, which is what the
# ruler below relies on.  Leave PROJECTION=perspective for the classic look.
PROJECTION=perspective # 'perspective' or 'orthographic'
PARALLEL_SCALE=""      # orthographic half-height in M; empty = fit to domain
CAMERA_TRACK=""        # optional zoom-track file: "time_ms scale [ruler]" rows

# How the 3D frames are labelled.  "mathtext" composites the clock, colour bar
# and ruler with Matplotlib after the render, so maths is typeset properly
# (log_10 rho, Y_e, M_sun) exactly as in the 2D plots.  "vtk" uses VTK's own
# text, which has no maths layout and can only print "log10 rho [g cm^-3]".
ANNOTATE=mathtext      # mathtext | vtk

# Annotation typeface.  sans = DejaVu Sans (the matplotlib default, standard
# for figures), serif = STIX (Times-alike), cm = Computer Modern (LaTeX
# papers); arial/times/courier are VTK built-ins used only by ANNOTATE=vtk.
FONT=sans

SCALEBAR=""            # "auto" (1/2/5 ruler) or a length in M; empty = off
                       # NOTE: setting this forces orthographic projection.
AXES=false             # true = draw an orientation triad
SMOOTH_ITERATIONS=20   # volume-preserving Taubin passes per contour (0 = off)

# ======================== Y_e isosurface settings ======
# Composition painted on density isosurfaces.  Volume rendering cannot show
# Y_e once the ejecta is optically thick -- every ray then leaves through tens
# of M of Y_e ~ 0.4 outflow and the frame saturates -- so this stays surfaces.
# Levels are log10(rho) in the active UNITS (cgs: 11 and 13 g/cm^3).
# Keep YE_CLIM fixed across a run so colors mean the same thing in every frame.
COLOR_LEVELS="11,13"  # density shells used for both T and Ye colouring
COLOR_OPACITY="0.18,0.65"  # reveal the inner shell through the outer one
YE_LEVELS="${COLOR_LEVELS}"
YE_OPACITY="${COLOR_OPACITY}"
YE_CLIM="0.05,0.45"
YE_CMAP="coolwarm"

# ======================== volume settings ==============
# Density volume rendering.  Its colour ramp paints tenuous material dark,
# which is why density works here and composition does not.
VOL_SCALAR=log_dens
VOL_CLIM="7,15"        # cgs log10(g/cm^3); peak is ~15 near merger
VOL_CMAP="inferno"
VOL_OPACITY=sigmoid    # linear | sqrt | sigmoid | geom
VOL_OPACITY_SCALE=1.0
VOL_SAMPLES=""         # resample to N^3 before rendering; empty = keep grid
VOL_MODE=auto          # auto | volume | layers  (auto falls back to layers
                       # where the GL stack cannot composite a true volume)
VOL_LAYERS=24          # more, lighter shells reduce fallback banding
# A 480^3 float32 field is ~422 MiB and VTK keeps its own copy, so volume
# frames are the memory-hungry ones.  Jobs are capped for them alone.
VOL_MAX_JOBS=2

# ======================== movie settings ================
MOVIE_FPS=10               # frames per second
MOVIE_CRF=18               # H.264 quality: 0 (lossless) - 51 (worst)
MOVIE_FAST_AFTER=100       # keep these first frames at normal speed
MOVIE_FAST_STRIDE=4        # after that, keep every Nth frame (1 = disabled)
# MAKE_MOVIES=false        # uncomment to skip movie generation

# ======================== dry run ======================
# DRY_RUN=true prints the render commands it would run, then exits.
# Useful for checking new presets before spending node time.
DRY_RUN=${DRY_RUN:-false}

# ======================== skip existing ================
# SKIP_EXISTING=false
SKIP_EXISTING=true

# ======================== run banner ===================
mkdir -p "${FIGPATH}"
echo "SIMPATH : ${SIMPATH}"
echo "FIGPATH : ${FIGPATH}"
echo "PYTHON  : ${PYTHON} ($(${PYTHON} --version 2>&1))"
echo "UNITS   : ${UNITS}   SCALARS: ${SCALARS[*]}"
echo ""

contains() { local e; for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done; return 1; }

# ---------------------------------------------------------------------------
# Flags shared by every PyVista renderer (isosurface / volume)
# ---------------------------------------------------------------------------
COMMON_3D_FLAGS=""
COMMON_3D_FLAGS+=" --interpolation ${INTERPOLATION}"
COMMON_3D_FLAGS+=" --anti-aliasing ${ANTI_ALIASING}"
[[ -n "${T0}" ]]              && COMMON_3D_FLAGS+=" --t0 ${T0}"
[[ "${PROJECTION}" == "orthographic" ]] && COMMON_3D_FLAGS+=" --projection orthographic"
[[ -n "${PARALLEL_SCALE}" ]]  && COMMON_3D_FLAGS+=" --parallel-scale ${PARALLEL_SCALE}"
[[ -n "${CAMERA_TRACK}" ]]    && COMMON_3D_FLAGS+=" --camera-track ${CAMERA_TRACK}"
[[ -n "${SCALEBAR}" ]]        && COMMON_3D_FLAGS+=" --scalebar ${SCALEBAR}"
[[ "${AXES}" == "true" ]]     && COMMON_3D_FLAGS+=" --axes"
[[ -n "${FONT}" ]]            && COMMON_3D_FLAGS+=" --font ${FONT}"
[[ -n "${ANNOTATE}" ]]        && COMMON_3D_FLAGS+=" --annotate ${ANNOTATE}"

if [[ -n "${CAMERA_TRACK}" && ! -f "${CAMERA_TRACK}" ]]; then
    echo "ERROR: CAMERA_TRACK file not found: ${CAMERA_TRACK}"
    exit 1
fi

# Volume frames hold several copies of the grid, so they get their own cap.
EFFECTIVE_JOBS=${JOBS}
if contains volume "${PLOTS[@]}" && (( JOBS > VOL_MAX_JOBS )); then
    echo "NOTE: volume rendering is memory-heavy; reducing JOBS ${JOBS} -> ${VOL_MAX_JOBS}"
    echo "      (raise VOL_MAX_JOBS if you know the node can take it)"
    EFFECTIVE_JOBS=${VOL_MAX_JOBS}
fi

# ---------------------------------------------------------------------------
# Discover snapshot files across all output-XXXX segments
# ---------------------------------------------------------------------------
mapfile -t ALL_FILES < <(
    find "${SIMPATH}"/output-*/bin -name 'bns.prim3D.*.bin' 2>/dev/null | sort
)

if [[ ${#ALL_FILES[@]} -eq 0 ]]; then
    echo "ERROR: no bns.prim3D.*.bin files found under ${SIMPATH}/output-*/bin/"
    exit 1
fi
echo "Found ${#ALL_FILES[@]} prim3D snapshots under ${SIMPATH}"

# ---------------------------------------------------------------------------
# Apply stride / start / end filters
# ---------------------------------------------------------------------------
SELECTED=()
for f in "${ALL_FILES[@]}"; do
    idx=$(basename "$f" .bin | grep -oP '\d+$' | sed 's/^0*//' )
    idx=${idx:-0}
    [[ -n "${START:-}" && $idx -lt $START ]] && continue
    [[ -n "${END:-}"   && $idx -gt $END   ]] && continue
    (( idx % STRIDE == 0 )) || continue
    SELECTED+=("$f")
done
echo "Selected ${#SELECTED[@]} snapshots (stride=${STRIDE}), JOBS=${EFFECTIVE_JOBS}"

# ---------------------------------------------------------------------------
# Build command list — one line per plot
# ---------------------------------------------------------------------------
CMDFILE=$(mktemp "${TMPDIR:-/tmp}/bns3d_cmds.XXXXXX")
TOTAL_QUEUED=0
SKIPPED=0

for FILE in "${SELECTED[@]}"; do
    TAG=$(basename "$FILE" .bin | grep -oP '\d+$')

    for SCALAR in "${SCALARS[@]}"; do
        CMAP="${ISO_CMAP[$SCALAR]:-inferno}"

        # --- orthoslice ---
        if contains orthoslice "${PLOTS[@]}"; then
            OUT="${FIGPATH}/ortho_${SCALAR}_${TAG}.png"
            if [[ "${SKIP_EXISTING}" == "true" && -f "$OUT" ]]; then
                (( SKIPPED++ ))
            else
                echo "${PYTHON} ${PLOT3D_DIR}/plot_orthoslices.py \
-i ${FILE} -o ${OUT} \
--scalar ${SCALAR} --resolution ${RESOLUTION} --crop ${CROP} \
--interpolation ${INTERPOLATION} --image-interpolation ${IMAGE_INTERPOLATION} \
--cmap ${CMAP} --dpi 200 --units ${UNITS}" >> "$CMDFILE"
                (( TOTAL_QUEUED++ ))
            fi
        fi

        # --- isosurface ---
        if contains isosurface "${PLOTS[@]}"; then
            OUT="${FIGPATH}/iso_${SCALAR}_${TAG}.png"
            if [[ "${SKIP_EXISTING}" == "true" && -f "$OUT" ]]; then
                (( SKIPPED++ ))
            else
                CMD="${PYTHON} ${PLOT3D_DIR}/plot_isosurfaces.py \
-i ${FILE} -o ${OUT} \
--scalar ${SCALAR} --resolution ${RESOLUTION} --crop ${CROP} \
--cmap ${CMAP} --off-screen \
--window-size ${WINDOW_SIZE} \
--camera-elevation ${CAMERA_ELEVATION} --camera-azimuth ${CAMERA_AZIMUTH} \
--view-angle ${VIEW_ANGLE} --background ${BACKGROUND} --show-scalar-bar \
--units ${UNITS}"
                [[ -n "${ISO_LEVELS[$SCALAR]:-}" ]]  && CMD+=" --levels=${ISO_LEVELS[$SCALAR]}"
                [[ -n "${ISO_OPACITY[$SCALAR]:-}" ]] && CMD+=" --opacity ${ISO_OPACITY[$SCALAR]}"
                [[ -n "${ISO_CLIM[$SCALAR]:-}" ]]    && CMD+=" --clim=${ISO_CLIM[$SCALAR]}"
                [[ -n "${STREAMLINES}" ]]            && CMD+=" --streamlines ${STREAMLINES} --streamline-seeds dual-density-max"
                CMD+=" --smooth-iterations ${SMOOTH_ITERATIONS}"
                CMD+="${COMMON_3D_FLAGS}"
                echo "$CMD" >> "$CMDFILE"
                (( TOTAL_QUEUED++ ))
            fi
        fi

        # --- isosurface colored by temperature ---
        if contains isosurface_colorT "${PLOTS[@]}" && [[ "$SCALAR" != "temperature" ]]; then
            OUT="${FIGPATH}/iso_${SCALAR}_colorT_${TAG}.png"
            if [[ "${SKIP_EXISTING}" == "true" && -f "$OUT" ]]; then
                (( SKIPPED++ ))
            else
                CMD="${PYTHON} ${PLOT3D_DIR}/plot_isosurfaces.py \
-i ${FILE} -o ${OUT} \
--scalar ${SCALAR} --color-by temperature \
--resolution ${RESOLUTION} --crop ${CROP} \
--cmap viridis --clim=0,50 --off-screen \
--window-size ${WINDOW_SIZE} \
--camera-elevation ${CAMERA_ELEVATION} --camera-azimuth ${CAMERA_AZIMUTH} \
--view-angle ${VIEW_ANGLE} --background ${BACKGROUND} --show-scalar-bar \
--units ${UNITS}"
                CMD+=" --levels=${COLOR_LEVELS} --opacity ${COLOR_OPACITY}"
                [[ -n "${STREAMLINES}" ]]            && CMD+=" --streamlines ${STREAMLINES} --streamline-seeds dual-density-max"
                CMD+=" --smooth-iterations ${SMOOTH_ITERATIONS}"
                CMD+="${COMMON_3D_FLAGS}"
                echo "$CMD" >> "$CMDFILE"
                (( TOTAL_QUEUED++ ))
            fi
        fi

        # --- density isosurface colored by Y_e ---
        # Y_e comes from 'ye' in coarse-grained HDF5 boxes and from the passive
        # scalar 's_00' in AthenaK prim3D output; the loader accepts either.
        if contains isosurface_colorYe "${PLOTS[@]}"; then
            OUT="${FIGPATH}/iso_${SCALAR}_colorYe_${TAG}.png"
            if [[ "${SKIP_EXISTING}" == "true" && -f "$OUT" ]]; then
                (( SKIPPED++ ))
            else
                CMD="${PYTHON} ${PLOT3D_DIR}/plot_isosurfaces.py \
-i ${FILE} -o ${OUT} \
--scalar ${SCALAR} --color-by ye \
--resolution ${RESOLUTION} --crop ${CROP} \
--levels=${YE_LEVELS} --opacity ${YE_OPACITY} \
--clim=${YE_CLIM} --cmap ${YE_CMAP} --off-screen \
--window-size ${WINDOW_SIZE} \
--camera-elevation ${CAMERA_ELEVATION} --camera-azimuth ${CAMERA_AZIMUTH} \
--view-angle ${VIEW_ANGLE} --background ${BACKGROUND} --show-scalar-bar \
--units ${UNITS}"
                CMD+=" --smooth-iterations ${SMOOTH_ITERATIONS}"
                CMD+="${COMMON_3D_FLAGS}"
                echo "$CMD" >> "$CMDFILE"
                (( TOTAL_QUEUED++ ))
            fi
        fi

    done

    # --- density volume rendering (one per snapshot, not per scalar) ---
    if contains volume "${PLOTS[@]}"; then
        OUT="${FIGPATH}/vol_${VOL_SCALAR}_${TAG}.png"
        if [[ "${SKIP_EXISTING}" == "true" && -f "$OUT" ]]; then
            (( SKIPPED++ ))
        else
            CMD="${PYTHON} ${PLOT3D_DIR}/plot_volume.py \
-i ${FILE} -o ${OUT} \
--scalar ${VOL_SCALAR} --resolution ${RESOLUTION} --crop ${CROP} \
--smooth-iterations ${SMOOTH_ITERATIONS} \
--clim=${VOL_CLIM} --cmap ${VOL_CMAP} \
--opacity ${VOL_OPACITY} --opacity-scale ${VOL_OPACITY_SCALE} \
--mode ${VOL_MODE} --layers ${VOL_LAYERS} --off-screen \
--window-size ${WINDOW_SIZE} \
--camera-elevation ${CAMERA_ELEVATION} --camera-azimuth ${CAMERA_AZIMUTH} \
--view-angle ${VIEW_ANGLE} --background ${BACKGROUND} --show-scalar-bar \
--units ${UNITS}"
            [[ -n "${VOL_SAMPLES}" ]] && CMD+=" --samples ${VOL_SAMPLES}"
            CMD+="${COMMON_3D_FLAGS}"
            echo "$CMD" >> "$CMDFILE"
            (( TOTAL_QUEUED++ ))
        fi
    fi
done

echo "Queued ${TOTAL_QUEUED} plots, skipped ${SKIPPED}"

# ---------------------------------------------------------------------------
# Execute — serial or parallel.  Uses xargs -P (GNU parallel not required).
# ---------------------------------------------------------------------------
if [[ $TOTAL_QUEUED -eq 0 ]]; then
    echo "Nothing to render."
    rm -f "$CMDFILE"
    exit 0
fi

if [[ "${DRY_RUN}" == "true" ]]; then
    echo ""
    echo "=== DRY RUN: ${TOTAL_QUEUED} command(s) ==="
    cat "$CMDFILE"
    rm -f "$CMDFILE"
    exit 0
fi

START_TIME=$SECONDS

if [[ ${EFFECTIVE_JOBS} -le 1 ]]; then
    DONE=0
    while IFS= read -r cmd; do
        (( DONE++ ))
        echo "[${DONE}/${TOTAL_QUEUED}] $(echo "$cmd" | grep -oP '(ortho|iso|vol)\S*\.png')"
        eval "$cmd" 2>&1 | grep -vaE "WARN\|"
    done < "$CMDFILE"
else
    echo "Running ${TOTAL_QUEUED} plots, ${EFFECTIVE_JOBS} at a time (xargs)..."
    # Each line of CMDFILE is one full command; run it with bash -c, JOBS in flight.
    xargs -P "${EFFECTIVE_JOBS}" -I CMD bash -c CMD < "$CMDFILE" 2>&1 | grep -vaE "WARN\|"
fi

ELAPSED=$(( SECONDS - START_TIME ))
rm -f "$CMDFILE"

echo ""
echo "=== Done: ${TOTAL_QUEUED} plots in ${ELAPSED}s (${EFFECTIVE_JOBS} jobs), ${SKIPPED} skipped ==="
echo "Output in ${FIGPATH}"

# ---------------------------------------------------------------------------
# Movie generation — create one MP4 per frame sequence
# ---------------------------------------------------------------------------
if [[ "${MAKE_MOVIES:-true}" != "true" ]]; then
    echo "Movie generation skipped (MAKE_MOVIES=false)"
    exit 0
fi

# Prefer the bundled ffmpeg (has libx264) from imageio-ffmpeg over the system one
_FFMPEG="$("${PYTHON}" -c "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())" 2>/dev/null)" \
    || _FFMPEG=""
if [[ -z "${_FFMPEG}" ]] && command -v ffmpeg &>/dev/null; then
    _FFMPEG="ffmpeg"
fi
if [[ -z "${_FFMPEG}" ]]; then
    echo "WARNING: ffmpeg not found — skipping movie generation."
    echo "  pip install --target <pythonpath_extra> imageio-ffmpeg  to enable."
    exit 0
fi

# Auto-detect encoder
_VCODEC=""
_encoders=$("${_FFMPEG}" -hide_banner -encoders 2>/dev/null)
for _c in libx264 h264_nvenc libx265 mpeg4; do
    if echo "${_encoders}" | grep -q " ${_c} "; then
        _VCODEC="${_c}"
        break
    fi
done
[[ -z "${_VCODEC}" ]] && _VCODEC="gif"

case "${_VCODEC}" in
    libx264|libx265|h264_nvenc)  _QFLAGS=(-crf "${MOVIE_CRF}") ;;
    gif)                         _QFLAGS=()                     ;;
    *)                           _QFLAGS=(-q:v 5)               ;;
esac

echo ""
echo "=== Making movies (encoder: ${_VCODEC}, fps: ${MOVIE_FPS}) ==="

if ! [[ "${MOVIE_FAST_AFTER}" =~ ^[0-9]+$ ]] ||
   ! [[ "${MOVIE_FAST_STRIDE}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MOVIE_FAST_AFTER must be >= 0 and MOVIE_FAST_STRIDE must be >= 1." >&2
    exit 1
fi

_MOVIE_TIMING_FILTER=""
if (( MOVIE_FAST_STRIDE > 1 )); then
    # Preserve frames [0, MOVIE_FAST_AFTER) and decimate the rest.  Resetting
    # PTS plays the selected frames continuously at MOVIE_FPS, so the tail is
    # MOVIE_FAST_STRIDE times faster and requires fewer encoded frames.
    _MOVIE_TIMING_FILTER="select='if(lt(n\\,${MOVIE_FAST_AFTER})\\,1\\,not(mod(n-${MOVIE_FAST_AFTER}\\,${MOVIE_FAST_STRIDE})))',setpts=N/(${MOVIE_FPS}*TB),fps=${MOVIE_FPS},"
    echo "    timing: first ${MOVIE_FAST_AFTER} frames normal, then ${MOVIE_FAST_STRIDE}x faster"
fi

cd "${FIGPATH}"

# Discover unique prefixes from iso_/ortho_/vol_<name>_NNNNN.png
shopt -s nullglob
_ALL_FRAMES=( iso_*_?????.png ortho_*_?????.png vol_*_?????.png )
shopt -u nullglob

if [[ ${#_ALL_FRAMES[@]} -eq 0 ]]; then
    echo "No frame sequences found in ${FIGPATH}"
    exit 0
fi

unset _seen; declare -A _seen
_PREFIXES=()
for f in "${_ALL_FRAMES[@]}"; do
    p="${f%_?????.png}"
    if [[ -z "${_seen[$p]+x}" ]]; then
        _seen[$p]=1
        _PREFIXES+=("$p")
    fi
done

MOVIEDIR="${FIGPATH}/movies"
mkdir -p "${MOVIEDIR}"

_MOVIE_COUNT=0
for prefix in "${_PREFIXES[@]}"; do
    shopt -s nullglob
    # The five-digit suffix makes the shell's lexical order numerical.  Keep
    # this as a real globbed array (rather than parsing `ls`) so unusual paths
    # cannot split frame names.
    FRAMES=( "${prefix}"_?????.png )
    shopt -u nullglob

    [[ ${#FRAMES[@]} -lt 2 ]] && continue

    [[ "${_VCODEC}" == "gif" ]] && _ext="gif" || _ext="mp4"
    outfile="${MOVIEDIR}/${prefix}.${_ext}"
    # Encode beside the destination and publish with one atomic rename.  An
    # interrupted batch then leaves the last good movie intact instead of an
    # MP4 containing only an unfinished mdat atom and no moov index.
    tmp_out="${MOVIEDIR}/.${prefix}.tmp.${BASHPID}.${_ext}"
    rm -f -- "${tmp_out}"

    echo "  ${prefix}: ${#FRAMES[@]} frames -> $(basename "${outfile}")"

    # image2's glob input assigns one timestamp per PNG at -framerate.  The
    # previous concat list contained no duration records, which can produce
    # invalid or collapsed DTS/PTS values for single-image files.
    input_glob="${FIGPATH}/${prefix}_?????.png"
    _encode_ok=false
    if [[ "${_VCODEC}" == "gif" ]]; then
        if "${_FFMPEG}" -y \
            -framerate "${MOVIE_FPS}" -pattern_type glob -i "${input_glob}" \
            -vf "${_MOVIE_TIMING_FILTER}scale=trunc(iw/2)*2:trunc(ih/2)*2:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse=dither=bayer" \
            -loop 0 \
            "${tmp_out}" \
            -loglevel warning -stats; then
            _encode_ok=true
        fi
    else
        if "${_FFMPEG}" -y \
            -framerate "${MOVIE_FPS}" -pattern_type glob -i "${input_glob}" \
            -vf "${_MOVIE_TIMING_FILTER}scale=trunc(iw/2)*2:trunc(ih/2)*2" \
            -c:v "${_VCODEC}" -pix_fmt yuv420p \
            "${_QFLAGS[@]}" \
            -movflags +faststart \
            "${tmp_out}" \
            -loglevel warning -stats; then
            _encode_ok=true
        fi
    fi

    # A successful encoder exit should already imply a finalized container;
    # open and decode one frame as a cheap guard before replacing outfile.
    if [[ "${_encode_ok}" == "true" ]] && \
       "${_FFMPEG}" -v error -i "${tmp_out}" -map 0:v:0 -frames:v 1 -f null -; then
        mv -f -- "${tmp_out}" "${outfile}"
        (( _MOVIE_COUNT++ ))
    else
        echo "WARNING: ${prefix} movie failed validation; keeping any existing ${outfile}" >&2
        rm -f -- "${tmp_out}"
    fi
done

echo ""
echo "=== ${_MOVIE_COUNT} movie(s) saved to ${MOVIEDIR} ==="
