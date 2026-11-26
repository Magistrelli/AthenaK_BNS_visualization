#!/usr/bin/env bash
# movie_tools.sh — shared movie encoding for make_movies.sh and BNS_compare.sh.
# -----------------------------------------------------------------------------
# SOURCE this file. Settings (set before calling, defaults below):
#   FPS, STRIDE (use every STRIDEth frame), CRF, VCODEC (empty = auto),
#   PACK_FRAMES (1 = after a movie is made, pack its frames into frames.tar).
# Functions:
#   _movie_setup                    find ffmpeg and the encoder; sets _EXT
#   _encode_movie OUT FRAME...      encode the frames, in the order given
#   _movie_folder DIR OUT           movie of all frames of a frame folder
#
# A frame folder holds the PNG frames of one series, and after its movie was
# made only frames.tar (all frames, exactly as rendered: extract them with
# `tar -xf frames.tar`) plus frames.json (see frame_cache.py).  New or changed
# frames rendered later lie next to frames.tar as PNGs; the next _movie_folder
# unpacks the rest, re-encodes the movie with all frames, and packs them again.
# A movie that is newer than its frames, with the same settings, is kept.

: "${FPS:=10}"
: "${STRIDE:=1}"
: "${CRF:=18}"            # H.264 quality: 0 (lossless) – 51 (worst); 18 ≈ visually lossless
: "${VCODEC:=}"
: "${PACK_FRAMES:=1}"

# Locate ffmpeg: in PATH, else the Lmod `ffmpeg` module (e.g. Aurora, Hunter),
# else the pip-installed imageio-ffmpeg binary (pythonpath_extra in config.ini).
# Then pick the encoder: libx264 > h264_nvenc > libx265 > libopenh264 > mpeg4,
# else GIF.
_movie_setup() {
    [[ -n "${FFMPEG:-}" ]] && return 0
    if ! command -v ffmpeg &>/dev/null; then
        if ! command -v module &>/dev/null; then
            for _init in "${MODULESHOME:-/usr/share/lmod/lmod}/init/bash" \
                         /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh; do
                [[ -r "${_init}" ]] && source "${_init}" && break
            done
        fi
        command -v module &>/dev/null && module load ffmpeg &>/dev/null
    fi
    if command -v ffmpeg &>/dev/null; then
        FFMPEG="ffmpeg"
    else
        local dir ppx py
        dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        ppx="$(awk -F= '/^[[:space:]]*#/{next}
            /^[[:space:]]*pythonpath_extra[[:space:]]*=/{
                gsub(/^[[:space:]]+|[[:space:]]+$/,"",$2); print $2; exit}' \
            "${CONFIG_FILE:-${dir}/config.ini}" 2>/dev/null)"
        command -v module &>/dev/null && module load python/3.12.12 &>/dev/null
        for py in python3.12 python3 python; do
            command -v "${py}" &>/dev/null || continue
            FFMPEG="$(PYTHONPATH="${ppx}:${PYTHONPATH:-}" "${py}" -c \
                "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())" 2>/dev/null)"
            [[ -n "${FFMPEG}" ]] && break
        done
    fi
    if [[ -z "${FFMPEG:-}" ]]; then
        echo "ERROR: ffmpeg not found."
        echo "  module load ffmpeg, conda install -c conda-forge ffmpeg, or"
        echo "  python3 -m pip install --target <pythonpath_extra> imageio-ffmpeg"
        return 1
    fi
    if [[ -z "${VCODEC}" ]]; then
        local encoders c
        encoders="$("${FFMPEG}" -hide_banner -encoders 2>/dev/null)"
        for c in libx264 h264_nvenc libx265 libopenh264 mpeg4; do
            grep -q " ${c} " <<< "${encoders}" && { VCODEC="${c}"; break; }
        done
        [[ -z "${VCODEC}" ]] && VCODEC="gif"
        [[ "${VCODEC}" == mpeg4 ]] && echo "WARNING: no H.264 encoder; mpeg4 movies may not play in browsers"
    fi
    case "${VCODEC}" in
        libx264|libx265|h264_nvenc)  QUALITY_FLAGS=(-crf "${CRF}") ;;
        libopenh264)                 QUALITY_FLAGS=(-b:v 8M)        ;;
        gif)                         QUALITY_FLAGS=()               ;;
        *)                           QUALITY_FLAGS=(-q:v 5)         ;;
    esac
    [[ "${VCODEC}" == "gif" ]] && _EXT="gif" || _EXT="mp4"
    echo "ffmpeg  : ${FFMPEG}  (encoder ${VCODEC})"
}

_encode_movie() {
    local outfile="$1"; shift
    local listfile f status
    listfile="$(mktemp "${TMPDIR:-/tmp}/mm_fflist.XXXXXX")"
    for f in "$@"; do printf "file '%s'\n" "${f}" >> "${listfile}"; done
    if [[ "${VCODEC}" == "gif" ]]; then
        # High-quality GIF via palette generation (single-pass filtergraph).
        "${FFMPEG}" -y -r "${FPS}" -f concat -safe 0 -i "${listfile}" \
            -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse=dither=bayer" \
            -loop 0 "${outfile}" -loglevel warning -stats
    else
        "${FFMPEG}" -y -r "${FPS}" -f concat -safe 0 -i "${listfile}" \
            -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
            -c:v "${VCODEC}" -pix_fmt yuv420p "${QUALITY_FLAGS[@]}" \
            -movflags +faststart -f mp4 "${outfile}" -loglevel warning -stats
    fi
    status=$?
    rm -f "${listfile}"
    return "${status}"
}

# Every STRIDEth of the (sorted) arguments, one per line.
_subsample() {
    local -a all=( "$@" )
    local j
    for (( j=0; j<${#all[@]}; j+=STRIDE )); do printf '%s\n' "${all[j]}"; done
}

_movie_folder() {
    local dir="${1%/}" out="$2"
    local name settings loose frames used n
    name="$(basename "${dir}")"
    settings="fps=${FPS} stride=${STRIDE} codec=${VCODEC} crf=${CRF} out=$(basename "${out}")"
    shopt -s nullglob
    loose=( "${dir}"/*.png )
    shopt -u nullglob

    if [[ ${#loose[@]} -eq 0 && ! -f "${dir}/frames.tar" ]]; then
        return 0                                    # no frames at all
    fi
    if [[ ${#loose[@]} -eq 0 && -f "${out}" && ! "${dir}/frames.tar" -nt "${out}" \
          && "$(cat "${dir}/movie.settings" 2>/dev/null)" == "${settings}" ]]; then
        echo "[${name}]  movie up to date"
        return 0
    fi

    # Unpack the frames packed before, without replacing newer re-rendered ones.
    if [[ -f "${dir}/frames.tar" ]]; then
        tar -xf "${dir}/frames.tar" -C "${dir}" --skip-old-files || return 1
    fi
    shopt -s nullglob
    frames=( "${dir}"/*.png )
    shopt -u nullglob
    IFS=$'\n' frames=( $(printf '%s\n' "${frames[@]}" | sort) ); unset IFS
    IFS=$'\n' used=( $(_subsample "${frames[@]}") ); unset IFS
    echo "[${name}]  ${#frames[@]} frames, ${#used[@]} used  ->  $(basename "${out}")"

    if ! _encode_movie "${out}.tmp.${_EXT}" "${used[@]}"; then
        rm -f "${out}.tmp.${_EXT}"
        echo "WARNING: no movie for ${name}; its frames are kept as PNGs"
        return 1
    fi
    mv "${out}.tmp.${_EXT}" "${out}"
    echo "${settings}" > "${dir}/movie.settings"

    [[ "${PACK_FRAMES}" == 1 ]] || return 0
    # Pack all frames into frames.tar; remove the PNGs only once the archive
    # is complete.
    n=${#frames[@]}
    if printf '%s\n' "${frames[@]##*/}" | tar -cf "${dir}/frames.tar.tmp" -C "${dir}" -T - \
       && [[ "$(tar -tf "${dir}/frames.tar.tmp" | wc -l)" -eq "${n}" ]]; then
        mv "${dir}/frames.tar.tmp" "${dir}/frames.tar"
        touch -r "${dir}/frames.tar" "${out}"       # movie as new as its frames
        rm -f "${frames[@]}"
        echo "  packed ${n} frames into ${name}/frames.tar"
    else
        rm -f "${dir}/frames.tar.tmp"
        echo "WARNING: could not pack the frames of ${name}; they are kept as PNGs"
    fi
}
