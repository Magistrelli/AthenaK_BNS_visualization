#!/usr/bin/env bash
# =============================================================================
# archive_run.sh  —  copy what later analyses need from an AthenaK run
#
# Usage:
#   bash archive_run.sh [--dry-run] [--radii "100 200 300 400"] [-j N] SRC DEST
#
# Copies the run folder SRC into DEST, keeping the AthenaK layout
# (DEST/output-XXXX/...), so every script in this repo works on the archive as
# on the run itself:
#   keep      everything small: parfiles, logs, *.hst, bns.co_*.txt,
#             waveforms/, and the run-level folders (analysis*/, Figs/, ...)
#   reduce    sph/: only the shells nearest RADII (sph_data_analysis/sph_extract.py)
#   compress  pdf_*/: every file gzipped (np.loadtxt reads .gz directly)
#   drop      bin/ (2D/3D dumps) and rst/ (restart files)
# A run without output-*/ (all outputs in one folder, e.g. data passed on by
# collaborators) goes to DEST/output-0000/, except its run-level folders.
#
# Nothing in SRC is changed.  Re-running updates DEST: files already archived
# and unchanged in SRC are skipped.  Check with --dry-run first: it lists what
# would be kept, reduced, compressed and dropped, with sizes, and copies nothing.
# Every real run appends what it did to DEST/ARCHIVE_INFO.txt.
# Needs only bash, coreutils, gzip and python3 (no packages).
# =============================================================================

set -euo pipefail
shopt -s nullglob dotglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
RADII="100 200 300 400"        # sph shells to keep [M_sun]
JOBS=4                         # parallel workers for sph and pdf files
DROP=(bin rst)                 # segment folders that are not archived
# Folders of a run without output-*/ that belong to the run, not to a segment.
RUN_LEVEL_RE='^(analysis.*|Figs|BATCH|compare.*|comparison.*)$'

usage() { echo "Usage: $0 [--dry-run] [--radii \"R1 R2 ...\"] [-j N] SRC DEST"; exit 1; }

DRY=0
while [[ $# -gt 0 && $1 == -* ]]; do
    case "$1" in
        -n|--dry-run) DRY=1; shift ;;
        --radii)      RADII="${2:?}"; shift 2 ;;
        -j|--jobs)    JOBS="${2:?}"; shift 2 ;;
        *) usage ;;
    esac
done
[[ $# -eq 2 ]] || usage
SRC="$(realpath "$1")"
DEST="$(realpath -m "$2")"
[[ -d "${SRC}" ]] || { echo "ERROR: run folder not found: ${SRC}"; exit 1; }
case "${DEST}/" in
    "${SRC}/"*) echo "ERROR: DEST must not be inside SRC"; exit 1 ;;
esac
case "${SRC}/" in
    "${DEST}/"*) echo "ERROR: SRC must not be inside DEST"; exit 1 ;;
esac

NKEEP=$(wc -w <<< "${RADII}")
EXTRACT="${SCRIPT_DIR}/sph_data_analysis/sph_extract.py"
LOG="${DEST}/ARCHIVE_INFO.txt"

size_of() { du -sb "$1" 2>/dev/null | cut -f1; }
human()   { numfmt --to=iec --suffix=B "$1"; }

TOTAL_IN=0; TOTAL_OUT=0
report() {   # ACTION SIZE_IN SIZE_OUT PATH  (SIZE_OUT empty = not known yet)
    local out="..."
    TOTAL_IN=$(( TOTAL_IN + $2 ))
    if [[ -n "$3" ]]; then
        TOTAL_OUT=$(( TOTAL_OUT + $3 )); out=$(human "$3")
    fi
    printf '%-9s %8s -> %8s  %s\n' "$1" "$(human "$2")" "${out}" "${4#"${SRC}"/}"
}

# Estimated size of sph/ after the reduction: the multi-radius files shrink to
# NKEEP of their nradii shells (read from the first file's header).
sph_estimate() {
    local vtk=("$1"/*.vtk) all vtk_size n
    all=$(size_of "$1")
    (( ${#vtk[@]} )) || { echo "${all}"; return; }
    vtk_size=$(find "$1" -maxdepth 1 -name '*.vtk' -printf '%s\n' | awk '{s += $1} END {print s + 0}')
    n=$(sed -n '2s/.*nradii=\([0-9]*\).*/\1/p' "${vtk[0]}" | head -c 9)
    if [[ -n "$n" && "$n" -gt "${NKEEP}" ]]; then
        echo $(( all - vtk_size + vtk_size * NKEEP / n ))
    else
        echo "${all}"
    fi
}

keep() {      # SRC_PATH DEST_DIR
    local s; s=$(size_of "$1")
    report keep "$s" "$s" "$1"
    (( DRY )) && return
    mkdir -p "$2"
    cp -a -u "$1" "$2/"
}

drop() {
    local s; s=$(size_of "$1")
    report drop "$s" 0 "$1"
}

reduce_sph() {   # SPH_DIR DEST_DIR
    local out="$2/sph" s
    if (( DRY )); then
        s=$(size_of "$1")
        report reduce "$s" "$(sph_estimate "$1")" "$1"
        return
    fi
    report reduce "$(size_of "$1")" "" "$1"
    mkdir -p "${out}"
    find "$1" -maxdepth 1 -type f ! -name '*.vtk' -exec cp -a -u {} "${out}/" \;
    if ! find "$1" -maxdepth 1 -type f -name '*.vtk' -print0 | sort -z |
            xargs -0 -r -n 50 -P "${JOBS}" python3 -B "${EXTRACT}" \
                --radii ${RADII} --out-dir "${out}" >> "${LOG}"; then
        echo "ERROR: sph extraction failed for $1 (see ERROR lines in ${LOG})"
        exit 1
    fi
}

compress_dir() {   # DIR DEST_DIR
    local out="$2/${1##*/}" s
    s=$(size_of "$1")
    if (( DRY )); then
        # text histograms, mostly zeros: gzip shrinks them about 20x
        report compress "$s" $(( s / 20 )) "$1"
        return
    fi
    report compress "$s" "" "$1"
    mkdir -p "${out}"
    find "$1" -mindepth 1 -maxdepth 1 ! -type f -exec echo "WARNING: not archived: {}" \;
    find "$1" -maxdepth 1 -type f -print0 |
        xargs -0 -r -n 100 -P "${JOBS}" sh -c '
            out=$1; shift
            for f; do
                o="${out}/${f##*/}.gz"
                [ -e "$o" ] && [ ! "$f" -nt "$o" ] && continue
                gzip -c -n "$f" > "$o.tmp" && touch -r "$f" "$o.tmp" && mv "$o.tmp" "$o"
            done' sh "${out}"
}

segment_entry() {   # PATH SEGMENT_DEST_DIR
    local name="${1##*/}"
    if [[ -d "$1" && " ${DROP[*]} " == *" ${name} "* ]]; then
        drop "$1"
    elif [[ -d "$1" && "${name}" == sph ]]; then
        reduce_sph "$1" "$2"
    elif [[ -d "$1" && "${name}" == pdf_* ]]; then
        compress_dir "$1" "$2"
    else
        keep "$1" "$2"
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
echo "============================================================"
echo "  Archive $( (( DRY )) && echo '(dry run) ')${SRC}"
echo "       -> ${DEST}"
echo "  sph shells kept: ${RADII} M_sun"
echo "============================================================"
if (( ! DRY )); then
    mkdir -p "${DEST}"
    {
        echo "# ---- $(date -Iseconds)  $(whoami)@$(hostname)"
        echo "# source : ${SRC}"
        echo "# radii  : ${RADII}"
        echo "# script : $(git -C "${SCRIPT_DIR}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    } >> "${LOG}"
fi

SEGMENTS=("${SRC}"/output-*/)
for entry in "${SRC}"/*; do
    name="${entry##*/}"
    if (( ${#SEGMENTS[@]} )); then
        if [[ -d "${entry}" && "${name}" == output-* ]]; then
            for item in "${entry}"/*; do
                segment_entry "${item}" "${DEST}/${name}"
            done
        else
            keep "${entry}" "${DEST}"
        fi
    elif [[ -d "${entry}" && "${name}" =~ ${RUN_LEVEL_RE} ]]; then
        keep "${entry}" "${DEST}"
    else
        segment_entry "${entry}" "${DEST}/output-0000"
    fi
done

echo "------------------------------------------------------------"
if (( DRY )); then
    echo "  total: $(human "${TOTAL_IN}") -> about $(human "${TOTAL_OUT}")  (dry run, nothing copied)"
else
    echo "  done: $(human "$(size_of "${DEST}")") in ${DEST}"
    echo "# done : $(date -Iseconds)" >> "${LOG}"
fi
