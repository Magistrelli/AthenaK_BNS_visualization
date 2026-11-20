#!/usr/bin/env bash
# load_config.sh — shared machine-config loader for the BNS visualization scripts.
# -----------------------------------------------------------------------------
# SOURCE this file (don't execute it). After sourcing you get:
#   * _cfg_get <key>   — read a value from config.ini
#   * _cfg_sims        — the simulations listed in config.ini, one "path|label"
#                        line each (one per `simpath` line)
#   * _for_each_sim    — rerun the calling script once per simulation
#   * the machine environment set up (Lmod made available + the config's
#     `env_setup` command run, e.g. `module load python`)
#
# The point: ALL machine-specific settings live in ONE file, config.ini —
# both the paths (simpath/figpath/eos_table/pythonpath_extra) AND the
# environment activation (env_setup). Move to a new machine → edit config.ini
# only; none of the *.sh scripts need touching.
#
# config.ini itself stays a plain INI file (parsed, not executed); only the
# single `env_setup` value is run, as your explicit environment command.

# Directory of THIS file = repo root, where config.ini lives beside it.
_LC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${_LC_DIR}/config.ini}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: ${CONFIG_FILE} not found."
    echo "       Copy config.example.ini to config.ini and edit it for your machine."
    exit 1
fi

# Values of one key in config.ini: the first match, or all of them with "all"
# as second argument. Comments (#/;) and [section] headers are ignored.
_cfg_values() {
    awk -v key="$1" -v all="${2:-}" '
        /^[[:space:]]*[#;]/ { next }              # comment line
        /^[[:space:]]*\[/   { next }              # section header
        {
            k = $0; sub(/=.*/, "", k)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", k)
            if (k == key && index($0, "=")) {
                v = substr($0, index($0, "=") + 1)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
                print v
                if (all != "all") exit
            }
        }
    ' "${CONFIG_FILE}"
}

# The simulations: one "path|label" line per `simpath` line of config.ini
# ("simpath = /path/to/run | label"); the label defaults to the folder name.
_cfg_sims() {
    local v path label
    while IFS= read -r v; do
        path="${v%%|*}"; label=""
        [[ "${v}" == *"|"* ]] && label="${v#*|}"
        path="${path%"${path##*[![:space:]]}"}"
        label="${label#"${label%%[![:space:]]*}"}"
        [[ -n "${path}" ]] && echo "${path}|${label:-$(basename "${path}")}"
    done < <(_cfg_values simpath all)
}

# Read one key from config.ini; the first "key = value" match wins.  For
# `simpath` this is the path of the simulation being processed (see
# _for_each_sim), or else of the first one listed.
_cfg_get() {
    if [[ "$1" == simpath ]]; then
        local sims
        mapfile -t sims < <(_cfg_sims)
        local sim="${_SIM_CURRENT:-${sims[0]:-}}"
        echo "${sim%%|*}"
    else
        _cfg_values "$1"
    fi
}

# Settings that belong to one simulation cannot be shared by several.
_cfg_check_multi_sim() {
    local key bad=0
    for key in figpath t_merger; do
        if [[ -n "$(_cfg_get "${key}")" ]]; then
            echo "ERROR: ${CONFIG_FILE} lists several simulations but sets" \
                 "'${key}', which belongs to a single one; leave it blank."
            bad=1
        fi
    done
    if [[ -n "${T_MERGER:-}" ]]; then
        echo "ERROR: T_MERGER is set, but ${CONFIG_FILE} lists several simulations."
        bad=1
    fi
    return "${bad}"
}

# For scripts that work on one simulation: when config.ini lists several, run
# the calling script once per simulation (one after the other, passing on its
# arguments), report the ones that failed, and exit.  With a single simulation
# it just returns.  Call it right after sourcing this file.
_for_each_sim() {
    [[ -n "${_SIM_CURRENT:-}" ]] && return 0
    local sims sim n=0 failed=()
    mapfile -t sims < <(_cfg_sims)
    (( ${#sims[@]} > 1 )) || return 0
    _cfg_check_multi_sim || exit 1
    for sim in "${sims[@]}"; do
        n=$(( n + 1 ))
        echo "################ simulation ${n}/${#sims[@]}: ${sim#*|}  (${sim%%|*})"
        _BNS_ENV_READY=1 _SIM_CURRENT="${sim}" bash "${BASH_SOURCE[1]}" "$@" \
            || failed+=("${sim#*|}")
    done
    if (( ${#failed[@]} )); then
        echo "FAILED for ${#failed[@]} of ${#sims[@]} simulations: ${failed[*]}"
        exit 1
    fi
    exit 0
}

# --- Machine environment setup ------------------------------------------------
# (Skipped in the per-simulation runs started by _for_each_sim / run_all.sh,
# which inherit the environment already set up by their parent.)
if [[ -z "${_BNS_ENV_READY:-}" ]]; then
    # Lmod's `module` is usually only a shell function in interactive/login shells,
    # so it is NOT available in a plain `bash script.sh`. Source Lmod's init first
    # so `module` (and thus env_setup) works in non-interactive shells too.
    if ! command -v module &>/dev/null; then
        for _init in "${MODULESHOME:-/usr/share/lmod/lmod}/init/bash" \
                     /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh \
                     /etc/profile.d/00-modulepath.sh; do
            [[ -r "${_init}" ]] && source "${_init}" && break
        done
    fi

    # Run the machine-specific environment command(s) from config.ini's env_setup.
    # (`|| ...` keeps a harmless module warning from aborting `set -e` scripts.)
    _ENV_SETUP="$(_cfg_get env_setup)"
    if [[ -n "${_ENV_SETUP}" ]]; then
        eval "${_ENV_SETUP}" \
            || echo "WARNING: env_setup returned non-zero: ${_ENV_SETUP}" >&2
    else
        echo "WARNING: 'env_setup' is empty in ${CONFIG_FILE};" \
             "relying on the current environment." >&2
    fi
fi
