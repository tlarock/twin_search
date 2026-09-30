#!/usr/bin/env bash
# Shared helpers for the reproduction scripts in this directory.
#
# Source this, do not run it:
#     source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
#
# What it provides:
#   REPO, BUILD, MAX_THREADS, SEED, DRY_RUN
#   need_binary <name>          - resolve and check a driver exists
#   start_manifest <file>       - begin a provenance record for this run
#   log_cell <fields...>        - append one TSV row to the manifest
#   run_guarded <secs> <mb> <cmd...>
#                               - run with wall-clock and peak-RSS limits,
#                                 killing the job rather than the machine
#   nCk <n> <k>                 - binomial coefficient
#   human <seconds>             - pretty-print a duration

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${TWIN_SEARCH_BUILD:-$REPO/build}"

# Cap parallelism by default. These drivers will happily take every core.
MAX_THREADS="${MAX_THREADS:-6}"

# Always seed. The published sampled data predates --seed and so cannot be
# reproduced; anything generated here should not repeat that mistake.
SEED="${SEED:-20260930}"

DRY_RUN="${DRY_RUN:-0}"

need_binary() {
    local name="$1"
    if [[ ! -x "$BUILD/$name" ]]; then
        echo "error: $BUILD/$name not found or not executable." >&2
        echo "  Build first, e.g.:  cmake --preset conan-release && cmake --build build -j $MAX_THREADS" >&2
        exit 1
    fi
    printf '%s\n' "$BUILD/$name"
}

nCk() {
    local n="$1" k="$2" num=1 den=1 i
    (( k > n - k )) && k=$(( n - k ))
    for (( i = 1; i <= k; ++i )); do
        num=$(( num * (n + 1 - i) ))
        den=$(( den * i ))
    done
    echo $(( num / den ))
}

human() {
    local s="${1%.*}"
    if (( s < 60 )); then printf '%ss' "$s"
    elif (( s < 3600 )); then printf '%dm%02ds' $(( s / 60 )) $(( s % 60 ))
    else printf '%dh%02dm' $(( s / 3600 )) $(( (s % 3600) / 60 )); fi
}

MANIFEST=""
start_manifest() {
    MANIFEST="$1"
    mkdir -p "$(dirname "$MANIFEST")"
    {
        echo "# twin_search reproduction run"
        echo "# date_utc      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "# git_commit    $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
        echo "# git_describe  $(git -C "$REPO" describe --always --dirty 2>/dev/null || echo unknown)"
        echo "# git_dirty     $(git -C "$REPO" status --porcelain 2>/dev/null | grep -q . && echo yes || echo no)"
        echo "# host          $(uname -srm)"
        echo "# max_threads   $MAX_THREADS"
        echo "# seed          $SEED"
        echo "# script        $(basename "${BASH_SOURCE[1]:-?}")"
        echo "# argv          $*"
        printf '#\n'
        printf 'status\tseconds\tpeak_rss_mb\toutput\tcommand\n'
    } > "$MANIFEST"
}

log_cell() {
    [[ -n "$MANIFEST" ]] && printf '%s\t%s\t%s\t%s\t%s\n' "$@" >> "$MANIFEST"
}

# run_guarded <max_seconds> <max_rss_mb> <command...>
#
# Runs the command in the background and polls it once a second, killing it if
# it exceeds either limit. macOS has no timeout(1) and no working `ulimit -v`,
# so this is done by hand. Sets RG_STATUS (ok|timeout|oom|failed), RG_SECONDS
# and RG_PEAK_MB.
run_guarded() {
    local max_s="$1" max_mb="$2"; shift 2
    local start peak=0 rss pid rc

    start=$(date +%s)
    "$@" >/dev/null 2>&1 &
    pid=$!

    RG_STATUS=ok
    while kill -0 "$pid" 2>/dev/null; do
        rss=$(ps -o rss= -p "$pid" 2>/dev/null | tr -d ' ')
        rss=${rss:-0}
        (( rss > peak )) && peak=$rss
        if (( max_mb > 0 && rss / 1024 > max_mb )); then
            kill -9 "$pid" 2>/dev/null; RG_STATUS=oom; break
        fi
        if (( max_s > 0 && $(date +%s) - start > max_s )); then
            kill -9 "$pid" 2>/dev/null; RG_STATUS=timeout; break
        fi
        sleep 1
    done
    wait "$pid" 2>/dev/null; rc=$?
    [[ "$RG_STATUS" == ok && $rc -ne 0 ]] && RG_STATUS=failed

    RG_SECONDS=$(( $(date +%s) - start ))
    RG_PEAK_MB=$(( peak / 1024 ))
    [[ "$RG_STATUS" == ok ]]
}
