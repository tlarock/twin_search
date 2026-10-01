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
#   guard_output <file>         - name the file the current cell is writing, so
#                                 an interrupt can tidy it rather than leaving
#                                 a torn row behind
#   rows_in <file>              - completed samples in a results file, 0 if absent
#   trim_torn_line <file>       - drop a final row cut off mid-write
#   clear_output                - the cell finished; nothing to clean up
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

# Begin a provenance block for this run.
#
# APPENDS rather than truncates. Runs are resumable, so a directory is often
# filled by several invocations; truncating would throw away the timings of
# every cell an earlier run completed and replace them with "skipped 0s". The
# file is therefore a log of runs, each introduced by its own "# ---" header.
start_manifest() {
    MANIFEST="$1"
    shift
    mkdir -p "$(dirname "$MANIFEST")"
    # Decide this BEFORE the block below starts writing: its first echo makes
    # the file non-empty, so testing -s inside the block always says "resuming".
    local resuming=0
    [[ -s "$MANIFEST" ]] && resuming=1
    {
        (( resuming )) && printf '\n'
        echo "# --------------------------------------------------------------"
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
        # Column header only on a fresh file; a resumed run just adds rows.
        (( resuming )) || printf 'status\tseconds\tpeak_rss_mb\toutput\tcommand\n'
    } >> "$MANIFEST"
}

log_cell() {
    [[ -n "$MANIFEST" ]] && printf '%s\t%s\t%s\t%s\t%s\n' "$@" >> "$MANIFEST"
}

# The output file the current cell is writing, if any.
#
# These drivers append as they go, so a run killed from outside - Ctrl-C, or a
# `pkill` aimed at the script - leaves a partial file behind. Because cells are
# skipped when their output already exists, that partial file would be silently
# accepted as finished by the next run. So an interrupt must delete it. The
# in-cell guards in run_guarded handle the limits they enforce themselves; this
# trap covers everything else.
CURRENT_OUTPUT=""
RG_PID=""

# Completed samples in a results file. The drivers write exactly one row per
# sample, so this is the sample count - which is what "is this cell done?" must
# ask. Testing -s only asks "did anything get written", which silently accepts
# a cell that stopped after one sample.
rows_in() {
    [[ -s "$1" ]] || { echo 0; return; }
    wc -l < "$1" | tr -d ' '
}

# Drop a final row that was cut off mid-write.
#
# Each row is written under a lock and terminated by std::endl, which flushes,
# so a file ending in a newline has only complete rows. If it does not end in
# one, the process died partway through writing the last row; everything before
# the last newline is intact. A file left with nothing at all is removed, so
# that -s and rows_in agree it is absent.
trim_torn_line() {
    local f="$1"
    [[ -s "$f" ]] || return 0
    if [[ -n "$(tail -c 1 "$f")" ]]; then
        sed '$d' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
    fi
    [[ -s "$f" ]] || rm -f "$f"
}

guard_output() { CURRENT_OUTPUT="$1"; }
clear_output()  { CURRENT_OUTPUT=""; }

_on_interrupt() {
    trap - INT TERM
    [[ -n "$RG_PID" ]] && kill -9 "$RG_PID" 2>/dev/null
    # KEEP what the cell managed. Every completed sample is already on disk;
    # only the row in flight is lost. Deleting the file threw away real work -
    # on 2026-09-30 a TERM at a harness time limit destroyed ~850 samples of
    # k=4 n=10 m=16 that had taken 70 minutes. Callers decide what to do with a
    # short cell by asking rows_in, not by testing existence.
    local kept=0
    if [[ -n "$CURRENT_OUTPUT" && -e "$CURRENT_OUTPUT" ]]; then
        trim_torn_line "$CURRENT_OUTPUT"
        kept=$(rows_in "$CURRENT_OUTPUT")
        echo >&2
        if (( kept > 0 )); then
            echo "interrupted: keeping $kept completed samples in $CURRENT_OUTPUT" >&2
        else
            echo "interrupted: no completed samples, removed $CURRENT_OUTPUT" >&2
        fi
    fi
    [[ -n "${MANIFEST:-}" ]] && printf 'interrupted\t0\t0\t%s\tkept=%s rows\n' \
        "$(basename "${CURRENT_OUTPUT:-none}")" "$kept" >> "$MANIFEST"
    exit 130
}
trap _on_interrupt INT TERM

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
    RG_PID=$pid

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

    RG_PID=""
    RG_SECONDS=$(( $(date +%s) - start ))
    RG_PEAK_MB=$(( peak / 1024 ))
    [[ "$RG_STATUS" == ok ]]
}
