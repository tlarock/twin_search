#!/usr/bin/env bash
# Exact Gram-mate pairs for each exhaustive projection file.
#
# This is the dataset behind paper Figure A.9. It is a POST-PROCESSING pass:
# it reads the CSVs written by exhaustive_nonuniform.sh and, for every pair of
# twins sharing a projection, tests whether they are exact Gram mates - equal
# projection matrix AND equal line-graph matrix, rather than merely isomorphic.
# For each input FILE it writes FILE.exact_mates_pairs alongside it.
#
#   ./scripts/exhaustive_nonuniform.sh          # produce the inputs first
#   ./scripts/exact_mates.sh [input_dir]
#
# The comparison is O(T^2) in the number of twins on a row, so cost is driven
# by the largest rows rather than by the number of rows. Guarded per file.
#
# Env: MAX_THREADS (6), MAX_SECONDS (3600), MAX_RSS_MB (12000), DRY_RUN=1

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

INPUT_DIR="${1:-$REPO/results/reproduce/exhaustive}"
MAX_SECONDS="${MAX_SECONDS:-3600}"
MAX_RSS_MB="${MAX_RSS_MB:-12000}"

BIN="$(need_binary exact_mates)" || exit 1

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "error: no such directory: $INPUT_DIR" >&2
    echo "  Run ./scripts/exhaustive_nonuniform.sh first." >&2
    exit 1
fi

shopt -s nullglob
inputs=( "$INPUT_DIR"/*_exhaustive_projections.csv )
if (( ${#inputs[@]} == 0 )); then
    echo "error: no *_exhaustive_projections.csv in $INPUT_DIR" >&2
    exit 1
fi

start_manifest "$INPUT_DIR/RUN_MANIFEST_exact_mates.tsv" "$INPUT_DIR"
echo "exact mates over ${#inputs[@]} file(s) in $INPUT_DIR"
echo

total_start=$(date +%s)
for in_file in "${inputs[@]}"; do
    out="$in_file.exact_mates_pairs"
    # NOTE: --output-path is the INPUT file for this driver, despite the name;
    # it writes <that path>.exact_mates_pairs.
    cmd=( "$BIN" --output-path "$in_file" --max-threads "$MAX_THREADS" )

    if [[ -s "$out" ]]; then
        echo "  $(basename "$in_file")  skipped (exists)"
        log_cell skipped 0 0 "$(basename "$out")" "${cmd[*]}"
        continue
    fi
    if [[ "$DRY_RUN" == 1 ]]; then echo "  DRY RUN: ${cmd[*]}"; continue; fi

    printf '  %-56s ' "$(basename "$in_file")"
    guard_output "$out"
    if run_guarded "$MAX_SECONDS" "$MAX_RSS_MB" "${cmd[@]}"; then
        printf '%8s  %5sMB  %s pair-rows\n' \
            "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$(wc -l < "$out" | tr -d ' ')"
        log_cell ok "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
    else
        # This is an exhaustive ENUMERATION, not a sample, so the two cases
        # differ: n rows of a sample are n valid independent draws, but n rows
        # of a truncated enumeration are a prefix, and nothing in the file says
        # how many rows there should have been. Keeping it under the real name
        # would let a later run skip it as finished. So it is preserved under
        # .partial - nothing is destroyed, and the cell correctly re-runs.
        trim_torn_line "$out"
        rows=$(rows_in "$out")
        if (( rows > 0 )); then
            mv "$out" "$out.partial-${rows}rows"
            printf '%8s  %5sMB  %s -- kept %s rows as .partial-%srows\n' \
                "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS" "$rows" "$rows"
        else
            rm -f "$out"
            printf '%8s  %5sMB  %s -- no rows written\n' \
                "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS"
        fi
        log_cell "$RG_STATUS" "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
    fi
    clear_output
done

echo
echo "total $(human $(( $(date +%s) - total_start )))"
