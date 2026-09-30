#!/usr/bin/env bash
# Sampled twin counts across the full density range, for k=3 and n=6..9.
#
# This is the dataset behind the paper Figure 7 LINE PANEL and Figure A.10,
# which read a different set from the heatmap: fewer nodes, but m swept across
# the whole feasible range rather than a fixed 6..16 window. The published set
# used 500 samples per cell.
#
#   ./scripts/sample_density.sh [n ...]        # default: 6 7 8 9
#
# COST WARNING. Unlike the heatmap, this gets genuinely expensive. Cost climbs
# steeply with m because the search tree grows with the number of hyperedges.
# Measured on this machine, 6 threads, 500 samples:
#
#     n=9 m=20    7s
#     n=9 m=40    over 10 minutes (not run to completion)
#
# So the upper m range for n=8 and n=9 is where the time goes. Every cell is
# guarded by MAX_SECONDS and skipped, with the skip recorded in the manifest,
# rather than being allowed to run for hours. Raise MAX_SECONDS deliberately if
# you want the dense tail; M_MAX lets you stop short of it instead.
#
# m runs 1..C(n,3)-1. m = C(n,3) is excluded because it is the complete
# 3-uniform hypergraph: a single hypergraph, so sampling it teaches nothing.
#
# Env: SAMPLES (500), M_MIN (1), M_MAX (0 = C(n,3)-1), SEED, MAX_THREADS (6),
#      MAX_SECONDS (300), MAX_RSS_MB (12000), OUTPUT_DIR, DRY_RUN=1

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

NS=( "$@" )
[[ $# -eq 0 ]] && NS=( 6 7 8 9 )
K=3
SAMPLES="${SAMPLES:-500}"
M_MIN="${M_MIN:-1}"
M_MAX="${M_MAX:-0}"
MAX_SECONDS="${MAX_SECONDS:-300}"
MAX_RSS_MB="${MAX_RSS_MB:-12000}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO/results/reproduce/density-$SAMPLES}"

BIN="$(need_binary count_twins_random)" || exit 1
mkdir -p "$OUTPUT_DIR"
start_manifest "$OUTPUT_DIR/RUN_MANIFEST.tsv" "n=${NS[*]} samples=$SAMPLES m_max=$M_MAX"

echo "density sweep k=$K, n=${NS[*]}, $SAMPLES samples, seed=$SEED"
echo "  per-cell limit ${MAX_SECONDS}s / ${MAX_RSS_MB}MB -- cells over budget are SKIPPED, not waited on"
echo "  -> $OUTPUT_DIR"
echo

total_start=$(date +%s)
for n in "${NS[@]}"; do
    top=$(( $(nCk "$n" "$K") - 1 ))
    (( M_MAX > 0 && M_MAX < top )) && top=$M_MAX
    echo "  n=$n: m=$M_MIN..$top"
    for (( m = M_MIN; m <= top; m++ )); do
        out="$OUTPUT_DIR/n-${n}_m-${m}_k-${K}_samples-${SAMPLES}.csv"
        cmd=( "$BIN" -n "$n" -m "$m" -k "$K" --samples "$SAMPLES"
              --min-k "$K" --max-k "$K" --seed "$SEED"
              --max-threads "$MAX_THREADS" --output-path "$OUTPUT_DIR/" )

        if [[ -s "$out" ]]; then
            log_cell skipped 0 0 "$(basename "$out")" "${cmd[*]}"
            continue
        fi
        if [[ "$DRY_RUN" == 1 ]]; then echo "    DRY RUN: ${cmd[*]}"; continue; fi

        printf '    m=%-3s ' "$m"
        guard_output "$out"
        if run_guarded "$MAX_SECONDS" "$MAX_RSS_MB" "${cmd[@]}"; then
            printf '%8s  %5sMB\n' "$(human "$RG_SECONDS")" "$RG_PEAK_MB"
            log_cell ok "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
        else
            printf '%8s  %5sMB  %s -- removing partial output\n' \
                "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS"
            rm -f "$out"
            log_cell "$RG_STATUS" "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
            # Cost is monotone in m within an n, so once one cell blows the
            # budget the rest of this n will too. Stop instead of grinding.
            if [[ "$RG_STATUS" == timeout || "$RG_STATUS" == oom ]]; then
                echo "    -- $RG_STATUS at m=$m; skipping the rest of n=$n (cost rises with m)"
                for (( mm = m + 1; mm <= top; mm++ )); do
                    log_cell "not-attempted" 0 0 \
                        "n-${n}_m-${mm}_k-${K}_samples-${SAMPLES}.csv" "-"
                done
                break
            fi
        fi
        clear_output
    done
done

echo
echo "total $(human $(( $(date +%s) - total_start )))   manifest: $OUTPUT_DIR/RUN_MANIFEST.tsv"
