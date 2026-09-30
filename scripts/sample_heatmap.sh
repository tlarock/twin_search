#!/usr/bin/env bash
# Sampled twin counts over the (n, m) grid, for both k=3 and k=4.
#
# This is the dataset behind the paper Figure 7 HEATMAP, whose cell value is
#
#     p = #{samples with total_filtered > 1} / num_samples
#
# For each cell, draw --samples k-uniform hypergraphs on n nodes with m
# hyperedges and count the twins of each one's projection. The search is
# k-uniform too (min_k = max_k = k), which the filenames record by carrying no
# min-k suffix.
#
#   ./scripts/sample_heatmap.sh [k ...]        # default: 3 4
#
# Unlike the published files, every run here is seeded, so the output is a
# deterministic function of (seed, sample index) and is reproducible
# independently of thread count. Line order within a file still varies under
# parallel writes - sort before diffing.
#
# Two cells cannot exist and are skipped with an explanation: k=4 n=6 m=15 and
# m=16. C(6,4) = 15, so m=16 asks for more distinct 4-subsets than exist, and
# m=15 is the complete 4-uniform hypergraph - a single hypergraph, so 1000
# draws give an effective sample size of 1 and p is identically 0.
#
# Env: N_MIN/N_MAX/M_MIN/M_MAX (6/16/6/16), SAMPLES (1000), SEED,
#      MAX_THREADS (6), MAX_SECONDS (3600), MAX_RSS_MB (12000),
#      OUTPUT_DIR, DRY_RUN=1

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

KS=( "${@:-3 4}" )
[[ $# -eq 0 ]] && KS=( 3 4 )
N_MIN="${N_MIN:-6}";  N_MAX="${N_MAX:-16}"
M_MIN="${M_MIN:-6}";  M_MAX="${M_MAX:-16}"
SAMPLES="${SAMPLES:-1000}"
MAX_SECONDS="${MAX_SECONDS:-3600}"
MAX_RSS_MB="${MAX_RSS_MB:-12000}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO/results/reproduce/heatmap-$SAMPLES}"

BIN="$(need_binary count_twins_random)" || exit 1
mkdir -p "$OUTPUT_DIR"
start_manifest "$OUTPUT_DIR/RUN_MANIFEST.tsv" "k=${KS[*]} n=$N_MIN..$N_MAX m=$M_MIN..$M_MAX samples=$SAMPLES"

echo "heatmap grid k=${KS[*]}, n=$N_MIN..$N_MAX, m=$M_MIN..$M_MAX, $SAMPLES samples, seed=$SEED"
echo "  -> $OUTPUT_DIR"
echo

total_start=$(date +%s)
for k in "${KS[@]}"; do
  for (( n = N_MIN; n <= N_MAX; n++ )); do
    max_m=$(nCk "$n" "$k")
    for (( m = M_MIN; m <= M_MAX; m++ )); do
        out="$OUTPUT_DIR/n-${n}_m-${m}_k-${k}_samples-${SAMPLES}.csv"

        if (( m > max_m )); then
            echo "  k=$k n=$n m=$m  impossible: C($n,$k)=$max_m < m"
            log_cell impossible 0 0 "$(basename "$out")" "-"
            continue
        fi
        if (( m == max_m )); then
            echo "  k=$k n=$n m=$m  degenerate: m == C($n,$k), the complete k-uniform hypergraph (one sample, p=0)"
            log_cell degenerate 0 0 "$(basename "$out")" "-"
            continue
        fi

        cmd=( "$BIN" -n "$n" -m "$m" -k "$k" --samples "$SAMPLES"
              --min-k "$k" --max-k "$k" --seed "$SEED"
              --max-threads "$MAX_THREADS" --output-path "$OUTPUT_DIR/" )

        if [[ -s "$out" ]]; then
            log_cell skipped 0 0 "$(basename "$out")" "${cmd[*]}"
            continue
        fi
        if [[ "$DRY_RUN" == 1 ]]; then echo "  DRY RUN: ${cmd[*]}"; continue; fi

        printf '  k=%s n=%-3s m=%-3s ' "$k" "$n" "$m"
        if run_guarded "$MAX_SECONDS" "$MAX_RSS_MB" "${cmd[@]}"; then
            printf '%8s  %5sMB\n' "$(human "$RG_SECONDS")" "$RG_PEAK_MB"
            log_cell ok "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
        else
            printf '%8s  %5sMB  %s -- removing partial output\n' \
                "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS"
            rm -f "$out"
            log_cell "$RG_STATUS" "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
        fi
    done
  done
done

echo
echo "total $(human $(( $(date +%s) - total_start )))   manifest: $OUTPUT_DIR/RUN_MANIFEST.tsv"
