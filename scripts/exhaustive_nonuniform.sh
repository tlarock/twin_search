#!/usr/bin/env bash
# Exhaustive non-uniform twin search over every n=6, k=3 projection.
#
# This is the dataset behind paper Figures 2-6 and A.9: for each m, enumerate
# the projections of every 3-uniform hypergraph on 6 nodes up to isomorphism,
# then find every hypergraph sharing each projection, allowing hyperedges of
# size 2..n in the search (min_k=2, which is what makes it "non-uniform").
#
#   ./scripts/exhaustive_nonuniform.sh [m_min] [m_max]
#
# Defaults to m=2..13. Cost is NOT monotone in m: the number of distinct
# projections peaks at m=10 (342) and falls away, while the search tree per
# projection keeps growing. Measured on this machine, 6 threads:
#
#     m=2..9   under 8s each, 13s total,   under 0.2 GB
#     m=10     94s     342 projections     1.2 GB
#     m=11     39s     304 projections     1.2 GB
#     m=12     71s     245 projections     1.9 GB
#     m=13    124s     159 projections     4.2 GB
#     m=14    395s      94 projections     8.7 GB
#
# Peak RSS roughly doubles per step above m=12, so m=15 is about 17 GB and
# m=16 about 34 GB. m>=15 is therefore opt-in: ask for it explicitly and watch
# the machine. Each cell is guarded and skipped rather than allowed to swap.
#
# Env: MAX_THREADS (6), MAX_SECONDS (1800), MAX_RSS_MB (12000),
#      OUTPUT_DIR, TWIN_SEARCH_BUILD, DRY_RUN=1

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

M_MIN="${1:-2}"
M_MAX="${2:-13}"
N=6
K=3
MIN_K=2
MAX_SECONDS="${MAX_SECONDS:-1800}"
MAX_RSS_MB="${MAX_RSS_MB:-12000}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO/results/reproduce/exhaustive}"

BIN="$(need_binary exhaustive_search_projections)" || exit 1
mkdir -p "$OUTPUT_DIR"
start_manifest "$OUTPUT_DIR/RUN_MANIFEST.tsv" "m=$M_MIN..$M_MAX"

echo "exhaustive n=$N k=$K min_k=$MIN_K, m=$M_MIN..$M_MAX -> $OUTPUT_DIR"
echo "limits per cell: ${MAX_SECONDS}s, ${MAX_RSS_MB}MB, threads=$MAX_THREADS"
echo

total_start=$(date +%s)
for (( m = M_MIN; m <= M_MAX; m++ )); do
    out="$OUTPUT_DIR/n-${N}_m-${m}_k-${K}_non-uniform_exhaustive_projections.csv"
    cmd=( "$BIN" -n "$N" -m "$m" -k "$K" --min-k "$MIN_K"
          --max-threads "$MAX_THREADS" --output-path "$OUTPUT_DIR/" )

    if [[ -s "$out" ]]; then
        echo "  m=$m  skipped (exists: $(basename "$out"))"
        log_cell skipped 0 0 "$(basename "$out")" "${cmd[*]}"
        continue
    fi
    if [[ "$DRY_RUN" == 1 ]]; then
        echo "  m=$m  DRY RUN: ${cmd[*]}"
        continue
    fi

    printf '  m=%-3s ' "$m"
    # The driver appends, so a killed cell must not leave a partial file behind
    # for the next run to mistake for a finished one.
    guard_output "$out"
    if run_guarded "$MAX_SECONDS" "$MAX_RSS_MB" "${cmd[@]}"; then
        printf '%8s  %5sMB  %s projections\n' \
            "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$(wc -l < "$out" | tr -d ' ')"
        log_cell ok "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
    else
        printf '%8s  %5sMB  %s -- removing partial output\n' \
            "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS"
        rm -f "$out"
        log_cell "$RG_STATUS" "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
    fi
    clear_output
done

echo
echo "total $(human $(( $(date +%s) - total_start )))   manifest: $OUTPUT_DIR/RUN_MANIFEST.tsv"
