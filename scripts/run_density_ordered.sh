#!/usr/bin/env bash
# Run the Figure 7 line-panel sweep in estimated-cost order, cheapest first.
#
# sample_density.sh walks m in numerical order for one k. That is the wrong
# order for a bounded run: cost over this sweep spans seven orders of
# magnitude and peaks in a band near m = C(n,k)/2, so walking m in order sinks
# the whole budget into the band and never reaches the cheap dense tail. It
# also hardcodes k=3, while the figure needs k=3 for n=6..9 AND k=4 for
# n=6..8 - three of the seven curves it cannot produce at all.
#
# This instead consumes the ordered list from estimate_density_cost.py, so a
# run that is stopped at any point has completed the most cells it could have.
#
#   scripts/estimate_density_cost.py <published_dir>   # EMIT_LIST=order.tsv
#   scripts/run_density_ordered.sh order.tsv
#
# Three bounds, because an estimate is not a promise:
#   EST_LIMIT           do not attempt a cell estimated above this. Cells in
#                       the expensive band are 10-1000x past any reachable
#                       budget, so attempting one only burns its cap and
#                       deletes the output.
#   CAP_FACTOR          per-cell wall limit, as a multiple of the estimate.
#                       The calibration spread is 1.36-4.37x and our draw is
#                       not their draw, so the cap must be loose or it kills
#                       cells that would have finished.
#   TOTAL_BUDGET_HOURS  stop STARTING cells after this. A cell already running
#                       is left to finish; interrupting it would delete it.
#
# Env: SAMPLES(500) SEED MAX_THREADS(6) MAX_RSS_MB(12000) OUTPUT_DIR DRY_RUN
#      EST_LIMIT(7200) CAP_FACTOR(3) CAP_MIN(1200) CAP_MAX(10800)
#      TOTAL_BUDGET_HOURS(11)

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

LIST="${1:?usage: run_density_ordered.sh <ordered_list.tsv>}"
[[ -r "$LIST" ]] || { echo "error: cannot read $LIST" >&2; exit 1; }

SAMPLES="${SAMPLES:-500}"
EST_LIMIT="${EST_LIMIT:-7200}"
CAP_FACTOR="${CAP_FACTOR:-3}"
CAP_MIN="${CAP_MIN:-1200}"
CAP_MAX="${CAP_MAX:-10800}"
TOTAL_BUDGET_HOURS="${TOTAL_BUDGET_HOURS:-11}"
MAX_RSS_MB="${MAX_RSS_MB:-12000}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO/results/reproduce/density-$SAMPLES}"

BIN="$(need_binary count_twins_random)" || exit 1
mkdir -p "$OUTPUT_DIR"
start_manifest "$OUTPUT_DIR/RUN_MANIFEST.tsv" \
    "list=$(basename "$LIST") samples=$SAMPLES est_limit=$EST_LIMIT budget=${TOTAL_BUDGET_HOURS}h threads=$MAX_THREADS"

budget_s=$(awk -v h="$TOTAL_BUDGET_HOURS" 'BEGIN{printf "%d", h*3600}')
start=$(date +%s)

echo "density sweep in estimated-cost order, $SAMPLES samples, seed=$SEED, ${MAX_THREADS} threads"
echo "  list       $LIST"
echo "  est limit  ${EST_LIMIT}s per cell   budget ${TOTAL_BUDGET_HOURS}h total"
echo "  -> $OUTPUT_DIR"
echo

n_ok=0; n_skip=0; n_over=0; n_fail=0
while IFS=$'\t' read -r k n m est cum note <&3; do
    [[ -z "${k:-}" || "$k" == \#* || "$k" == "k" ]] && continue

    elapsed=$(( $(date +%s) - start ))
    if (( elapsed > budget_s )); then
        echo "  budget of ${TOTAL_BUDGET_HOURS}h reached after $(human $elapsed) -- stopping"
        break
    fi

    out="$OUTPUT_DIR/n-${n}_m-${m}_k-${k}_samples-${SAMPLES}.csv"
    cmd=( "$BIN" -n "$n" -m "$m" -k "$k" --samples "$SAMPLES"
          --min-k "$k" --max-k "$k" --seed "$SEED"
          --max-threads "$MAX_THREADS" --output-path "$OUTPUT_DIR/" )

    if [[ -s "$out" ]]; then
        log_cell skipped 0 0 "$(basename "$out")" "${cmd[*]}"
        n_skip=$((n_skip+1)); continue
    fi

    # est is -1 where no published cell existed to estimate from. Those are
    # m=1 and the near-complete cells, structurally trivial, so they run with
    # the floor cap rather than being treated as unbounded.
    if awk -v e="$est" -v l="$EST_LIMIT" 'BEGIN{exit !(e > l)}'; then
        echo "  k=$k n=$n m=$m  skipped: estimated $(human "${est%.*}") > limit"
        log_cell too-expensive 0 0 "$(basename "$out")" "${cmd[*]}"
        n_over=$((n_over+1)); continue
    fi
    cap=$(awk -v e="$est" -v f="$CAP_FACTOR" -v lo="$CAP_MIN" -v hi="$CAP_MAX" \
        'BEGIN{c=(e<0?lo:e*f); if(c<lo)c=lo; if(c>hi)c=hi; printf "%d", c}')

    if [[ "$DRY_RUN" == 1 ]]; then
        printf '  DRY RUN k=%s n=%-2s m=%-3s est=%-8s cap=%ss\n' \
            "$k" "$n" "$m" "$(human "${est%.*}")" "$cap"
        continue
    fi

    printf '  [%3d] k=%s n=%-2s m=%-3s est %-7s ' \
        "$((n_ok+n_fail+1))" "$k" "$n" "$m" "$(human "${est%.*}")"
    guard_output "$out"
    if run_guarded "$cap" "$MAX_RSS_MB" "${cmd[@]}"; then
        printf 'took %8s  %5sMB\n' "$(human "$RG_SECONDS")" "$RG_PEAK_MB"
        log_cell ok "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
        n_ok=$((n_ok+1))
    else
        printf 'took %8s  %5sMB  %s -- removing partial output\n' \
            "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS"
        rm -f "$out"
        log_cell "$RG_STATUS" "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
        n_fail=$((n_fail+1))
    fi
    clear_output
done 3< "$LIST"

echo
echo "done in $(human $(( $(date +%s) - start )))  ok=$n_ok failed=$n_fail "\
"already-present=$n_skip over-est-limit=$n_over"
echo "manifest: $OUTPUT_DIR/RUN_MANIFEST.tsv"
