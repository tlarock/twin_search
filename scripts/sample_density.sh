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
# COST WARNING. Unlike the heatmap, this gets genuinely expensive - and the
# cost is NOT monotone in m. It peaks near m = C(n,3)/2 and then collapses, by
# seven orders of magnitude, as the hypergraph approaches complete and the twin
# set becomes forced. From the published samples-500 data (mean ms per sample):
#
#     n=8, k=3:  m=12 -> 0.7       m=25 -> 7,482 (peak)   m=55 -> 0.0
#     n=9, k=3:  m=20 -> 123       m=43 -> 6,938,925      m=80 -> 0.1
#
# The expense is therefore a BAND in the middle, not a tail. Cells over
# MAX_SECONDS are skipped and recorded; the script keeps going rather than
# abandoning the rest of an n, because everything past the peak is cheap.
# Raise MAX_SECONDS if you want the middle band; M_MAX stops short of it.
#
# These runtimes come from the old code, before the isomorphism pre-filter sped
# post-processing up 7-15x, so they overstate current cost - but the shape is
# the point, not the absolute numbers. This is an NP-hard search with a
# heavy-tailed cost distribution: individual samples in the middle band really
# can take hours.
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

        have=$(rows_in "$out")
        if (( have >= SAMPLES )); then
            log_cell skipped 0 0 "$(basename "$out")" "${cmd[*]}"
            continue
        fi
        if (( have > 0 )); then
            echo "    m=$m  present but short ($have/$SAMPLES) -- left as is"
            log_cell short "$have" 0 "$(basename "$out")" "${cmd[*]}"
            continue
        fi
        if [[ "$DRY_RUN" == 1 ]]; then echo "    DRY RUN: ${cmd[*]}"; continue; fi

        printf '    m=%-3s ' "$m"
        guard_output "$out"
        if run_guarded "$MAX_SECONDS" "$MAX_RSS_MB" "${cmd[@]}"; then
            printf '%8s  %5sMB\n' "$(human "$RG_SECONDS")" "$RG_PEAK_MB"
            log_cell ok "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
        else
            # Keep whatever samples completed; only the row in flight is lost.
            trim_torn_line "$out"
            have=$(rows_in "$out")
            if (( have > 0 )); then
                printf '%8s  %5sMB  %s at %s/%s samples -- KEPT\n' \
                    "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS" "$have" "$SAMPLES"
                log_cell partial "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
            else
                printf '%8s  %5sMB  %s with 0 samples\n' \
                    "$(human "$RG_SECONDS")" "$RG_PEAK_MB" "$RG_STATUS"
                rm -f "$out"
                log_cell "$RG_STATUS" "$RG_SECONDS" "$RG_PEAK_MB" "$(basename "$out")" "${cmd[*]}"
            fi
            # Deliberately NOT breaking out of the m loop on a failure.
            #
            # An earlier version did, assuming cost rises monotonically with m.
            # It does not, over this script's range. Measured from the
            # published samples-500 data (runtime field, mean ms per sample):
            #
            #   n=9, k=3:  m=20 -> 123 ms,  m=43 -> 6,938,925 ms (peak),
            #              m=60 -> 57,102,  m=80 -> 0.1,  m=83 -> 0.8
            #   n=8, k=3:  peak at m=25 (7,482 ms), 0.0 ms by m=55
            #
            # Cost peaks near m = C(n,k)/2 and then collapses by seven orders
            # of magnitude, because as the hypergraph approaches complete the
            # twin set becomes forced - at m == C(n,k) it is exactly one
            # hypergraph. Breaking at the peak would skip the whole dense tail,
            # which contains the CHEAPEST cells in the sweep.
            #
            # The heatmap script never sees this: its window stops at m=16,
            # which is only 23% of C(8,4), so it only ever observes the rising
            # limb. MAX_SECONDS bounds the cost of continuing.
        fi
        clear_output
    done
done

echo
echo "total $(human $(( $(date +%s) - total_start )))   manifest: $OUTPUT_DIR/RUN_MANIFEST.tsv"
