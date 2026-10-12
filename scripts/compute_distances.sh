#!/usr/bin/env bash
# Pairwise twin-distance matrices for the exhaustive cells (Figure 6).
#
#   ./scripts/compute_distances.sh --distance hyperNS,jaccard --cores 6
#   ./scripts/compute_distances.sh --distance hyperPD --cores 6
#
# Options
#   --distance  comma list of hyperNS, hyperPD, jaccard   (default: all three)
#   --cores     worker processes                          (default: 6)
#   --m-min/--m-max  cell range                           (default: 2..13)
#   --input     directory of *_exhaustive_projections.csv.zst
#   --output    where the matrices go
#   --force     recompute cells that already have output
#   --dry-run   print the plan and stop
#
# WHY m=2..13 AND NOT 14..20: the matrix for one projection is T(T-1)/2
# distances. m=13's largest row has 40,696 twins (828M pairs); m=19's has
# 1,500,229, which is 1.1e12 pairs - about 4.5 TB as float32 for a single row.
# The published set stops at 13 and so does this.
#
# COST, measured on this machine (not estimated):
#     per pair   hyperNS 5.7 us   hyperPD 33.1 us   jaccard 9.3 us
#     12,165,654,831 pairs over m=2..13
#                hyperNS 19.3 core-h   hyperPD 112.0   jaccard 31.4
#   at 6 cores:  hyperNS  ~3.2 h       hyperPD ~18.7 h  jaccard  ~5.2 h
#
# The inner script's Pool parallelises over PROJECTIONS, and one projection's
# O(T^2) loop runs on a single core. So wall time has a FLOOR set by the
# largest projection - hyperPD m=13 is 7.6 h on one core - and adding cores
# past ~22 buys nothing. 6 cores is a sensible laptop setting, not a
# compromise that costs much.
#
# hyperNS WILL NOT match the published matrices. feature_vec went from 45 dims
# to 21 upstream; the vendored header records differences up to 0.20 on m=4.
# That is the agreed update to the current measure, not a reproduction
# failure. hyperPD and jaccard should reproduce.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GRAM=/Users/tl2574/git/gram_mates
PY="$GRAM/.venv/bin/python"
DRIVER="$GRAM/results/python/parallel_distances.py"

DISTS=hyperNS,hyperPD,jaccard
CORES=6
M_MIN=2
M_MAX=13
IN="$REPO/results/reproduce/exhaustive"
OUT="$REPO/results/reproduce/distance-matrices"
FORCE=0
DRY=0

while [ $# -gt 0 ]; do
    case "$1" in
        --distance) DISTS="$2"; shift 2 ;;
        --cores)    CORES="$2"; shift 2 ;;
        --m-min)    M_MIN="$2"; shift 2 ;;
        --m-max)    M_MAX="$2"; shift 2 ;;
        --input)    IN="$2";    shift 2 ;;
        --output)   OUT="$2";   shift 2 ;;
        --force)    FORCE=1;    shift ;;
        --dry-run)  DRY=1;      shift ;;
        -h|--help)  sed -n '2,40p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

for f in "$PY" "$DRIVER"; do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done
"$PY" -c 'import xgi, scipy, numpy, tqdm' 2>/dev/null || {
    echo "the gram_mates venv is missing a dependency (xgi/scipy/numpy/tqdm)" >&2; exit 1; }

mkdir -p "$OUT"
LOG="$OUT/compute_distances.log"

say() { printf '%s %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"; }

IFS=',' read -r -a DIST_ARR <<< "$DISTS"
for d in "${DIST_ARR[@]}"; do
    case "$d" in hyperNS|hyperPD|jaccard) ;;
        *) echo "unknown distance: $d (want hyperNS, hyperPD or jaccard)" >&2; exit 2 ;;
    esac
done

say "=== compute_distances: ${DISTS}, m=${M_MIN}..${M_MAX}, ${CORES} cores ==="
say "    input  $IN"
say "    output $OUT"

TOTAL_START=$(date +%s)
for d in "${DIST_ARR[@]}"; do
  for m in $(seq "$M_MIN" "$M_MAX"); do
    BASE="n-6_m-${m}_k-3_non-uniform_exhaustive_projections"
    SRC="$IN/$BASE.csv.zst"
    DST="$OUT/${BASE}_${d}.csv.gz"
    [ -e "$SRC" ] || { say "  m=$m $d: no enumeration at $SRC - skipping"; continue; }
    if [ -s "$DST" ] && [ "$FORCE" -eq 0 ]; then
        say "  m=$m $d: already present ($(wc -c < "$DST" | tr -d ' ') bytes) - skipping"
        continue
    fi
    if [ "$DRY" -eq 1 ]; then say "  DRY RUN would compute m=$m $d"; continue; fi

    # The driver reads a PLAIN csv and writes its per-projection temp files
    # next to that input, so stage into a work dir rather than letting it
    # litter results/reproduce/exhaustive - a stray .csv there is exactly what
    # sync-della.sh refuses to pull and what resolve_path would prefer over
    # the real archive.
    WORK=$(mktemp -d "${TMPDIR:-/tmp}/dist-${BASE}-${d}.XXXXXX")
    trap 'rm -rf "$WORK"' EXIT
    zstd -dqf "$SRC" -o "$WORK/$BASE.csv" || { say "  m=$m $d: decompress FAILED"; rm -rf "$WORK"; continue; }

    START=$(date +%s)
    say "  m=$m $d: starting ($(wc -l < "$WORK/$BASE.csv" | tr -d ' ') projections)"
    # nice: this runs for hours on a laptop and should yield to everything else.
    ( cd "$GRAM/results/python" && \
      nice -n 10 "$PY" "$DRIVER" -n 6 -k 3 -m "$m" \
        --results-dir "$WORK/" --output-dir "$WORK/" \
        "--${d}" --nprocs "$CORES" ) >> "$LOG" 2>&1
    rc=$?
    if [ "$rc" -ne 0 ] || [ ! -s "$WORK/${BASE}_${d}.csv.gz" ]; then
        say "  m=$m $d: FAILED rc=$rc after $(( $(date +%s) - START ))s"
        rm -rf "$WORK"; trap - EXIT; continue
    fi
    mv "$WORK/${BASE}_${d}.csv.gz" "$DST"
    say "  m=$m $d: ok in $(( $(date +%s) - START ))s -> $(wc -c < "$DST" | tr -d ' ') bytes"
    rm -rf "$WORK"; trap - EXIT
  done
done
say "=== total $(( ($(date +%s) - TOTAL_START) / 60 )) min ==="
