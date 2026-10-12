#!/usr/bin/env python3
"""Build the Figure 6 distance-stats tables from the .f32.zst matrices.

A faster, equivalent replacement for calling distance_stats.write_new_stats
once per (distance, match_type). `get_stat_dist` stays untouched and is the
reference this verifies against (--verify).

Three things make it cheap enough to run on a laptop:

  1. ONE PASS PER DISTANCE, not one per (distance, match_type). Each
     projection's row is read once and feeds all four match types. 12 passes
     over 26 GiB become 3.
  2. The filter is a vectorised gather, not a Python double loop over every
     pair. The reference walks all N(N-1)/2 pairs testing `i in indices and
     j in indices` - about 110 billion iterations across the full set. The
     pairs wanted are exactly the pairs WITHIN the selected set, so their
     condensed indices can be computed directly.
  3. A light row parser for the enumeration. `read_data` builds every twin as
     tuples (or xgi objects) and throws them away; all that is needed per
     projection is total_twins, max_width, num_cliques and the number of
     DISTINCT size distributions.

Everything is streamed in chunks: a single m=13 projection is 828M values,
and `np.std(x, dtype=np.float64)` on a float32 array of that size allocates a
6.6 GB float64 temporary for `x - mean`.
"""
import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'analysis', 'python'))
from distance_stats import (open_text, read_refined_sets, iter_pairwise_rows,  # noqa: E402
                            stats_path, write_stats_file, get_N)

DISTANCES = ('hyperNS', 'hyperPD', 'jaccard')
MATCH_TYPES = ('all', 'm', 'diag', 'both')
# The published tables are CUMULATIVE and nested: m-2-10 is a prefix of
# m-2-13. Only the longest is computed; the rest are exact slices on the `m`
# column. (`update_stats` exists for the same reason - this was always how the
# tables extended.)
M_HI_RANGES = (8, 9, 10, 11, 12, 13)
CHUNK = 1 << 22            # values per streamed block: 32 MB as float64
GATHER_LIMIT = 50_000_000  # build a pair-index array in one shot below this


# ------------------------------------------------------- light row parser ---
def read_projection_meta(path):
    """Per-projection scalars, without parsing a single twin.

    Row layout: <csv scalars> / <size-dist> / <unfilt,filt> / ... / <twins>.
    `num_size_dists` counts DISTINCT size distributions, because read_data
    keys its dict on the distribution and merges repeats - so it is not the
    number of middle pairs.
    """
    out = {"total_twins": [], "max_width": [], "num_cliques": [], "num_size_dists": []}
    with open_text(path) as fin:
        for line in fin:
            parts = line.rstrip("\n").split("/")
            scal = parts[0].split(",")
            total_twins = 0
            seen = set()
            idx = 1
            while idx < len(parts) - 1:
                seen.add(parts[idx])
                total_twins += int(parts[idx + 1].split(",")[0])
                idx += 2
            out["total_twins"].append(total_twins)
            out["max_width"].append(int(scal[3]))
            out["num_cliques"].append(int(scal[5]))
            out["num_size_dists"].append(len(seen))
    return out


# --------------------------------------------------------- pair selection ---
def pair_indices(sel, N, lo, hi):
    """Condensed indices of the pairs within `sel`, for sel rows [lo, hi).

    The condensed index of (a, b), a < b, in an upper triangle over N items is
        a*N - a*(a+1)/2 + (b - a - 1)
    """
    out = []
    for t in range(lo, hi):
        a = int(sel[t])
        b = sel[t + 1:]
        if b.size:
            out.append(a * N - (a * (a + 1)) // 2 + (b - a - 1))
    return np.concatenate(out) if out else np.empty(0, dtype=np.int64)


def iter_selected(row, sel, N):
    """Yield float64 blocks of the pairwise distances within `sel`."""
    s = sel.size
    if s < 2:
        return
    if s * (s - 1) // 2 <= GATHER_LIMIT:
        yield row[pair_indices(sel, N, 0, s - 1)].astype(np.float64)
        return
    # Rare: a large projection where nearly every twin matches. Walk the
    # selected rows in blocks so no index array gets out of hand.
    step = max(1, GATHER_LIMIT // max(s, 1))
    for lo in range(0, s - 1, step):
        p = pair_indices(sel, N, lo, min(lo + step, s - 1))
        if p.size:
            yield row[p].astype(np.float64)


def iter_all(row):
    for lo in range(0, row.size, CHUNK):
        yield row[lo:lo + CHUNK].astype(np.float64)


def mean_std(make_iter, count):
    """Population mean and std, two-pass, in float64 and bounded memory.

    Two passes rather than sum/sumsq: the one-pass identity cancels
    catastrophically when the variance is small against the mean, and these
    distances sit in a narrow band. This is the same algorithm numpy uses, so
    it agrees with the reference to rounding.
    """
    total = 0.0
    for block in make_iter():
        total += float(block.sum())
    mean = total / count
    ss = 0.0
    for block in make_iter():
        d = block - mean
        ss += float(np.dot(d, d))
    return mean, math.sqrt(ss / count)


def _empty_stats():
    return {"min": [], "max": [], "mean": [], "std": [], "entropy": [],
            "max_width": [], "num_cliques": [], "num_size_dists": [], "m": []}


# ---------------------------------------------------------------- driver ----
def collect(dist, args, log=None):
    """One pass over every projection, filling all four match types at once."""
    n, k = args.n, args.k
    stats = {mt: _empty_stats() for mt in MATCH_TYPES}
    t0 = time.time()

    for m in range(args.m_min, args.m_max + 1):
        base = f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections"
        meta = read_projection_meta(os.path.join(args.exhaustive_dir, base + ".csv"))
        row_lengths = [t * (t - 1) // 2 for t in meta["total_twins"]]
        sel_sets = {
            mt: read_refined_sets(args.indices_dir + os.sep, f"{base}_matching_{mt}_indices.csv")
            for mt in ("m", "diag")
        }

        for proj_idx, row in enumerate(
                iter_pairwise_rows(args.dist_dir + os.sep,
                                   base + f"_{dist}.csv.gz", row_lengths)):
            N = get_N(row.shape[0])
            sel = {"m": np.asarray(sel_sets["m"][proj_idx], dtype=np.int64),
                   "diag": np.asarray(sel_sets["diag"][proj_idx], dtype=np.int64)}
            sel["both"] = np.intersect1d(sel["m"], sel["diag"])

            for mt in MATCH_TYPES:
                if mt == "all":
                    count = row.shape[0]
                    make_iter = lambda r=row: iter_all(r)
                else:
                    s_idx = sel[mt]
                    count = s_idx.size * (s_idx.size - 1) // 2
                    make_iter = lambda r=row, q=s_idx, nn=N: iter_selected(r, q, nn)
                if count == 0:
                    continue
                if count == 1:
                    # The reference special-cases this and stores a plain int 0
                    # for std. write_stats_file keeps int and float apart, and
                    # 69 of the 92 published tables carry that int.
                    mean, std = float(next(make_iter())[0]), 0
                else:
                    mean, std = mean_std(make_iter, count)
                col = stats[mt]
                col["mean"].append(mean)
                col["std"].append(std)
                col["max_width"].append(meta["max_width"][proj_idx])
                col["num_cliques"].append(meta["num_cliques"][proj_idx])
                col["num_size_dists"].append(meta["num_size_dists"][proj_idx])
                col["m"].append(m)
            del row
        if log:
            log(f"  {dist:<8} m={m:<3} {len(meta['total_twins']):>4} projections"
                f"  {time.time() - t0:>6.0f}s")
    return stats


def write_tables(dist, stats, args):
    """The published ranges are cumulative and nested, so slice, do not recompute."""
    written = []
    for mt in MATCH_TYPES:
        full = stats[mt]
        for m_hi in M_HI_RANGES:
            if m_hi > args.m_max or m_hi < args.m_min:
                continue
            keep = [i for i, mm in enumerate(full["m"]) if mm <= m_hi]
            sliced = _empty_stats()
            for c in ("mean", "std", "max_width", "num_cliques", "num_size_dists", "m"):
                sliced[c] = [full[c][i] for i in keep]
            path = stats_path(args.out_dir + os.sep, args.n, args.k,
                              args.m_min, m_hi, dist, mt) + ".csv.zst"
            write_stats_file(path, sliced)
            written.append(path)
    return written


def verify(dist, args, log=print):
    """Compare every column against distance_stats.get_stat_dist."""
    from distance_stats import get_stat_dist, read_data
    m_vals = list(range(args.m_min, args.m_max + 1))
    output_datas = {
        m: read_data(os.path.join(args.exhaustive_dir,
                                  f"n-{args.n}_m-{m}_k-{args.k}"
                                  f"_non-uniform_exhaustive_projections.csv"),
                     xgi_hypergraphs=False)
        for m in m_vals
    }
    mine = collect(dist, args)
    worst = 0.0
    for mt in MATCH_TYPES:
        ref = get_stat_dist(args.n, args.k, output_datas, dist, mt, m_vals,
                            args.indices_dir + os.sep, dist_dir=args.dist_dir + os.sep)
        got = mine[mt]
        assert len(ref["mean"]) == len(got["mean"]), \
            f"{mt}: {len(ref['mean'])} reference rows vs {len(got['mean'])}"
        for c in ("max_width", "num_cliques", "num_size_dists", "m"):
            assert list(ref[c]) == list(got[c]), f"{mt}: column {c} differs"
        for c in ("mean", "std"):
            for a, b in zip(ref[c], got[c]):
                worst = max(worst, abs(float(a) - float(b)))
        # the int-0 std convention has to survive too
        ref_ints = [i for i, v in enumerate(ref["std"]) if isinstance(v, int)]
        got_ints = [i for i, v in enumerate(got["std"]) if isinstance(v, int)]
        assert ref_ints == got_ints, f"{mt}: int-0 std rows differ"
        log(f"  {dist:<8} {mt:<5} {len(ref['mean']):>6} rows  identical columns, "
            f"{len(ref_ints)} int-0 std")
    log(f"  {dist}: largest |difference| in mean or std = {worst:.3e}")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exhaustive-dir', required=True)
    ap.add_argument('--indices-dir', default=None, help='defaults to --exhaustive-dir')
    ap.add_argument('--dist-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--distances', nargs='*', default=list(DISTANCES))
    ap.add_argument('--n', type=int, default=6)
    ap.add_argument('--k', type=int, default=3)
    ap.add_argument('--m-min', type=int, default=2)
    ap.add_argument('--m-max', type=int, default=13)
    ap.add_argument('--cores', type=int, default=3,
                    help='processes, one per distance; 3 is the natural cap')
    ap.add_argument('--verify', action='store_true',
                    help='compare against get_stat_dist instead of writing')
    args = ap.parse_args()
    if args.indices_dir is None:
        args.indices_dir = args.exhaustive_dir
    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)

    if args.verify:
        worst = 0.0
        for d in args.distances:
            worst = max(worst, verify(d, args))
        print(f"\nlargest |difference| across all distances: {worst:.3e}")
        return 0

    t0 = time.time()
    if args.cores > 1 and len(args.distances) > 1:
        from multiprocessing import Pool
        with Pool(min(args.cores, len(args.distances))) as p:
            out = p.starmap(_worker, [(d, args) for d in args.distances])
        written = [f for sub in out for f in sub]
    else:
        written = []
        for d in args.distances:
            written += write_tables(d, collect(d, args, log=print), args)
    print(f"\n{len(written)} table(s) in {(time.time() - t0) / 60:.1f} min -> {args.out_dir}")
    return 0


def _worker(dist, args):
    return write_tables(dist, collect(dist, args, log=lambda t: print(t, flush=True)), args)


if __name__ == '__main__':
    sys.exit(main())
