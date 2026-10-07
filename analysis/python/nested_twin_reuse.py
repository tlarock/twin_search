#!/usr/bin/env python3
"""How much of cell m+1's twin set is already answered by cell m's?

THE STRUCTURE THIS EXPLOITS
---------------------------
sample_uniform_random draws hyperedges from one std::uniform_int_distribution
and inserts them into a std::set until it holds m distinct ones. The generator
is seeded from (seed, sample index, attempt) - m is NOT in the stream - so for
a fixed index the draw sequence is the same whatever m is, and the hyperedge
SET for m is a prefix-set of the one for m+1:

    H_{m+1} = H_m + {e}          e = the (m+1)-th distinct draw

(Both must cover all n nodes, or the generator restarts; at these densities it
never does. The script checks rather than assumes.)

Projections therefore nest too, P_{m+1} = P_m + delta(e), and the twin sets
decompose EXACTLY:

    T(P_{m+1})  =  { G' + {e} : G' in T(P_m), e not in G' }
                 U { G in T(P_{m+1}) : e not in G }

The first part is a bijection with a set we have already computed and written
to disk. The second is a search for twins avoiding one named hyperedge.

WHAT THIS MEASURES
------------------
1. The bijection itself, as an equality that must hold exactly:

       |{G in T(P_{m+1}) : e in G}|  ==  |{G' in T(P_m) : e not in G'}|

   on real output. If this fails the idea is wrong, not merely unprofitable.

2. The fraction of T(P_{m+1}) the first part accounts for. That fraction is
   the work a DP would skip, and it is the CEILING on any speedup - the
   remaining search still has to run.

3. The same fraction weighted by each sample's measured runtime, because the
   cells are heavy-tailed and an unweighted mean over 450 samples says almost
   nothing about where the hours go.

WHAT IT DOES NOT MEASURE
------------------------
Whether the constrained search in the second part is actually cheaper in
proportion to the twins it no longer returns. Backtracking cost is tree nodes,
not leaves, and forbidding one clique prunes branches without shrinking the
rest of the tree. Treat the weighted fraction as an upper bound and measure the
constrained search directly before believing it.

Usage:
    nested_twin_reuse.py --results <density-500 dir> --fp <fingerprint dir> \
                         --k 3 --n 9 --m 47 [--limit N]

    --fp holds fp-m<M>.txt produced by
        count_twins_random ... --seed <SEED> --dry-run | grep '^FP '
    RUN THE DRY RUN ON LINUX. Sampling is libstdc++-specific; see docker/.
"""
import os
import sys


def die(msg):
    sys.exit("error: " + msg)


def parse_args(argv):
    a = {"--results": None, "--fp": None, "--k": "3", "--n": "9",
         "--m": None, "--limit": "0"}
    i = 0
    while i < len(argv):
        if argv[i] in a and i + 1 < len(argv):
            a[argv[i]] = argv[i + 1]
            i += 2
        else:
            die("unrecognised argument %r\n%s" % (argv[i], __doc__))
    # No defaults for the paths on purpose. A wrong-but-plausible default turns
    # "I cannot find your data" into a confident wrong answer.
    for need in ("--results", "--fp", "--m"):
        if not a[need]:
            die("%s is required\n%s" % (need, __doc__))
    return a


def read_fingerprints(path, n):
    """-> {index: tuple of upper-triangle entries}, in (i<j) row-major order."""
    want = n * (n - 1) // 2
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("FP "):
                continue
            i, _w, _e, _c, proj = line[3:].rstrip("\n").split(",", 4)
            vals = tuple(int(x) for x in proj.split("-"))
            if len(vals) != want:
                die("%s: index %s has %d projection entries, expected %d"
                    % (path, i, len(vals), want))
            out[int(i)] = vals
    if not out:
        die("no 'FP ' lines in %s" % path)
    return out


def pair_list(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def new_hyperedge(small, big, n):
    """The hyperedge present in H_{m+1} but not H_m, from the two projections.

    Adding one k=3 hyperedge {a,b,c} raises exactly the three pairs it spans by
    one. Anything else means the two samples are not nested - which happens if
    the generator restarted for one m and not the other - so return None and
    let the caller skip the index rather than silently comparing unrelated
    hypergraphs.
    """
    raised = [p for p, (s, b) in zip(pair_list(n), zip(small, big)) if b - s == 1]
    if len(raised) != 3 or any(b - s not in (0, 1) for s, b in zip(small, big)):
        return None
    nodes = sorted({u for p in raised for u in p})
    if len(nodes) != 3:
        return None
    # The three raised pairs must be exactly the triangle on those nodes.
    a, b, c = nodes
    if set(raised) != {(a, b), (a, c), (b, c)}:
        return None
    return (a, b, c)


def split_row(line):
    """-> (n, runtime_ms, num_twins, twins_field) for one result row.

    Layout is scalars/(size-dist/counts)*/twins. There is normally one
    (size-dist, counts) pair, because min_k == max_k == k forces every twin to
    the same size distribution, but the writer permits several.
    """
    parts = line.rstrip("\n").split("/")
    scal = parts[0].split(",")
    n, runtime = int(scal[0]), int(scal[2])
    idx, total = 1, 0
    while idx + 1 < len(parts) and ":" in parts[idx] and "," in parts[idx + 1]:
        total += int(parts[idx + 1].split(",")[0])
        idx += 2
    twins = parts[idx] if idx < len(parts) else ""
    return n, runtime, total, twins


def projection_key_of_row(twins_field, n):
    """The sample's projection, rebuilt from any one of its twins.

    Every twin has the row's projection by definition, and twins share node
    labels, so the first one is enough.
    """
    if not twins_field:
        return None
    mat = [[0] * n for _ in range(n)]
    for he in twins_field.split(";", 1)[0].split("|"):
        nodes = [int(x) for x in he.split(":")]
        for x in range(len(nodes)):
            for y in range(x + 1, len(nodes)):
                u, v = sorted((nodes[x], nodes[y]))
                mat[u][v] += 1
    return tuple(mat[i][j] for i, j in pair_list(n))


def index_rows(path, n):
    """-> {projection: byte offset of a row with that projection}.

    One pass, storing offsets rather than text: these files run to hundreds of
    megabytes and only a few rows are ever needed in full.
    """
    out = {}
    with open(path, "rb") as fh:
        off = fh.tell()
        for raw in fh:
            line = raw.decode()
            _n, _rt, _tot, twins = split_row(line)
            key = projection_key_of_row(twins, n)
            if key is not None and key not in out:
                out[key] = off
            off += len(raw)
    return out


def row_at(path, off):
    with open(path, "rb") as fh:
        fh.seek(off)
        return fh.readline().decode()


def main(argv):
    a = parse_args(argv)
    k, n, m = int(a["--k"]), int(a["--n"]), int(a["--m"])
    limit = int(a["--limit"])
    res, fpdir = a["--results"], a["--fp"]

    paths, fps = {}, {}
    for mm in (m, m + 1):
        csv = os.path.join(res, "n-%d_m-%d_k-%d_samples-500.csv" % (n, mm, k))
        fp = os.path.join(fpdir, "fp-m%d.txt" % mm)
        for f in (csv, fp):
            if not os.path.exists(f):
                die("missing %s" % f)
        paths[mm] = csv
        fps[mm] = read_fingerprints(fp, n)

    print("cell k=%d n=%d, m=%d -> m=%d" % (k, n, m, m + 1))
    print("  %s" % os.path.basename(paths[m]))
    print("  %s" % os.path.basename(paths[m + 1]))

    idx = {mm: index_rows(paths[mm], n) for mm in (m, m + 1)}
    print("  distinct projections on disk: m=%d %d rows, m=%d %d rows\n"
          % (m, len(idx[m]), m + 1, len(idx[m + 1])))

    shared, not_nested = [], 0
    for i in sorted(set(fps[m]) & set(fps[m + 1])):
        e = new_hyperedge(fps[m][i], fps[m + 1][i], n)
        if e is None:
            not_nested += 1
            continue
        if fps[m][i] in idx[m] and fps[m + 1][i] in idx[m + 1]:
            shared.append((i, e))
    if limit:
        shared = shared[:limit]
    print("indices with BOTH cells computed and nested: %d" % len(shared))
    if not_nested:
        print("indices where the hypergraphs are NOT nested: %d" % not_nested)
    if not shared:
        die("no index has both cells on disk; nothing to measure")

    tot_big = tot_reuse = 0
    rt_total = rt_saved = 0
    mismatches = []
    per_sample = []

    for i, e in shared:
        he = "%d:%d:%d" % e
        needle = "|" + he + "|"

        line_b = row_at(paths[m + 1], idx[m + 1][fps[m + 1][i]])
        _n, rt_b, n_b, tw_b = split_row(line_b)
        with_e = sum(1 for t in tw_b.split(";") if needle in "|" + t + "|")

        line_s = row_at(paths[m], idx[m][fps[m][i]])
        _n, _rt, n_s, tw_s = split_row(line_s)
        without_e = sum(1 for t in tw_s.split(";") if needle not in "|" + t + "|")

        if with_e != without_e:
            mismatches.append((i, he, with_e, without_e, n_b, n_s))

        tot_big += n_b
        tot_reuse += with_e
        rt_total += rt_b
        rt_saved += rt_b * (with_e / n_b) if n_b else 0
        per_sample.append((rt_b, with_e / n_b if n_b else 0.0, n_b, i))

    print("\n=== 1. the bijection, checked exactly ===")
    print("  |{G in T(P_{m+1}) : e in G}| == |{G' in T(P_m) : e not in G'}|")
    if mismatches:
        print("  FAILED on %d of %d samples" % (len(mismatches), len(shared)))
        for row in mismatches[:5]:
            print("    i=%d e=%s  with_e=%d  without_e=%d  (|T_big|=%d |T_small|=%d)" % row)
    else:
        print("  HOLDS on all %d samples (%d twins at m=%d checked)"
              % (len(shared), tot_big, m + 1))

    print("\n=== 2. share of T(P_{m+1}) that the previous cell already answers ===")
    print("  twins at m=%d       : %d" % (m + 1, tot_big))
    print("  of those, contain e : %d  (%.1f%% pooled)"
          % (tot_reuse, 100.0 * tot_reuse / tot_big if tot_big else 0))
    fr = sorted(f for _rt, f, _nb, _i in per_sample)
    print("  per-sample fraction : median %.1f%%  min %.1f%%  max %.1f%%"
          % (100 * fr[len(fr) // 2], 100 * fr[0], 100 * fr[-1]))

    print("\n=== 3. weighted by measured runtime (where the hours actually are) ===")
    print("  total search time in these samples : %.1f s" % (rt_total / 1000.0))
    print("  CEILING on time a DP could remove  : %.1f s  (%.1f%%)"
          % (rt_saved / 1000.0, 100.0 * rt_saved / rt_total if rt_total else 0))
    print("  (assumes cost falls in proportion to twins not returned - an")
    print("   upper bound, since backtracking pays for tree nodes, not leaves)")

    per_sample.sort(reverse=True)
    print("\n  ten most expensive samples:")
    print("    %-7s %10s %10s %8s" % ("index", "runtime_s", "twins", "reusable"))
    for rt, f, nb, i in per_sample[:10]:
        print("    %-7d %10.1f %10d %7.1f%%" % (i, rt / 1000.0, nb, 100 * f))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
