#!/usr/bin/env python3
"""Regenerate the `_matching_m_indices` / `_matching_diag_indices` companions.

These feed the refinement-set-size panels and, through `get_stat_dist`, the
pair filtering behind Figure 6. The original producer is in neither repo - the
deposit traces the published files to `updated_increasing_density.tar.gz`, old
working code - so this is a rewrite against the published output, which exists
for all 19 exhaustive cells and is diffed in --verify.

THE FORMAT STORES INDEX LISTS, NOT COUNTS. One line per projection, in the
source file's row order; each line is the comma-separated indices of the twins
in THAT row's twin list which satisfy the criterion. Lines are frequently a
single number, which is what makes "these are counts" such an easy misreading:
the refinement-size plot takes len() of each line, but get_stat_dist uses the
actual indices to select rows and columns of the pairwise distance matrix.
Emit counts and Figure 6 silently gets the wrong pairs.

THE INDICES ARE RELATIVE TO ONE FILE'S TWIN ORDER, and twins are discovered in
parallel, so that order differs between runs of the same cell. Indices
generated here will NOT equal the published ones numerically even though they
select the same twins. They are only meaningful alongside the enumeration they
were produced from - never mix a generation's indices with another's
enumeration or distance matrix.

Two criteria, both verified against all 1,895 published projections:

  matching_m    the twin has exactly m hyperedges, the same count as the
                k-uniform hypergraph that generated the projection.

  matching_diag the twin's projection DIAGONAL matches the generator's. For
                incidence matrix B the projection is B B^T, whose (i,i) entry
                is the number of hyperedges containing i - the node's degree.
                The search zeroes the diagonal, which is exactly why twins may
                differ there and why this is a strict refinement.
                The generator's degrees come from the projection itself:
                sum_v A[u][v] = (k-1)*deg(u) for a k-uniform generator, so
                deg(u) = rowsum(u)/(k-1). Any twin rebuilds A, since all twins
                of a row share the projection by definition.

    matching_indices.py <exhaustive.csv[.zst]> [-o OUTDIR]
    matching_indices.py --verify <published_dir> [m_min] [m_max]
"""
import os
import re
import subprocess
import sys

NAME_RE = re.compile(r"n-(\d+)_m-(\d+)_k-(\d+)_")


def read(path):
    if path.endswith(".zst"):
        return subprocess.run(["zstd", "-dcq", path],
                              capture_output=True, check=True).stdout.decode()
    if path.endswith(".gz"):
        return subprocess.run(["gzip", "-dc", path],
                              capture_output=True, check=True).stdout.decode()
    with open(path) as fh:
        return fh.read()


def twin_edges(section):
    return [[tuple(int(x) for x in e.split(":")) for e in t.split("|") if e]
            for t in section.split(";")] if section else []


def iter_twins(section):
    """Yield one twin's edge list at a time.

    m=20 is a SINGLE row holding 3,882,273 twins. Materialising them all as
    nested tuples costs ~10 GB of Python objects for no reason: both criteria
    are per-twin and need nothing but the twin itself, and only the FIRST twin
    is needed twice (to rebuild the projection). So stream.
    """
    if not section:
        return
    start = 0
    n = len(section)
    while start < n:
        end = section.find(";", start)
        if end == -1:
            end = n
        chunk = section[start:end]
        if chunk:
            yield [tuple(int(x) for x in e.split(":"))
                   for e in chunk.split("|") if e]
        start = end + 1


def degrees(edges, n):
    d = [0] * n
    for e in edges:
        for v in e:
            d[v] += 1
    return tuple(d)


def generator_degrees(twin, n, k):
    """Degree vector of the k-uniform hypergraph that produced this projection."""
    A = [[0] * n for _ in range(n)]
    for e in twin:
        for a in range(len(e)):
            for b in range(a + 1, len(e)):
                u, v = e[a], e[b]
                A[u][v] += 1
                A[v][u] += 1
    out = []
    for u in range(n):
        s = sum(A[u])
        # Integrality is a free check on k: a non-k-uniform generator, or the
        # wrong k parsed from the filename, shows up here rather than silently
        # producing a plausible wrong answer.
        if s % (k - 1):
            raise ValueError(f"row sum {s} not divisible by k-1={k-1}; "
                             f"is k={k} right for this file?")
        out.append(s // (k - 1))
    return tuple(out)


def indices_for(path, k=None):
    """-> (matching_m_lines, matching_diag_lines), one string per projection."""
    mo = NAME_RE.search(os.path.basename(path))
    if k is None:
        if not mo:
            raise ValueError(f"cannot read k from {os.path.basename(path)}; "
                             f"pass it explicitly")
        k = int(mo.group(3))
    out_m, out_d = [], []
    for line in read(path).splitlines():
        if not line.strip():
            continue
        parts = line.split("/")
        scal = parts[0].split(",")
        n, m = int(scal[0]), int(scal[1])
        section = parts[-1]
        if not section:
            out_m.append("")
            out_d.append("")
            continue
        gen = None
        hit_m, hit_d = [], []
        for i, t in enumerate(iter_twins(section)):
            if gen is None:
                gen = generator_degrees(t, n, k)
            if len(t) == m:
                hit_m.append(i)
            if degrees(t, n) == gen:
                hit_d.append(i)
        out_m.append(",".join(str(i) for i in hit_m))
        out_d.append(",".join(str(i) for i in hit_d))
    return out_m, out_d


def write_for(path, outdir=None):
    base = os.path.basename(path)
    for suffix in (".zst", ".gz"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    stem = base[: -len(".csv")] if base.endswith(".csv") else base
    d = outdir or os.path.dirname(path)
    mm, md = indices_for(path)
    written = []
    for kind, lines in (("m", mm), ("diag", md)):
        p = os.path.join(d, f"{stem}_matching_{kind}_indices.csv")
        with open(p, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        written.append((p, len(lines)))
    return written


def verify(pubdir, m_min, m_max):
    bad = checked = 0
    print(f"{'m':>3} {'projs':>6} {'matching_m':>12} {'matching_diag':>14}")
    for m in range(m_min, m_max + 1):
        stem = f"n-6_m-{m}_k-3_non-uniform_exhaustive_projections"
        src = os.path.join(pubdir, stem + ".csv.zst")
        if not os.path.exists(src):
            continue
        mm, md = indices_for(src)
        res = []
        for kind, got in (("m", mm), ("diag", md)):
            ref = os.path.join(pubdir, f"{stem}_matching_{kind}_indices.csv.zst")
            if not os.path.exists(ref):
                res.append("no reference")
                continue
            want = read(ref).splitlines()
            if want == got:
                res.append(f"{len(got)}/{len(got)} ok")
            else:
                n_bad = sum(1 for a, b in zip(want, got) if a != b) \
                        + abs(len(want) - len(got))
                res.append(f"DIFFERS on {n_bad}")
                bad += 1
        checked += 1
        print(f"{m:3d} {len(mm):6d} {res[0]:>12} {res[1]:>14}")
    print(f"\nchecked {checked} cells, {bad} criterion/cell pairs differing")
    return 1 if bad else 0


def main():
    argv = sys.argv[1:]
    if not argv:
        sys.exit(__doc__)
    if argv[0] == "--verify":
        if len(argv) < 2:
            sys.exit(__doc__)
        lo = int(argv[2]) if len(argv) > 2 else 2
        hi = int(argv[3]) if len(argv) > 3 else 20
        return verify(argv[1], lo, hi)
    outdir = None
    if "-o" in argv:
        i = argv.index("-o")
        outdir = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    for src in argv:
        for p, n in write_for(src, outdir):
            print(f"{n:6d} lines  {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
