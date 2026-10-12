"""Tests for projection_utils.get_exact_max_width.

Runs standalone (`python3 test_projection_utils.py`) or under pytest; it needs
only numpy, which projection_utils already requires.

The function computes the paper's worst case search tree size,
prod_e binom(|eta_e|, w_e). The version used to produce the published figures
had two bugs, both of which pushed the value DOWN:

  1. it skipped k == 2, dropping the 2-clique {i,j} from eta_e;
  2. it seeded with `if max_width == 0`, re-seeding the product instead of
     multiplying whenever a factor was 0.

These tests pin both, and pin the k-uniform case where nothing should change.
"""
import numpy as np

from projection_utils import get_projection, get_cliques, get_cliques_nx, get_exact_max_width


def triangle():
    """One hyperedge {0,1,2}: a triangle with all weights 1."""
    return get_projection(3, [[0, 1, 2]])


def test_two_clique_is_counted_once():
    """Worked example small enough to verify by hand.

    Cliques of the triangle are the three 2-cliques and the one 3-clique, so
    each edge {u,v} has eta_e = {{u,v}, {0,1,2}} and |eta_e| = 2. The product
    is binom(2,1)^3 = 8.

    The old code skipped k == 2, giving |eta_e| = 1 and binom(1,1)^3 = 1.
    """
    assert get_exact_max_width(triangle()) == 8


def test_zero_factor_zeroes_the_product():
    """An edge that cannot be satisfied makes the whole product 0.

    Here the weight of edge (0,1) is 3 but only two cliques contain it, so
    binom(2,3) == 0. The old re-seeding logic discarded that factor and
    returned a positive number instead.
    """
    proj = np.array([[0., 3., 1.],
                     [3., 0., 1.],
                     [1., 1., 0.]])
    assert get_exact_max_width(proj) == 0


def test_min_k_above_two_excludes_the_two_clique():
    """With min_k = 3 a 2-hyperedge is not a legal way to satisfy an edge, so
    eta_e holds only the triangle and the product is binom(1,1)^3 = 1."""
    assert get_exact_max_width(triangle(), min_k=3) == 1


def test_both_clique_enumerations_agree():
    """get_cliques and get_cliques_nx must give the same answer, since callers
    pass either one."""
    hypergraphs = [
        [[0, 1, 2]],
        [[0, 1, 2], [1, 2, 3]],
        [[0, 1, 2], [1, 2, 3], [0, 2, 3], [0, 1, 3]],
        [[0, 1, 2], [2, 3, 4], [0, 4, 5]],
    ]
    for hg in hypergraphs:
        n = max(v for he in hg for v in he) + 1
        proj = get_projection(n, hg)
        assert (get_exact_max_width(proj, cliques=get_cliques(proj))
                == get_exact_max_width(proj, cliques=get_cliques_nx(proj))), hg


def test_matches_independent_brute_force():
    """Cross-check against a deliberately naive implementation: enumerate every
    vertex subset, keep the cliques, and count the ones containing each edge."""
    from itertools import combinations
    from math import comb

    def brute(proj, min_k=2):
        n = proj.shape[0]
        cliques = []
        for size in range(min_k, n + 1):
            for cand in combinations(range(n), size):
                if all(proj[a, b] > 0 for a, b in combinations(cand, 2)):
                    cliques.append(set(cand))
        out = 1
        for i in range(n):
            for j in range(i + 1, n):
                if proj[i, j] > 0:
                    deg = sum(1 for c in cliques if i in c and j in c)
                    out *= comb(deg, int(proj[i, j]))
        return out

    hypergraphs = [
        [[0, 1, 2]],
        [[0, 1, 2], [1, 2, 3]],
        [[0, 1, 2], [1, 2, 3], [0, 2, 3], [0, 1, 3]],
        [[0, 1, 2], [2, 3, 4], [0, 4, 5]],
        [[0, 1], [1, 2], [0, 2], [0, 1, 2]],
    ]
    for hg in hypergraphs:
        n = max(v for he in hg for v in he) + 1
        proj = get_projection(n, hg)
        for min_k in (2, 3):
            assert get_exact_max_width(proj, min_k=min_k) == brute(proj, min_k), (hg, min_k)


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL {name}: {exc}")
    print(f"\n{failures} failure(s)")
    raise SystemExit(1 if failures else 0)
