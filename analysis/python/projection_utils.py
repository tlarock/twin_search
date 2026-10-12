import numpy as np
import networkx as nx
from math import comb

def get_projection(n, hypergraph, diagonal=False):
    proj = np.zeros((n,n))
    for he in hypergraph:
        for i in range(len(he)):
            u = he[i]
            for j in range(i+1, len(he)):
                v = he[j]
                proj[u,v] += 1
                proj[v,u] += 1
            if diagonal:
                proj[u,u] += 1
    return proj


def get_cliques_nx(proj_in):
    proj = np.array(proj_in)
    proj[proj > 0] = 1
    cliques = {2: set()}
    G = nx.from_numpy_array(proj)
    for c in nx.clique.enumerate_all_cliques(G):
        if len(c) > 1:
            clique_len = len(c)
            if clique_len in cliques:
                cliques[clique_len].add(tuple(sorted(c)))
            else:
                cliques[clique_len] = set([tuple(sorted(c))])
    return cliques

def get_cliques(proj_in):
    proj = np.array(proj_in)
    # binarize the projection
    proj[proj > 0] = 1
    cliques = {2: set()}
    # get the edges
    for i in range(proj.shape[0]):
        for j in range(i+1, proj.shape[1]):
            if proj[i,j] > 0:
                cliques[2].add(tuple(sorted([i,j])))
    
    k = 3
    while len(cliques[k-1]) > 0:
        cliques[k] = set()
        for c in cliques[k-1]:
            # sum the rows
            rowsum = proj[list(c), :].sum(axis=0)
            for u in range(rowsum.shape[0]):
                if rowsum[u] == k-1:
                    new_clique = set(c)
                    new_clique.add(u)
                    cliques[len(new_clique)].add(tuple(sorted(new_clique)))
        k += 1
    return cliques

def get_exact_max_width(proj, cliques = None, min_k = 2):
    """Worst case search tree size: prod_e binom(|eta_e|, w_e).

    |eta_e| is the number of DISTINCT cliques that could satisfy edge e, i.e.
    every clique of size >= min_k containing both endpoints. When min_k <= 2
    that INCLUDES the 2-clique {i,j} itself, because a 2-hyperedge is a
    legitimate way to satisfy the edge in a non-uniform search.

    Two corrections relative to the version this replaces, which was used to
    produce the published figures:

      1. It skipped k == 2 outright (`if k == 2: continue`), so |eta_e| was one
         too small for every edge whenever min_k <= 2. The plotted worst case
         tree sizes were correspondingly too LOW - a median of roughly 1 to 2.6
         orders of magnitude over the n=6, k=3 data.
      2. It started at 0 and used `if max_width == 0: max_width = ...`, which
         re-seeds the product rather than multiplying whenever a factor is 0.
         That silently discarded unsatisfiable edges instead of making the whole
         product 0, which is why sparse projections came out as 1 (10^0) rather
         than their true value. math.comb already returns 0 for deg < w, so
         multiplying from 1 gives the right answer with no special case.

    Note this is the same pair of mistakes, in the opposite direction, to the
    C++ compute_log_width_product, which used boost::degree and so counted the
    2-clique's self-loop TWICE.
    """
    if cliques is None:
        # Get the cliques
        cliques = get_cliques_nx(proj)
    smallest = max(2, min_k)
    max_width = 1
    # For each edge in the projection
    for i in range(proj.shape[0]):
        for j in range(i+1, proj.shape[1]):
            if proj[i,j] > 0:
                # |eta_e|: distinct cliques of size >= smallest containing i and j
                deg = sum(1
                          for k, kcliques in cliques.items() if k >= smallest
                          for c in kcliques
                          if i in c and j in c)
                max_width *= comb(deg, int(proj[i,j]))

    return max_width