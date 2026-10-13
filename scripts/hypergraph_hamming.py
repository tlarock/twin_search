"""Hypergraph Hamming and symmetric-difference distances.

Helper for parallel_distances.py. Self-contained: numpy and xgi only.
"""
import xgi
import numpy as np

def get_incidence_mat(n, hypergraph):
    incidence_mat = np.zeros((n, len(hypergraph)))
    for he_col, he in enumerate(hypergraph):
        for u in he:
            incidence_mat[u-1, he_col] = 1.0

    return incidence_mat

# Function that takes an incidence matrix and a hyperedge
# size and returns the corresponding tensor
# NOTE: This function does not construct a symmetric tensor,
# instead each hyperedge is counted exaclty once
def get_k_tensor(incidence_mat, k):
    # Get the columns with sum == k
    cols = incidence_mat[:, np.where(incidence_mat.sum(axis=0) == k)[0]]
    # Construct the tensor
    t = np.zeros([incidence_mat.shape[0]]*k)
    for col_idx in range(cols.shape[1]):
        c = cols[:, col_idx]
        t[tuple(np.where(c>0)[0])] = 1.0

    return t

# Compute a normalization factor as the total number
# of possible nodes that could be substituted based
# on the number of nodes in hyperedges of size k
def normalization(A, B, k):
    total_size_k = np.where(A.sum(axis=0) == k)[0].shape[0] + np.where(B.sum(axis=0) == k)[0].shape[0]
    return  total_size_k * k

def hypergraph_hamming_dist(n, HA, HB, normalized=False):
    # Need hypergraphs as lists-of-lists
    hypergraph_A = [sorted(list(HA.edges.members(he_idx))) for he_idx in HA.edges]
    hypergraph_B = [sorted(list(HB.edges.members(he_idx))) for he_idx in HB.edges]
    A = get_incidence_mat(n, hypergraph_A)
    B = get_incidence_mat(n, hypergraph_B)
    max_k = max([len(he) for he in hypergraph_A + hypergraph_B])
    tensor_hamming = 0.0
    for k in range(2, max_k+1):
        TA = get_k_tensor(A, k)
        TB = get_k_tensor(B, k)
        if not normalized:
            tensor_hamming += abs(TA - TB).sum()
        else:
            diff = abs(TA - TB).sum()
            norm = normalization(A, B, k)
            if diff > 0:
                tensor_hamming += diff / norm

    return tensor_hamming

def symmetric_difference(n, HA, HB, normalized=False):
    hypergraph_A = set([tuple(sorted(list(HA.edges.members(he_idx)))) for he_idx in HA.edges])
    hypergraph_B = set([tuple(sorted(list(HB.edges.members(he_idx)))) for he_idx in HB.edges])
    d = len(hypergraph_A.union(hypergraph_B) - hypergraph_A.intersection(hypergraph_B))
    if normalized:
        d /= len(hypergraph_A.union(hypergraph_B))
    return d
