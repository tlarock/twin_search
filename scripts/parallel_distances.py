"""Pairwise twin-distance worker: one (n, m, k) cell, one measure, in parallel.

Driven by scripts/compute_distances.sh, which is how the distance matrices
behind Figure 6 were produced. Runnable directly too - see --help.

Computes hyperNS, hyperPD, jaccard and symmetric-difference distances between
every pair of twins in a projection's twin set, writing one gzipped row per
projection.

The hyperNS and hyperPD measures come from the third-party higher-order
dissimilarity package of Agostinelli, Mancastroppa & Barrat
(arXiv:2503.16959), which is NOT vendored here - see the import below for how
to point at it. The exact commit used for the published matrices is recorded
in the data archive at provenance/HOR_COMMIT.
"""
import gzip
import numpy as np
from scipy.spatial.distance import canberra
from scipy.stats import entropy
import xgi

from hypergraph_hamming import *

import os
import sys
# The distance functions live in a sibling checkout. That repo was
# reorganised after the published runs - its HEAD commit is literally
# "Delete hypergraph_distances.py" - and the module now sits under code/, so
# the original single path no longer imports. Try both, and say which are
# missing rather than failing with a bare ModuleNotFoundError.
_HOR = os.environ.get("HOR_PATH") or os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..",
    "Hor_dissimilarity_measures")
for _cand in (os.path.join(_HOR, "code"), _HOR):
    if os.path.isdir(_cand):
        sys.path.append(_cand)
try:
    from hypergraph_distances import feature_vec, hyper_portrait_divergence, hyperedge_portrait
except ModuleNotFoundError:
    sys.exit(
        "hypergraph_distances not found.\n\n"
        "This needs the third-party higher-order dissimilarity measures of\n"
        "Agostinelli, Mancastroppa & Barrat (arXiv:2503.16959), which are NOT\n"
        "vendored here. Clone that repository beside this one, or point\n"
        "HOR_PATH at your checkout:\n\n"
        "    HOR_PATH=/path/to/Hor_dissimilarity_measures " + sys.argv[0] + " ...\n\n"
        "The exact commit the published matrices were computed with is in the\n"
        "data archive at provenance/HOR_COMMIT.\n\n"
        f"looked in: {_HOR}")

def parse_hypergraphs_string(hypergraphs_str, xgi_hypergraphs=False):
    hypergraphs = []
    # Split up hypergraphs
    hypergraphs_strs = hypergraphs_str.strip().split(";")
    for hg_str in hypergraphs_strs:
        hypergraph = []
        # split on hyperedges
        for he_str in hg_str.split("|"):
            # parse hyperedge
            he = tuple(map(int, he_str.split(":")))
            hypergraph.append(he)
        if xgi_hypergraphs:
            hypergraphs.append(xgi.Hypergraph(hypergraph))
        else:
            hypergraphs.append(hypergraph)

    return hypergraphs

def read_hypergraphs(filepath, xgi_hypergraphs=True):
    all_hypergraphs = []
    with open(filepath, "r") as fin:
        for line in fin:
            # Split on slash
            split_vec = line.strip().split("/")

            # Hypergraphs are the last string
            hypergraphs = parse_hypergraphs_string(split_vec[-1], xgi_hypergraphs=xgi_hypergraphs)

            # Add to all hypergraphs
            all_hypergraphs.append(hypergraphs)

    return all_hypergraphs

def read_one_twinset(filepath, row_idx, xgi_hypergraphs=True):
    with open(filepath, "r") as fin:
        curr_row = 0
        for line in fin:
            if curr_row == row_idx:
                # Split on slash
                split_vec = line.strip().split("/")

                # Hypergraphs are the last string
                hypergraphs = parse_hypergraphs_string(split_vec[-1], xgi_hypergraphs=xgi_hypergraphs)
                return hypergraphs
            curr_row += 1

    return curr_row

def count_rows(path):
    """Projections already written to a distance file. 0 if it is not there."""
    if not os.path.exists(path):
        return 0
    n = 0
    with gzip.open(path, "rt") as fin:
        for _ in fin:
            n += 1
    return n


def compute_jaccard(A, B):
    Aset = set([tuple(sorted(A.edges.members(he_idx))) for he_idx in A.edges])
    Bset = set([tuple(sorted(B.edges.members(he_idx))) for he_idx in B.edges])
    return len(Aset.intersection(Bset)) / (len(Aset.union(Bset)))

def get_gzip_str(input_str):
    return gzip.compress(input_str.encode("utf-8"))

def compute_dist(dist, hypergraphs, feature_vectors, hyper_portraits, i, j):
    if dist == "hyperNS":
        val = canberra(feature_vectors[i], feature_vectors[j]) / len(feature_vectors[i])
    elif dist == "hyperPD":
        val = hyper_portrait_divergence(hyper_portraits[i], hyper_portraits[j])
    elif dist == "jaccard":
        val = compute_jaccard(hypergraphs[i], hypergraphs[j])
    elif dist == "symmetricDiff":
        val = symmetric_difference(n, hypergraphs[i], hypergraphs[j])

    return val

def run_one_proj(dist, n, filepath, proj_idx):
    hypergraphs = read_one_twinset(filepath, proj_idx)
    num_twins = len(hypergraphs)
    # Pre-compute feature vectors/portraits if needed
    feature_vectors = dict()
    hyper_portraits = dict()
    if dist == "hyperNS" or dist == "hyperPD":
        for i in range(num_twins):
            if dist == "hyperNS":
                feature_vectors[i] = feature_vec(hypergraphs[i])
            if dist == "hyperPD":
                hyper_portraits[i] = hyperedge_portrait(hypergraphs[i])

    # Open a temporary output file
    with gzip.open(filepath + f"_{dist}_{proj_idx}.tmp.gz", "wt") as tmp_out:
        for i in range(num_twins):
            for j in range(i+1, num_twins):
                val = compute_dist(dist, hypergraphs, feature_vectors, hyper_portraits, i, j)
                if j < num_twins-1 or (j == num_twins-1 and (i+1) < num_twins-1):
                    tmp_out.write(str(val) + ",")
                elif j == num_twins-1 and (i+1) == num_twins-1:
                    tmp_out.write(str(val))

if __name__ == "__main__":
    from multiprocessing import Pool
    from functools import partial
    import tqdm
    import warnings
    warnings.filterwarnings("ignore")
    print("Warning: This script suppresses all warnings.")

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="Order k of input file.",
                        type=int, default=3)
    parser.add_argument("-n", help="Number of nodes for input file.",
                        type=int, default=6)
    parser.add_argument("-m", help="Number of hyperedges for input file.",
                        type=int, default=2)
    parser.add_argument("--results-dir", help="Location of input file.",
                        type=str, default="../updated_increasing_density/")
    parser.add_argument("--output-dir", help="Location to write output files. Defaults to input file location.",
                        type=str, default="../updated_increasing_density/")
    parser.add_argument("--hyperNS", help="Compute hyper NetSimile",
                        action="store_true")
    parser.add_argument("--jaccard", help="Compute Jaccard index.",
                       action="store_true")
    parser.add_argument("--hyperPD", help="Compute hyper Portrait Divergence.",
                        action="store_true")
    parser.add_argument("--symmetricDiff", help="Compute symmetric difference distance.",
                        action="store_true")
    parser.add_argument("--nprocs", help="Number of processors. Default 8.",
                        type=int, default=8)
    parser.add_argument("--resume", help="Keep the existing output file and compute "
                        "only the projections missing from its end. Without this the "
                        "output is truncated and everything is recomputed.",
                        action="store_true")
    parser.add_argument("--only", help="Compute only these projection indices and "
                        "append them, in the order given. For filling a known gap.",
                        nargs="+", type=int, default=None)

    args = parser.parse_args()
    k = args.k
    n = args.n
    m = args.m
    nprocs = args.nprocs
    results_dir = args.results_dir
    output_dir = args.output_dir
    dist = None
    if args.hyperNS:
        dist = "hyperNS"
    elif args.hyperPD:
        dist = "hyperPD"
    elif args.jaccard:
        dist = "jaccard"
    elif args.symmetricDiff:
        dist = "symmetricDiff"
    else:
        print("Exactly one distance must be specified. Exiting.")
        sys.exit(1)

    if sum([args.hyperNS, args.hyperPD, args.jaccard, args.symmetricDiff]) > 1:
        print("Exactly one distance must be specified. Exiting.")
        sys.exit(1)

    exhaustive_filename = f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections.csv"
    dist_filename = f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections_{dist}.csv.gz"

    print(f"Starting pairwise distance computation for k={k}, n={n}, m={m}.", flush=True)
    print(f"Input file: {results_dir + exhaustive_filename}", flush=True)

    # Get the number of lines in the input file
    num_projections = read_one_twinset(results_dir + exhaustive_filename, -1)

    # Each projection is one line, written in order, so an interrupted run
    # leaves a valid PREFIX rather than a corrupt file - which is how several
    # of the published distance files came to be short. That makes the missing
    # projections exactly the trailing ones, so they can be appended.
    if args.only is not None:
        todo = list(args.only)
        append = True
        print(f"Computing only projections {todo} and appending.", flush=True)
    elif args.resume:
        have = count_rows(output_dir + dist_filename)
        todo = list(range(have, num_projections))
        append = True
        print(f"Resuming: {have} of {num_projections} projections already present, "
              f"computing {len(todo)}.", flush=True)
        if not todo:
            print("Nothing to do.", flush=True)
            sys.exit(0)
    else:
        todo = list(range(num_projections))
        append = False
        print("Overwriting distance file at:", flush=True)
        print(output_dir + dist_filename, flush=True)
        with gzip.open(output_dir + dist_filename, "wt"):
            pass

    # Create a partial function so that the only input to run_one_proj is the
    # projection index/row in the input file
    pfunct = partial(run_one_proj, dist, n, results_dir + exhaustive_filename)
    args_list = todo
    print("Number of projections to compute: " + str(len(args_list)))
    # Build the result separately so we can have nice TQDM progress bar
    with Pool(nprocs) as p:
        # Our partial function doesn't return anything, it just writes the
        # pairwise distances for each projection to a temporary file
        list(tqdm.tqdm(p.imap_unordered(pfunct, args_list), total=len(args_list)))

    import os
    # Append, never truncate: a fresh run has already emptied the file above.
    with gzip.open(output_dir + dist_filename, "at") as fout:
        # Read all of the small temp files in order
        for proj_idx in args_list:
            # Open input file
            with gzip.open(results_dir + exhaustive_filename + f"_{dist}_{proj_idx}.tmp.gz", "rt") as tmp_in:
                for line in tmp_in:
                    # Concat directly to output file
                    fout.write(line + "\n")
            # Close and delete input file
            os.remove(results_dir + exhaustive_filename + f"_{dist}_{proj_idx}.tmp.gz")
