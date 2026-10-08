import pickle
import gzip
import io
import os
import numpy as np


# The reproducibility dataset stores its CSVs compressed, so a caller asking
# for "...projections.csv" has to be handed "...projections.csv.zst". Keeping
# that fallback here means the notebooks only ever change a results_path.
_COMPRESSED_SUFFIXES = (".zst", ".gz")


def resolve_path(filepath):
    """-> the path that actually holds `filepath`'s data, or None."""
    filepath = str(filepath)
    if os.path.exists(filepath):
        return filepath
    for suffix in _COMPRESSED_SUFFIXES:
        if os.path.exists(filepath + suffix):
            return filepath + suffix
    return None


def data_exists(filepath):
    """Path.exists() for callers that should not care about compression.

    Several notebooks skip (n, m, k) cells that were never run, so they ask
    before reading. A plain .exists() answers False for every file in the
    compressed dataset, which would silently blank those cells rather than
    fail - so use this instead.
    """
    return resolve_path(filepath) is not None


def _open_zstd(path):
    """Open a .zst file as text.

    zstd is stdlib from Python 3.14 (`compression.zstd`). The notebook venvs
    here are 3.10, which needs the `zstandard` package - hence the fallback,
    and the explicit error rather than a bare ImportError if neither is there.
    """
    try:
        from compression import zstd
        return zstd.open(path, "rt")
    except ImportError:
        pass
    try:
        import zstandard
    except ImportError:
        raise ImportError(
            f"reading {path} needs zstd support: either Python >= 3.14 "
            f"(stdlib compression.zstd) or `pip install zstandard`"
        ) from None
    import io
    fh = open(path, "rb")
    reader = zstandard.ZstdDecompressor().stream_reader(fh)
    return io.TextIOWrapper(reader, encoding="utf-8")


# The stats tables get_stat_dist produces. min/max/entropy are created and
# never filled, so they carry no data and a CSV cannot hold a zero-length
# column next to populated ones; read_stats restores them as empty lists.
STATS_COLUMNS = ("m", "mean", "std", "max_width", "num_cliques", "num_size_dists")
STATS_EMPTY_COLUMNS = ("min", "max", "entropy")


def read_stats_file(path):
    """Load a stats table, from CSV if there is one, else from the pickle.

    Named ...__file because figure-7 uses `read_stats` and `write_stats` as
    boolean flags; a function of either name would be shadowed by them.

    Pickle ties the payload to a protocol version and, for numpy scalars, to
    the unpickling numpy; the CSV form written by
    validation/convert_stats_pickles.py does not. Pass either path - a
    ".pickle" is resolved to its ".csv.zst" when that exists, so callers can
    keep naming the pickle.

    Values come back with their original Python types. A token with no decimal
    point or exponent is an int, which matters because `std` holds a plain
    int 0 wherever a projection had a single pair.
    """
    if path.endswith(".pickle"):
        for alt in (path[:-len(".pickle")] + ".csv.zst",
                    path[:-len(".pickle")] + ".csv"):
            if os.path.exists(alt):
                path = alt
                break
        else:
            with open(path, "rb") as fin:
                return pickle.load(fin)

    stats = {c: [] for c in STATS_COLUMNS}
    stats.update({c: [] for c in STATS_EMPTY_COLUMNS})
    with open_text(path) as fin:
        header = next(fin).strip().split(",")
        for line in fin:
            line = line.strip()
            if not line:
                continue
            for col, tok in zip(header, line.split(",")):
                stats[col].append(
                    float(tok) if ("." in tok or "e" in tok or "E" in tok
                                   or "inf" in tok or "nan" in tok) else int(tok))
    return stats


def stats_is_int(v):
    """True for a Python or numpy integer.

    numpy.float64 subclasses float, so an isinstance check against float would
    call an integer column a float. The distinction matters: `std` holds a
    plain int 0 wherever a projection had a single pair.
    """
    return (isinstance(v, int) and not isinstance(v, bool)) or \
        type(v).__name__.startswith("int")


def write_stats_file(path, stats, level=19):
    """Write a stats table as CSV, zstd-compressed if `path` ends in .zst.

    Replaces pickling. Values go out through repr(), which for a float64 is the
    shortest string that reads back to identical bits, so this is lossless for
    everything except the numpy-ness of the scalars - and shedding that is the
    point, since the archive should not need numpy to be read.
    """
    for c in STATS_EMPTY_COLUMNS:
        if len(stats.get(c, [])) != 0:
            raise ValueError(
                f"column {c!r} has {len(stats[c])} values, but the CSV layout "
                f"drops it because it is empty in every published table. Add it "
                f"to STATS_COLUMNS before writing data into it."
            )
    n = len(stats[STATS_COLUMNS[0]])
    if any(len(stats[c]) != n for c in STATS_COLUMNS):
        raise ValueError(f"columns have unequal lengths: "
                         f"{ {c: len(stats[c]) for c in STATS_COLUMNS} }")

    buf = io.StringIO()
    buf.write(",".join(STATS_COLUMNS) + "\n")
    for row in zip(*(stats[c] for c in STATS_COLUMNS)):
        buf.write(",".join(repr(int(v)) if stats_is_int(v) else repr(float(v))
                           for v in row) + "\n")
    blob = buf.getvalue().encode()

    if path.endswith(".zst"):
        import subprocess
        blob = subprocess.run(["zstd", f"-{level}", "-q", "-c"],
                              input=blob, capture_output=True, check=True).stdout
        with open(path, "wb") as fout:
            fout.write(blob)
    else:
        with open(path, "wb") as fout:
            fout.write(blob)


def stats_path(results_dir, n, k, m_lo, m_hi, dist, match_type):
    """The stats table for one (m range, distance, match type), sans extension."""
    return (results_dir + f"n-{n}_m-{m_lo}-{m_hi}_k-{k}"
            f"_non-uniform_exhaustive_projections_{dist}_{match_type}_stats")


def open_binary(path):
    """Open a possibly-zstd-compressed file for reading as BYTES."""
    if path.endswith(".zst"):
        try:
            from compression import zstd
            return zstd.open(path, "rb")
        except ImportError:
            pass
        try:
            import zstandard
        except ImportError:
            raise ImportError(
                f"reading {path} needs zstd support: either Python >= 3.14 "
                f"(stdlib compression.zstd) or `pip install zstandard`"
            ) from None
        return zstandard.ZstdDecompressor().stream_reader(open(path, "rb"))
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def iter_pairwise_rows(results_dir, dist_filename, row_lengths):
    """Yield each projection's condensed pairwise-distance vector in turn.

    Prefers the float32 binary form written by
    validation/convert_distance_matrices.py, falling back to the original
    gzipped decimal text when it is not there, so both layouts read the same.

    The binary form stores no row delimiters: row r is the next
    total_twins[r] * (total_twins[r] - 1) / 2 values, which the caller knows
    from the exhaustive CSV it has already loaded. Rows are read one at a time
    because a single row at m=13 reaches 828 million values (3.3 GB), so the
    file must never be materialised whole.
    """
    stem = dist_filename[:-len(".csv.gz")] if dist_filename.endswith(".csv.gz") \
        else dist_filename
    binary = results_dir + stem + ".f32.zst"
    if os.path.exists(binary):
        with open_binary(binary) as fin:
            for want in row_lengths:
                need = want * 4
                buf = bytearray()
                while len(buf) < need:
                    block = fin.read(need - len(buf))
                    if not block:
                        break
                    buf += block
                if len(buf) < need:
                    if not buf:
                        # The file ended cleanly on a row boundary: it simply
                        # holds fewer projections than the exhaustive CSV. Some
                        # PUBLISHED distance files are short like this - m=8
                        # hyperNS stops one projection early - so stopping here
                        # is what reproduces the original reader, which just ran
                        # out of lines. Anything already yielded was complete.
                        return
                    raise EOFError(
                        f"{binary} ends mid-projection: wanted {want} values, "
                        f"got {len(buf) // 4}. A file that merely holds fewer "
                        f"projections ends on a row boundary, so this one is "
                        f"damaged rather than short."
                    )
                yield np.frombuffer(bytes(buf), dtype="<f4")
        return
    with gzip.open(results_dir + dist_filename, "rb") as fin:
        for line in fin:
            yield np.array(list(map(float, line.decode().strip().split(","))))


def open_text(filepath):
    """Open a results file for reading as text, decompressing if needed."""
    path = resolve_path(filepath)
    if path is None:
        raise FileNotFoundError(filepath)
    if path.endswith(".zst"):
        return _open_zstd(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


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

def read_data(filepath, old_file_format=False, xgi_hypergraphs=True):
    """Read an exhaustive-search or sampled CSV from a path.

    .zst and .gz are handled transparently, so the same call works against the
    reproducibility dataset and against a loose directory of plain CSVs.

    `old_file_format` selects the pre-TwinSearch row layout used by the
    increasing_density samples-500 files that Figure 7's line panel and Figure
    A.10 read. It lived only as a duplicated cell inside those two notebooks;
    it is here so the dataset ships with one reader rather than three.
    """
    with open_text(filepath) as fin:
        if old_file_format:
            return read_data_old_lines(fin)
        return read_data_lines(fin, xgi_hypergraphs=xgi_hypergraphs)


def read_sampled_cell(filepath, expected_samples, old_file_format=False,
                      xgi_hypergraphs=False, min_samples=None):
    """Read one sampled (n, m, k) cell, refusing to return a censored one.

    Returns (data, info). `info` always has path/rows/expected/status, where
    status is one of:

        "ok"       rows >= min_samples; `data` is the parsed cell
        "short"    the file exists but holds fewer rows; `data` is None
        "missing"  nothing on disk;                      `data` is None

    A cell cut short by a wall-clock limit is CENSORED, not merely small. The
    samples that finished are the fast ones, and runtime correlates with the
    twin count (rho 0.67-0.81 in this data), so a mean taken over "whatever
    rows happen to be present" is biased LOW - by up to 39% in the cells
    measured here. Five cells in the published k=3 n=9 sweep hold a single
    degenerate 19-byte record carrying no twin data at all.

    That is the failure this exists to prevent: the previous notebook code
    called num_iso_classes.mean() on whatever it found, so a censored cell
    produced an ordinary-looking point with no visual cue that it was wrong.
    The only reason the five degenerate cells did not plot as zeros is that a
    later `mean_isos_m[mean_isos_m == 0] = np.nan` happened to catch them,
    which is luck rather than a check.

    `min_samples` defaults to `expected_samples`, i.e. accept only complete
    cells. Lower it deliberately to admit reduced-N cells that were RUN TO
    COMPLETION - those are unbiased, just wider - and plot their N so the
    reader can see which points rest on fewer samples. Never lower it to admit
    a cell that was stopped early; that is the biased case.
    """
    if min_samples is None:
        min_samples = expected_samples
    info = {"path": str(filepath), "rows": 0,
            "expected": expected_samples, "status": "missing"}
    if not data_exists(filepath):
        return None, info
    data = read_data(filepath, old_file_format=old_file_format,
                     xgi_hypergraphs=xgi_hypergraphs)
    info["rows"] = len(data["n"])
    if info["rows"] < min_samples:
        info["status"] = "short"
        return None, info
    info["status"] = "ok"
    return data, info


def report_excluded_cells(infos, label=""):
    """Print the cells a figure dropped, so exclusions are visible not silent.

    Takes the `info` dicts from read_sampled_cell. Missing cells are usually
    expected (m*k < n cannot span n nodes, and m = C(n,k) is a single
    hypergraph); SHORT cells are the ones worth looking at, because each is a
    point someone intended to have and does not.
    """
    short = [i for i in infos if i["status"] == "short"]
    missing = [i for i in infos if i["status"] == "missing"]
    head = f"[{label}] " if label else ""
    print(f"{head}{sum(1 for i in infos if i['status'] == 'ok')} cells used, "
          f"{len(short)} censored and excluded, {len(missing)} absent")
    for i in sorted(short, key=lambda d: -d["rows"]):
        print(f"    CENSORED {os.path.basename(i['path'])}: "
              f"{i['rows']}/{i['expected']} rows")
    return short, missing


def _empty_exhaustive_data():
    return {
        "n": [],
        "m": [],
        "runtime": [],
        "max_width":  [],
        "num_edges":  [],
        "num_cliques": [],
        "total_twins": [],
        "total_filtered": [],
        "total_mates": [],
        "dist_results_dict": [],
        "hypergraphs": []
    }


def read_data_old_lines(lines):
    """Parse the old sampled-results layout: seven ints, then the hypergraph.

    The layout is fixed by the writer this data came from - count_mates_random
    before commit 6cd48c5 ("Standardized output format"), at
    src/gram_mates/count_mates_random.cpp:125:

        n, m, num_unfiltered, num_filtered, runtime, max_log_width,
        num_cliques, <hyperedges separated by "|">

    An earlier version of this function unpacked only six fields and called the
    sixth `num_cliques`, so it labelled max_log_width as num_cliques, dropped
    the real num_cliques, and reported max_width as absent when these rows do
    carry it. The notebook cells it replaced had the identical bug. Only
    total_twins and total_filtered were ever read downstream, so no published
    figure was affected - but it was a trap for anyone reaching for the rest.

    `runtime` is in milliseconds. `num_edges` and `total_mates` really are
    absent from this layout and come back as -1 rather than 0, so a caller that
    plots them gets something obviously wrong instead of a plausible zero.
    """
    exhaustive_data = _empty_exhaustive_data()
    for line in lines:
        split_vec = line.strip().split(",")
        (n, m, total_twins, total_filtered,
         runtime, max_log_width, num_cliques) = list(map(int, split_vec[0:7]))
        exhaustive_data["n"].append(n)
        exhaustive_data["m"].append(m)
        exhaustive_data["runtime"].append(runtime)
        exhaustive_data["max_width"].append(max_log_width)
        exhaustive_data["num_edges"].append(-1)
        exhaustive_data["num_cliques"].append(num_cliques)
        exhaustive_data["total_twins"].append(total_twins)
        exhaustive_data["total_filtered"].append(total_filtered)
        exhaustive_data["total_mates"].append(-1)
        exhaustive_data["dist_results_dict"].append(dict())
        exhaustive_data["hypergraphs"].append([])
    return exhaustive_data


def read_data_lines(lines, xgi_hypergraphs=True):
    """Same as read_data, but over any iterable of lines.

    Split out so the data can be read straight from a tar member, a gzip
    stream or anything else line-iterable without first writing it to disk.
    See archive_io.py, which uses this to read the results tarballs without
    unpacking them.
    """
    exhaustive_data = {
        "n": [],
        "m": [],
        "runtime": [],
        "max_width":  [],
        "num_edges":  [],
        "num_cliques": [],
        "total_twins": [],
        "total_filtered": [],
        "total_mates": [],
        "dist_results_dict": [],
        "hypergraphs": []
    }
    for line in lines:
        # Split on slash
        split_vec = line.strip().split("/")

        # The first 7 columns are just CSV. Rows written after the phase-timing
        # change carry three more (ms_traversal, ms_mates, ms_iso), so slice
        # rather than unpack the whole list - this reads both formats.
        data = list(map(int, split_vec[0].strip().split(",")))
        (n, m, runtime, log_max_width, num_edges, num_cliques, total_mates) = data[:7]
        (ms_traversal, ms_mates, ms_iso) = data[7:10] if len(data) >= 10 else (-1, -1, -1)

        # The last column is the hypergraphs
        if len(split_vec[-1]) > 0:
            hypergraphs = parse_hypergraphs_string(split_vec[-1], xgi_hypergraphs=xgi_hypergraphs)
        else:
            hypergraphs = []

        # The middle are the size distribution statistics
        idx = 1
        res_dict = dict()
        total_twins = 0
        total_filtered = 0
        while idx < len(split_vec)-1:
            # Get the size distribution as a string
            size_dist_str = split_vec[idx]
            # Get list of strings representing size distribution
            dist_pair_strs = size_dist_str.split(",")

            # Parse the size distribution
            dist = dict()
            for dist_pair_str in dist_pair_strs:
                knum = dist_pair_str.split(":")
                if len(knum) > 0:
                    (k, num) = list(map(int, knum))
                    dist[k] = num

            # Now get the stats
            stats_strs = split_vec[idx+1].split(",")
            (num_unfilt, num_filt) = list(map(int, stats_strs))
            hashable_dist = tuple(sorted(dist.items()))
            if not hashable_dist in res_dict:
                res_dict[hashable_dist] = {
                    "unfiltered":num_unfilt,
                    "filtered":num_filt
                }
            else:
                res_dict[hashable_dist]["unfiltered"] += num_unfilt
                res_dict[hashable_dist]["filtered"] += num_filt

            total_twins += num_unfilt
            total_filtered += num_filt
            idx += 2

        exhaustive_data["n"].append(n)
        exhaustive_data["m"].append(m)
        exhaustive_data["runtime"].append(runtime)
        exhaustive_data["max_width"].append(log_max_width)
        exhaustive_data["num_edges"].append(num_edges)
        exhaustive_data["num_cliques"].append(num_cliques)
        exhaustive_data["total_twins"].append(total_twins)
        exhaustive_data["total_filtered"].append(total_filtered)
        exhaustive_data["total_mates"].append(total_mates)
        exhaustive_data["dist_results_dict"].append(dict(res_dict))
        exhaustive_data["hypergraphs"].append(hypergraphs)

    return exhaustive_data

def get_stat_dist(n, k, output_datas, distance_name, match_type, m_vals, results_dir,
                  dist_dir=None):
    stats = {
            "min": [],
            "max": [],
            "mean": [],
            "std": [],
            "entropy": [],
            "max_width": [],
            "num_cliques": [],
            "num_size_dists": [],
            "m": []
    }

    for m in m_vals:
        print(f"m: {m}")
        dist_filename = f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections_{distance_name}.csv.gz"
        if match_type in ["m", "diag"]:
            matching_indices = read_refined_sets(results_dir, f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections_matching_{match_type}_indices.csv")
        elif match_type == "both":
            m_indices = read_refined_sets(results_dir, f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections_matching_m_indices.csv")
            diag_indices = read_refined_sets(results_dir, f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections_matching_diag_indices.csv")
            matching_indices = []
            for proj_idx in range(len(m_indices)):
                matching_indices.append(sorted(list(set(m_indices[proj_idx]).intersection(diag_indices[proj_idx]))))

        # One projection at a time (bad for cpu, good for memory). Row lengths
        # are N(N-1)/2 from total_twins, which is what lets the float32 form
        # dispense with a row index entirely.
        row_lengths = [t * (t - 1) // 2 for t in output_datas[m]["total_twins"]]
        # The distance matrices are far bigger than everything else, so the
        # reproducibility dataset keeps them in their own directory. Default to
        # results_dir so a flat directory of published files still works.
        matrices_dir = dist_dir if dist_dir is not None else results_dir
        if True:
            # each row holds the upper triangular of a distance matrix for a single projection
            for proj_idx, all_pairwise_dists in enumerate(
                    iter_pairwise_rows(matrices_dir, dist_filename, row_lengths)):
                if match_type == "all":
                    # If we are comparing all, we don't need to do any filtering
                    pairwise_dists = all_pairwise_dists
                else:
                    # Otherwise, we need to filter to only the indices we are interested in
                    indices = set(matching_indices[proj_idx])
                    l = all_pairwise_dists.shape[0]
                    num_hypergraphs = get_N(l)
                    pairwise_dists = []
                    p_idx = 0
                    for i in range(num_hypergraphs):
                        for j in range(i+1, num_hypergraphs):
                            if i in indices and j in indices:
                                pairwise_dists.append(all_pairwise_dists[p_idx])

                            p_idx += 1

                    pairwise_dists = np.array(pairwise_dists)

                if len(pairwise_dists) > 0:
                    if len(pairwise_dists) == 1:
                        stats["mean"].append(pairwise_dists[0])
                        stats["std"].append(0)
                    else:
                        stats["mean"].append(pairwise_dists.mean())
                        stats["std"].append(pairwise_dists.std())

                    stats["max_width"].append(output_datas[m]["max_width"][proj_idx])
                    stats["num_cliques"].append(output_datas[m]["num_cliques"][proj_idx])
                    stats["num_size_dists"].append(len(output_datas[m]["dist_results_dict"][proj_idx]))
                    stats["m"].append(m)
                else:
                    # Defaults if there were no distances?
                    pass

                del pairwise_dists, all_pairwise_dists

    return stats

# Silly function to recover the original array dimension
# there is probably a closed form way...
def get_N(l):
    x = 0
    N = 1
    while x < l:
        x += N
        N +=1
    return N

# Read a pairwise distance file and construct a vector of dictionaries
# with sorted keys corresponding to a triangular of the distance matrix
def read_pairwise_distances(results_dir, dist_filename):
    distance_dicts = []
    with gzip.open(results_dir + dist_filename, "rb") as fin:
        # Each row is one projection
        for line in fin:
            # Read the line of pairwise distances for the twins associated to this projection
            pairwise_dists = np.array(list(map(float, line.decode().strip().split(","))))
            l = pairwise_dists.shape[0]
            num_hypergraphs = get_N(l)
            p_idx = 0
            ddict = dict()
            for i in range(num_hypergraphs):
                for j in range(i+1, num_hypergraphs):
                    ddict[(i,j)] = pairwise_dists[p_idx]
                    p_idx += 1
            distance_dicts.append(dict(ddict))

    return distance_dicts

def read_refined_sets(results_dir, set_filename):
    refined_sets = []
    with open_text(results_dir + set_filename) as fin:
        for line in fin:
            refined_sets.append(sorted(list(map(int, line.strip().split(",")))))
    return refined_sets

def write_new_stats(k, n, m_vals, dist, match_type, output_datas, results_dir,
                    dist_dir=None):
    """Compute the plotting stats from the big distance files and write them out.

    Massively speeds up plotting: the distance matrices are hundreds of GB, the
    stats are a megabyte.

    Writes CSV rather than a pickle. The old .pickle form tied the payload to a
    protocol version and to the unpickling numpy, which is a poor property for
    something meant to be re-read years later; read_stats_file loads either.
    """
    stats = get_stat_dist(n, k, output_datas, dist, match_type, m_vals,
                          results_dir, dist_dir=dist_dir)
    write_stats_file(
        stats_path(results_dir, n, k, min(m_vals), max(m_vals), dist, match_type)
        + ".csv.zst", stats)


def update_stats(k, n, old_max, m_vals, dist, match_type, output_datas, results_dir,
                 dist_dir=None):
    """Extend an existing table for min(m_vals)..old_max with old_max+1..max(m_vals).

    Reads whatever is there - CSV or the original pickle - and writes CSV.
    """
    stats = read_stats_file(
        stats_path(results_dir, n, k, min(m_vals), old_max, dist, match_type)
        + ".pickle")

    ms_to_update = list(range(old_max + 1, max(m_vals) + 1))
    new_stats = get_stat_dist(n, k, output_datas, dist, match_type, ms_to_update,
                              results_dir, dist_dir=dist_dir)

    for key in stats:
        stats[key].extend(new_stats[key])

    write_stats_file(
        stats_path(results_dir, n, k, min(m_vals), max(m_vals), dist, match_type)
        + ".csv.zst", stats)


# The old names, kept so anything outside this repo still runs. They now write
# CSV, not pickles, so the names are wrong - prefer the ones above.
write_new_pickle = write_new_stats
update_pickle = update_stats

def read_output_datas(k, n, m_vals, results_dir):
    output_datas = dict()
    for m in m_vals:
        filename = f"n-{n}_m-{m}_k-{k}_non-uniform_exhaustive_projections.csv"
        output_datas[m] = read_data(results_dir + filename, xgi_hypergraphs=False)

    return output_datas

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="Order k of input file.",
                        type=int, default=3)
    parser.add_argument("-n", help="Number of nodes for input file.",
                        type=int, default=6)
    parser.add_argument("--min-m", help="Minimum m of hyperedges for input file.",
                        type=int, default=2)
    parser.add_argument("--max-m", help="Maximum m of hyperedges for input file.",
                        type = int, default=13)
    parser.add_argument("--results-dir", help="Location of input file.",
                        type=str, default="../updated_increasing_density/")
    parser.add_argument("--match-type", help="Match type.",
                        type=str, default="all")
    parser.add_argument("--distance", help="Distance.",
                        type=str, default="jaccard")
    args = parser.parse_args()

    results_dir = args.results_dir
    k = args.k
    n = args.n
    min_m = args.min_m
    max_m = args.max_m
    m_vals = list(range(min_m, max_m+1))

    output_datas = read_output_datas(k, n, m_vals, results_dir)

    match_type = args.match_type
    dist = args.distance

    print(f"Writing new pickle for distance {dist}, match type {match_type}, m in {min_m}...{max_m}")
    write_new_stats(k, n, m_vals, dist, match_type, output_datas, results_dir)
