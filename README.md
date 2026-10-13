# twin_search: C++ library for enumerating sets of twins from hypergraph node co-occurrence adjacency matrices

Code for enumerating all hypergraphs that correspond to a weighted node co-occurence matrix, called the set of twins. The method was described in the following pre-print:

* Timothy LaRock & Renaud Lambiotte. Exploring the Non-uniqueness of Node Co-occurrence Matrices of Hypergraphs. June 2025. [arxiv:2506.01479](https://arxiv.org/abs/2506.01479). Accepted in Journal of Complex Networks, 2026.

# Zenodo data archive

The calculations are computationally intensive. An archive of the data used to
generate the figures in the paper is on Zenodo:

**DOI:** [10.5281/zenodo.23321964](https://doi.org/10.5281/zenodo.23321964) — the concept DOI, which always resolves to the latest version.

A single 841 MiB archive, `twin-search-reproducibility-core.tar`, reproduces
**all seven figures**. It holds the exhaustive enumerations, the sampled grids,
and the per-projection distance statistics that Figure 6 plots, with the
commit, seed and thread count each cell was produced with recorded alongside
it. It carries its own `README.md` describing the layout and file formats, and
a `MANIFEST.tsv` of SHA-256 checksums so the contents can be verified after
download.

The raw pairwise distance matrices that Figure 6's statistics were derived from
(29 GB compressed, 103 GB uncompressed) are **not** deposited: they are 97% of
the data by volume, and no notebook reads them — the figure is plotted from the
precomputed tables, which are in the archive. Contact us if you want them.

The notebooks in `analysis/notebooks/` read the archive directly and reproduce
each figure from it.

# Commit provenance

Code for reproducibility is stable since
[`7a44c3e`](https://github.com/tlarock/twin_search/commit/7a44c3e), the commit
that generated every dataset in the archive; each cell records that hash in its
own metadata. Notebooks for plotting the figures were added in
[`a913c49`](https://github.com/tlarock/twin_search/commit/a913c49), which
touched no compiled source. Later commits do not affect the generated data.

# Requirements
We assume C++20 or later and require `gcc`/`g++` 13 or higher.
 
The code relies heavily on [uxlfoundation/oneTBB](https://github.com/uxlfoundation/oneTBB) for parallel functionality. As `cmake` and `TBB` do not always play nice together, we recommend installing this requirement via `Conan2` (see below), rather than linking to your own install. If you choose to install it standalone, you need to point `cmake` to the `FindTBB.cmake` file provided in the `cmake/` directory of your `TBB` installation.

Similarly, the code requires `Boost`, specifically [Boost graph](https://www.boost.org/doc/libs/master/libs/graph/doc/index.html) and [uBLAS](https://www.boost.org/doc/libs/1_88_0/libs/numeric/ublas/doc/). As above, it is recommended to use `Conan2` to manage `Boost`, rather than linking to a separate installation.

The code expects that you have the [morrisfranken/argparse](https://github.com/morrisfranken/argparse) argument parsing repository cloned at the same level as this repository (e.g., in `../`). It is header-only, so no further installation is required.

# Compilation
Assuming the above requirements are in place, the recommended way to compile the code is to use the [Conan2](https://docs.conan.io/2/) package management system. We provide a `conanfile.txt` at the top level of this repository that specifies the dependencies.

Conan2 can be installed via python. You will need to set up a default conan profile and point it to an appropriate compiler (one that supports at least C++20 and implements `std::format`), which may require modifying `~/.conan2/profiles/default`.

With a profile configured, you can then run the following from the top-level directory of this repository to install requirements and generate makefiles:

`conan install . --output-folder=build --build=missing`

The output of the previous command should provide (1) a path to the appropriate `cmake` for compilation and (2) a `cmake` command you can run to generate the `cmake` build files. Enter the build directory (`cd build`) then identify and run the `cmake` command that usually looks something like:

`/path/to/conan/cmake/bin/cmake .. -G "Unix Makefiles" -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake  -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_BUILD_TYPE=Release`

Finally, you can run:

`/path/to/conan/cmake/build/bin/cmake --build .`

from within the `./build/` directory to actually compile the project.

If you opt not to use `Conan2`, the next best way to compile the code is to configure a `cmake` toolchain for your system. We provide an example toolchain file that should provide a starting point. Once configured, the steps are to create and enter a build directory (`mkdir build; cd build`), generate build files with cmake (something like `cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/your/toolchain.cmake ..`), the finally build (`cmake --build .` from within `build/` directory).  

Note: The main language feature that causes problems seems to be `std::format`. It is implemented in the latest versions of the major compilers, but is sometimes still finnicky. It should be relatively straightforward to remove the usages of `std::format` from the code if you can't use a complier that supports it. Feel free to open an issue or reach out to discuss if needed.

# Package Structure

```
.
|-- ./CMakeLists.txt
|-- ./conanfile.txt                 conan dependency manifest
|-- ./ExampleToolchain.cmake        starting point if you build without conan
|-- ./LICENSE
|-- ./README.md
|-- ./include/twin_search/
|   |-- combinations.hpp
|   |-- factor_graph.hpp
|   |-- hypergraph.hpp
|   |-- projected_graph.hpp
|   |-- random_hypergraph_generators.hpp
|   |-- search_checkpoint.hpp
|   |-- twin_search.hpp
|   |-- utils.hpp
|-- ./src/
|   |-- CMakeLists.txt
|   |-- ./src/twin_search/
|       |-- CMakeLists.txt
|       |-- twin_search.cpp                     the search itself
|       |-- hypergraph.cpp
|       |-- projected_graph.cpp
|       |-- factor_graph.cpp
|       |-- random_hypergraph_generators.cpp
|       |-- search_checkpoint.cpp
|       |-- utils.cpp
|       |-- count_twins_random.cpp              driver: sampled grids
|       |-- exhaustive_search_projections.cpp   driver: exhaustive enumeration
|       |-- exact_mates.cpp                     driver: exact mate pairs
|       |-- real_data_twins.cpp                 driver: empirical hypergraphs
|       |-- twin_frontier_probe.cpp             driver: search frontier size
|-- ./tests/
|   |-- CMakeLists.txt
|   |-- test_chung_lu_hypergraph.cpp
|   |-- test_combinations.cpp
|   |-- test_factor_graph.cpp
|   |-- test_hypergraph.cpp
|   |-- test_hypergraphs.cpp                    shared fixtures
|   |-- test_projected_graph.cpp
|   |-- test_random_hypergraph_generators.cpp
|   |-- test_twin_search.cpp
|   |-- test_uniform_hypergraph.cpp
|   |-- test_utils.cpp
|-- ./scripts/                      regenerate the datasets from the binaries
|   |-- README.md
|   |-- common.sh
|   |-- exhaustive_nonuniform.sh    exhaustive enumeration over m
|   |-- exact_mates.sh
|   |-- sample_heatmap.sh           sampled (n, m) grids
|   |-- sample_density.sh
|   |-- run_density_ordered.sh
|   |-- estimate_density_cost.py
|   |-- compute_distances.sh        pairwise twin-distance matrices
|   |-- convert_distance_matrices.py
|   |-- build_distance_stats.py
|   |-- matching_indices.py
|   |-- compress_results.py
|   |-- build_manifest.py           checksum manifest for the deposit
|   |-- verify_f32_matrices.py
|   |-- gen_golden.cpp
|   |-- exhaustive_row_counts.tsv
|-- ./analysis/
|   |-- ./analysis/notebooks/       one notebook per figure (2-7, A10)
|   |-- ./analysis/python/
|       |-- distance_stats.py       readers for the result files
|       |-- projection_utils.py
|       |-- compare_density_sweep.py
|       |-- nested_twin_reuse.py
|       |-- test_projection_utils.py
|-- ./docker/                       container used for the cluster runs
    |-- Dockerfile
    |-- run.sh
```

Each run under `scripts/` writes a `RUN_MANIFEST.tsv` recording the commit,
seed and exact commands it used. Generated results are not tracked in git; see
the Zenodo archive above.
