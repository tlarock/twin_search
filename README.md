# gram_mates
Experimental C++ library for constructing gram mates from hypergraph projection adjacency matrices.

_Definition_: Given a weighted node co-occurence matrix corresponding to a hypergraph, find all isomorphism classes of hypergraphs that correspond to the input matrix. Pairs of non-isomorphic hypergraphs that have the same node co-occurence matrix are called gram mates.

# Requirements
The code is written assuming C++20 or later and requires `gcc`/`g++` 13 or higher.
 
The code relies heavily on [uxlfoundation/oneTBB](https://github.com/uxlfoundation/oneTBB) for parallel functionality. As cmake and TBB do not always play nice together, I recommend installing this requirement via Conan (see below), rather than linking to your own install. If you choose to install it standalone, you need to point cmake to the `FindTBB.cmake` file provided in the `cmake/` directory of your installation.

We also use Boost, specifically [Boost graph](https://www.boost.org/doc/libs/master/libs/graph/doc/index.html) and [uBLAS](https://www.boost.org/doc/libs/1_88_0/libs/numeric/ublas/doc/). As above, it is recommended to use Conan2 to manage Boost, rather than linking to a separate installation.

The code expects that you have two repositories cloned at the same level as this repository (e.g., in ../). Both are header-only, so no further installation is required. The repositories are:
* [morrisfranken/argparse](https://github.com/morrisfranken/argparse)
* [mraggi/discreture](https://github.com/mraggi/discreture/tree/master)


# Compilation
Assuming the above requirements are in place, the recommended way to compile the code is to use the [Conan2](https://docs.conan.io/2/) package management system. We provide a `conanfile.txt` at the top level of this repository.

Conan2 can be installed via python. You will need to set up a default conan profile and point it to an appropriate compiler (one that supports at least C++20 and implements `std::format`), which may require modifying `~/.conan2/profiles/default`.

With a profile configured, you can then run `conan install . --output-folder=build --build=missing` to install requirements and generate makefiles. The output of the previous command will also provide a `cmake` command you can run to generate the `cmake` build, as well as the path to the conan-installed `cmake` that you can finally use to run `/path/to/conan/cmake/build/bin/cmake --build .` from within the ./build/ directory to actually compile the project.

If you opt not to use Conan2, the next best way to compile the code is to configure a cmake toolchain for your system, then install with cmake. We provide a few example toolchain files in this repository that should provide starting points. Once configured, the steps are to create and enter a build directory (`mkdir build; cd build`), generate build files with cmake (something like `cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/your/toolchain.cmake ..`), the finally build (`cmake --build .` from within `build/` directory).  


Note: The main language feature that causes problems seems to be `std::format`. It is implemented in the latest versions of the major compilers, but is sometimes still finnicky. 

# Package Structure

```
.
|-- ./CMakeLists.txt
|-- ./LICENSE
|-- ./include
|   |-- ./include/gram_mates/
|       |-- hypergraph.hpp
|       |-- projected_graph.hpp
|       |-- factor_graph.hpp
|       |-- gram_mates.hpp
|       |-- random_hypergraph_generators.hpp
|-- ./src
|   |-- ./src/CMakeLists.txt
|   |-- ./src/gram_mates/
|       |-- CMakeLists.txt
|       |-- hypergraph.cpp
|       |-- count_mates_random.cpp
|       |-- factor_graph.cpp
|       |-- projected_graph.cpp
|       |-- gram_mates.cpp
|       |-- random_hypergraph_generators.cpp
|       |-- test_parallel.cpp
|       |-- test_hypergraphs.cpp
|       |-- test_filter_isomorphic.cpp
|-- ./README.md
```
