# twin_search: C++ library for enumerating sets of twins from hypergraph node co-occurrence adjacency matrices

Code for enumerating all hypergraphs that correspond to a weighted node co-occurence matrix.

# Requirements
We assume C++20 or later and require `gcc`/`g++` 13 or higher.
 
The code relies heavily on [uxlfoundation/oneTBB](https://github.com/uxlfoundation/oneTBB) for parallel functionality. As `cmake` and `TBB` do not always play nice together, I recommend installing this requirement via `Conan2` (see below), rather than linking to your own install. If you choose to install it standalone, you need to point `cmake` to the `FindTBB.cmake` file provided in the `cmake/` directory of your `TBB` installation.

Similarly, the code requires `Boost`, specifically [Boost graph](https://www.boost.org/doc/libs/master/libs/graph/doc/index.html) and [uBLAS](https://www.boost.org/doc/libs/1_88_0/libs/numeric/ublas/doc/). As above, it is recommended to use `Conan2` to manage `Boost`, rather than linking to a separate installation.

The code expects that you have two repositories cloned at the same level as this repository (e.g., in ../). Both are header-only, so no further installation is required. The repositories are:
* [morrisfranken/argparse](https://github.com/morrisfranken/argparse)
* [mraggi/discreture](https://github.com/mraggi/discreture/tree/master)


# Compilation
Assuming the above requirements are in place, the recommended way to compile the code is to use the [Conan2](https://docs.conan.io/2/) package management system. We provide a `conanfile.txt` at the top level of this repository.

Conan2 can be installed via python. You will need to set up a default conan profile and point it to an appropriate compiler (one that supports at least C++20 and implements `std::format`), which may require modifying `~/.conan2/profiles/default`.

With a profile configured, you can then run
`conan install . --output-folder=build --build=missing`
to install requirements and generate makefiles. The output of the previous command will also provide a `cmake` command you can run to generate the `cmake` build, as well as the path to the conan-installed `cmake` that you can finally use to run
`/path/to/conan/cmake/build/bin/cmake --build .`
from within the `./build/` directory to actually compile the project.

If you opt not to use `Conan2`, the next best way to compile the code is to configure a `cmake` toolchain for your system. We provide an example toolchain file that should provide a starting point. Once configured, the steps are to create and enter a build directory (`mkdir build; cd build`), generate build files with cmake (something like `cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/your/toolchain.cmake ..`), the finally build (`cmake --build .` from within `build/` directory).  

Note: The main language feature that causes problems seems to be `std::format`. It is implemented in the latest versions of the major compilers, but is sometimes still finnicky. 

# Package Structure

```
.
|-- ./CMakeLists.txt
|-- ./LICENSE
|-- ./include
|   |-- ./include/twin_search/
|       |-- hypergraph.hpp
|       |-- projected_graph.hpp
|       |-- factor_graph.hpp
|       |-- twin_search.hpp
|       |-- random_hypergraph_generators.hpp
|-- ./src
|   |-- ./src/CMakeLists.txt
|   |-- ./src/twin_search/
|       |-- CMakeLists.txt
|       |-- hypergraph.cpp
|       |-- count_mates_random.cpp
|       |-- factor_graph.cpp
|       |-- projected_graph.cpp
|       |-- twin_search.cpp
|       |-- random_hypergraph_generators.cpp
|-- ./tests
|   |-- ./tests/CMakeLists.txt
|   |-- ./tests/test_factor_graph.cpp
|   |-- ./tests/test_hypergraph.cpp
|   |-- ./tests/test_hypergraphs.cpp
|   |-- ./tests/test_random_hypergraph_generators.cpp
|   |-- ./tests/test_twin_search.cpp
|   |-- ./tests/test_uniform_hypergraph.cpp
|   |-- ./tests/test_utils.cpp
|-- ./README.md
```
