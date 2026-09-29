cmake_minimum_required(VERSION 3.10)

cmake_policy(SET CMP0135 NEW)

project(
    TwinSearch
    VERSION 0.0
    )

# Choose C++ standard libs and max optimization
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Build type Release")
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Build type Debug")
elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    message(STATUS "Build type RelWithDebugInfo")
endif()

# Declare the C++ standard as a CMake property rather than as an explicit -std
# flag in CMAKE_CXX_FLAGS_*. A toolchain file can set CMAKE_CXX_STANDARD itself
# - conan's sets it from the profile's compiler.cppstd - and CMake appends the
# matching -std flag AFTER CMAKE_CXX_FLAGS_*, so a hand-written -std there is
# silently overridden by the last flag on the command line. With the default
# conan profile (compiler.cppstd=gnu17) that produced
# "-std=c++23 ... -std=gnu++17", building this C++20/23 code as C++17.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# IMPORTANT: none of these define NDEBUG, and that is deliberate.
#
# Boost uBLAS bounds-checks every matrix access through BOOST_UBLAS_CHECK
# unless NDEBUG is defined, and that check is the last line of defence against
# an out-of-range write to a projection matrix. Adding -DNDEBUG - which is what
# CMake's own stock Release flags do, and is the obvious thing to "tidy up"
# here - silently turns every one of those checks off across the whole
# codebase, converting a loud bad_index exception into silent heap corruption.
# It was verified that an out-of-range node id throws at -O3 as written and
# corrupts memory at -O3 -DNDEBUG.
#
# If you want NDEBUG for speed, measure it first and add explicit bounds checks
# at the matrix accesses that matter, rather than relying on uBLAS.
#
# NOTE: these assignments also REPLACE the flags CMake would otherwise supply
# for each build type, so anything added here is the complete set. That is why
# sanitizer builds have to pass their flags through CMAKE_CXX_FLAGS with a
# build type this file does not name; see docker/run.sh.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -Wextra -Wpedantic -g")
# NOTE: -O0 here is almost certainly a mistake for a "with debug info" build,
# but changing it would silently change the performance of anyone using it.
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0 -Wall -Wextra -Wpedantic -g")

# Sets current directory to receive the executables
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR})

# Find TBB for parallel algorithms and
# Boost for graphs/ublas
# NOTE: linked in src/TwinSearch/CMakeLists
find_package(TBB REQUIRED)
find_package(Boost CONFIG 1.79.0 REQUIRED graph)

# Get googletests
include(FetchContent)
FetchContent_Declare(
  gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.15.2
)

FetchContent_MakeAvailable(gtest)

# Some boost settings recommended by someone on the internet
# Not sure if I need them but they are not causing problems
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

include_directories(${Boost_INCLUDE_DIRS})
include_directories("./include/twin_search/")

# NOTE: This is not managed by conan so must be included directly
include_directories("../argparse/include")

# Add the src directory and subsequent CMakeLists
add_subdirectory(src)

# Add the test drirectory
add_subdirectory(tests)
