cmake_minimum_required(VERSION 3.10)

cmake_policy(SET CMP0135 NEW)

project(
    GramMates
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

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++23 -O3 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++23 -O0 -Wall -Wextra -Wpedantic -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-std=c++23 -O0 -Wall -Wextra -Wpedantic -g")

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

# NOTE: These are not managed by conan so must be included directly
include_directories("../argparse/include")
include_directories("../discreture/include")

# Add the src directory and subsequent CMakeLists
add_subdirectory(src)

# Add the test drirectory
add_subdirectory(tests)
