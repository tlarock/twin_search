# construct a cmake object library including everything except TwinSearch.cpp
# because it does not properly link oneTBB if I include it
# This allows object files to be created only once, without compiling into a library on disk
# TODO: For some reason, compiling twin_search.cpp into the objlib fails to link. It would be good to fix this if possible.
add_library(objlib OBJECT projected_graph.cpp hypergraph.cpp factor_graph.cpp random_hypergraph_generators.cpp utils.cpp)

add_executable(count_twins_random count_twins_random.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

add_executable(exhaustive_search_projections exhaustive_search_projections.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

add_executable(real_data_twins real_data_twins.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

add_executable(exact_mates exact_mates.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

# TBB::tbbmalloc_proxy replaces global new/delete with oneTBB's scalable
# allocator. It is a link-time change only; no code calls it directly.
#
# OFF by default because it is a measured trade, not a free win, and the
# direction of the trade depends on the platform. Measured on cell
# n=9 m=23 k=3 (macOS, 500 samples) and n=9 m=28 k=3 (Linux, 500 samples,
# 11 threads, all three builds back to back on one della node):
#
#   macOS / Apple clang : 9-10% SLOWER at 2, 4 and 6 threads.
#   Linux / gcc, glibc  : 2.8% slower, but 18% lower peak RSS
#                         (10.43 GB -> 8.54 GB).
#
# On macOS the proxy works by replacing the malloc zone; on Linux it uses
# symbol interposition, which is the mechanism it is actually designed for.
# The memory saving is most likely glibc's per-thread arenas retaining freed
# memory, which tbbmalloc's unified pools avoid.
#
# Turn this ON for memory-bound runs - the large-m cells are limited by peak
# RSS, not by CPU, so trading 2.8% of wall time for 18% of memory is usually
# the right call there. Leave it OFF for development and for macOS.
option(TWIN_SEARCH_SCALABLE_ALLOC
       "Link oneTBB's scalable allocator (tbbmalloc_proxy) into the drivers"
       OFF)

set(TWIN_SEARCH_DRIVERS
    count_twins_random
    exhaustive_search_projections
    real_data_twins
    exact_mates)

if(TWIN_SEARCH_SCALABLE_ALLOC AND NOT TARGET TBB::tbbmalloc_proxy)
    message(FATAL_ERROR
        "TWIN_SEARCH_SCALABLE_ALLOC=ON but the target TBB::tbbmalloc_proxy does "
        "not exist. The oneTBB package in use was built without its proxy "
        "component. Check that conan resolved onetbb with the proxy enabled - "
        "`grep COMPONENT_NAMES <build>/TBB-*-data.cmake` should list "
        "TBB::tbbmalloc_proxy alongside TBB::tbb.")
endif()

if(TWIN_SEARCH_SCALABLE_ALLOC)
    message(STATUS "Scalable allocator ON (tbbmalloc_proxy linked into drivers)")
else()
    message(STATUS "Scalable allocator OFF (system allocator; -DTWIN_SEARCH_SCALABLE_ALLOC=ON to enable)")
endif()

# Linked in a loop rather than one call per driver so that adding a fifth
# driver cannot silently pick up a different allocator from the other four.
foreach(driver IN LISTS TWIN_SEARCH_DRIVERS)
    target_link_libraries(${driver} PRIVATE TBB::tbb ${Boost_LIBRARIES})
    if(TWIN_SEARCH_SCALABLE_ALLOC)
        target_link_libraries(${driver} PRIVATE TBB::tbbmalloc_proxy)
    endif()
endforeach()
