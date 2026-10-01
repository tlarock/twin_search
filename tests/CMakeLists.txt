# NOTE: This relies on objlib being constructed in ../src/twin_search/CMakeLists.txt
# enable_testing() is called in the top-level CMakeLists; see the note there.

add_executable(test_utils test_utils.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_utils PRIVATE GTest::gtest_main)

add_executable(test_hypergraph test_hypergraph.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_hypergraph PRIVATE GTest::gtest_main)

add_executable(test_projected_graph test_projected_graph.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_projected_graph PRIVATE GTest::gtest_main)

add_executable(test_random_hypergraph_generators test_random_hypergraph_generators.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_random_hypergraph_generators PRIVATE GTest::gtest_main)

add_executable(test_factor_graph test_factor_graph.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_factor_graph PRIVATE GTest::gtest_main)

# Pins combinations_of() against reference data generated from
# discreture::combinations (see scripts/gen_golden.cpp), so the thread-safe
# replacement used by TwinSearch::get_combinations cannot change search results.
add_executable(test_combinations test_combinations.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_combinations PRIVATE GTest::gtest_main)


add_executable(test_twin_search test_twin_search.cpp ../src/twin_search/twin_search.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_twin_search PRIVATE TBB::tbb ${Boost_LIBRARIES})
target_link_libraries(test_twin_search PRIVATE GTest::gtest_main)

include(GoogleTest)

# Register every test binary with CTest.
#
# gtest_discover_tests() takes ONE target followed by keyword options, so the
# previous call - which listed six comma-separated names after the target -
# silently registered only test_hypergraph and discarded the rest as unparsed
# arguments. 49 of the 72 tests were invisible to `ctest` as a result, and
# test_projected_graph and test_random_hypergraph_generators were never in the
# list to begin with.
#
# Driving it from one list rather than seven hand-written calls means adding a
# test binary above and forgetting to register it here is not possible: the
# same name is used for both.
set(TWIN_SEARCH_TEST_TARGETS
    test_utils
    test_hypergraph
    test_projected_graph
    test_random_hypergraph_generators
    test_factor_graph
    test_combinations
    test_twin_search
)

foreach(test_target IN LISTS TWIN_SEARCH_TEST_TARGETS)
    if(NOT TARGET ${test_target})
        message(FATAL_ERROR
            "${test_target} is listed in TWIN_SEARCH_TEST_TARGETS but no such "
            "target exists. Add the add_executable() above, or remove it here.")
    endif()
    gtest_discover_tests(${test_target})
endforeach()
