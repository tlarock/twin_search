# NOTE: This relies on objlib being constructed in ../src/twin_search/CMakeLists.txt
enable_testing()

add_executable(test_utils test_utils.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_utils PRIVATE GTest::gtest_main)

add_executable(test_hypergraph test_hypergraph.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_hypergraph PRIVATE GTest::gtest_main)

add_executable(test_random_hypergraph_generators test_random_hypergraph_generators.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_random_hypergraph_generators PRIVATE GTest::gtest_main)

add_executable(test_factor_graph test_factor_graph.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_factor_graph PRIVATE GTest::gtest_main)


add_executable(test_twin_search test_twin_search.cpp ../src/twin_search/twin_search.cpp $<TARGET_OBJECTS:objlib>)
target_link_libraries(test_twin_search PRIVATE TBB::tbb ${Boost_LIBRARIES})
target_link_libraries(test_twin_search PRIVATE GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_hypergraph uniform_hypergraph, test_factor_graph, test_twin_search, test_utils)
