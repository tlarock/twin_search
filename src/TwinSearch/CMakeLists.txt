# construct a cmake object library including everything except TwinSearch.cpp
# because it does not properly link oneTBB if I include it
# This allows object files to be created only once, without compiling into a library on disk
# TODO: For some reason, compiling twin_search.cpp into the objlib fails to link. It would be good to fix this if possible.
add_library(objlib OBJECT projected_graph.cpp hypergraph.cpp factor_graph.cpp random_hypergraph_generators.cpp utils.cpp)

add_executable(count_twins_random count_twins_random.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

add_executable(exhaustive_search_projections exhaustive_search_projections.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

add_executable(real_data_mates real_data_twins.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

add_executable(exact_mates exact_mates.cpp twin_search.cpp $<TARGET_OBJECTS:objlib>)

target_link_libraries(count_twins_random PRIVATE TBB::tbb ${Boost_LIBRARIES})
target_link_libraries(exhaustive_search_projections PRIVATE TBB::tbb ${Boost_LIBRARIES})
target_link_libraries(real_data_mates PRIVATE TBB::tbb ${Boost_LIBRARIES})
target_link_libraries(exact_mates PRIVATE TBB::tbb ${Boost_LIBRARIES})
