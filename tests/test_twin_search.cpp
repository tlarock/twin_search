#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/assignment.hpp>
#include "projected_graph.hpp"
#include "factor_graph.hpp"
#include "utils.hpp"
#include "twin_search.hpp"
#include "test_hypergraphs.cpp"
#include <gtest/gtest.h>

namespace ublas=boost::numeric::ublas;

TEST(TwinSearchTest, IsomorphicCountTest) {
    Hypergraph H = h5();
    ProjectedGraph proj(H);
    
    // Test sequential 
    TwinSearch twins(proj, 3, 3, true, false, true, false);
    EXPECT_EQ(twins.twins.size(), 3);
    EXPECT_EQ(twins.filtered_twins.size(), 2);

    // Test parallel
    twins = TwinSearch(proj, 3, 3, true, true, true, false);
    EXPECT_EQ(twins.twins.size(), 3);
    EXPECT_EQ(twins.filtered_twins.size(), 2);

    // Second hypergraph
    H = h4();
    proj = ProjectedGraph(H);
    twins = TwinSearch(proj, 3, 3, true, false, true, false);
    EXPECT_EQ(twins.twins.size(), 1);
    EXPECT_EQ(twins.filtered_twins.size(), 1);

    twins = TwinSearch(proj, 3, 3, true, true, true, false);
    EXPECT_EQ(twins.twins.size(), 1);
    EXPECT_EQ(twins.filtered_twins.size(), 1);
}

TEST(TwinSearchTest, TestParallelIsomorphicFilter) {
    // Should keep only 1
    std::vector<Hypergraph> hypergraphs = {h1(), h1(), h1(), h1(), h1(), h1(), h1(), h1(), h1()};
    std::vector<UndirectedGraph> bipartites;
    for (Hypergraph h : hypergraphs) {
        bipartites.push_back(h.get_bipartite());
        EXPECT_EQ(h.n + h.m, boost::num_vertices(bipartites.back()));
    }

    std::vector<int> to_filter = TwinSearch::run_iso_tests_parallel(bipartites);
    int sum = 0;
    for(int i : to_filter)
        sum += i;
    EXPECT_EQ(sum, hypergraphs.size()-1);

    // Should remove 2 copies of h1, 1 copy of h2, 1 copy of h3
    hypergraphs = {h1(), h1(), h2(), h3(), h3(), h4(), h5(), h1(), h2()};
    bipartites = std::vector<UndirectedGraph>(0);
    for (Hypergraph h : hypergraphs) {
        bipartites.push_back(h.get_bipartite());
        EXPECT_EQ(h.n + h.m, boost::num_vertices(bipartites.back()));
    }

    to_filter = TwinSearch::run_iso_tests_parallel(bipartites);
    sum = 0;
    for(int i : to_filter)
        sum += i;
    EXPECT_EQ(sum, 4);
}

TEST(TwinSearchTest, TestMultiGraphIsomorphisms) {
    // Does the vf2 implementation correctly distinguish multigraphs?
    UndirectedGraph g1;
    boost::add_edge(0, 1, g1);
    boost::add_edge(1, 2, g1);
    boost::add_edge(2, 0, g1);

    UndirectedGraph g2;
    boost::add_edge(0, 1, g2);
    boost::add_edge(1, 2, g2);
    boost::add_edge(2, 0, g2);

    std::vector<UndirectedGraph> gs {g1, g2};
    std::vector<int> to_filter = TwinSearch::run_iso_tests_parallel(gs);
    int sum = 0;
    for (int i : to_filter)
        sum += i;
    EXPECT_EQ(sum, 1);

    // Add a multi-edge to g1 and test that they are not isomorphic
    boost::add_edge(0, 1, gs[0]);
    to_filter = TwinSearch::run_iso_tests_parallel(gs);
    sum = 0;
    for (int i : to_filter)
        sum += i;
    EXPECT_EQ(sum, 0);

    // Now add the same multi-edge to g2 and test that they are isomorphic
    boost::add_edge(0, 1, gs[1]);
    to_filter = TwinSearch::run_iso_tests_parallel(gs);
    sum = 0;
    for (int i : to_filter)
        sum += i;
    EXPECT_EQ(sum, 1);

    // Note that it will NOT reliably deal with self-loops
    boost::add_edge(0, 0, gs[0]);
    boost::add_edge(0, 0, gs[1]);
    to_filter = TwinSearch::run_iso_tests_parallel(gs);
    sum = 0;
    for (int i : to_filter)
        sum += i;
    // NOTE: This actually *should* be 1, since I've added the same self-loop
    // to a pair of isomorphic graphs. This is a bug in boost::vf2 using
    // undirectedS graph data structure.
    EXPECT_EQ(sum, 0);

}


TEST(TwinSearchTest, TestNonUniform) {
    // Kind of a useless test, but leaving here is harmless
    // and can be a baseline for a better test in the future
    Hypergraph h = h5();
    int min_k = 3;
    int max_k = h.n;
    ProjectedGraph proj(h);
    TwinSearch twins(proj, min_k, max_k, true, true, true, false);
    EXPECT_GE(twins.twins.size(), 0);
    EXPECT_GE(twins.filtered_twins.size(), 0);

    twins = TwinSearch(proj, min_k, max_k, true, false, true, false);
    EXPECT_GE(twins.twins.size(), 0);
    EXPECT_GE(twins.filtered_twins.size(), 0);

    // more useful test: A triangle should have two twins and
    // those twins should be non-isomorphic, since 1 is the
    // 3-hyperedge and the other is 3 2-hyperedges.
    h = h9();
    min_k = 2;
    max_k = h.n;

    proj = ProjectedGraph(h);
    twins = TwinSearch(proj, min_k, max_k, true, true, true, false);
    EXPECT_GE(twins.twins.size(), 2);
    EXPECT_GE(twins.filtered_twins.size(), 2);

    twins = TwinSearch(proj, min_k, max_k, true, false, true, false);
    EXPECT_GE(twins.twins.size(), 2);
    EXPECT_GE(twins.filtered_twins.size(), 2);
}

TEST(TwinSearchTest, TestLineGraphEquiv) {
    // very simple tests to check whether the static
    // function equivalent_lg works as intended
    std::vector<ublas::matrix<int> > line_graphs;
    ublas::matrix<int> m1(10, 10, 0); 
    ublas::matrix<int> m2(10, 10, 0);
    line_graphs.push_back(m1);
    line_graphs.push_back(m2);

    // Empty matrices of same dimension are equivalent
    EXPECT_EQ(TwinSearch::equivalent_lg(line_graphs, 0, 1), 1);

    // Same matrix with 1 entry, equivalent
    line_graphs[0](1, 1) = 1;
    line_graphs[1](1, 1) = 1;
    EXPECT_EQ(TwinSearch::equivalent_lg(line_graphs, 0, 1), 1);

    // Modify 1 entry, not equivalent
    line_graphs[0](1, 1) = 2;
    EXPECT_EQ(TwinSearch::equivalent_lg(line_graphs, 0, 1), 0);

    // Same modification to m2, equivalent
    line_graphs[1](1, 1) = 2;
    EXPECT_EQ(TwinSearch::equivalent_lg(line_graphs, 0, 1), 1);

    // Diffrent dimensions, always unequal
    ublas::matrix<int> m3(5, 5);
    m3(1, 1) = 2;
    EXPECT_EQ(TwinSearch::equivalent_lg(line_graphs, 0, 2), 0);
    EXPECT_EQ(TwinSearch::equivalent_lg(line_graphs, 1, 2), 0);
}


TEST(TwinSearchTest, TestGramMateExample) {
    // Test that sequential and parallel implementations give the same answer
    Hypergraph h = GM();
    ProjectedGraph proj(h);
    TwinSearch twins(proj, 2, h.n, true, true, true, false);
    TwinSearch seq_twins(proj, 2, h.n, true, false, true, false);

    EXPECT_EQ(twins.mates.size(), seq_twins.mates.size());
    EXPECT_EQ(twins.twins.size(), seq_twins.twins.size());
    EXPECT_EQ(twins.filtered_twins.size(), seq_twins.filtered_twins.size());
    
    // Check that equivalent_lg correctly identifies the known mate pair
    ublas::matrix<int> A(6,6);
    ublas::matrix<int> B(6,6);
    A <<= 1,1,0,0,0,0,
            1,1,1,1,0,0,
            1,0,1,1,1,0,
            0,1,1,1,1,0,
            0,0,1,0,1,1,
            0,0,0,1,1,1;

    B <<= 0,0,0,0,1,1,
           0,0,1,1,1,1,
           1,0,1,1,1,0,
           0,1,1,1,1,0,
           1,1,1,0,0,0,
           1,1,0,1,0,0;

    std::vector<ublas::matrix<int> > line_graphs {ublas::prod(A, ublas::trans(A)), ublas::prod(B, ublas::trans(B))};
    EXPECT_EQ(TwinSearch::equivalent_lg(line_graphs, 0, 1), 1);

    // Check that run_mates_tests correctly finds the known mate pair
    UndirectedGraph lgA;
    UndirectedGraph lgB;
    for (std::size_t r = 0; r < line_graphs[0].size1(); ++r) {
        for (std::size_t c = 0; c < line_graphs[0].size2(); ++c) {
            if (r == c)
                continue;

            if (line_graphs[0](r, c) > 0) {
                for (int v = 0; v < line_graphs[0](r,c); ++v) {
                    boost::add_edge(r, c, lgA);
                }
            }
            if (line_graphs[1](r, c) > 0) {
                for (int v = 0; v < line_graphs[1](r,c); ++v) {
                    boost::add_edge(r, c, lgB);
                }
            }
        }
    }

    std::vector<UndirectedGraph> lgs {lgA, lgB};
    std::vector<std::vector<int> > mates_pairs = TwinSearch::run_mates_tests_parallel(lgs);
    EXPECT_EQ(mates_pairs.size(), 1);

    // Check that run_iso_tests also correctly filters the known mate pair
    std::vector<int> to_filter = TwinSearch::run_iso_tests_parallel(lgs);
    int sum = 0;
    for (int v : to_filter)
        sum += v;
    EXPECT_EQ(sum, 1);

    // Check that run_mates_tests_parallel finds the right pairs
    std::vector<Hypergraph> hypergraphs = {h1(), h1(), h1(), h1(), h1(), h1(), h1(), h1(), h1()};
    lgs.clear();
    for (Hypergraph h : hypergraphs) {
        lgs.push_back(h.get_line_graph());
    }
    mates_pairs = TwinSearch::run_mates_tests_parallel(lgs);
    EXPECT_EQ(mates_pairs.size(), 36);

    // Should pair up the equals
    hypergraphs = {h1(), h1(), h2(), h3(), h3(), h4(), h5(), h1(), h2()};
    lgs = std::vector<UndirectedGraph>(0);
    for (Hypergraph h : hypergraphs) {
        lgs.push_back(h.get_line_graph());
    }
    mates_pairs = TwinSearch::run_mates_tests_parallel(lgs);
    EXPECT_EQ(mates_pairs.size(), 5);
}

// With use_diagonal false the constructor is supposed to ignore the diagonal by
// zeroing it. It guarded that on `proj_diag_sum > 1`, so a diagonal summing to
// exactly 1 survived; matsum(proj_rem) could then never reach 0 and the search
// reported no twins at all. Silently returning an empty twin set is the worst
// available failure mode, so this pins every diagonal sum around the boundary.
TEST(TwinSearchTest, NonZeroDiagonalIsIgnoredWhenUseDiagonalIsFalse) {
    std::vector<std::vector<int> > he = {{0,1,2},{1,2,3}};
    Hypergraph h(he);

    std::size_t expected = 0;
    for (int stray = 0; stray <= 3; stray++) {
        ProjectedGraph proj(h);              // zero diagonal
        if (stray > 0)
            proj.proj_mat(0, 0) = stray;

        TwinSearch gm(proj, 2, 4, false, false, true, false);   // use_diagonal false
        gm.search(false);

        // the diagonal must have been zeroed regardless of what it summed to
        for (std::size_t u = 0; u < gm.proj.proj_mat.size1(); ++u)
            EXPECT_EQ(gm.proj.proj_mat(u, u), 0) << "stray=" << stray << " u=" << u;

        if (stray == 0) {
            expected = gm.twins.size();
            EXPECT_GT(expected, 0u) << "fixture should produce twins";
        } else {
            EXPECT_EQ(gm.twins.size(), expected)
                << "diagonal summing to " << stray << " changed the twin count";
        }
    }
}

// The worst case search tree size is prod_e binom(|eta_e|, w_e). It must be
// computed over the SET of cliques that could satisfy each edge. It previously
// used boost::degree, which counts an edge-node's self-loop twice when
// min_k <= 2, so every published value was inflated.
//
// Worked example, small enough to check by hand: the single hyperedge {0,1,2}
// projects to a triangle with all weights 1. Cliques are the three 2-cliques
// and the one 3-clique. Each edge-node {u,v} therefore has
//     eta_e = { the 2-clique {u,v} itself, the triangle {0,1,2} },  |eta_e| = 2
// so the product is binom(2,1)^3 = 8 and its log10 is log10(8).
// With the old degree it was binom(3,1)^3 = 27.
TEST(TwinSearchTest, WidthProductUsesSetSizeNotBoostDegree) {
    std::vector<std::vector<int> > single_triangle = {{0, 1, 2}};
    Hypergraph h(single_triangle);
    ProjectedGraph proj(h);
    FactorGraph fact(proj, 2, 3);

    EXPECT_DOUBLE_EQ(static_cast<double> (TwinSearch::compute_width_product(proj, fact)), 8.0);
    EXPECT_NEAR(TwinSearch::compute_log_width_product(proj, fact), std::log10(8.0), 1e-12);

    // Guard against a silent regression to the old behaviour.
    EXPECT_NE(static_cast<double> (TwinSearch::compute_width_product(proj, fact)), 27.0);
}

// With min_k > 2 there is no self-loop, so the corrected width must be
// identical to what the degree-based version produced. This is what keeps the
// k-uniform results in the paper unchanged.
TEST(TwinSearchTest, WidthProductUnchangedForUniformSearches) {
    std::vector<std::vector<int> > hes = {{0,1,2},{1,2,3},{0,2,3},{0,1,3}};
    Hypergraph h(hes);
    ProjectedGraph proj(h);
    FactorGraph fact(proj, 3, 4);

    // No edge-node carries a self-loop here, so degree == |eta_e| throughout
    // and the product below is exactly what the old code computed.
    long double expected = 1.0;
    for (int eid = 0; eid < fact.num_edge_nodes; eid++) {
        ASSERT_EQ(fact.neighborhood_size(eid), fact.node_degree(eid));
        std::vector<int> e = fact.node_map.at(eid);
        expected *= binom(fact.node_degree(eid), proj.proj_mat(e[0], e[1]));
    }
    EXPECT_DOUBLE_EQ(static_cast<double> (TwinSearch::compute_width_product(proj, fact)),
                     static_cast<double> (expected));
}
