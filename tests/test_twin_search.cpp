#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/assignment.hpp>
#include "projected_graph.hpp"
#include "factor_graph.hpp"
#include "utils.hpp"
#include <filesystem>
#include <fstream>
#include "twin_search.hpp"
#include "search_checkpoint.hpp"
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

// --- checkpointing: drain, serialise, reload, resume ---
//
// The gate for the whole idea is one property: a search that was interrupted
// and resumed must produce the SAME answer as one that ran straight through.
// Not a similar twin count - the same twins, the same isomorphism classes, the
// same mate pairs.
//
// n=7 m=14 was chosen by search rather than guessed: 20 twins, 11 after the
// isomorphism filter and 9 mate pairs, so all three phases of parallel_search
// actually produce something. A fixture whose filter removed nothing, or whose
// mate set was empty, would pass while testing a third of the work.
namespace {
std::vector<std::vector<int> > resume_fixture() {
    return {{0,1,2},{0,1,5},{0,2,3},{0,2,4},{0,3,4},{0,3,6},{0,4,6},
            {0,5,6},{1,3,6},{1,4,5},{2,3,4},{2,3,5},{2,4,6},{2,5,6}};
}
std::vector<std::vector<int> > canon_twins(std::vector<std::vector<int> > v) {
    for (std::vector<int> &g : v) std::sort(g.begin(), g.end());
    std::sort(v.begin(), v.end());
    return v;
}
// cnode vectors -> hyperedges -> cnode vectors, the round trip a checkpoint makes.
std::vector<std::vector<int> > to_cnodes(
        const std::vector<std::vector<std::vector<int> > > &graphs, TwinSearch &ts) {
    std::vector<std::vector<int> > out;
    for (const std::vector<std::vector<int> > &g : graphs) {
        std::vector<int> cn;
        for (const std::vector<int> &he : g) cn.push_back(ts.cnode_for(he));
        out.push_back(cn);
    }
    return out;
}
}  // namespace

TEST(TwinSearchTest, DrainedSearchResumesThroughAFileToTheSameAnswer) {
    Hypergraph h(resume_fixture());
    ProjectedGraph proj(h);

    TwinSearch full(proj, 3, 3, true, true, false, false);
    ASSERT_TRUE(full.feasible);
    full.parallel_search(true);
    ASSERT_EQ(full.twins.size(), 20u);
    ASSERT_EQ(full.filtered_twins.size(), 11u);
    ASSERT_EQ(full.mates.size(), 9u);

    const std::string path =
        (std::filesystem::temp_directory_path() / "twin_ckpt_test").string();

    // Several budgets so the drain lands at different depths, including right
    // at the root. Per-thread, so the exact node is not deterministic - which
    // is the point: the resume has to be correct wherever it stopped.
    int drained_at_least_once = 0;
    for (std::uint64_t budget : {std::uint64_t(1), std::uint64_t(3),
                                 std::uint64_t(10), std::uint64_t(40)}) {
        TwinSearch part(proj, 3, 3, true, true, false, false);
        part.arm_drain_nodes(budget);
        part.parallel_search(true);
        if (!part.drained) continue;        // finished inside the budget
        drained_at_least_once++;

        // A drained search must NOT have run the two pairwise phases: they
        // range over the whole twin set and would be wrong from a prefix.
        EXPECT_TRUE(part.mates.empty());
        EXPECT_TRUE(part.filtered_twins.empty());

        SearchCheckpoint cp;
        cp.n = 7; cp.m = 14; cp.k = 3; cp.min_k = 3; cp.max_k = 3;
        cp.commit = "test";
        cp.proj = twin_projection_key(proj.proj_mat);
        for (const std::vector<int> &g : part.twins)
            cp.twins.push_back(part.inflate_cnodes(g));
        for (const std::vector<int> &g : part.frontier)
            cp.frontier.push_back(part.inflate_cnodes(g));

        std::string err;
        ASSERT_TRUE(write_checkpoint(path, cp, err)) << err;
        SearchCheckpoint back;
        ASSERT_TRUE(read_checkpoint(path, back, err)) << err;
        EXPECT_EQ(back.proj, cp.proj);
        EXPECT_EQ(back.twins.size(), cp.twins.size());
        EXPECT_EQ(back.frontier.size(), cp.frontier.size());

        TwinSearch res(proj, 3, 3, true, true, false, false);
        ASSERT_TRUE(res.parallel_search_from(to_cnodes(back.frontier, res),
                                             to_cnodes(back.twins, res), true));
        EXPECT_EQ(canon_twins(res.twins), canon_twins(full.twins))
            << "budget " << budget;
        EXPECT_EQ(res.filtered_twins.size(), full.filtered_twins.size());
        EXPECT_EQ(res.mates.size(), full.mates.size());
    }
    EXPECT_GT(drained_at_least_once, 0);
    std::filesystem::remove(path);
}

TEST(TwinSearchTest, TruncatedCheckpointIsRejectedRatherThanPartlyRead) {
    // The reason the format has an end sentinel. A file cut short by a kill
    // mid-write would otherwise parse as a shorter but perfectly valid
    // checkpoint, and the resume would silently drop part of the search tree -
    // producing too few twins with nothing to indicate why.
    Hypergraph h(resume_fixture());
    ProjectedGraph proj(h);
    TwinSearch part(proj, 3, 3, true, true, false, false);
    part.arm_drain_nodes(3);
    part.parallel_search(true);
    ASSERT_TRUE(part.drained);

    SearchCheckpoint cp;
    cp.n = 7; cp.m = 14; cp.k = 3; cp.min_k = 3; cp.max_k = 3;
    cp.proj = twin_projection_key(proj.proj_mat);
    for (const std::vector<int> &g : part.frontier)
        cp.frontier.push_back(part.inflate_cnodes(g));
    const std::string path =
        (std::filesystem::temp_directory_path() / "twin_ckpt_trunc").string();
    std::string err;
    ASSERT_TRUE(write_checkpoint(path, cp, err)) << err;

    std::vector<std::string> lines;
    { std::ifstream is(path); std::string l; while (std::getline(is, l)) lines.push_back(l); }
    ASSERT_GT(lines.size(), 3u);
    { std::ofstream os(path, std::ios::trunc);
      for (std::size_t i = 0; i + 1 < lines.size(); i++) os << lines[i] << '\n'; }

    SearchCheckpoint back;
    EXPECT_FALSE(read_checkpoint(path, back, err));
    EXPECT_NE(err.find("#end"), std::string::npos) << "got: " << err;
    std::filesystem::remove(path);
}

TEST(TwinSearchTest, FrontierFromTheWrongProjectionIsRefused) {
    // residual_for returns false when a parked hypergraph is not a valid prefix
    // of this projection. Without that, a checkpoint from a different sample
    // would search a tree that has nothing to do with the projection in hand
    // and return a confident wrong answer.
    Hypergraph a(resume_fixture());
    std::vector<std::vector<int> > other = resume_fixture();
    other[0] = {1, 2, 3};                       // perturb one hyperedge
    Hypergraph b(other);
    ProjectedGraph pa(a), pb(b);

    TwinSearch sa(pa, 3, 3, true, true, false, false);
    sa.arm_drain_nodes(5);
    sa.parallel_search(true);
    ASSERT_TRUE(sa.drained);
    ASSERT_FALSE(sa.frontier.empty());

    TwinSearch sb(pb, 3, 3, true, true, false, false);
    std::vector<std::vector<int> > mapped;
    bool mappable = true;
    for (const std::vector<int> &g : sa.frontier) {
        std::vector<int> cn;
        for (const std::vector<int> &he : sa.inflate_cnodes(g)) {
            const int c = sb.cnode_for(he);
            if (c < 0) { mappable = false; break; }
            cn.push_back(c);
        }
        if (!mappable) break;
        mapped.push_back(cn);
    }
    // Either a hyperedge has no clique in the other projection, or it does and
    // the residual goes negative. Both are refusals; neither may be a search.
    if (mappable)
        EXPECT_FALSE(sb.parallel_search_from(mapped, {}, true));
}
