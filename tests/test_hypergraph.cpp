#include <sstream>
#include <thread>
#include <regex>
#include <functional>
#include <string>
#include "hypergraph.hpp"
#include "test_hypergraphs.cpp"
#include <gtest/gtest.h>


TEST(HypergraphTest, NodeAndEdgeCounts) {
    // Check that hypergraphs have the correct number
    // of nodes and edges
    Hypergraph h = h1();
    EXPECT_EQ(h.n, 6);
    EXPECT_EQ(h.m, 3);
    int max_size = 0;
    for(const auto& [he_idx, s] : h.hyperedge_sizes) {
        if (s > max_size) {
            max_size = s;
        }
    }

    EXPECT_EQ(max_size, 4);

    h = h2();
    EXPECT_EQ(h.n, 10);
    EXPECT_EQ(h.m, 7);
    for(const auto& [he_idx, he] : h.hyperedges) {
        EXPECT_EQ(h.hyperedge_sizes[he_idx], 3);
        EXPECT_EQ(he.size(), 3);
    }

    h = h3();
    EXPECT_EQ(h.n, 10);
    EXPECT_EQ(h.m, 8);
    for(const auto& [he_idx, he] : h.hyperedges) {
        EXPECT_EQ(h.hyperedge_sizes[he_idx], 3);
        EXPECT_EQ(he.size(), 3);
    }

    h = h4();
    EXPECT_EQ(h.n, 10);
    EXPECT_EQ(h.m, 6);
    for(const auto& [he_idx, he] : h.hyperedges) {
        EXPECT_EQ(h.hyperedge_sizes[he_idx], 3);
        EXPECT_EQ(he.size(), 3);
    }

    h = h5();
    EXPECT_EQ(h.n, 10);
    EXPECT_EQ(h.m, 10);
    for(const auto& [he_idx, he] : h.hyperedges) {
        EXPECT_EQ(h.hyperedge_sizes[he_idx], 3);
        EXPECT_EQ(he.size(), 3);
    }
}

TEST(HypergraphTest, DefaultConstructor) {
    // Check that hypergraphs have the correct number
    // of nodes and edges
    Hypergraph h;
    EXPECT_EQ(h.n, 0);
    EXPECT_EQ(h.m, 0);
    EXPECT_EQ(h.hyperedges.size(), 0);
    EXPECT_EQ(h.hyperedge_sizes.size(), 0);
    EXPECT_EQ(h.node_memberships.size(), 0);
}

TEST(HypergraphTest, NConstructor) {
    // Check that hypergraphs have the correct number
    // of nodes and edges
    int n = 7;
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 2};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 2, 3};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {1, 2, 4, 5};
    input_hyperedges.insert(he3);
    Hypergraph h(input_hyperedges, n);

    EXPECT_EQ(h.n, n);
    EXPECT_EQ(h.m, input_hyperedges.size());
    EXPECT_EQ(h.hyperedges.size(), h.m);
    EXPECT_EQ(h.hyperedge_sizes.size(), h.m);
    EXPECT_EQ(h.node_memberships.size(), n);

    // Test for singleton node handling
    int num_singletons = 0;
    int hyperedges_with_node = 0;
    for(const auto& [node, membs] : h.node_memberships) {
        hyperedges_with_node = 0; 
        for(const auto& [he_idx, he] : h.hyperedges) {
            if(std::find(he.begin(), he.end(),  node) != he.end()) {
                hyperedges_with_node += 1;
            }
        }
        if(membs.size() == 0) {
            num_singletons += 1;
            EXPECT_EQ(hyperedges_with_node, 0);
        } else {
            EXPECT_EQ(hyperedges_with_node, h.node_memberships[node].size());
        }
    }

    EXPECT_EQ(num_singletons, n-6);
}

TEST(HypergraphTest, NMConstructor) {
    // Check that hypergraphs have the correct number
    // of nodes and edges
    int n = 7;
    int m = 10;
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 2};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 2, 3};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {1, 2, 4, 5};
    input_hyperedges.insert(he3);
    Hypergraph h(input_hyperedges, n, m);

    EXPECT_EQ(h.n, n);
    EXPECT_EQ(h.m, m);
    EXPECT_EQ(h.hyperedges.size(), m);
    EXPECT_EQ(h.hyperedge_sizes.size(), m);
    EXPECT_EQ(h.node_memberships.size(), n);

    // Test for singleton node handling
    int num_singletons = 0;
    int hyperedges_with_node = 0;
    for(const auto& [node, membs] : h.node_memberships) {
        hyperedges_with_node = 0; 
        for(const auto& [he_idx, he] : h.hyperedges) {
            if(std::find(he.begin(), he.end(),  node) != he.end()) {
                hyperedges_with_node += 1;
            }
        }
        if(membs.size() == 0) {
            num_singletons += 1;
            EXPECT_EQ(hyperedges_with_node, 0);
        } else {
            EXPECT_EQ(hyperedges_with_node, h.node_memberships[node].size());
        }
    }

    EXPECT_EQ(num_singletons, n-6);

    // Check for empty hyperedge handling
    int num_empty_edges = 0;
    for(const auto& [he_idx, he] : h.hyperedges) {
        if (he.size() == 0) {
            num_empty_edges += 1;
            EXPECT_EQ(h.hyperedge_sizes[he_idx], 0);
        } else {
            EXPECT_EQ(h.hyperedge_sizes[he_idx], h.hyperedges[he_idx].size());
        }
    }
    EXPECT_EQ(num_empty_edges, h.hyperedges.size() - 3);
}

TEST(HypergraphTest, RemapTest) {
    // Check that re-mapping the nodes of a hypergraph whose input is not
    // 0...n-1 works correctly
    std::vector<Hypergraph> hypergraphs = all_hypergraphs();
    for (Hypergraph h : hypergraphs) {
        for (const auto& [u, umembs] : h.node_memberships)
            EXPECT_LT(u, h.n); // u must be less than n

        for (int u = 0; u < h.n; u++)
            EXPECT_TRUE(h.node_memberships.contains(u));
    }
}

// A node repeated inside a hyperedge is not a valid simple hypergraph, and the
// search is not defined for one. Hypergraph drops the repeats at construction
// so that downstream invariants hold - in particular so that
// ProjectedGraph::num_edges keeps matching the number of edge-nodes
// FactorGraph builds. Without this, the parallel search looks up an edge-node
// id that is absent from FactorGraph::node_map, and std::map::operator[]
// inserts it concurrently from every worker thread.
// Captures whatever the constructor wrote to stderr, so the user-facing
// warning is tested rather than assumed.
static std::string capture_stderr(const std::function<void()> &fn) {
    std::ostringstream buf;
    std::streambuf *old = std::cerr.rdbuf(buf.rdbuf());
    fn();
    std::cerr.rdbuf(old);
    return buf.str();
}

TEST(HypergraphTest, RepeatedNodesInHyperedgeAreDropped) {
    std::vector<std::vector<int> > input = {{0,1,2},{1,2,3},{0,3,3},{0,2,3}};
    Hypergraph h(input);

    EXPECT_EQ(h.m, 4);
    EXPECT_EQ(h.n, 4);

    // {0,3,3} must be stored as {0,3}
    EXPECT_EQ(h.hyperedges[2], (std::vector<int>{0,3}));
    EXPECT_EQ(h.hyperedge_sizes[2], 2);

    // node 3 must be recorded as a member of hyperedge 2 exactly once
    int count = 0;
    for (int he_idx : h.node_memberships[3])
        if (he_idx == 2) count++;
    EXPECT_EQ(count, 1);
}

TEST(HypergraphTest, RepeatedNodesDroppedInSizedConstructors) {
    std::vector<std::vector<int> > input = {{0,1,1},{1,2,2,2}};

    Hypergraph h2(input, 3);
    EXPECT_EQ(h2.hyperedges[0], (std::vector<int>{0,1}));
    EXPECT_EQ(h2.hyperedges[1], (std::vector<int>{1,2}));
    EXPECT_EQ(h2.hyperedge_sizes[0], 2);
    EXPECT_EQ(h2.hyperedge_sizes[1], 2);

    Hypergraph h3(input, 3, 2);
    EXPECT_EQ(h3.hyperedges[0], (std::vector<int>{0,1}));
    EXPECT_EQ(h3.hyperedges[1], (std::vector<int>{1,2}));
}

// A hyperedge of all-identical nodes collapses to a singleton, which
// contributes nothing off-diagonal. It must not be counted as an edge.
TEST(HypergraphTest, AllRepeatedNodesCollapseToSingleton) {
    std::vector<std::vector<int> > input = {{0,1},{2,2,2}};
    Hypergraph h(input);
    EXPECT_EQ(h.hyperedges[1], (std::vector<int>{2}));
    EXPECT_EQ(h.hyperedge_sizes[1], 1);
}

// An exact repeat of an earlier hyperedge violates simplicity too: the
// projection would require the same clique twice, while the factor graph holds
// one clique-node per clique.
TEST(HypergraphTest, DuplicateHyperedgesAreDropped) {
    std::vector<std::vector<int> > input = {{0,1,2},{1,2,3},{0,1,2},{0,2,3}};
    Hypergraph h(input);

    EXPECT_EQ(h.m, 3);
    EXPECT_EQ(h.hyperedges.size(), 3u);
    EXPECT_EQ(h.hyperedges[0], (std::vector<int>{0,1,2}));
    EXPECT_EQ(h.hyperedges[1], (std::vector<int>{1,2,3}));
    EXPECT_EQ(h.hyperedges[2], (std::vector<int>{0,2,3}));

    // node 0 belongs to the surviving {0,1,2} once, not twice
    int count = 0;
    for (int he_idx : h.node_memberships[0])
        if (he_idx == 0) count++;
    EXPECT_EQ(count, 1);
}

// De-duplicating nodes can itself create a duplicate hyperedge: {0,1,1}
// becomes {0,1}. The duplicate check must therefore run after node repair.
TEST(HypergraphTest, NodeRepairCanExposeADuplicateHyperedge) {
    std::vector<std::vector<int> > input = {{0,1},{0,1,1}};
    Hypergraph h(input);
    EXPECT_EQ(h.m, 1);
    EXPECT_EQ(h.hyperedges[0], (std::vector<int>{0,1}));
}

TEST(HypergraphTest, WarnsWhenRepeatedNodesAreRemoved) {
    std::vector<std::vector<int> > input = {{0,1,2},{0,3,3}};
    std::string msg = capture_stderr([&]{ Hypergraph h(input); });
    EXPECT_NE(msg.find("has been modified"), std::string::npos) << msg;
    EXPECT_NE(msg.find("repeated node"), std::string::npos) << msg;
    EXPECT_NE(msg.find("1 hyperedge(s) contained"), std::string::npos) << msg;
}

TEST(HypergraphTest, WarnsWhenDuplicateHyperedgesAreDropped) {
    std::vector<std::vector<int> > input = {{0,1,2},{0,1,2}};
    std::string msg = capture_stderr([&]{ Hypergraph h(input); });
    EXPECT_NE(msg.find("has been modified"), std::string::npos) << msg;
    EXPECT_NE(msg.find("duplicated an earlier hyperedge"), std::string::npos) << msg;
}

TEST(HypergraphTest, ReportsBothRepairsTogether) {
    std::vector<std::vector<int> > input = {{0,1,2},{0,1,2},{0,3,3},{1,1,2}};
    std::string msg = capture_stderr([&]{ Hypergraph h(input); });
    EXPECT_NE(msg.find("2 hyperedge(s) contained a repeated node"), std::string::npos) << msg;
    EXPECT_NE(msg.find("1 hyperedge(s) duplicated"), std::string::npos) << msg;
}

// The common case must stay silent, or the warning is worthless.
TEST(HypergraphTest, SaysNothingForWellFormedInput) {
    std::vector<std::vector<int> > input = {{0,1,2},{1,2,3},{0,2,3}};
    std::string msg = capture_stderr([&]{ Hypergraph h(input); });
    EXPECT_EQ(msg, "") << msg;
}

// The repair report is emitted from constructors that run inside TBB loops, so
// several threads can report at once. It goes through std::osyncstream, which
// buffers each report and emits it in one piece; this pins that, because a
// report interleaved with another is worse than no report at all - it looks
// like corrupted data rather than a warning.
TEST(HypergraphTest, RepairReportIsNotInterleavedAcrossThreads) {
    const int kThreads = 16;
    std::ostringstream buf;
    std::streambuf *old = std::cerr.rdbuf(buf.rdbuf());
    {
        std::vector<std::thread> threads;
        for (int t = 0; t < kThreads; t++) {
            threads.emplace_back([]{
                std::vector<std::vector<int> > bad = {{0,1,2},{0,1,2},{0,3,3}};
                Hypergraph h(bad);
            });
        }
        for (auto &th : threads) th.join();
    }
    std::cerr.rdbuf(old);

    const std::regex header("^Warning: input hypergraph was not simple and has been modified\\.$");
    const std::regex detail("^  - [0-9]+ hyperedge\\(s\\) .*\\.$");
    const std::regex footer("^  Results below describe the modified hypergraph, not the input as given\\.$");

    int headers = 0, footers = 0, unmatched = 0, total = 0;
    std::istringstream in(buf.str());
    for (std::string line; std::getline(in, line); ) {
        if (line.empty()) continue;
        total++;
        if (std::regex_match(line, header)) headers++;
        else if (std::regex_match(line, footer)) footers++;
        else if (!std::regex_match(line, detail)) unmatched++;
    }

    EXPECT_EQ(headers, kThreads);
    EXPECT_EQ(footers, kThreads);
    // A spliced report shows up as a line matching none of the three shapes.
    EXPECT_EQ(unmatched, 0) << "interleaved output:\n" << buf.str();
    EXPECT_EQ(total, kThreads * 4);   // header + 2 details + footer, per thread
}
