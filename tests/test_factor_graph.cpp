#include "hypergraph.hpp"
#include "projected_graph.hpp"
#include "factor_graph.hpp"
#include "test_hypergraphs.cpp"
#include <gtest/gtest.h>


TEST(FactorGraphTest, NodeAndEdgeCounts) {
    // Check that hypergraphs have the correct number
    // of nodes and edges
    std::vector<Hypergraph> hgs = all_hypergraphs();
    int min_k = 100;
    int max_k = 0;
    int clique_count;
    for(Hypergraph &h : hgs) {
        min_k = 100;
        max_k = 0;
        for (const auto &[he_idx, he] : h.hyperedges) {
            if (he.size() < static_cast<std::size_t> (min_k))
                min_k = static_cast<std::size_t> (he.size());
            if (he.size() > static_cast<std::size_t> (max_k))
                max_k = static_cast<std::size_t> (he.size());
        }
        ProjectedGraph proj(h);
        FactorGraph fact(proj, min_k, max_k);
        EXPECT_EQ(fact.num_edge_nodes, proj.num_edges);
        CliqueMap cliques = proj.compute_cliques(min_k, max_k);
        clique_count = 0;
        for (const auto &[k, kcliques] : cliques) {
            if (k > 2 && k <= max_k) {
                clique_count += kcliques.size();
            }
        }
            
        EXPECT_EQ(fact.num_clique_nodes, clique_count);
        EXPECT_EQ(boost::num_vertices(fact.get_graph()), fact.num_edge_nodes+fact.num_clique_nodes);
        for(int i = 0; i < fact.num_edge_nodes; i++)
            EXPECT_EQ(fact.node_map[i].size(), 2);

        for(int i = fact.num_edge_nodes; i < fact.num_edge_nodes + fact.num_clique_nodes; i++) {
            EXPECT_GE(fact.node_map[i].size(), min_k);
            EXPECT_LE(fact.node_map[i].size(), max_k);
        }

    }
}

TEST(FactorGraphTest, MinK2) {
    Hypergraph H = h8();
    ProjectedGraph proj(H);
    int min_k = 2;
    int max_k = H.n;
    FactorGraph fact(proj, min_k, max_k);
    for (const auto& [node_id, edge] : fact.node_map) {
        auto boost_neighbors = boost::make_iterator_range(boost::adjacent_vertices(node_id, fact.get_graph()));
        std::vector<int> neighbors = fact.get_vertex_neighbors(node_id);
        if (2 == edge.size())
            EXPECT_EQ(neighbors.size(), boost_neighbors.size()-1);
        else
            EXPECT_EQ(neighbors.size(), boost_neighbors.size());

    }
    
    min_k = 3;
    fact = FactorGraph(proj, min_k, max_k);
    for (const auto& [node_id, edge] : fact.node_map) {
        auto boost_neighbors = boost::make_iterator_range(boost::adjacent_vertices(node_id, fact.get_graph()));
        std::vector<int> neighbors = fact.get_vertex_neighbors(node_id);
        EXPECT_EQ(neighbors.size(), boost_neighbors.size());
    }
}

// FactorGraph::num_edge_nodes comes from ProjectedGraph::num_edges, but the
// edge-node vertices are counted by a separate loop. The two must agree, or
// the search is handed edge-node ids that are not in node_map. This is the
// invariant that a repeated node inside a hyperedge used to break.
TEST(FactorGraphTest, EdgeNodeCountMatchesProjectionEdgeCount) {
    std::vector<std::vector<int> > with_repeat = {{0,1,2},{1,2,3},{0,3,3},{0,2,3}};
    Hypergraph h(with_repeat);
    ProjectedGraph proj(h);
    FactorGraph fact(proj, 2, 4);

    int true_edges = 0;
    for (std::size_t i = 0; i + 1 < proj.proj_mat.size1(); i++)
        for (std::size_t j = i + 1; j < proj.proj_mat.size2(); j++)
            if (proj.proj_mat(i, j) > 0) true_edges++;

    EXPECT_EQ(fact.num_edge_nodes, true_edges);

    // Every edge-node id the search will iterate over must already be present,
    // so that .at() never throws and operator[] would never have inserted.
    for (int eid = 0; eid < fact.num_edge_nodes; eid++) {
        ASSERT_TRUE(fact.node_map.contains(eid)) << "missing edge-node " << eid;
        EXPECT_EQ(fact.node_map.at(eid).size(), 2u);
    }
}
