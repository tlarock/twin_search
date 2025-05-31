#ifndef FACTOR_GRAPH_H
#define FACTOR_GRAPH_H
#include <vector>
#include <map>
#include <set>

#include <boost/graph/adjacency_list.hpp>

#include "projected_graph.hpp"
#include "utils.hpp"

// Implements a factor graph, an undirected bipartite graph where one set of nodes
// represents edges and the other represents cliques in which those edges participate
class FactorGraph {
    private:
        // a boost::graph representing the hypergraph with int node ids
		UndirectedGraph g;

    public: 
        // a mapping from node_id in g to node, where each node is a vector
        std::map<int, std::vector<int> > node_map;

        // a mapping from node to node_id in g
        std::map<std::vector<int>, int> rev_node_map;

        // The number of edge and clique nodes
        // NOTE: Can be used to iterate, since edge-nodes
        // run from {0,...,num_edge_nodes-1} and clique-nodes
        // run from {num_edge_nodes-1,...,num_edge_nodes+num_clique_nodes-1}
        int num_edge_nodes;
        int num_clique_nodes;
        int min_k;
        int max_k;

        // Constructor that takes a ProjectedGraph object and constructs
        // the FactorGraph associated to it with cliques of sizes
        // k in {min_k,...,max_k}
		FactorGraph(ProjectedGraph &, int min_k_, int max_k_);
        FactorGraph();

        UndirectedGraph get_graph();
        // convenience function wrapping boost::degree
        int node_degree(int node_id);
        std::vector<int> get_vertex_neighbors(int node_id);
};

#endif
