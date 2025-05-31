#include "factor_graph.hpp"

// FactorGraph constructor from ProjectedGraph. Uses proj
// to compute all cliques between min_k and max_k, then constructs 
// the factor graph from those cliques.
//
// Note: Using an encapsulation DAG or a g-trie in place of
// the clique map object would be more efficient, since that
// structure would already contain the factor relationships.
// This verion requires a doubly nested loop over each clique
// that could be partially avoided in the future by just moving
// the logic from ProjectedGraph.compute_cliques here.
FactorGraph::FactorGraph(ProjectedGraph & proj, int min_k_, int max_k_){
    min_k = min_k_;
    max_k = max_k_;
    // Get cliques map from proj
    std::map<int, std::set<std::vector<int> > > cliques = proj.compute_cliques(min_k, max_k);

    // edge-node ids will run from 0,...,n-1
	int enode_idx = 0;
	// pre-load edge-nodes
	for (std::size_t i = 0; i < proj.proj_mat.size1()-1; i++) {
		for (std::size_t j = i+1; j < proj.proj_mat.size2(); j++) {
			if (proj.proj_mat(i,j) > 0 || proj.proj_mat(j,i) > 0) {
					enode_idx = boost::add_vertex(g);
					std::vector<int> e {static_cast <int> (i), static_cast <int> (j)};
					node_map[enode_idx] = e;
					rev_node_map[e] = enode_idx;
			}
		}
	}
	if (enode_idx != proj.num_edges-1) {
        std::cout << "Something is wrong. enode_idx: " << enode_idx << " != proj.num_edges-1: " << proj.num_edges-1 << std::endl;
	}
    num_edge_nodes = proj.num_edges;

    // clique-node ids will run from n,...,n+number of clique-nodes-1
    int cnode_idx = num_edge_nodes;
    std::vector<int> e(2);
    for (const auto& [k, kcliques] : cliques) {
        if (k >= min_k) {
            for (std::vector<int> clique : kcliques) {
                // Add clique-node to maps
                if (!rev_node_map.contains(clique)) {
                    node_map[cnode_idx] = clique;
                    rev_node_map[clique] = cnode_idx;
                    cnode_idx += 1;
                }
                // Get all edges from cliuqe
                for (std::size_t i = 0; i < clique.size()-1; i++) {
                    for (std::size_t j = i+1; j < clique.size(); j++) {
                        e[0] = clique[i];
                        e[1] = clique[j];
                        // Add factor graph edge
                        boost::add_edge(rev_node_map[e], rev_node_map[clique], g);
                    }
                }
            }
        }
    }
    num_clique_nodes = boost::num_vertices(g) - num_edge_nodes;
}

// minimal default constructor
FactorGraph::FactorGraph() {}

// Getter function for the graph to make
// it more difficult to "accidentally" modify.
UndirectedGraph FactorGraph::get_graph() { return g; }

// Convenience function to get degree of a node without having
// to access fact.g directly
int FactorGraph::node_degree(int node_id) {
	if (node_id < static_cast<int> (boost::num_vertices(g))) {
    	return boost::degree(node_id, g);
    } else {
        std::cout << "WARNING: Tried to get degree of node_id: " << node_id << " which is larger than number of vertices: " << boost::num_vertices(g) << std::endl;
		return 0;
    }
}

// Function that gets the set of unique neighbors of a vertex in the factor graph.
// NOTE: This is needed to deal with self-loops, which are used for factor
// graphs with minimum hyperedge size of 2 (e.g., edges).
//
// This is because boost includes the self-loop twice in the return of
// adjacent_vertices, so that the degree count is correct, but I don't
// want this behavior and can't find a way to turn it off.
// TODO: It may be more memory efficient to avoid using both a set and a vector.
std::vector<int> FactorGraph::get_vertex_neighbors(int node_id) {
    if (node_id < static_cast<int> (boost::num_vertices(g))) {
        std::vector<int> ne_vect;
        auto neighbors = boost::make_iterator_range(boost::adjacent_vertices(node_id, g));
        if (min_k <= 2) {
            std::set<int> ne_set;
            for (int u : neighbors)
                ne_set.insert(u);
            ne_vect = std::vector<int>(ne_set.size());
            std::copy(ne_set.begin(), ne_set.end(), ne_vect.begin());
        } else {
            ne_vect = std::vector<int>(neighbors.size());
            std::copy(neighbors.begin(), neighbors.end(), ne_vect.begin());
        }

    	return ne_vect;
    } else {
        std::cout << "WARNING: Tried to get neighbors of node_id: " << node_id << " which is larger than number of vertices: " << boost::num_vertices(g) << std::endl;
		return std::vector<int>(0);
    }
}

