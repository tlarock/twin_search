#include "factor_graph.hpp"

#include <stdexcept>
#include <format>

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
	int enodes_added = 0;
	// pre-load edge-nodes
	// NOTE: guard the loop bound; size1() == 0 would underflow size_t here.
	for (std::size_t i = 0; i + 1 < proj.proj_mat.size1(); i++) {
		for (std::size_t j = i+1; j < proj.proj_mat.size2(); j++) {
			if (proj.proj_mat(i,j) > 0 || proj.proj_mat(j,i) > 0) {
					enode_idx = boost::add_vertex(g);
					std::vector<int> e {static_cast <int> (i), static_cast <int> (j)};
					node_map[enode_idx] = e;
					rev_node_map[e] = enode_idx;
					enodes_added += 1;
			}
		}
	}
	// num_edge_nodes is taken from proj.num_edges but the vertices above are
	// counted independently, so the two must agree exactly. If they do not,
	// default_edge_execution_order() will hand the search edge-node ids that
	// are absent from node_map, and the parallel search would then corrupt
	// node_map looking them up. Fail here instead: a mismatch means the input
	// projection is inconsistent, not that the caller should carry on.
	//
	// NOTE: the previous check compared proj.num_edges-1 against enode_idx,
	// the last vertex id assigned, which is 0 when no vertices were added at
	// all - so it silently passed the case it most needed to catch.
	if (enodes_added != proj.num_edges) {
        throw std::runtime_error(std::format(
            "FactorGraph: projection is inconsistent - built {} edge-nodes but "
            "ProjectedGraph::num_edges is {}. This usually means a hyperedge "
            "contained a repeated node, which Hypergraph should have removed.",
            enodes_added, proj.num_edges));
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
// Raw boost::degree.
//
// WARNING: this double-counts the self-loop that an edge-node carries when
// min_k <= 2, so it is one larger than |eta_e| for every edge-node in that
// case. It is kept as a plain accessor; callers that mean |eta_e| must use
// neighborhood_size() instead.
int FactorGraph::node_degree(int node_id) {
	if (node_id < static_cast<int> (boost::num_vertices(g))) {
    	return boost::degree(node_id, g);
    } else {
        diagnostic() << "WARNING: Tried to get degree of node_id: " << node_id << " which is larger than number of vertices: " << boost::num_vertices(g) << std::endl;
		return 0;
    }
}

// The size of eta_e: the number of DISTINCT neighbours of node_id in the
// factor graph. For an edge-node e this is the number of candidate cliques
// that could satisfy e, which is what the paper's worst-case search tree
// size, prod_e binom(|eta_e|, w_e), is defined over.
//
// This is deliberately not boost::degree. When min_k <= 2 the 2-clique {u,v}
// gets no clique-node of its own - the pair is already an edge-node, so the
// constructor calls add_edge(e, e) and creates a SELF-LOOP - and boost counts
// a self-loop twice so that the degree sum stays even. The self-loop is a
// genuine member of eta_e (a 2-hyperedge is a legitimate way to satisfy the
// edge when min_k <= 2), but it must be counted once, not twice.
//
// Using boost::degree here inflated the worst-case tree size by one degree for
// every edge-node; see tests/test_factor_graph.cpp.
int FactorGraph::neighborhood_size(int node_id) {
    if (node_id >= static_cast<int> (boost::num_vertices(g))) {
        diagnostic() << "WARNING: Tried to get neighborhood size of node_id: " << node_id << " which is larger than number of vertices: " << boost::num_vertices(g) << std::endl;
        return 0;
    }
    std::set<int> unique_neighbors;
    for (int u : boost::make_iterator_range(boost::adjacent_vertices(node_id, g)))
        unique_neighbors.insert(u);
    return static_cast<int> (unique_neighbors.size());
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
        diagnostic() << "WARNING: Tried to get neighbors of node_id: " << node_id << " which is larger than number of vertices: " << boost::num_vertices(g) << std::endl;
		return std::vector<int>(0);
    }
}

