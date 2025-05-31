#include "projected_graph.hpp"

namespace ublas=boost::numeric::ublas;

// minimal default constructor
ProjectedGraph::ProjectedGraph()
{
    ProjectedGraph::proj_mat = ProjMatT(0,0,0);
    ProjectedGraph::num_edges = 0;
    ProjectedGraph::num_multi_edges = 0;
}

// Constructor that creates a ProjectedGraph from a Hypergraph
// TODO: Currently can't deal with missing/singleton nodes in h (since
// Hypergraph does not deal with them)
ProjectedGraph::ProjectedGraph(Hypergraph h, bool zero_diagonal)
{
	// initialize projection matrix
    ProjectedGraph::proj_mat = ProjMatT(h.n, h.n, 0);
   
    // initialize counter for unique undirected edges
    ProjectedGraph::num_edges = 0;

    // initialize counter for total undirected multi-edges
    ProjectedGraph::num_multi_edges = 0;

    // loop over hyperedges map
    for (const auto& [he_idx, he] : h.hyperedges)
    {
        // Add each pairwise edge
        for (std::size_t i = 0; i < he.size(); i++)
        {
            if (!zero_diagonal)
                proj_mat(he[i], he[i]) += 1;

            for (std::size_t j = i+1; j < he.size(); j++)
            {
                // Count edges and multi-edges
                if (proj_mat(he[i], he[j]) < 1)
                    num_edges += 1;

                num_multi_edges += 1;

                // Increment weight in projection adjacency matrix
                proj_mat(he[i], he[j]) += 1;
                proj_mat(he[j], he[i]) += 1;
            }
        }
    }
}


// Compute all cliques of sizes min_k,...,max_k in this ProjectedGraph
// adjacency matrix. Edge weights are ignored in this function, any
// entry with value > 0 is considered an edge. Cliques are guaranteed
// to be sorted by node ids.
CliqueMap ProjectedGraph::compute_cliques(int min_k, int max_k)
{
    // The output data structure, a map from clique size k
    // to a set of vectors representing sorted cliques
	CliqueMap cliques;

    // Declaring a vector for storing common neighbors
    std::vector<int> common_ne;

    // Start by gathering the edges
    int k = 2;
    std::vector<int> curr_nodes (k);
    for (std::size_t i = 0; i < proj_mat.size1()-1; i++) {
        for (std::size_t j = i+1; j < proj_mat.size2(); j++) {
            if (proj_mat(i, j) > 0) {
                curr_nodes[0] = i;
                curr_nodes[1] = j;
                sort(curr_nodes.begin(), curr_nodes.end());
                cliques[k].insert(curr_nodes);
            }
        }
    }

    // For k from 2...max_k-1, loop over each clique
    // of size k and find all cliques of size k+1 by
    // searching the adjacency matrix for common neighbors
    while (k < max_k) {
        std::vector<int> clique (k+1);
        for (std::vector<int> curr_nodes : cliques[k]) {
            common_ne = common_neighbors(curr_nodes);
            
            // If there are no common neighbors, ignore
            if (common_ne.size() < 1) {continue;}

            // For each common neighbor, add a clique
            for (int w : common_ne) {
                // Fill the clique vector with the current items
                for (int i = 0; i < k; i++) { clique[i] = curr_nodes[i]; }
                clique[k] = w;

                // Sort the clique to avoid duplicates
                sort(clique.begin(), clique.end());

                // add to the clique set
                cliques[k+1].insert(clique);
            }
        }
        // if cliques[k+1] is empty, we can break here
        if (cliques[k+1].size() == 0)
            break;

        k++;
    }

    // if min_k is larger than edges, remove all smaller cliques 
    if (min_k > 2) {
        for (int k = 2; k < min_k; k++) {
            cliques.erase(k);
        }
    }

	return cliques;
}

// Given a group of nodes, find neighbors they all have in common by summing up
// their adjacency matrix rows into a vector. Entries in this vector with value
// the same as the size of the group of nodes are common neighbors.
std::vector<int> ProjectedGraph::common_neighbors(std::vector<int> nodes)
{
    int num_nodes = static_cast<int> (nodes.size());
    int n = proj_mat.size1();
    ublas::vector<int> sum_vect (n);
    std::fill(sum_vect.begin(),sum_vect.end(),0);
    for (int node : nodes)
    {
        // Get binarized projection matrix row
        // Add them element-wise
        for (int j = 0; j < n; ++j) {
            if (j == node)
                continue;

            if (proj_mat(node, j) > 0) {
                sum_vect(j) += 1;
            }
        }
    }

    // Entries with size(nodes) are common neighbors
    std::vector<int> out;
    for (int i = 0; i < n; i++) {
        if (sum_vect(i) == num_nodes) {
            out.push_back(i);
        }
    }

    return out;
}

UndirectedGraph ProjectedGraph::get_boost_graph() {
    UndirectedGraph projg;
    for (std::size_t r = 0; r < proj_mat.size1(); ++r) {
        for (std::size_t c = 0; c < proj_mat.size2(); ++c) {
            if (r == c)
                continue;

            if (proj_mat(r,c) > 0) {
                for (int v = 0; v < proj_mat(r,c); ++v)
                    boost::add_edge(r, c, projg);
            }
        }
    }
    return projg;
}
