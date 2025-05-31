#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

// Including to get UndirectedGraph for get_bipartite
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <utils.hpp>

namespace ublas=boost::numeric::ublas;
// A very simple hypergraph class that stores the basic structure of
// the hypergraph as well as some maps for convenient traversal
class Hypergraph {
	public:
        // number of nodes, number of hyperedges
		int n, m;

        // map from hyperedge id to hyperedge vector
        // hyperedge_id is in {0,...,m} with no guaranteed ordering
		std::map<int, std::vector<int> > hyperedges;

        // map from hyperedge id to the size of the hyperedge
		std::map<int, int> hyperedge_sizes;

        // map from node_id (in {0,...,n}) to hyperedge_ids containing node
		std::map<int, std::vector<int> > node_memberships;

        // constructors for either set or vector of vectors as input
        Hypergraph();
        template<typename T>
        Hypergraph(const T&);
        template<typename T>
        Hypergraph(const T &input_hyperedges, int n_);
        template<typename T>
        Hypergraph(const T &input_hyperedges, int n_, int m_);
        
        UndirectedGraph get_bipartite();
        ublas::matrix<int> get_incidence_matrix();
        UndirectedGraph get_line_graph();
        ublas::matrix<int> get_lg_mat();
        // method for printing to console
        void pretty_print();
};

// Constructor for a hypergraph using template input.
//
// The template type must be an iterable containing vector<int>s of hyperedges,
// or at least a container that is sortable.
// Generally the expected types are:
//      std::set<vector<int> > or std::vector<vector<int> >
//
// This version does not accept the number of nodes or hyperedges as input.
// Instead, it infers the number of nodes based on the input, using the number
// of unique nodes as n (NOT the node id).
//
// NOTE: Does NOT support std::map<int, vector<int> >. May want to add this
// possibility later, but for now must be iterable of iterable.
//
// NOTE: If nodes are not on 0,...,n-1, automatically re-maps from input ids.
template<typename T>
Hypergraph::Hypergraph(const T &input_hyperedges) {
    Hypergraph::m = input_hyperedges.size();
    // Use default initializers if empty
    if (m > 0) {
        std::set<int> nodes;
        int he_idx = 0;
        //for (std::vector<int> he : input_hyperedges)
        for(const auto& he : input_hyperedges)
        {
            Hypergraph::hyperedges[he_idx] = he;
            sort(hyperedges[he_idx].begin(), hyperedges[he_idx].end());
            Hypergraph::hyperedge_sizes[he_idx] = he.size();
            for (int node_id : he)
            {
                Hypergraph::node_memberships[node_id].push_back(he_idx);
                nodes.insert(node_id);
            }
            he_idx += 1;
        }

        Hypergraph::n = static_cast<int> (nodes.size());

        // check if nodes are 0,...,n, remap if not
        int max_node_id = *nodes.rbegin();
        if (max_node_id != n-1) {
            std::vector<int> nodes_vect(nodes.begin(), nodes.end());
            sort(nodes_vect.begin(), nodes_vect.end());
            std::map<int, int> node_map;
            for (int node_id = 0; node_id < static_cast <int> (nodes_vect.size()); node_id++) {
                node_map[nodes_vect[node_id]] = node_id;
            }

            for(const auto& [old_id, new_id] : node_map) {
                // Replace old edge id in hyperedges with new id
                for(int edge_id : node_memberships[old_id]) {
                    std::replace(hyperedges[edge_id].begin(), hyperedges[edge_id].end(), old_id, new_id);
                    sort(hyperedges[edge_id].begin(), hyperedges[edge_id].end());
                }

                if (!node_memberships.contains(new_id)) {
                    // move node memberships to new_id location
                    node_memberships[new_id] = std::vector<int>(node_memberships[old_id].size());
                    std::copy(node_memberships[old_id].begin(), node_memberships[old_id].end(), node_memberships[new_id].begin());
                    // Remove old id from node memberships
                    node_memberships.erase(old_id);
                } else {
                    // Need to swap
                    std::vector<int> tmp(node_memberships[new_id].size());
                    std::copy(node_memberships[new_id].begin(), node_memberships[new_id].end(), tmp.begin());
                    std::copy(node_memberships[old_id].begin(), node_memberships[old_id].end(), node_memberships[new_id].begin());
                    node_memberships[old_id] = std::vector<int>(tmp.size());
                    std::copy(tmp.begin(), tmp.end(), node_memberships[old_id].begin());
                }
            }
        }
    }
}


// Constructor for a hypergraph using template input.
//
// The template type must be an iterable containing vector<int>s of hyperedges,
// or at least a container that is sortable.
// Generally the expected types are:
//      std::set<vector<int> > or std::vector<vector<int> >
//
// This version accepts the number of nodes as input and supports singleton
// nodes in the sense that empty vectors are instantiated in node_memberships
// for all nodes from 0,...,n-1.
//
// NOTE: Does NOT support std::map<int, vector<int> >. May want to add this
// possibility later, but for now must be iterable of iterable.
//
// NOTE: If nodes are not on 0,...,n-1, automatically re-maps from input ids.
template<typename T>
Hypergraph::Hypergraph(const T &input_hyperedges, int n_) {
    Hypergraph::n = n_;
    // Gaurantee entries in node_memberships
    for(int u = 0; u < n; u++) {
        node_memberships[u] = std::vector<int>(0);
    }

    Hypergraph::m = input_hyperedges.size();
    if (m > 0) {
        std::set<int> nodes;
        int he_idx = 0;
        //for (std::vector<int> he : input_hyperedges)
        for(const auto& he : input_hyperedges) {
            Hypergraph::hyperedges[he_idx] = he;
            sort(hyperedges[he_idx].begin(), hyperedges[he_idx].end());
            Hypergraph::hyperedge_sizes[he_idx] = he.size();
            for (int node_id : he) {
                Hypergraph::node_memberships[node_id].push_back(he_idx);
                nodes.insert(node_id);
            }
            he_idx += 1;
        }

        if (static_cast<int> (nodes.size()) > Hypergraph::n) {
            std::cout << "Number of nodes in the input " << nodes.size() << " is larger than input value for n " << n_ << ". Value of Hypergraph.n is actual number of nodes." << std::endl;
            Hypergraph::n = static_cast<int> (nodes.size());
        }
    }
}

// Constructor for a hypergraph using template input.
//
// The template type must be an iterable containing vector<int>s of hyperedges,
// or at least a container that is sortable.
// Generally the expected types are:
//      std::set<vector<int> > or std::vector<vector<int> >
//
// This version accepts the number of nodes AND number of hyperedges as input
// and supports singleton nodes and empty hyperedges in the sense that empty
// vectors are instantiated in node_memberships for all nodes from 0,...,n-1
// and in hyperedges for edge_id 0,...,m-1
//
// NOTE: Does NOT support std::map<int, vector<int> >. May want to add this
// possibility later, but for now must be iterable of iterable.
//
// NOTE: If nodes are not on 0,...,n-1, automatically re-maps from input ids.
template<typename T>
Hypergraph::Hypergraph(const T &input_hyperedges, int n_, int m_) {
    Hypergraph::n = n_;
    Hypergraph::m = m_;

    // Gaurantee entries in node_memberships
    for(int u = 0; u < n; u++) {
        node_memberships[u] = std::vector<int>(0);
    }

    // Gaurantee entries in node_memberships
    for(int eid = 0; eid < m; eid++) {
        Hypergraph::hyperedges[eid] = std::vector<int>(0);
        Hypergraph::hyperedge_sizes[eid] = 0;
    }

    // Use default initializers if empty
    if (input_hyperedges.size() > 0) {
        std::set<int> nodes;
        int he_idx = 0;
        //for (std::vector<int> he : input_hyperedges)
        for(const auto& he : input_hyperedges) {
            Hypergraph::hyperedges[he_idx] = he;
            sort(hyperedges[he_idx].begin(), hyperedges[he_idx].end());
            Hypergraph::hyperedge_sizes[he_idx] = he.size();
            for (int node_id : he) {
                Hypergraph::node_memberships[node_id].push_back(he_idx);
                nodes.insert(node_id);
            }
            he_idx += 1;
        }


        if (static_cast<int> (nodes.size()) > n_) {
            std::cout << "Number of nodes in the input " << nodes.size() << " is larger than input value for n " << n_ << ". Value of Hypergraph.n is actual number of nodes." << std::endl;
            Hypergraph::n = static_cast<int> (nodes.size());
        }

        if (static_cast<int> (hyperedges.size()) > m_) {
            std::cout << "Number of hyperedges in the input " << Hypergraph::hyperedges.size() << " is larger than input value for m " << m_ << ". Value of Hypergraph.m is actual number of hyperedges." << std::endl;
            Hypergraph::m = static_cast<int> (hyperedges.size());
        } 
    }
}

#endif
