#ifndef PROJECTED_GRAPH_H
#define PROJECTED_GRAPH_H
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <discreture.hpp>
#include "hypergraph.hpp"
#include "utils.hpp"

// Forward declaring FactorGraph so that ProjectedGraph is aware of its
// existence without creating a dependency loop between the two classes
class FactorGraph;

// generic matrix type to facilitate swapping across the codebase without
// changing multiple files and lines
// NOTE: Benchmarked with sparse matrices and speed was degraded
typedef boost::numeric::ublas::matrix<int> ProjMatT;

// Class implementing a pairwise node co-occurrence projection of a hypergraph as a weighted
// adjacency matrix where entry proj_mat(u,v) stores the number of hyperedges
// of which both u and v are members. Also includes functionality for computing
// all cliques in the projection.
//
// TODO: ProjectedHypergraph would be a clearer name.
class ProjectedGraph {
	public: 
        // the matrix storing the projection
        ProjMatT proj_mat;
        
        // number of unique edges, multiedges (equivalently: sum(proj_mat) / 2)
        int num_edges, num_multi_edges;

        ProjectedGraph();
        // Constructor the defaults to 0s on the diagonal
		ProjectedGraph(Hypergraph h) : ProjectedGraph(h, true) {};
        ProjectedGraph(Hypergraph, bool zero_diagonal);

        // function for computing cliques of size min_k,...,max_k with defaults
		CliqueMap compute_cliques(int min_k = 2, int max_k = 3);

        // function for computing common neighbors of group of nodes
        std::vector<int> common_neighbors(std::vector<int> nodes);
        UndirectedGraph get_boost_graph();

};

#endif
