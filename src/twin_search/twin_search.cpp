#include "twin_search.hpp"
#include <float.h>
#include <algorithm>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <oneapi/tbb.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>


// Including here because not needed in hpp
#include "factor_graph.hpp"
#include "utils.hpp"

namespace ublas=boost::numeric::ublas;

// A mutex for writing to cout from within this file
// Note: does not prevent threads in other translation
// units from writing to cout at the same time!
// TODO: Switch to wrapping std::cout in std::basic_osyncstream
// instead of this ad-hoc non-solution!
tbb::spin_mutex COUT_MUTEX;

// A simple struct to store a partial hypergraph, its projection
// remainder, and the index into the execution order of the current
// edge to be satisfied.
// NOTE: ProjMatT defined in projected_graph.hpp
struct TwinSearch::StackItem {
    std::vector<int> hypergraph;
    ProjMatT proj_rem;
    std::size_t edge_execution_index;
};

// Searches for all hypergraphs that correspond to this projected adjacency
// matrix, constraining the minimum and maximum hyperedge size. Also compares
// each pair of hypergraphs to check whether their line graphs are equivalent,
// which would make them Gram Mates.
//
// If filter_isomorphic is true, we also run an isomorphism test on every
// pair of hypergraphs (via their bipartite representations using
// boost::vf2_graph_iso) and fill filtered_twins with the index of 1 representative of
// each isomorphism class found.
TwinSearch::TwinSearch(ProjectedGraph proj_, int min_k, int max_k, bool filter_isomorphic, bool parallel, bool run_search, bool use_diag_)
{
    TwinSearch::proj = proj_;
    TwinSearch::use_diagonal = use_diag_;

    // Validate use_diagonal setting
    int proj_diag_sum = 0;
    for (std::size_t u = 0; u < proj.proj_mat.size1(); ++u)
        proj_diag_sum += proj.proj_mat(u,u);

    if (use_diagonal && proj_diag_sum < 1) {
        std::cout << "Warning: use_diagonal set to true, but sum of diagonal is 0." << std::endl;
    } else if (!use_diagonal && proj_diag_sum > 1) {
        std::cout << "Warning: use_diagonal set to false, but sum of diagonal is greater than 0. Setting proj.proj_mat(i,i) entries to 0." << std::endl;
        for (std::size_t u = 0; u < proj.proj_mat.size1(); ++u)
            proj.proj_mat(u,u) = 0;
    }

    // Get a FactorGraph of the projection
    TwinSearch::fact = FactorGraph(proj, min_k, max_k);

    // Initialize containers for twins
    TwinSearch::twins = std::vector<std::vector<int> >(0);
    TwinSearch::filtered_twins = std::vector<int>(0);
    TwinSearch::mates = std::vector<std::vector<int> >(0);

    // Detrmine if a search is feasible based on input
    TwinSearch::feasible = test_feasibility();

    if ( run_search && feasible) {
        // Actually run the twins search
        if (parallel)
            parallel_search(filter_isomorphic);
        else
            search(filter_isomorphic);
    } else if (!feasible){
        {
            tbb::spin_mutex::scoped_lock lock(COUT_MUTEX);
            std::cout << "Warning: encountered infeasible search due to an edge-node with degree smaller than its weight." << std::endl;
        }
    }
}

// TODO: Second constructor that takes a hyperedge size distribution 

bool TwinSearch::test_feasibility() {
    int degree;
    int weight;
    std::vector<int> edge;
    for(int eid = 0; eid < fact.num_edge_nodes; eid++) {
	    degree = fact.node_degree(eid);
        edge = fact.node_map[eid];
        weight = proj.proj_mat(edge[0], edge[1]);
        // Degree must be larger than or equal to weight
        if (degree == 0 || degree < weight)
            return false;
    }
    return true;
}

bool TwinSearch::equivalent_lg(std::vector<ublas::matrix<int> > &line_graphs, const int i, const int j) {
    // Check that the dimensions agree, otherwise they cannot be equivalent
    if ( !( (line_graphs[i].size1() == line_graphs[j].size1()) && (line_graphs[i].size2() == line_graphs[j].size2()) ) )
        return false;

    // Check if the matrices are equivalent
    for (std::size_t r = 0; r < line_graphs[i].size1(); r++) {
        for (std::size_t c = r; c < line_graphs[i].size2(); c++) {
            if ( line_graphs[i](r, c) != line_graphs[j](r, c) ) {
                return false;
            }
        }
    }

    return true;
}

void TwinSearch::process_item(std::vector<StackItem> &stack, StackItem &s, std::vector<UndirectedGraph> &bipartites, std::vector<UndirectedGraph> &line_graphs, bool filter_isomorphic) {
    // check if the sum of the modified projection is 0
    int proj_rem_sum = matsum(s.proj_rem);
    if (proj_rem_sum < 1) {
        // Add every twin to twins
        twins.push_back(s.hypergraph);

        // Compute the bipartite representation and the incidence matrix for
        // s.hypergraph to use for comparisons
        compute_bipartite_and_linegraph(bipartites, line_graphs, -1, s.hypergraph);

        // If required, decide whether to put this hypergraph into
        // filtered_twins or not.
        // NOTE: When filter_isomorphic is false, we are wasting some
        // time/space by always constructing the bipartite representation
        // However, the speed savings by not doing the isomorhpism
        // tests is going to be much larger than the cost of computing the
        // bipartite graphs, and when we do want to do the filtering it would be a
        // waste to compute them separately, so I think this is a reasonable trade-off.
        if (filter_isomorphic) {
            // if the filtered vector is empty *OR*
            // the current hypergraph is not isomorphic to anything in bipartites
            if (filtered_twins.empty() || !is_isomorphic(bipartites, bipartites.size()-1)) {
                // Add the index of this twin to filtered_twins
                filtered_twins.push_back(twins.size()-1);
            }
        }
    }
    else {
        // if the projection is not yet satisfied, choose another
        // edge and add combinations to the stack

        // If edge_execution_index is already exhausted, do nothing
        if (s.edge_execution_index >= edge_execution_order.size())
            return;

        // Pop the next edge of the execution order vector
        int enode_id = edge_execution_order[s.edge_execution_index];	
        std::vector<int> e = fact.node_map[enode_id];
        // skip to the next unsatisfied edge
        while (s.proj_rem(e[0], e[1]) < 1) {
            s.edge_execution_index++;
            if (s.edge_execution_index >= edge_execution_order.size())
                return;
            enode_id = edge_execution_order[s.edge_execution_index];
            e = fact.node_map[enode_id];
        }
        std::vector<std::vector<int> > comb_vect = get_combinations(enode_id, s.proj_rem(e[0], e[1]), s);
        for (std::vector<int> comb : comb_vect)
            add_to_stack(s, comb, stack);
    }
}

void TwinSearch::search(bool filter_isomorphic) {
    if (!feasible) {
        {
            tbb::spin_mutex::scoped_lock lock(COUT_MUTEX);
            std::cout << "search() was called on infeasible projection. Returning without running search." << std::endl;
        }
    }
    // Initialize container for bipartite representations
    std::vector<UndirectedGraph> bipartites;
    std::vector<UndirectedGraph> line_graphs;


    // Initialize a stack representation using the StackItem struct
    std::vector<StackItem> stack;

    // Get a vector of edge ids ordered by their constraint values
    // such that deterministic edges are solved first.
    edge_execution_order = default_edge_execution_order();

    // The first stackitem is always an empty hypergraph and the
    // ProjectedGraph.proj_mat matrix from the input
    StackItem s(std::vector<int> (0), proj.proj_mat, 0);
    process_item(stack, s, bipartites, line_graphs, filter_isomorphic);

    // Stores the sum of StackItem.proj_rem to check
    // whether we have satisfied every edge
    while (!stack.empty()) {
        // pop an item off the stack
        s = stack.back();
        stack.pop_back();
        process_item(stack, s, bipartites, line_graphs, filter_isomorphic);
    }

    mates = run_mates_tests_parallel(line_graphs);

}

// Takes a StackItem and returns the first unsatisfied edge based
// on a row-major iteration over the projection.
// NOTE: Have replaced this with execution_order functionality, although
// this function performs similarly to default_edge_execution_order.
int TwinSearch::choose_random_edge(const TwinSearch::StackItem &s) {
    int enode_id;
    // loop over non-zero entries in s.proj_rem
    for (std::size_t i = 0; i < s.proj_rem.size1()-1; i++) {
        for (std::size_t j = i+1; j < s.proj_rem.size2(); j++) {
            if (s.proj_rem(i, j) > 0) {
                // get the enode_id corresponding to the entry
                std::vector<int> e {static_cast<int> (i), static_cast<int> (j)};
                enode_id = fact.rev_node_map.at(e);
                return enode_id;
            }
        }
    }

    // TODO: This should return a bad index (such as a negative number)
    // so that it can be handled appropriately, rather than just trying
    // to use the 0th index when it is not really positive.
    return 0;
}

// utility function that accepts a vector v and returns a vector
// of indices into v sorted in ascending order.
// NOTE: Copied from StackOverflow:
// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
// NOTE: Used only with compute_edge_execution_order function. If needed
// elsewhere, should move to utils.cpp.
template <typename T>
std::vector<std::size_t> sort_indices(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// NOTE: Unused function, since default_edge_execution_order performed
// better in benchmarking. Keeping this here for future research.
//
// Uses TwinSearch::fact and TwinSearch::proj to compute an execution
// ordering for the edges based on increasing constraint value with
// ties broken arbitrarily (preferring fewer index swaps by using stable_sort).
std::vector<std::size_t> TwinSearch::compute_edge_execution_order() {
    std::vector<int> constraints;
    double constraint;
    int i, j;
    std::vector<int> neighbors_vect;
    std::vector<int> e; 
    // loop over non-zero entries in s.proj_rem
    // loop over edge ids
    for (int enode_id = 0; enode_id < fact.num_edge_nodes; enode_id++) {
        e = fact.node_map[enode_id];
        i = e[0];
        j = e[1];
        neighbors_vect = fact.get_vertex_neighbors(enode_id);
        constraint = binom(neighbors_vect.size(), proj.proj_mat(i, j));
        constraints.push_back(constraint);
    }

    std::vector<std::size_t> execution_order = sort_indices(constraints);
    
    return execution_order;
}

// Returns a "random" (really default) ordering for execution of
// the edges, just their indices. Originally designed to compare
// against the function compute_edge_execution_order() defined above,
// it turned out to be faster in benchmarks to just use this order.
std::vector<std::size_t> TwinSearch::default_edge_execution_order() {
    std::vector<std::size_t> indices(fact.num_edge_nodes);
    std::iota(indices.begin(), indices.end(), 0);
    return indices;
}

// Add an item to the stack, first checking whether it is a viable candidate
// based on the remaining entries in curr.proj_rem. 
void TwinSearch::add_to_stack(const TwinSearch::StackItem &curr, const std::vector<int> &comb, std::vector<TwinSearch::StackItem> &stack) {
    // Copy curr.proj_rem to a new ublas::matrix
    ProjMatT new_proj_rem(curr.proj_rem);

    // Check whether the addition of any of the cliques in comb would over-use
    // any edge (making proj_rem(edge) negative). If it would, then this
    // combination is not viable so we will return without adding it to the stack.
    std::vector<int> new_hypergraph;
    std::vector<int> clique;
    for (int cnode_id : comb) {
        clique = fact.node_map.at(cnode_id);
        for (std::size_t i = 0; i < clique.size(); i++) {
            for (std::size_t j = i+1; j < clique.size(); j++) {
                new_proj_rem(clique[i], clique[j]) -= 1;
                new_proj_rem(clique[j], clique[i]) -= 1;
                if ( new_proj_rem(clique[i], clique[j]) < 0 ) {
                    // add nothing to the stack and return
                    return;
                }
            }

            // If diagonals are non-zero, decrement
            // TODO: This is not exactly correct. I think we need a flag here,
            // otherwise we can't tell whether proj_rem(0,0) == 0 is just
            // because the matrix was 0-diagonal or because this entry should
            // be disallowed.
            if (use_diagonal) {
                if ( new_proj_rem(clique[i], clique[i]) > 0 ) {
                    new_proj_rem(clique[i], clique[i]) -= 1;
                } else {
                    // add nothing to the stack and return
                    return;
                }
            }
        }
        new_hypergraph.push_back(cnode_id);
    }

    // Add existing edges to the new hypergraph
    for (int cnode_id : curr.hypergraph) {
        new_hypergraph.push_back(cnode_id);
    }

    // Add a StackItem representing the new hypergraph and proj_rem
    TwinSearch::StackItem to_add(new_hypergraph, new_proj_rem, curr.edge_execution_index+1);

    stack.push_back(to_add);
}


// Uses the discretures library to construct combinations(cnode_neighbors, weight) for enode_id.
std::vector<std::vector<int> > TwinSearch::get_combinations(int enode_id, int weight, const TwinSearch::StackItem &s) {
    std::vector<int> neighbors_vect = get_filtered_neighbors(s, enode_id);
    std::vector<std::vector<int> > comb_vect(0);
    if (neighbors_vect.size() > 0) {
        // NOTE: discreture uses some rvalue magic that I don't fully understand,
        // so I just push the combinations into an std::vector.
        auto combs = discreture::combinations(neighbors_vect, weight);
        for (auto&& comb : combs) {
            std::vector<int> new_comb;
            for (int c : comb) {
                new_comb.push_back(c);
            }
            comb_vect.push_back(new_comb);
        }
    }

    return comb_vect;
}

// Callback function/struct for vf2_sub_graph_iso
// NOTE: Found via SO.
// TODO: Add proper tests of this callback
template <typename Graph1, typename Graph2>
struct my_callback {
    my_callback(const Graph1& graph1, const Graph2& graph2)
      : graph1_(graph1), graph2_(graph2) {}

    template <typename CorrespondenceMap1To2,
              typename CorrespondenceMap2To1>
    bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) const {
      return true;
    }

    private:
        const Graph1& graph1_;
        const Graph2& graph2_;
};

// Takes a vector of bipartite twin representations and a candidate twin ID,
// then compares the candidate twin to those found in filtered_twins.
//
//  NOTE: boost::vf2_graph_iso can not handle self-loops with undirectedS
//  graphs. Proceed with caution.
//
// Returns true if isomorphic to an existing twin isomorphism class, false otherwise.
bool TwinSearch::is_isomorphic(const std::vector<UndirectedGraph> &bipartites, const int cand_idx) {
    for (int i : filtered_twins) {
        if (i != cand_idx) {
            my_callback<UndirectedGraph, UndirectedGraph> my_callback(bipartites[i], bipartites[cand_idx]);
            if ( boost::vf2_graph_iso(bipartites[i], bipartites[cand_idx], my_callback) )
                return true;
        }
    }
    return false;
}

// Compute an UndirectedGraph corresponding to the bipartite representation and
// a ublas::matrix<int> corresponding to the line graph of the input
// hypergraph. If a non-negative value is given for i, the new structures will
// be placed in the vectors at position i. Otherwise, they will be appended to
// the end.
// NOTE: For the purposes of computing isomorphisms, it does not matter if the
// identifiers in the bipartite graph match the real nodes or factor graph
// identifiers, so I will not bother with vertexpropertymaps and so on to
// preserve ids between fact and bipartite.
//
// TODO: I could simplify and standardize this using member functions of 
// hypergraph objects, but it would require constructing those objects,
// which will add a multiplier to how this function works.
void TwinSearch::compute_bipartite_and_linegraph(std::vector<UndirectedGraph> &bipartites, std::vector<UndirectedGraph> &line_graphs, const int i, std::vector<int> &hypergraph) { 
    // Initialize objects for new structures
    UndirectedGraph bipartite;
    ublas::matrix<int> incidence_matrix(hypergraph.size(), proj.proj_mat.size1(), 0);

    // id_map will keep track of which ids have
    // already been given new identities
    std::map<int, int> id_map;
    std::vector<int> e;
    int bp_node_id = 0;
    int incidence_row = 0;
    for (int cnode_id : hypergraph) {
        // add cnode_id to the map
        id_map[cnode_id] = bp_node_id;
        bp_node_id++;
        // Loop over the nodes in the clique and add them to
        // the bipartite and incidence matrix representations
        for (int u : fact.node_map.at(cnode_id)) {
            // nodes are already in 0,..,n-1, which works for
            // the incidence matrix as is
            incidence_matrix(incidence_row, u) = 1;

            // For the bipartite graph we need to map since cnode_ids
            // will overlap with 0,...,n-1
            if ( id_map.find(u) == id_map.end() ) {
                id_map[u] = bp_node_id;
                bp_node_id++;
            }
            boost::add_edge(id_map[u], id_map[cnode_id], bipartite);
        }
        incidence_row++;
    }

    UndirectedGraph line_graph;
    ublas::matrix<int> lg_mat = ublas::prod(incidence_matrix, ublas::trans(incidence_matrix));
    for (std::size_t r = 0; r < lg_mat.size1(); ++r) {
        for (std::size_t c = 0; c < lg_mat.size2(); ++c) {
            if (r == c)
                continue;

            if (lg_mat(r, c) > 0) {
                for (int v = 0; v < lg_mat(r,c); ++v)
                    boost::add_edge(r, c, line_graph);
            }
        }
    }

    if (i < 0) { 
        line_graphs.push_back(line_graph);
        bipartites.push_back(bipartite);
    } else {
        bipartites[i] = bipartite;
        line_graphs[i] = line_graph;
    }
}

// Take a stackitem and an enode_id and get a vector of cnode_ids corresponding
// to the neighbors of enode_id that do not appear in s.hypergraph
std::vector<int> TwinSearch::get_filtered_neighbors(const StackItem &s, int enode_id) {
    // Get an iterator over the clique-neighbors of enode_id
    std::vector<int> neighbors = fact.get_vertex_neighbors(enode_id);
    // Remove clique-neighbors that are already in s.hypergraph
    std::vector<int> neighbors_vect;
    for (int cnode_id : neighbors) {
        if ( std::find(s.hypergraph.begin(), s.hypergraph.end(), cnode_id) == s.hypergraph.end() ) {
            neighbors_vect.push_back(cnode_id);
        }
    }

    return neighbors_vect;
}

// Convenience function to inflate a hypergraph from a vector of clique-node id
// ints to a vector of vectors of ints representing the actual hyperedges
std::vector<std::vector<int> > TwinSearch::inflate_cnodes(const std::vector<int> &cnode_ids) {
    std::vector<std::vector<int> > inflated_hypergraph(cnode_ids.size());
    std::vector<int> inflated_hyperedge;
    // For each hyperedge
    for (std::size_t i = 0; i < cnode_ids.size(); i++) {
        // Inflate the hyperedge from the factor graph clique-node
        inflated_hyperedge = std::vector<int>(fact.node_map.at(cnode_ids[i]).size());
        for (std::size_t j = 0; j < inflated_hyperedge.size(); j++) {
            inflated_hyperedge[j] = fact.node_map.at(cnode_ids[i])[j];
        }

        // Add the inflated hyperedge to the inflated hypergraph
        for (std::size_t j = 0; j < inflated_hyperedge.size(); j++)
            inflated_hypergraph[i].push_back(inflated_hyperedge[j]);
    }
    return inflated_hypergraph;
}

// prints a container of containers of hypegraphs to the console
void TwinSearch::print_twins(const std::vector<std::vector<int> > &twins){
    for (std::vector<int> mate : twins) {
        for (int cnode_id : mate) {
            std::vector<int> he = fact.node_map.at(cnode_id);
            for (std::size_t i = 0; i < he.size(); i++) {
                if (i < he.size()-1) 
                    std::cout << he[i] << ",";
                else
                    std::cout << he[i] << " ";
            }
        }
        std::cout << std::endl;
    }
}

// START PARALLEL FUNCTIONS
void TwinSearch::parallel_search(bool filter_isomorphic) {
    if (!feasible) {
        {
            tbb::spin_mutex::scoped_lock lock(COUT_MUTEX);
            std::cout << "parallel_search() was called on infeasible projection. Returning without running search." << std::endl;
        }
    }
    tbb::concurrent_vector<std::vector<int> > concurrent_twins;
    
    // Initialize a stack representation using the StackItem struct
    // Note: In the parallel version this is not really implementing
    // a stack, it is actually a queue.
    std::vector<StackItem> stack;

    // Get a vector of edge ids ordered by their constraint values
    // such that deterministic edges are solved first.
    edge_execution_order = default_edge_execution_order();

    // The first stackitem is always an empty hypergraph and the
    // ProjectedGraph.proj_mat matrix from the input
    StackItem s(std::vector<int> (0), ProjMatT(proj.proj_mat), 0);
    std::vector<StackItem> init_stack(0);
    parallel_process_item(s, concurrent_twins, init_stack);
    for(StackItem new_item : init_stack) {
        stack.push_back(StackItem(new_item));
    }
    tbb::parallel_for_each(stack.begin(), stack.end(),
            [&](StackItem &s, tbb::feeder<StackItem>& feeder) {
        std::vector<StackItem> tmp_stack(0);
        parallel_process_item(s, concurrent_twins, tmp_stack);
        if(!tmp_stack.empty()) {
            for(std::size_t i = 0; i < tmp_stack.size(); i++) {
                feeder.add(StackItem(std::vector<int>(tmp_stack[i].hypergraph), ProjMatT(tmp_stack[i].proj_rem), tmp_stack[i].edge_execution_index));
            }
        }
    });

    for(std::vector<int> hg : concurrent_twins) {
        twins.push_back(hg);
    }


    // Initialize container for bipartite and line graph representations
    std::vector<UndirectedGraph> bipartites(twins.size());
    std::vector<UndirectedGraph> line_graphs(twins.size());
    // Construct in parallel
    tbb::parallel_for(tbb::blocked_range<int>(0, twins.size()),
                       [&](tbb::blocked_range<int> r) {
        for (int i=r.begin(); i<r.end(); ++i) {
            compute_bipartite_and_linegraph(bipartites, line_graphs, i, twins[i]);
        }
    });


    // Check for mates
    mates = run_mates_tests_parallel(line_graphs);

    if (filter_isomorphic) {
        std::vector<int> to_filter = TwinSearch::run_iso_tests_parallel(bipartites);
        for(std::size_t i = 0; i < twins.size(); i++) {
            if (to_filter[i] == 0)
                filtered_twins.push_back(i);
        }
    }
}

void TwinSearch::parallel_process_item(StackItem &s, tbb::concurrent_vector<std::vector<int> > &concurrent_twins, std::vector<StackItem> &tmp_stack) {
    //std::vector<StackItem> tmp_stack;
    // check if the sum of the modified projection is 0
    int proj_rem_sum = matsum(s.proj_rem);
    if (proj_rem_sum < 1) {
        // Add every twin to twins 
        // safe with concurrent vector
        concurrent_twins.push_back(s.hypergraph);
    } else {
        // if the projection is not yet satisfied, choose another
        // edge and add combinations to the stack

        // Get the next edge from the execution order vector
        if (s.edge_execution_index >= edge_execution_order.size())
            return;

        int enode_id = edge_execution_order[s.edge_execution_index];
        std::vector<int> e = fact.node_map[enode_id];
        // skip to the next unsatisfied edge
        while (s.proj_rem(e[0], e[1]) < 1) {
            s.edge_execution_index++;
            if (s.edge_execution_index >= edge_execution_order.size())
                return;
            enode_id = edge_execution_order[s.edge_execution_index];
            e = fact.node_map[enode_id];
        }
        std::vector<std::vector<int> > comb_vect = get_combinations(enode_id, s.proj_rem(e[0], e[1]), s);
        for (std::vector<int> comb : comb_vect) {
            add_to_stack(s, comb, tmp_stack);
        }
    }
}

// Static function that accepts a vector of UndirectedGraphs and runs a parallel loop
// that implements a pairwise comparison, adding 1 to the to_filter concurrent vector
// at position j if the bipartite graph in that position is isomorphic to a graph at some
// position i < j.
//
//  NOTE: boost::vf2_graph_iso can not handle self-loops with undirectedS
//  graphs. Proceed with caution.
//
std::vector<int> TwinSearch::run_iso_tests_parallel(std::vector<UndirectedGraph> &bipartites) {
    std::vector<int> to_filter(bipartites.size());
    if (bipartites.size() < 2)
        return to_filter;

    for(std::size_t i = 0; i < bipartites.size()-1; i++) {
            if(to_filter[i] > 0) {
                continue;
            }
        tbb::parallel_for(std::size_t(i+1), bipartites.size(), [&](std::size_t j){
            if (to_filter[j] < 1) {
                my_callback<UndirectedGraph, UndirectedGraph> mc(bipartites[i], bipartites[j]);
                if ( boost::vf2_graph_iso(bipartites[i], bipartites[j], mc) ) {
                    to_filter[j] += 1;
                }
            }
        });
    }

    return to_filter;
}

// Static function that accepts a vector of boost graphs representing line
// graphs and runs a parallel loop that does a pairwise comparison and fills
// the mates vector with pairs of indices pointing to any mates found
//
// NOTE: This is not exactly a gram mates test, since it is about isomorphism
// rather than equivalence. However, this will also catch any equivalent line
// graphs since they are trivially isomorphic.
//
// NOTE: boost::vf2_graph_iso can not handle self-loops with undirectedS
// graphs. Proceed with caution.
//
std::vector<std::vector<int> > TwinSearch::run_mates_tests_parallel(std::vector<UndirectedGraph> &line_graphs) {
    tbb::concurrent_vector<std::vector<int> > mate_pairs;
    if (line_graphs.size() < 2)
        return std::vector<std::vector<int> >(0);

    for(std::size_t i = 0; i < line_graphs.size()-1; i++) {
        tbb::parallel_for(std::size_t(i+1), line_graphs.size(), [&](std::size_t j){
                thread_local my_callback<UndirectedGraph, UndirectedGraph> mc(line_graphs[i], line_graphs[j]);
                if ( boost::vf2_graph_iso(line_graphs[i], line_graphs[j], mc) ) {
                    // If line graphs are isomorphic, i and j are a pair of mates
                    mate_pairs.push_back( std::vector<int> {static_cast<int> (i), static_cast<int> (j)});
                }
        });
    }

    // put in an std vector for return
    std::vector<std::vector<int> > ret(mate_pairs.size());
    for (std::size_t i = 0; i < mate_pairs.size(); i++)
        ret[i] = mate_pairs[i];

    return ret;
}

// Static function that takes a factor graph and a projection (presumably the
// one that constructed the factor graph, but this is NOT tested - results in
// the case where the projection comes from elsewhere are undefined) and computes
// the maximum width of the search over this factor graph.
long double TwinSearch::compute_width_product(ProjectedGraph &proj, FactorGraph &fact) {
    // TODO I don't have a way of checking for overflow here
    long double prod = 0.0;
    int degree;
    int weight;
    std::vector<int> edge;
    for(int eid = 0; eid < fact.num_edge_nodes; eid++) {
        // Get num neighbors
        degree = fact.node_degree(eid);
        edge = fact.node_map[eid];
        weight = proj.proj_mat(edge[0], edge[1]);
        if (prod < 1)
            prod = binom(degree, weight);
        else
            prod *= binom(degree, weight);
    }

    return prod;
}

// Static function that takes a factor graph and a projection (presumably the
// one that constructed the factor graph, but this is NOT tested - results in
// the case where the projection comes from elsewhere are undefined) and computes
// the maximum width of the search over this factor graph.
double TwinSearch::compute_log_width_product(ProjectedGraph &proj, FactorGraph &fact) {
    double logsum = 0.0;
    int degree;
    int weight;
    std::vector<int> edge;
    for(int eid = 0; eid < fact.num_edge_nodes; eid++) {
        // Get num neighbors
        degree = fact.node_degree(eid);
        edge = fact.node_map[eid];
        weight = proj.proj_mat(edge[0], edge[1]);
        // TODO: There could be overflow/imprecision in binom
        logsum += std::log10(binom(degree, weight));
    }

    return logsum;
}
