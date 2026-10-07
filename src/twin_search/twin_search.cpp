#include <array>
#include <atomic>
#include <chrono>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include "twin_search.hpp"
#include <float.h>
#include <algorithm>
#include <limits>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <oneapi/tbb.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>


// Including here because not needed in hpp
#include "factor_graph.hpp"
#include "utils.hpp"

namespace ublas=boost::numeric::ublas;

// A simple struct to store a partial hypergraph, its projection
// remainder, and the index into the execution order of the current
// edge to be satisfied.
// NOTE: ProjMatT defined in projected_graph.hpp
struct TwinSearch::StackItem {
    std::vector<int> hypergraph;
    ProjMatT proj_rem;
    std::size_t edge_execution_index;
};

// Cheap isomorphism invariant: (|V|, |E|, sorted degree sequence). Isomorphic
// graphs necessarily agree on all three, so testing this before calling
// boost::vf2_graph_iso can only skip pairs that vf2 would have rejected - the
// set of reported pairs is unchanged. This matters because the pairwise loops
// below are O(T^2) in the number of twins and vf2 is expensive per call (it
// re-sorts vertices by multiplicity every time).
struct TwinSearch::GraphFingerprint {
    std::size_t num_vertices = 0;
    std::size_t num_edges = 0;
    std::vector<std::size_t> degrees;   // ascending
};

bool TwinSearch::fingerprints_match(const GraphFingerprint &a, const GraphFingerprint &b) {
    return a.num_vertices == b.num_vertices
        && a.num_edges == b.num_edges
        && a.degrees == b.degrees;
}

TwinSearch::GraphFingerprint TwinSearch::compute_fingerprint(const UndirectedGraph &g) {
    GraphFingerprint fp;
    fp.num_vertices = boost::num_vertices(g);
    fp.num_edges = boost::num_edges(g);
    fp.degrees.reserve(fp.num_vertices);
    for (std::size_t v = 0; v < fp.num_vertices; ++v)
        fp.degrees.push_back(boost::degree(v, g));
    std::sort(fp.degrees.begin(), fp.degrees.end());
    return fp;
}

// Deliberately serial: this is O(T * V log V) against the O(T^2) vf2 work it
// saves, so parallelising it buys nothing measurable, and it keeps one more
// nested TBB region out of a call path that is already nested two deep.
std::vector<TwinSearch::GraphFingerprint> TwinSearch::compute_fingerprints(const std::vector<UndirectedGraph> &graphs) {
    std::vector<GraphFingerprint> fps;
    fps.reserve(graphs.size());
    for (const UndirectedGraph &g : graphs)
        fps.push_back(compute_fingerprint(g));
    return fps;
}

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
        diagnostic() << "Warning: use_diagonal set to true, but sum of diagonal is 0." << std::endl;
    } else if (!use_diagonal && proj_diag_sum > 0) {
        // NOTE: this was `> 1`, which let a diagonal summing to exactly 1
        // through unzeroed. matsum(proj_rem) then never reached 0, so the
        // search reported no twins at all rather than ignoring the diagonal as
        // requested. The message always said "greater than 0"; the condition
        // did not.
        diagnostic() << "Warning: use_diagonal set to false, but sum of diagonal is greater than 0. Setting proj.proj_mat(i,i) entries to 0." << std::endl;
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
        diagnostic() << "Warning: encountered infeasible search due to an edge-node with degree smaller than its weight." << std::endl;
    }
}

// TODO: Second constructor that takes a hyperedge size distribution 

// NOTE: this uses neighborhood_size, not node_degree. It previously used the
// raw boost degree, which double-counts an edge-node's self-loop when
// min_k <= 2, so the test was off by one in the PERMISSIVE direction: an edge
// needing exactly one more hyperedge than it had candidate cliques was still
// declared feasible.
bool TwinSearch::test_feasibility() {
    int degree;
    int weight;
    std::vector<int> edge;
    for(int eid = 0; eid < fact.num_edge_nodes; eid++) {
	    degree = fact.neighborhood_size(eid);
        edge = fact.node_map.at(eid);
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

// State for the cooperative drain.
struct TwinSearch::DrainState {
    std::atomic<bool> drain{false};
    bool deadline_set = false;
    std::chrono::steady_clock::time_point deadline;
    std::uint64_t node_budget = 0;      // per thread; 0 = unused
    tbb::enumerable_thread_specific<std::uint64_t> local_count;
    tbb::concurrent_vector<std::vector<int> > parked;
};

// File-scope rather than a TwinSearch member: one process runs many samples
// concurrently and a signal applies to all of them.
static std::atomic<bool> g_drain_requested{false};
void twin_search_request_global_drain()   { g_drain_requested.store(true, std::memory_order_relaxed); }
bool twin_search_global_drain_requested() { return g_drain_requested.load(std::memory_order_relaxed); }
void twin_search_clear_global_drain()     { g_drain_requested.store(false, std::memory_order_relaxed); }

void TwinSearch::arm_drain(double seconds) {
    drain_state = std::make_shared<DrainState>();
    if (seconds > 0.0) {
        drain_state->deadline_set = true;
        drain_state->deadline = std::chrono::steady_clock::now()
            + std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                  std::chrono::duration<double>(seconds));
    }
}

void TwinSearch::arm_drain_nodes(std::uint64_t nodes_per_thread) {
    drain_state = std::make_shared<DrainState>();
    drain_state->node_budget = nodes_per_thread;
}

// Returns true if this node was parked instead of expanded.
bool TwinSearch::drain_check(const StackItem &s) {
    DrainState &d = *drain_state;
    std::uint64_t &c = d.local_count.local();
    ++c;
    if (!d.drain.load(std::memory_order_relaxed)) {
        bool fire = false;
        if (d.node_budget > 0) {
            fire = (c > d.node_budget);
        } else if ((c & 1023u) == 0u) {
            // Poll once per 1024 nodes per thread. steady_clock::now() is a
            // vDSO call - cheap, but this is the hottest loop in the program,
            // and neither trigger needs millisecond precision.
            fire = twin_search_global_drain_requested()
                || (d.deadline_set && std::chrono::steady_clock::now() >= d.deadline);
        }
        if (fire)
            d.drain.store(true, std::memory_order_relaxed);
    }
    if (!d.drain.load(std::memory_order_relaxed))
        return false;
    // Park rather than expand. Items already queued in the feeder are still
    // handed to this function once each and parked here too, so the loop
    // drains without having to reach inside TBB's feeder.
    d.parked.push_back(s.hypergraph);
    return true;
}

void TwinSearch::process_item(std::vector<StackItem> &stack, StackItem &s, std::vector<UndirectedGraph> &bipartites, std::vector<GraphFingerprint> &fingerprints, std::vector<UndirectedGraph> &line_graphs, bool filter_isomorphic) {
    // check if the sum of the modified projection is 0
    int proj_rem_sum = matsum(s.proj_rem);
    if (proj_rem_sum < 1) {
        // Add every twin to twins
        twins.push_back(s.hypergraph);

        // Compute the bipartite representation and the incidence matrix for
        // s.hypergraph to use for comparisons
        compute_bipartite_and_linegraph(bipartites, line_graphs, -1, s.hypergraph);

        // Kept index-aligned with bipartites unconditionally, so the two can
        // never drift apart.
        fingerprints.push_back(compute_fingerprint(bipartites.back()));

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
            if (filtered_twins.empty() || !is_isomorphic(bipartites, fingerprints, bipartites.size()-1)) {
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
        std::vector<int> e = fact.node_map.at(enode_id);
        // skip to the next unsatisfied edge
        while (s.proj_rem(e[0], e[1]) < 1) {
            s.edge_execution_index++;
            if (s.edge_execution_index >= edge_execution_order.size())
                return;
            enode_id = edge_execution_order[s.edge_execution_index];
            e = fact.node_map.at(enode_id);
        }
        // Streamed rather than materialised: the old form built a vector of
        // every combination and then copied each one again by value.
        std::vector<int> neighbors_vect = get_filtered_neighbors(s, enode_id);
        for_each_combination(neighbors_vect, s.proj_rem(e[0], e[1]),
            [&](const std::vector<int> &comb) { add_to_stack(s, comb, stack); });
    }
}

void TwinSearch::search(bool filter_isomorphic) {
    if (!feasible) {
        diagnostic() << "search() was called on infeasible projection. Returning without running search." << std::endl;
        // The message said so but the return was missing, so an infeasible
        // projection ran the whole search anyway. It could only ever produce an
        // empty twin set - test_feasibility fails when some edge-node has fewer
        // clique-neighbours than its weight, which no partial hypergraph can
        // satisfy - so this is wasted work rather than a wrong answer, but the
        // containers are already empty and the message was a lie.
        return;
    }
    // Initialize container for bipartite representations
    std::vector<UndirectedGraph> bipartites;
    std::vector<GraphFingerprint> fingerprints;
    std::vector<UndirectedGraph> line_graphs;


    // Initialize a stack representation using the StackItem struct
    std::vector<StackItem> stack;

    // Get a vector of edge ids ordered by their constraint values
    // such that deterministic edges are solved first.
    edge_execution_order = default_edge_execution_order();

    // The first stackitem is always an empty hypergraph and the
    // ProjectedGraph.proj_mat matrix from the input
    StackItem s(std::vector<int> (0), proj.proj_mat, 0);
    process_item(stack, s, bipartites, fingerprints, line_graphs, filter_isomorphic);

    // Stores the sum of StackItem.proj_rem to check
    // whether we have satisfied every edge
    while (!stack.empty()) {
        // pop an item off the stack
        s = stack.back();
        stack.pop_back();
        process_item(stack, s, bipartites, fingerprints, line_graphs, filter_isomorphic);
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
        e = fact.node_map.at(enode_id);
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
    // Validate the candidate BEFORE copying anything.
    //
    // This function used to copy curr.proj_rem into a new ublas::matrix as its
    // very first act, decrement entries in the copy, and bail out as soon as
    // one went negative. In a search that prunes heavily most candidates are
    // rejected, so the common path was: allocate a matrix, touch a few entries,
    // discover the candidate is infeasible, free the matrix. That makes the
    // allocator the contended resource under parallel search, which does not
    // scale with threads however many are idle.
    //
    // The accepted set is unchanged. The old code rejected on the first entry
    // to go negative; this accumulates the same decrements and rejects if any
    // entry would finish below zero. For an entry starting at v and decremented
    // d times both accept exactly when v - d >= 0. That includes the diagonal:
    // "must be > 0 before each of d decrements" is the same condition.
    //
    // deltas is thread_local to keep its capacity across calls, which is the
    // point - a plain local would allocate on first push_back every call. It
    // carries no state between calls: it is cleared on entry, nothing escapes,
    // and add_to_stack is a leaf, so it is never re-entered on one thread.
    static thread_local std::vector<std::array<int, 3> > deltas;   // row, col, count
    deltas.clear();

    // Returns the running decrement count for (r, c) including this one.
    auto bump = [](std::vector<std::array<int, 3> > &ds, int r, int c) -> int {
        for (std::array<int, 3> &d : ds)
            if (d[0] == r && d[1] == c)
                return ++d[2];
        ds.push_back({r, c, 1});
        return 1;
    };

    for (int cnode_id : comb) {
        // By reference: node_map.at() returned by value, copying a vector per
        // clique per candidate.
        const std::vector<int> &clique = fact.node_map.at(cnode_id);
        for (std::size_t i = 0; i < clique.size(); i++) {
            for (std::size_t j = i+1; j < clique.size(); j++) {
                if (curr.proj_rem(clique[i], clique[j]) - bump(deltas, clique[i], clique[j]) < 0)
                    return;
            }

            // If diagonals are non-zero, decrement
            // TODO: This is not exactly correct. I think we need a flag here,
            // otherwise we can't tell whether proj_rem(0,0) == 0 is just
            // because the matrix was 0-diagonal or because this entry should
            // be disallowed.
            if (use_diagonal) {
                if (curr.proj_rem(clique[i], clique[i]) - bump(deltas, clique[i], clique[i]) < 0)
                    return;
            }
        }
    }

    // Viable, so now pay for the copy.
    ProjMatT new_proj_rem(curr.proj_rem);
    for (const std::array<int, 3> &d : deltas) {
        new_proj_rem(d[0], d[1]) -= d[2];
        if (d[0] != d[1])
            new_proj_rem(d[1], d[0]) -= d[2];   // proj_rem is symmetric
    }

    std::vector<int> new_hypergraph;
    new_hypergraph.reserve(comb.size() + curr.hypergraph.size());
    for (int cnode_id : comb)
        new_hypergraph.push_back(cnode_id);
    for (int cnode_id : curr.hypergraph)
        new_hypergraph.push_back(cnode_id);

    stack.push_back(TwinSearch::StackItem(new_hypergraph, new_proj_rem, curr.edge_execution_index+1));
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
bool TwinSearch::is_isomorphic(const std::vector<UndirectedGraph> &bipartites, const std::vector<GraphFingerprint> &fingerprints, const int cand_idx) {
    for (int i : filtered_twins) {
        if (i != cand_idx && fingerprints_match(fingerprints[i], fingerprints[cand_idx])) {
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
// Bipartite representation of one hypergraph. Same construction as in
// compute_bipartite_and_linegraph, split out so a caller that needs only this
// does not pay for a line graph it will never read.
UndirectedGraph TwinSearch::compute_bipartite(const std::vector<int> &hypergraph) {
    UndirectedGraph bipartite;
    std::map<int, int> id_map;
    int bp_node_id = 0;
    for (int cnode_id : hypergraph) {
        id_map[cnode_id] = bp_node_id;
        bp_node_id++;
        for (int u : fact.node_map.at(cnode_id)) {
            // cnode_ids overlap 0..n-1, so node ids are remapped
            if ( id_map.find(u) == id_map.end() ) {
                id_map[u] = bp_node_id;
                bp_node_id++;
            }
            boost::add_edge(id_map[u], id_map[cnode_id], bipartite);
        }
    }
    return bipartite;
}

// Line graph of one hypergraph, via its incidence matrix.
UndirectedGraph TwinSearch::compute_linegraph(const std::vector<int> &hypergraph) {
    ublas::matrix<int> incidence_matrix(hypergraph.size(), proj.proj_mat.size1(), 0);
    int incidence_row = 0;
    for (int cnode_id : hypergraph) {
        for (int u : fact.node_map.at(cnode_id))
            incidence_matrix(incidence_row, u) = 1;
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
    return line_graph;
}

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
    // Named rather than per-line, so the whole dump is emitted as one unit.
    SyncStream out(std::cout);
    for (std::vector<int> mate : twins) {
        for (int cnode_id : mate) {
            std::vector<int> he = fact.node_map.at(cnode_id);
            for (std::size_t i = 0; i < he.size(); i++) {
                if (i < he.size()-1) 
                    out << he[i] << ",";
                else
                    out << he[i] << " ";
            }
        }
        out << "\n";
    }
}

// START PARALLEL FUNCTIONS
// TEMPORARY instrumentation for the memory experiment. getrusage reports the
// PEAK so far, so comparing the value at successive phase boundaries shows
// which phase actually raised the high-water mark. Enabled only when
// TWIN_RSS_TRACE is set, and writes to stderr so it cannot pollute output.
#include <sys/resource.h>
static void rss_mark(const char *phase, std::size_t n) {
    if (!std::getenv("TWIN_RSS_TRACE")) return;
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    const double mb = ru.ru_maxrss / 1048576.0;   // bytes on macOS
#else
    const double mb = ru.ru_maxrss / 1024.0;      // kilobytes on Linux
#endif
    std::cerr << "RSS_TRACE " << phase << " peak_mb=" << mb
              << " twins=" << n << std::endl;
}

int TwinSearch::cnode_for(const std::vector<int> &hyperedge) const {
    std::vector<int> key = hyperedge;
    std::sort(key.begin(), key.end());
    auto it = fact.rev_node_map.find(key);
    return it == fact.rev_node_map.end() ? -1 : it->second;
}

bool TwinSearch::residual_for(const std::vector<int> &hypergraph, ProjMatT &out) const {
    out = ProjMatT(proj.proj_mat);
    for (int cnode_id : hypergraph) {
        auto it = fact.node_map.find(cnode_id);
        if (it == fact.node_map.end())
            return false;
        const std::vector<int> &clique = it->second;
        for (std::size_t i = 0; i < clique.size(); i++) {
            for (std::size_t j = i + 1; j < clique.size(); j++) {
                if (--out(clique[i], clique[j]) < 0)
                    return false;
                out(clique[j], clique[i]) = out(clique[i], clique[j]);
            }
            if (use_diagonal && --out(clique[i], clique[i]) < 0)
                return false;
        }
    }
    return true;
}

bool TwinSearch::parallel_search_from(const std::vector<std::vector<int> > &frontier_cnodes,
                                      const std::vector<std::vector<int> > &twins_so_far,
                                      bool filter_isomorphic) {
    if (!feasible) {
        diagnostic() << "parallel_search_from called on infeasible projection." << std::endl;
        return false;
    }
    edge_execution_order = default_edge_execution_order();
    std::vector<StackItem> seeds;
    seeds.reserve(frontier_cnodes.size());
    for (const std::vector<int> &hg : frontier_cnodes) {
        ProjMatT r;
        if (!residual_for(hg, r))
            return false;     // not a valid prefix of this projection
        // Index 0 rather than a stored position: every edge before the real one
        // is already at zero, so process_item's skip loop lands on exactly the
        // edge this item stopped at. One fewer thing the format has to get
        // right, and one fewer thing that can go stale.
        seeds.push_back(StackItem(hg, r, 0));
    }
    parallel_search_impl(seeds, twins_so_far, filter_isomorphic);
    return true;
}

void TwinSearch::parallel_search(bool filter_isomorphic) {
    if (!feasible) {
        diagnostic() << "parallel_search() was called on infeasible projection. Returning without running search." << std::endl;
        // The message said so but the return was missing, so an infeasible
        // projection ran the whole search anyway. It could only ever produce an
        // empty twin set - test_feasibility fails when some edge-node has fewer
        // clique-neighbours than its weight, which no partial hypergraph can
        // satisfy - so this is wasted work rather than a wrong answer, but the
        // containers are already empty and the message was a lie.
        return;
    }
    // Get a vector of edge ids ordered by their constraint values
    // such that deterministic edges are solved first.
    edge_execution_order = default_edge_execution_order();

    // The first stackitem is always an empty hypergraph and the
    // ProjectedGraph.proj_mat matrix from the input.
    std::vector<StackItem> seeds;
    seeds.push_back(StackItem(std::vector<int> (0), ProjMatT(proj.proj_mat), 0));
    parallel_search_impl(seeds, std::vector<std::vector<int> >(), filter_isomorphic);
}

// Shared body of parallel_search and parallel_search_from. The only difference
// between them is where the stack starts and whether any twins are already
// known, so the three phases below - traversal, mates, isomorphism filter -
// cannot drift apart between a fresh run and a resumed one.
void TwinSearch::parallel_search_impl(std::vector<StackItem> &seeds,
                                      const std::vector<std::vector<int> > &initial_twins,
                                      bool filter_isomorphic) {
    std::chrono::high_resolution_clock phase_clock;
    auto phase_start = phase_clock.now();
    auto phase_ms = [&]() {
        const auto now = phase_clock.now();
        const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - phase_start).count();
        phase_start = now;
        return static_cast<std::int64_t>(ms);
    };

    tbb::concurrent_vector<std::vector<int> > concurrent_twins;
    for (const std::vector<int> &t : initial_twins)
        concurrent_twins.push_back(t);

    // Note: In the parallel version this is not really implementing
    // a stack, it is actually a queue.
    std::vector<StackItem> stack = seeds;
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

    // Move rather than copy. These vectors are the bulk of the twin storage and
    // concurrent_twins is dead afterwards, so copying held two full sets alive
    // at once for no reason.
    ms_traversal = phase_ms();
    rss_mark("after_traversal", concurrent_twins.size());
    twins.reserve(concurrent_twins.size());
    for (std::vector<int> &hg : concurrent_twins)
        twins.push_back(std::move(hg));
    concurrent_twins.clear();

    if (drain_state && drain_state->drain.load(std::memory_order_relaxed)) {
        drained = true;
        frontier.assign(drain_state->parked.begin(), drain_state->parked.end());
        drain_state.reset();
        return;
    }

    // Line graphs and bipartites are each read by exactly ONE phase, and the
    // phases are sequential: the mates test never looks at a bipartite, and the
    // isomorphism filter never looks at a line graph. Building both up front
    // therefore held two full-size graph vectors alive simultaneously, which on
    // these cells is where the memory goes - a mid-curve sample can produce
    // ~460k twins, and each graph heap-allocates per vertex.
    //
    // Build, use, free, then build the next. Peak becomes
    //   twins + max(line_graphs, bipartites)
    // instead of
    //   twins + line_graphs + bipartites.
    // The cost is recomputing each hypergraph's incidence matrix twice, which
    // is trivial next to the graph it feeds.
    {
        std::vector<UndirectedGraph> line_graphs(twins.size());
        tbb::parallel_for(tbb::blocked_range<int>(0, twins.size()),
                           [&](tbb::blocked_range<int> r) {
            for (int i=r.begin(); i<r.end(); ++i)
                line_graphs[i] = compute_linegraph(twins[i]);
        });
        rss_mark("line_graphs_built", twins.size());
        mates = run_mates_tests_parallel(line_graphs, &mates_stats, skip_vf2_for_measurement);
        rss_mark("after_mates", twins.size());
    }
    ms_mates = phase_ms();
    rss_mark("line_graphs_freed", twins.size());

    if (filter_isomorphic) {
        // Only built when the filter actually runs. Previously these were
        // allocated unconditionally even with filtering off.
        std::vector<UndirectedGraph> bipartites(twins.size());
        tbb::parallel_for(tbb::blocked_range<int>(0, twins.size()),
                           [&](tbb::blocked_range<int> r) {
            for (int i=r.begin(); i<r.end(); ++i)
                bipartites[i] = compute_bipartite(twins[i]);
        });
        rss_mark("bipartites_built", twins.size());
        std::vector<int> to_filter = TwinSearch::run_iso_tests_parallel(bipartites, &iso_stats, skip_vf2_for_measurement);
        rss_mark("after_iso", twins.size());
        for(std::size_t i = 0; i < twins.size(); i++) {
            if (to_filter[i] == 0)
                filtered_twins.push_back(i);
        }
        ms_iso = phase_ms();
    }
}

void TwinSearch::parallel_process_item(StackItem &s, tbb::concurrent_vector<std::vector<int> > &concurrent_twins, std::vector<StackItem> &tmp_stack) {
    if (drain_state && drain_check(s))
        return;
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
        std::vector<int> e = fact.node_map.at(enode_id);
        // skip to the next unsatisfied edge
        while (s.proj_rem(e[0], e[1]) < 1) {
            s.edge_execution_index++;
            if (s.edge_execution_index >= edge_execution_order.size())
                return;
            enode_id = edge_execution_order[s.edge_execution_index];
            e = fact.node_map.at(enode_id);
        }
        // Streamed rather than materialised: see the note in process_item.
        std::vector<int> neighbors_vect = get_filtered_neighbors(s, enode_id);
        for_each_combination(neighbors_vect, s.proj_rem(e[0], e[1]),
            [&](const std::vector<int> &comb) { add_to_stack(s, comb, tmp_stack); });
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
std::vector<int> TwinSearch::run_iso_tests_parallel(std::vector<UndirectedGraph> &bipartites,
                                                   PairStats *stats, bool skip_vf2) {
    std::vector<int> to_filter(bipartites.size());
    if (bipartites.size() < 2)
        return to_filter;

    std::vector<GraphFingerprint> fps = compute_fingerprints(bipartites);
    // Per-thread, summed at the end. A shared counter incremented T^2/2 times
    // would be measuring its own contention.
    tbb::enumerable_thread_specific<PairStats> local;

    // One parallel_for over j, rather than a serial loop over i each spawning a
    // parallel_for over j. The old shape put an implicit barrier after every i -
    // N barriers for N graphs - with the work per inner loop shrinking to
    // nothing as i grew, so the late iterations paid full task-spawn and
    // barrier cost to do almost nothing. The fingerprint pre-filter made that
    // worse, not better: it made the typical inner iteration so cheap that the
    // overhead dominated it.
    //
    // Parallelising over j instead also removes the write sharing: to_filter[j]
    // is now touched only by the task that owns j.
    //
    // Same result. The old loop skipped any i already filtered, so it marked j
    // exactly when some UNFILTERED i < j was isomorphic to it; this marks j when
    // ANY i < j is. Those agree because isomorphism is transitive: let i0 be the
    // smallest index isomorphic to j. If i0 were itself filtered there would be
    // an i' < i0 isomorphic to i0 and hence to j, contradicting minimality. So
    // i0 is unfiltered and the old loop marked j at i = i0.
    tbb::parallel_for(std::size_t(1), bipartites.size(), [&](std::size_t j){
        PairStats &ls = local.local();
        for (std::size_t i = 0; i < j; i++) {
            ls.pairs++;
            if (!fingerprints_match(fps[i], fps[j]))
                continue;
            ls.fp_match++;
            if (skip_vf2)
                continue;
            my_callback<UndirectedGraph, UndirectedGraph> mc(bipartites[i], bipartites[j]);
            if ( boost::vf2_graph_iso(bipartites[i], bipartites[j], mc) ) {
                ls.vf2_true++;
                to_filter[j] = 1;
                return;              // one witness is enough
            }
        }
    });

    if (stats) {
        *stats = PairStats();
        for (const PairStats &ls : local) {
            stats->pairs += ls.pairs;
            stats->fp_match += ls.fp_match;
            stats->vf2_true += ls.vf2_true;
        }
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
std::vector<std::vector<int> > TwinSearch::run_mates_tests_parallel(std::vector<UndirectedGraph> &line_graphs,
                                                   PairStats *stats, bool skip_vf2) {
    tbb::concurrent_vector<std::vector<int> > mate_pairs;
    if (line_graphs.size() < 2)
        return std::vector<std::vector<int> >(0);

    std::vector<GraphFingerprint> fps = compute_fingerprints(line_graphs);
    tbb::enumerable_thread_specific<PairStats> local;

    // One parallel_for over j rather than a barrier per i; see the note in
    // run_iso_tests_parallel. Every pair i < j is still tested, and unlike the
    // isomorphism filter there is no early exit, because every mate pair is
    // wanted rather than one witness.
    tbb::parallel_for(std::size_t(1), line_graphs.size(), [&](std::size_t j){
        PairStats &ls = local.local();
        for (std::size_t i = 0; i < j; i++) {
                ls.pairs++;
                if (!fingerprints_match(fps[i], fps[j]))
                    continue;
                ls.fp_match++;
                if (skip_vf2)
                    continue;
                // NOTE: mc must be a plain local. It was previously
                // thread_local, which constructs it once per thread and then
                // leaves it holding references to whichever two graphs that
                // thread happened to see first. Harmless only because the
                // callback never reads them.
                my_callback<UndirectedGraph, UndirectedGraph> mc(line_graphs[i], line_graphs[j]);
                if ( boost::vf2_graph_iso(line_graphs[i], line_graphs[j], mc) ) {
                    ls.vf2_true++;
                    // If line graphs are isomorphic, i and j are a pair of mates
                    mate_pairs.push_back( std::vector<int> {static_cast<int> (i), static_cast<int> (j)});
                }
        }
    });

    if (stats) {
        *stats = PairStats();
        for (const PairStats &ls : local) {
            stats->pairs += ls.pairs;
            stats->fp_match += ls.fp_match;
            stats->vf2_true += ls.vf2_true;
        }
    }

    // put in an std vector for return
    std::vector<std::vector<int> > ret(mate_pairs.size());
    for (std::size_t i = 0; i < mate_pairs.size(); i++)
        ret[i] = mate_pairs[i];

    // Sorted so the result does not depend on thread scheduling. It never did
    // before either - a concurrent_vector filled from a parallel_for is in
    // completion order - but the old shape at least grouped pairs by i, and
    // callers that write these out deserve a stable order.
    std::sort(ret.begin(), ret.end());

    return ret;
}

// Static function that takes a factor graph and a projection (presumably the
// one that constructed the factor graph, but this is NOT tested - results in
// the case where the projection comes from elsewhere are undefined) and computes
// the product of maximum widths of the TwinSearch tree over this factor graph.
long double TwinSearch::compute_width_product(ProjectedGraph &proj, FactorGraph &fact) {
    // TODO I don't have a way of checking for overflow here
    // NOTE: starts at 1 and always multiplies. It previously started at 0 and
    // used `if (prod < 1) prod = ...` to seed itself on the first edge, which
    // also re-seeded rather than multiplied whenever a factor was 0 - so an
    // unsatisfiable edge was silently discarded instead of zeroing the product.
    // For a feasible projection every factor is at least 1, so this is
    // equivalent there.
    long double prod = 1.0;
    int degree;
    int weight;
    std::vector<int> edge;
    for(int eid = 0; eid < fact.num_edge_nodes; eid++) {
        // Get num neighbors
        degree = fact.neighborhood_size(eid);
        edge = fact.node_map.at(eid);
        weight = proj.proj_mat(edge[0], edge[1]);
        prod *= binom(degree, weight);
    }

    return prod;
}

// Static function that takes a factor graph and a projection (presumably the
// one that constructed the factor graph, but this is NOT tested - results in
// the case where the projection comes from elsewhere are undefined) and computes
// the product of maximum widths of the TwinSearch tree over this factor graph.
// This computation is carried out in log space in an attempt to avoid
// overflow to the exten possible for large search trees. However, overflow is
// still possible, and not controlled, so results should be checked.
double TwinSearch::compute_log_width_product(ProjectedGraph &proj, FactorGraph &fact) {
    double logsum = 0.0;
    int degree;
    int weight;
    std::vector<int> edge;
    for(int eid = 0; eid < fact.num_edge_nodes; eid++) {
        // Get num neighbors
        degree = fact.neighborhood_size(eid);
        edge = fact.node_map.at(eid);
        weight = proj.proj_mat(edge[0], edge[1]);
        const double width = binom(degree, weight);

        // An edge with more required hyperedges than available clique-neighbours
        // cannot be satisfied, so the whole product is 0 and its log is -inf.
        // Returning immediately keeps that explicit: callers store this in an
        // int, and converting an infinity to an integer type is undefined, so
        // they must check feasibility first rather than discover it here.
        if (width <= 0.0)
            return -std::numeric_limits<double>::infinity();

        // TODO: There could be overflow/imprecision in binom
        logsum += std::log10(width);
    }

    return logsum;
}
