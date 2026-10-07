#ifndef TWIN_SEARCH_H
#define TWIN_SEARCH_H
#include <memory>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <oneapi/tbb.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "projected_graph.hpp"
#include "factor_graph.hpp"
#include "combinations.hpp"

namespace ublas=boost::numeric::ublas;

// Class implementing a search for Gram Mates based on a ProjectedGraph object. 
class TwinSearch {
	public: 
		ProjectedGraph proj;
        FactorGraph fact;

        // a vector of vectors containing all hypergraphs that
        // have proj as their node co-occurrence projection
        std::vector<std::vector<int> > twins;

        // twins with only one representative of
        // each isomorphism class
        std::vector<int> filtered_twins;

        // vector of pairs of indices into twins such
        // that the hypergraphs at each pair are Gram Mates
        std::vector<std::vector<int> > mates;

        // if false, the desired search is impossible
        bool feasible;

        // --- cooperative drain, for checkpointing and work-splitting ---
        //
        // A straggler fails because a single sample cannot be SPLIT, not
        // because 24h is too little compute: resume_cell.sh resumes at sample
        // granularity, so a sample needing 30h never finishes on the short
        // partition however often it is resubmitted.
        //
        // arm_drain(t) makes parallel_search stop expanding after t seconds and
        // PARK every outstanding node instead. What comes back is the frontier
        // of the search tree - a set of independent subtrees whose union is
        // exactly the remaining work - plus the twins found so far. Resuming
        // means re-seeding the stack from the frontier; splitting means handing
        // subtrees to different jobs.
        //
        // When this fires the twin set is PARTIAL, so parallel_search returns
        // without running the mate tests or the isomorphism filter: both range
        // over the whole set and would be wrong computed from a prefix. Check
        // `drained` before touching `mates` or `filtered_twins`.
        void arm_drain(double seconds);
        bool drained = false;
        // Parked partial hypergraphs, as clique-node id vectors. proj_rem is
        // deliberately NOT kept: it is recomputable from the projection and the
        // chosen cliques, and it is several times the size.
        std::vector<std::vector<int> > frontier;

        // if true, only return twins with diagonal entries that match the
        // input projection. Otherwise, ignore the diagonal.
        bool use_diagonal;

        // Main constructor
        TwinSearch(ProjectedGraph proj, int min_k, int max_k, bool filter_isomorphic, bool parallel, bool run_search, bool use_diag);

        // Constructor with parallel, run_search defaulted to true and use_diag defaulted to false
        TwinSearch(ProjectedGraph proj, int min_k, int max_k, bool filter_isomorphic) : TwinSearch(proj, min_k, max_k, filter_isomorphic, true, true, false) {};

        // Default constructor
        TwinSearch() {};

        // Function signatures
        bool test_feasibility();
        void print_twins(const std::vector<std::vector<int> > &twins);
        void search(bool);
        void parallel_search(bool);
        // Declaring this function static for easier testing access and potential multi-use
        static std::vector<int> run_iso_tests_parallel(std::vector<UndirectedGraph> &bipartites);
        static std::vector<std::vector<int> > run_mates_tests_parallel(std::vector<UndirectedGraph > &line_graphs);
        static long double compute_width_product(ProjectedGraph &proj, FactorGraph &fact);
        static double compute_log_width_product(ProjectedGraph &proj, FactorGraph &fact);
        std::vector<std::vector<int> > inflate_cnodes(const std::vector<int> &cnode_ids);
        // NOTE: This does not need to be public, but I wanted
        // to test it explicitly
        static bool equivalent_lg(std::vector<ublas::matrix<int> > &line_graphs, const int i, const int j);


    private:
        // A simple struct to store a partial hypergraph
        // and its projection remainder to use with the stack.
        struct StackItem;

        // Cheap isomorphism invariant used to avoid calling
        // boost::vf2_graph_iso on pairs that cannot possibly be isomorphic.
        struct GraphFingerprint;
        static GraphFingerprint compute_fingerprint(const UndirectedGraph &g);
        static std::vector<GraphFingerprint> compute_fingerprints(const std::vector<UndirectedGraph> &graphs);
        static bool fingerprints_match(const GraphFingerprint &a, const GraphFingerprint &b);
        //struct ParaReturn;
        std::vector<std::size_t> edge_execution_order;
        // Held by shared_ptr so TwinSearch stays copy-assignable:
        // count_twins_random does `twins = TwinSearch(...)`, and both
        // std::atomic and TBB's enumerable_thread_specific would delete that.
        struct DrainState;
        std::shared_ptr<DrainState> drain_state;
        bool drain_check(const StackItem &s);
        void process_item(std::vector<StackItem> &stack, StackItem &s, std::vector<UndirectedGraph> &bipartites, std::vector<GraphFingerprint> &fingerprints, std::vector<UndirectedGraph > &line_graphs, bool filter_isomorphic);
        void parallel_process_item(StackItem &s, tbb::concurrent_vector<std::vector<int> >&, std::vector<StackItem> &tmp_stack);
        void compute_bipartite_and_linegraph(std::vector<UndirectedGraph> &bipartites, std::vector<UndirectedGraph > &line_graphs, const int i, std::vector<int> &hypergraph);
        // Single-structure variants. parallel_search needs the two in separate
        // phases, so holding both at once is pure peak memory; each rebuilds
        // the small incidence matrix, which is far cheaper than the graph it
        // is used to make.
        UndirectedGraph compute_bipartite(const std::vector<int> &hypergraph);
        UndirectedGraph compute_linegraph(const std::vector<int> &hypergraph);
        bool is_isomorphic(const std::vector<UndirectedGraph> &bipartites, const std::vector<GraphFingerprint> &fingerprints, const int cand_idx);
        void add_to_stack(const StackItem &curr, const std::vector<int> &comb, std::vector<StackItem> &stack);
        std::vector<int> get_filtered_neighbors(const StackItem &s, int enode_id);
        int choose_edge(const StackItem &s);
        int choose_random_edge(const StackItem &s);
        std::vector<std::size_t> compute_edge_execution_order();
        std::vector<std::size_t> default_edge_execution_order();
};

#endif
