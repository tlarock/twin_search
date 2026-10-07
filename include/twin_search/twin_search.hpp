#ifndef TWIN_SEARCH_H
#define TWIN_SEARCH_H
#include <cstdint>
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
        // Three ways to ask for a drain; all land in the same mechanism.
        //
        //   arm_drain(seconds)       wall-clock budget. <= 0 means no deadline,
        //                            i.e. arm for the signal only.
        //   arm_drain_nodes(n)       after n nodes PER THREAD. Deterministic at
        //                            one thread, which is what the tests use -
        //                            a wall-clock trigger would make them flaky.
        //   request_global_drain()   from a signal handler. Free function below;
        //                            an armed search polls it.
        //
        // An unarmed search never checks anything, so none of this costs
        // production runs. See docs/reuse-experiments.md section 5.
        void arm_drain(double seconds);
        void arm_drain_nodes(std::uint64_t nodes_per_thread);
        bool drained = false;
        // Parked partial hypergraphs, as clique-node id vectors. proj_rem is
        // deliberately NOT kept: it is recomputable from the projection and the
        // chosen cliques, and it is several times the size.
        std::vector<std::vector<int> > frontier;

        // Re-seed from a drained frontier and continue. Both arguments are
        // clique-node id vectors, as `frontier` and `twins` are produced.
        // Returns false if the frontier does not fit this projection - a
        // corrupt or mismatched checkpoint - rather than searching garbage.
        bool parallel_search_from(const std::vector<std::vector<int> > &frontier_cnodes,
                                  const std::vector<std::vector<int> > &twins_so_far,
                                  bool filter_isomorphic);

        // Clique-node id for a hyperedge, or -1 if the projection has no such
        // clique. Checkpoints store hyperedges rather than ids, because ids
        // depend on FactorGraph's construction order and a format should not
        // rest on that.
        int cnode_for(const std::vector<int> &hyperedge) const;

        // if true, only return twins with diagonal entries that match the
        // input projection. Otherwise, ignore the diagonal.
        bool use_diagonal;

        // Main constructor
        TwinSearch(ProjectedGraph proj, int min_k, int max_k, bool filter_isomorphic, bool parallel, bool run_search, bool use_diag);

        // Constructor with parallel, run_search defaulted to true and use_diag defaulted to false
        TwinSearch(ProjectedGraph proj, int min_k, int max_k, bool filter_isomorphic) : TwinSearch(proj, min_k, max_k, filter_isomorphic, true, true, false) {};

        // Default constructor
        TwinSearch() {};

        // How the two O(T^2) phases actually spend themselves.
        //
        // Neither phase does T^2/2 vf2 calls: a (|V|, |E|, degree-sequence)
        // pre-filter runs first and vf2 only sees pairs that survive it. Which
        // of the two dominates decides how much a canonical-form rewrite would
        // buy, and it is not inferable from the totals.
        struct PairStats {
            std::uint64_t pairs = 0;      // (i, j) examined
            std::uint64_t fp_match = 0;   // survived the pre-filter -> vf2 called
            std::uint64_t vf2_true = 0;   // vf2 said isomorphic
        };
        PairStats iso_stats, mates_stats;

        // Run the pre-filter but never call vf2. The RESULTS are wrong; the
        // point is that timing a run with and without it attributes the phase
        // between the cheap scan and the expensive confirmations, without
        // putting a clock inside a loop that executes T^2/2 times.
        bool skip_vf2_for_measurement = false;

        // Wall time of each phase of parallel_search, in milliseconds. -1 if
        // the phase did not run.
        std::int64_t ms_traversal = -1;
        std::int64_t ms_mates = -1;
        std::int64_t ms_iso = -1;

        // Function signatures
        bool test_feasibility();
        void print_twins(const std::vector<std::vector<int> > &twins);
        void search(bool);
        void parallel_search(bool);
        // Declaring this function static for easier testing access and potential multi-use
        static std::vector<int> run_iso_tests_parallel(std::vector<UndirectedGraph> &bipartites,
                                                      PairStats *stats = nullptr,
                                                      bool skip_vf2 = false);
        static std::vector<std::vector<int> > run_mates_tests_parallel(std::vector<UndirectedGraph > &line_graphs,
                                                      PairStats *stats = nullptr,
                                                      bool skip_vf2 = false);
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

        // Cheap isomorphism invariants used to avoid calling
        // boost::vf2_graph_iso on pairs that cannot possibly be isomorphic.
        struct GraphFingerprint;
        // 1-WL / colour-refinement hash. See the definition in the .cpp for
        // why it is sound and why the degree sequence alone is not enough.
        static std::uint64_t colour_refinement_hash(const UndirectedGraph &g);
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
        // proj.proj_mat minus the deltas of every clique in `hypergraph`.
        // False if any entry would go negative, which means the partial
        // hypergraph is not a valid prefix for this projection.
        bool residual_for(const std::vector<int> &hypergraph, ProjMatT &out) const;
        void parallel_search_impl(std::vector<StackItem> &seeds,
                                  const std::vector<std::vector<int> > &initial_twins,
                                  bool filter_isomorphic);
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

// Drain request from a signal handler.
//
// A handler must do essentially nothing, so it sets this flag and returns; the
// search notices on its next poll. Separate from TwinSearch because one process
// runs many samples concurrently and a signal applies to all of them.
void twin_search_request_global_drain();
bool twin_search_global_drain_requested();
void twin_search_clear_global_drain();

#endif
