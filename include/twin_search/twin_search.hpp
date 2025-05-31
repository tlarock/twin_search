#ifndef TWIN_SEARCH_H
#define TWIN_SEARCH_H
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <oneapi/tbb.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <discreture.hpp>
#include "projected_graph.hpp"
#include "factor_graph.hpp"

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
        //struct ParaReturn;
        std::vector<std::size_t> edge_execution_order;
        void process_item(std::vector<StackItem> &stack, StackItem &s, std::vector<UndirectedGraph> &bipartites, std::vector<UndirectedGraph > &line_graphs, bool filter_isomorphic);
        void parallel_process_item(StackItem &s, tbb::concurrent_vector<std::vector<int> >&, std::vector<StackItem> &tmp_stack);
        void compute_bipartite_and_linegraph(std::vector<UndirectedGraph> &bipartites, std::vector<UndirectedGraph > &line_graphs, const int i, std::vector<int> &hypergraph);
        bool is_isomorphic(const std::vector<UndirectedGraph> &bipartites, const int cand_idx);
        std::vector<std::vector<int> > get_combinations(int enode_id, int weight, const StackItem &s);
        void add_to_stack(const StackItem &curr, const std::vector<int> &comb, std::vector<StackItem> &stack);
        std::vector<int> get_filtered_neighbors(const StackItem &s, int enode_id);
        int choose_edge(const StackItem &s);
        int choose_random_edge(const StackItem &s);
        std::vector<std::size_t> compute_edge_execution_order();
        std::vector<std::size_t> default_edge_execution_order();
};

#endif
