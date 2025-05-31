#include <iostream>
#include <fstream>
#include <chrono>
#include <iterator>
#include <oneapi/tbb.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>
#include <discreture.hpp>

#include "argparse/argparse.hpp"

#include "hypergraph.hpp"
#include "factor_graph.hpp"
#include "projected_graph.hpp"
#include "twin_search.hpp"
#include "random_hypergraph_generators.hpp"
#include "utils.hpp"


// This file implements an exhaustive search over k-uniform random hypergraphs
// where the total space of possible hypergraphs can be enumerated. The first
// goal is to implement a search for n=6 and m=16, which has 4845 hypergraphs.

// Mutex for writing to file in one_sample_write
tbb::spin_mutex WriteMutex;

// struct to read and store command line arguments with argparse
struct MyArgs : public argparse::Args {
    int &n = kwarg("n, num-nodes", "Number of nodes.").set_default(6);
    int &m = kwarg("m", "Number of hyperedges").set_default(16);
    int &k = kwarg("k", "Hyperedge dimension").set_default(3);
    int &min_k = kwarg("min-k", "Mininmum hyperedge dimension in search. Default 2").set_default(2);
    int &max_k = kwarg("max-k", "Maximum hyperedge dimension in search. Default -1 sets to n.").set_default(-1);
    bool &use_diagonal = flag("use-diagonal", "If given, use the diagonal entries when finding the unique projections up to isomorphism.");
    std::string &filepath = kwarg("output-path", "Output filepath. Default is results/increasing_density/").set_default("results/increasing_density/");
    int &max_threads = kwarg("max-threads, tbb-max-allowed-parallelism", "Value passed to TBB to limit number of worker threads. If <= 0, ignored.").set_default(0);
};

void write_all_twins(TwinSearch &twins, std::ofstream &outfile) {
    // write all of the twins
    std::size_t printed_twins = 0;
    for (std::size_t twin_id = 0; twin_id < twins.twins.size(); twin_id++) {
        std::vector<std::vector<int > > inflated = twins.inflate_cnodes(twins.twins[twin_id]);
        int printed_hyperedges = 0;
        for (std::vector<int> he : inflated) {
            for (std::size_t i = 0; i < he.size(); i++) {
                outfile << he[i];
                if (i < he.size()-1)
                    outfile <<  ":";
                else if (i == he.size()-1 && printed_hyperedges < static_cast<int> (inflated.size()-1))
                    outfile << "|";
            }
            printed_hyperedges += 1;
        }
        if (printed_twins < twins.twins.size()-1)
            outfile << ";";
        printed_twins += 1;
    }
}

bool one_sample_write(unsigned int n, unsigned int m, unsigned int min_k, unsigned int max_k, ProjectedGraph &proj, std::string filename, bool use_diag) {
    int num_cliques = 0;
    int num_edges = 0;
    std::int64_t runtime = 0;
    int max_log_width = 0;
    std::chrono::high_resolution_clock time;
    if (n != static_cast <unsigned int> (proj.proj_mat.size1()))
        std::cout << "input n: " << n << " does not match size of proj: " << proj.proj_mat.size1() << std::endl;

    // A map from a size distribution represented as a vector with entries
    // corresponding to min_k,...,max_k pointing to a 2-entry vector consisting
    // of the number of twins and filtered_twins
    std::map<std::vector<int>, std::vector<int> > twins_counts;
    int num_mate_pairs = 0;

    // Compute twins and mates
    TwinSearch twins;
    std::chrono::time_point start = time.now();
    twins = TwinSearch(proj, min_k, max_k, true, true, true, use_diag);
    const auto end = time.now();
    auto dur = end - start;
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur);
    runtime = elapsed_ms.count();

    // Use TwinSearch.fact to compute max width
    max_log_width = std::round(twins.compute_log_width_product(proj, twins.fact));

    // Loop over the filtered twins and count the number that correspond to
    // each unique hyperedge size distribution
    std::vector<int> size_dist(max_k-min_k+1);
    std::vector<std::vector<int> > mate;
    for (int i : twins.filtered_twins) {
        mate = twins.inflate_cnodes(twins.twins[i]);
        // Get size distribution
        for (std::size_t he_idx = 0; he_idx < mate.size(); he_idx++) {
            size_dist[mate[he_idx].size()-min_k] += 1;
        }

        if (!twins_counts.contains(size_dist)) {
            twins_counts[size_dist] = std::vector<int> {0, 1};
        } else {
            twins_counts[size_dist][1] += 1;
        }

        // Reset size distribution to 0
        for (std::size_t j = 0; j < size_dist.size(); j++) { size_dist[j] = 0; }
    }

    // Loop over all unfiltered twins and count the number that correspond to
    // each unique hyperedge size distribution
    for (std::size_t i = 0; i < twins.twins.size(); i++) {
        mate = twins.inflate_cnodes(twins.twins[i]);
        // Get size distribution
        for (std::size_t he_idx = 0; he_idx < mate.size(); he_idx++) {
            size_dist[mate[he_idx].size()-min_k] += 1;
        }

        if (!twins_counts.contains(size_dist)) {
            twins_counts[size_dist] = std::vector<int> {1, 0};
        } else {
            twins_counts[size_dist][0] += 1;
        }

        // Reset size distribution to 0
        for (std::size_t j = 0; j < size_dist.size(); j++) { size_dist[j] = 0; }
    }
    num_cliques = twins.fact.num_clique_nodes;
    num_edges = twins.fact.num_edge_nodes;
    num_mate_pairs = static_cast<int> (twins.mates.size());

    // Create local scope for mutex
    {
        // Get mutex lock
        tbb::spin_mutex::scoped_lock lock(WriteMutex);
        std::ofstream outfile;
        outfile.open(filename, std::ios::out | std::ios::app);

        // Write to output file
        if (outfile.is_open()) {
            outfile << n << "," << m << "," << runtime << "," << max_log_width << "," << num_edges << "," << num_cliques << "," << num_mate_pairs << "/";
            for (const auto& [size_dist, stats_vect] : twins_counts) {
                // Write the size_distribution/
                std::string pairs_str = "";
                for (std::size_t i = 0; i < size_dist.size(); i++) {
                    if (size_dist[i] > 0)
                        pairs_str += std::to_string(i+min_k) + ":" + std::to_string(size_dist[i]) + ",";
                }
                pairs_str.pop_back();
                // Write num_unf,num_filt/
                outfile << pairs_str << "/";
                outfile << std::to_string(stats_vect[0]) << "," << std::to_string(stats_vect[1]) << "/";
            }
            write_all_twins(twins, outfile); 
            outfile << std::endl;
            outfile.close();
        } else {
            std::cout << "Output file not open." << std::endl;
        }
    }

    return true;
}

std::string construct_filename(int n,  int m, int k) {
	return "n-" + std::to_string(n) +
		"_m-" + std::to_string(m) +
		"_k-" + std::to_string(k);
}


std::map<int, std::vector<int> > get_hyperedges(std::vector<int> nodes, int k) {
    auto combs = discreture::combinations(nodes, k);
    std::map<int, std::vector<int> > all_hyperedges;
    int he_idx = 0;
    for (auto&& comb : combs) {
        std::vector<int> new_comb;
        for (int c : comb)
            new_comb.push_back(c);
        sort(new_comb.begin(), new_comb.end());
        all_hyperedges[he_idx] = std::vector<int>(new_comb);
        he_idx += 1;
    }
    return all_hyperedges;
}

bool projections_equal(const ProjectedGraph &proj1, const ProjectedGraph &proj2) {
    for (std::size_t j = 0; j < proj1.proj_mat.size1(); j++) {
        for (std::size_t w = 0; w < proj1.proj_mat.size2(); w++) {
            if (proj1.proj_mat(j, w) - proj2.proj_mat(j, w) != 0)
                return false;
        }
    }
    return true;
}

std::vector<ProjectedGraph> get_unique_projections(std::map<int, std::vector<int> > &all_hyperedges, int n, int m, bool ignore_with_singletons, bool use_diag) {
    std::vector<int> he_ids;
    for(std::size_t i = 0; i < all_hyperedges.size(); i++)
        he_ids.push_back(static_cast <int> (i));

    // Get all combinations of m hyperedges by id
    // concretize into a vector
    auto disc_combs = discreture::combinations(he_ids, m);
    std::vector<std::vector<std::vector<int> > > combs;
    for (std::size_t i = 0; i < disc_combs.size(); ++i) {
        std::vector<std::vector<int> > new_hg;
        for (auto&& c : disc_combs[i]) {
            new_hg.push_back(all_hyperedges[c]);
        }
        combs.push_back(new_hg);
    }

    // Fill all_projections with every projection
    std::vector<ProjectedGraph> all_projections(combs.size());
    std::vector<UndirectedGraph> all_projections_boost(combs.size());
    // Loop over combinations
    std::vector<int> to_filter(all_projections.size());
    tbb::parallel_for(static_cast<std::size_t>(0), static_cast<std::size_t>(combs.size()), [&](std::size_t proj_idx){
        Hypergraph h(combs[proj_idx]);
        if (!ignore_with_singletons || h.n == n) {
            // Get the projection
            all_projections[proj_idx] = ProjectedGraph(combs[proj_idx], !use_diag);
            all_projections_boost[proj_idx] = all_projections[proj_idx].get_boost_graph();
        } else {
            to_filter[proj_idx] = 1;
        }
    });
    
    std::vector<int> iso_filter = TwinSearch::run_iso_tests_parallel(all_projections_boost);
    // Filter out everything in to_filter
    std::vector<ProjectedGraph> out;
    for (std::size_t i = 0; i < to_filter.size(); i++) {
        if (to_filter[i] < 1 && iso_filter[i] < 1) {
            out.push_back(all_projections[i]);
        }
    }

    return out;
}

int main(int argc, char *argv[]) {
	// Set cout to be unbuffered
	std::cout.setf(std::ios::unitbuf);
    // Read arguments
    auto args = argparse::parse<MyArgs>(argc, argv);
    const unsigned int n = args.n;
    const unsigned int m = args.m;
    const unsigned int k = args.k;
    const int min_k = args.min_k;
    const bool use_diag = args.use_diagonal;

    int max_k = args.max_k;
    if (max_k <= 0)
        max_k = n;

    const int max_threads = args.max_threads;

    // set up number of threads for TBB
    int default_threads = oneapi::tbb::global_control::active_value(oneapi::tbb::global_control::max_allowed_parallelism);
    int num_threads = default_threads;
    if (max_threads > default_threads) {
        std::cout << "Argument max_threads > default_threads. Ignoring.";
    } else if ( (max_threads > 0) && ( max_threads < default_threads) ) {
        num_threads = max_threads;
    }
    oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, num_threads);
    int thread_max = oneapi::tbb::global_control::active_value(oneapi::tbb::global_control::max_allowed_parallelism);
    if (thread_max != num_threads) {
        std::cout << "Potential issue with tbb max_allowed_parallelism. thread_max: " << thread_max << " num_threads: " << num_threads << std::endl;
    }
 
    // Construct output filename
    std::string filepath = args.filepath;
    std::string filename = filepath + construct_filename(n, m, k);

    if (min_k != max_k)
        filename += "_non-uniform";

    if (min_k > 2 && min_k != max_k)
        filename += "_min-k-" + std::to_string(min_k);

    if (use_diag)
        filename += "_with-diag";

    filename += "_exhaustive_projections.csv";

    // Construct the set of nodes
    std::cout << "n: " << n << " m: " << m << std::endl;
    std::vector<int> nodes;
    for(unsigned int i = 0; i < n; i++)
        nodes.push_back(i);

    // Construct all possible hypergraphs with discreture
    // First, get all possible hyperedges
    std::map<int, std::vector<int> > all_hyperedges = get_hyperedges(nodes, k);
    std::cout << "Num hyperedges: " << all_hyperedges.size() << std::endl;
    std::cout << "Total hypergraphs: " << binom(static_cast <double> (all_hyperedges.size()), m) << std::endl;

    // Second, get all of the hypergraphs and their projections
    std::vector<ProjectedGraph> all_projections = get_unique_projections(all_hyperedges, n, m, true, use_diag);

    std::cout << "Number of unique projections: " << all_projections.size() << std::endl;

    // touch/empty the output file
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    if (!outfile.is_open()) {
        std::cout << "Couldn't open file: " + filepath << ". Exiting." << std::endl;
        return 0;
    }
    outfile.close();

    std::cout << "Running mate search on each projection." << std::endl;
    tbb::parallel_for(size_t(0), all_projections.size(),
            [&](size_t i) {
            one_sample_write(n, m, min_k, max_k, all_projections[i], filename, use_diag);
    });
}
