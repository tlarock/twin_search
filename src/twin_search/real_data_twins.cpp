#include <iostream>
#include <fstream>
#include <chrono>
#include <iterator>
#include <oneapi/tbb.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>
#include <boost/algorithm/string.hpp>
#include "argparse/argparse.hpp"

#include "hypergraph.hpp"
#include "factor_graph.hpp"
#include "projected_graph.hpp"
#include "twin_search.hpp"
#include "random_hypergraph_generators.hpp"
#include "utils.hpp"




// This file implements a search over random hypergraphs of various sizes
// to find projection matrices that have gram mates. The search is run in
// parallel over all arguments and samples using tbb::parallel_for_each.

// Mutex for writing to file in one_sample_write
tbb::spin_mutex WriteMutex;

// struct to read and store command line arguments with argparse
struct MyArgs : public argparse::Args {
    int &min_k = kwarg("min-k", "Minimum allowed hyperedge size for mate search").set_default(2);
    int &max_k = kwarg("max-k", "Maximum allowed hyperedge size for mate search").set_default(0);
    std::string &input_filepath = kwarg("input-file", "Input filepath. No default.");
    std::string &filepath = kwarg("output-path", "Output filepath. Default is results/real_data/").set_default("results/real_data/");
    int &max_threads = kwarg("max-threads, tbb-max-allowed-parallelism", "Value passed to TBB to limit number of worker threads. If <= 0, ignored.").set_default(0);
};

// Takes a filename string and reads in a hypergraph. Filename
// is assumed to point to a text file, comma-delimited with
// a sorted hyperedge on each line.
Hypergraph read_data(std::string filename) {
    std::ifstream infile;
    infile.open(filename, std::ios::in);
    std::string line;
    std::vector<std::vector<int> > hyperedges;
    std::vector<std::string> strs;
    std::vector<int> he;
    while (std::getline(infile, line)) {
        boost::split(strs, line, boost::is_any_of(","));
        he = std::vector<int>(strs.size());
        for (size_t i = 0; i < strs.size(); i++) {
            he[i] = std::stoi( strs[i] );
        }

        // Add hyperedge to hyperedges
        hyperedges.push_back(std::vector<int>(he));
    }


    infile.close();

    return Hypergraph(hyperedges);
}

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

// Write a function analogous to one_sample_write from count_twins_random
void run_search(std::string filename, const Hypergraph &H, const int min_k, const int max_k) {
    std::map<std::vector<int>, std::vector<int> > twins_counts;
    int num_cliques = 0;
    int num_edges = 0;
    int num_mate_pairs = 0;
    std::int64_t runtime = 0;
    int max_log_width = 0; 
    std::chrono::high_resolution_clock time;

    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    if (!outfile.is_open()) {
        std::cout << "Couldn't open file: " + filename << ". Exiting." << std::endl;
        return;
    } else {
        std::cout << "Output will be written to: " << filename << std::endl;
    }

    // Construct the projected graph
    std::cout << "Constructing projection...";
    ProjectedGraph proj(H);
    std::cout << "Done." << std::endl;

    // Run the twins search
    std::cout << "Creating TwinSearch object with min_k " << min_k << " and max_k " << max_k << "...";
    TwinSearch twins(proj, min_k, max_k, true, true, false, false);
    std::cout << "Done." << std::endl;

    max_log_width = std::round(twins.compute_log_width_product(proj, twins.fact));
    std::cout << "Running mate search on dataset. Exponent of the width product: " << max_log_width << std::endl;
    std::chrono::time_point start = time.now();
    twins.parallel_search(true);
    const auto end = time.now();
    auto dur = end - start;
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur);
    runtime = elapsed_ms.count();
    std::cout << "Done with search. ";
    // Loop over the filtered twins and initialize the output data
    // structure for each size distribution we've encountered
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

    // Loop over all unfiltered twins and count the number of
    // unfiltered twins and pairs of mates
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

    std::cout << "Number of twins: " << twins.twins.size() << " filtered twins: " << twins.filtered_twins.size() << std::endl;

    // Write to output file
    if (outfile.is_open()) {
        outfile << H.n << "," << H.m << "," << runtime << "," << max_log_width << "," << num_edges << "," << num_cliques << "," << num_mate_pairs << "/";
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

    } else {
        std::cout << "Output file not open." << std::endl;
    }

    outfile.close();
}


int main(int argc, char *argv[]) {
	// Set cout to be unbuffered
	std::cout.setf(std::ios::unitbuf);
    // Read arguments
    auto args = argparse::parse<MyArgs>(argc, argv);
    const int min_k_in = args.min_k;
    const int max_k_in = args.max_k;
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

    // Read the input file
    std::cout << "Reading input...";
    Hypergraph H = read_data(args.input_filepath);
    std::cout << "Done." << std::endl;

    // Set maximum k if needed
    int max_k = max_k_in;
    if (max_k_in <= 0) {
        max_k = H.n;
    }
    run_search(args.input_filepath + ".twins", H, min_k_in, max_k);
}
