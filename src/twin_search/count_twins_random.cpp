#include <iostream>
#include <fstream>
#include <chrono>
#include <iterator>
#include <oneapi/tbb.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>

#include "argparse/argparse.hpp"

#include "hypergraph.hpp"
#include "factor_graph.hpp"
#include "projected_graph.hpp"
#include "twin_search.hpp"
#include "random_hypergraph_generators.hpp"
#include "utils.hpp"


// This file implements a search over random hypergraphs of various sizes,
// computes their projections, then finds the set of twins. The search is run in
// parallel over all arguments and samples using tbb::parallel_for_each.

// Mutex for writing to file in one_sample_write
tbb::spin_mutex WriteMutex;

// struct to read and store command line arguments with argparse
struct MyArgs : public argparse::Args {
    bool &config_model = flag("config-model", "If given, sample from configuration model. Otherwise sample from uniform model. Note: If given, -m is ignored and -g is used. If not given, -m is used and -g is ignored.");
    int &n = kwarg("n, num-nodes", "Number of nodes.").set_default(15);
    int &m = kwarg("m", "Number of hyperedges").set_default(15);
    double &gamma = kwarg("g", "Value of gamma for degree distribution.").set_default(3.0);
    int &k = kwarg("k", "Hyperedge dimension").set_default(3);
    int &num_samples = kwarg("samples", "Number of samples").set_default(100);
    int &min_k = kwarg("min-k", "Minimum allowed hyperedge size for mate search. If <=0 (default), automatically set to k.").set_default(0);
    int &max_k = kwarg("max-k", "Maximum allowed hyperedge size for mate search. If <=0 (default), automatically set to k.").set_default(0);
    bool &no_max_k = flag("no-max-k", "Set the maximum k to the maximum integer size. --max-k ignored if given.");
    std::string &filepath = kwarg("output-path", "Output filepath. Default is results/increasing_density/").set_default("results/increasing_density/");
    int &max_threads = kwarg("max-threads, tbb-max-allowed-parallelism", "Value passed to TBB to limit number of worker threads. If <= 0, ignored.").set_default(0);
    bool &sequential_samples = flag("sequential-samples", "If given, do not parallelize over samples. NOTE: parallelization of TwinSearch.*_search controlled by --sequential-twins.");
    bool &sequential_twins = flag("sequential-twins", "If given, run TwinSearch.search rather than TwinSearch.parallel_search. NOTE: parallelization of samples still controlled by --sequential-samples and --max-threads.");
    bool &append = flag("append", "If given, append to the appropriate file if it already exists. Useful for backfilling failed samplings.");
    int &start_sample = kwarg("start-sample", "If --append given, start with this index, meaning the true number of samples will be num_samples - start_sample. If --append not given, just leads to fewer samples in a mislabled file..").set_default(0);
    int &width_limit_exponent = kwarg("width-limit-exp", "If >0, no sample with product of binomials > 10^width-limit-exp will be searched.").set_default(0);
    bool &width_limit_auto = flag("width-limit-auto", "If given, automatically choose a width-limit by computing the width-product of --samples examples and choosing the median. Incompatible with width-limit-exp, which is ignored if this is given.");
    int &width_limit_samples = kwarg("width-limit-samples", "If >0 and --width-limit-auto is given, auto width limit will be set to the median of this number of samples").set_default(0);
    bool &no_max_rejections = flag("infinite-rejections", "If given with width_limit > 0, there will be no limit on the number of samples rejected. Warning: Could lead to infinite loops. No effect if width_limit_exp <= 0.");
};


// Convenience struct to store the parameters of a given search. Helpful
// for passing arguments to the std::for_each function.
struct Params {
    bool config_model;
    int i;
    int n;
    int m;
    double gamma;
    int k;
    int min_k;
    int max_k;
    bool filter_isomorphic;
    bool sequential;
    int width_limit_exp;
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

bool one_sample_write(Params &p, std::ofstream &outfile) {
    // A map from a size distribution represented as a vector with entries
    // corresponding to min_k,...,max_k pointing to a 2-entry vector consisting
    // of the number of twins and filtered_twins
    std::map<std::vector<int>, std::vector<int> > twins_counts;
    int num_cliques = 0;
    int num_edges = 0;
    int num_mate_pairs = 0;
    std::int64_t runtime = 0;
    int max_log_width = 0; 
    std::chrono::high_resolution_clock time;

    // Sample a hypergraph
    Hypergraph h;
    if (!p.config_model)
        h = sample_uniform_random(p.n, p.m, p.k);
    else
        h = uniform_hypergraph_configuration_model(p.n, p.gamma, p.k, p.n / 2);

    // Compute twins and mates
    TwinSearch twins;
    if (h.m > 0) {
        // Construct a projected graph object
        ProjectedGraph proj(h);
        
        // Construct the TwinSearch object
        twins = TwinSearch(proj, p.min_k, p.max_k, p.filter_isomorphic, !p.sequential, false, false);

        if (twins.feasible) {
            // Use TwinSearch.fact to compute max width
            max_log_width = std::round(twins.compute_log_width_product(proj, twins.fact));
            std::chrono::time_point start = time.now();
            if (p.width_limit_exp <= 0 || max_log_width <= p.width_limit_exp) {
                if (!p.sequential)
                    twins.parallel_search(p.filter_isomorphic);
                else
                    twins.search(p.filter_isomorphic);
            } else {
                return false;
            }
            const auto end = time.now();
            auto dur = end - start;
            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur);
            runtime = elapsed_ms.count();
            // Loop over the filtered twins and initialize the output data
            // structure for each size distribution we've encountered
            std::vector<int> size_dist(p.max_k-p.min_k+1);
            std::vector<std::vector<int> > mate;
            for (int i : twins.filtered_twins) {
                mate = twins.inflate_cnodes(twins.twins[i]);
                // Get size distribution
                for (std::size_t he_idx = 0; he_idx < mate.size(); he_idx++) {
                    size_dist[mate[he_idx].size()-p.min_k] += 1;
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
                    size_dist[mate[he_idx].size()-p.min_k] += 1;
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
        } else {
            {
                tbb::spin_mutex::scoped_lock lock(WriteMutex);
                std::cout << "WARNING: Infeasible sample detected. Check min-k and max-k. Writing 0s for this sample." << std::endl;
            }
        }
    }
    // Create local scope for mutex
    {
        // Get mutex lock
        tbb::spin_mutex::scoped_lock lock(WriteMutex);
        // Write to output file
        if (outfile.is_open()) {
            outfile << p.n << "," << p.m << "," << runtime << "," << max_log_width << "," << num_edges << "," << num_cliques << "," << num_mate_pairs << "/";
            for (const auto& [size_dist, stats_vect] : twins_counts) {
                // Write the size_distribution/
                std::string pairs_str = "";
                for (std::size_t i = 0; i < size_dist.size(); i++) {
                    if (size_dist[i] > 0)
                        pairs_str += std::to_string(i+p.min_k) + ":" + std::to_string(size_dist[i]) + ",";
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
    }

    return true;
}

int compute_width_limit(int n, int m, double gamma, int k, int min_k, int max_k, int num_samples, bool config_model) {
    std::vector<int> sampled_width_limits(num_samples);
    Hypergraph h;
    TwinSearch twins;
    ProjectedGraph proj;
    // Sample --samples hypergraphs and compute their width limits
    for(int i = 0; i < num_samples; i++) {
        if (!config_model)
            h = sample_uniform_random(n, m, k);
        else
            h = uniform_hypergraph_configuration_model(n, gamma, k, n / 2);
        
        // Construct a projected graph object
        proj = ProjectedGraph(h);
        
        // Construct the TwinSearch object
        twins = TwinSearch(proj, min_k, max_k, false, false, false, false);

        if (twins.feasible) {
            // Use TwinSearch.fact to compute max width
            sampled_width_limits[i] = std::round(twins.compute_log_width_product(proj, twins.fact));
        } else {
            std::cout << "WARNING: Infeasible sample detected. Check min-k and max-k. Setting a sampled_width_limit to 0." << std::endl;
            sampled_width_limits[i] = 0;
        }
    }

    // Set the width-limit to the median
    auto middle_idx = sampled_width_limits.size() / 2;
    std::nth_element(sampled_width_limits.begin(), sampled_width_limits.begin() + middle_idx, sampled_width_limits.end());
    return sampled_width_limits[middle_idx];
}

std::string construct_uniform_filename(int n,  int m, int k, int num_samples, int min_k, int max_k) {
    if (min_k == k && max_k == k) {
	    return "n-" + std::to_string(n) +
            "_m-" + std::to_string(m) +
            "_k-" + std::to_string(k) +
            "_samples-" + std::to_string(num_samples);
    } else {
        if (max_k <= n) {
            return "n-" + std::to_string(n) +
                "_m-" + std::to_string(m) +
                "_k-" + std::to_string(k) +
                "_min-k-" + std::to_string(min_k) +
                "_max-k-" + std::to_string(max_k) +
                "_samples-" + std::to_string(num_samples);
        } else {
            return "n-" + std::to_string(n) +
                "_m-" + std::to_string(m) +
                "_k-" + std::to_string(k) +
                "_min-k-" + std::to_string(min_k) +
                "_max-k-" + std::to_string(n) +
                "_samples-" + std::to_string(num_samples);
        }
    }
}

std::string construct_config_filename(int n, double gamma, int k, int num_samples, int min_k, int max_k) {
    if (min_k == k && max_k == k) {
	    return "n-" + std::to_string(n) +
            "_k-" + std::to_string(k) +
            "_gamma-" + std::format("{:.2f}", gamma) +
            "_samples-" + std::to_string(num_samples);
    } else {
        return "n-" + std::to_string(n) +
            "_k-" + std::to_string(k) +
            "_min-k-" + std::to_string(min_k) +
            "_max-k-" + std::to_string(max_k) +
            "_gamma-" + std::format("{:.2f}", gamma) +
            "_samples-" + std::to_string(num_samples);
    }
}


int main(int argc, char *argv[]) {
	// Set cout to be unbuffered
	std::cout.setf(std::ios::unitbuf);
    // Read arguments
    auto args = argparse::parse<MyArgs>(argc, argv);
    const bool config_model = args.config_model;
    const int n = args.n;
    const double gamma = args.gamma;
    const int m = args.m;
    const int k = args.k;
    const int num_samples = args.num_samples;
    int min_k = args.min_k;
    int max_k = args.max_k;
    const bool no_max_k = args.no_max_k;
    const int max_threads = args.max_threads;
    const bool sequential = args.sequential_twins;
    const bool sequential_samples = args.sequential_samples;
    const bool append = args.append;
    const int start_sample = args.start_sample;
    // Note: not const because modified if auto_width_limit is true
    int width_limit_exponent = args.width_limit_exponent;
    const bool filter_isomorphic = true;
    const bool auto_width_limit = args.width_limit_auto;
    const bool no_max_rejections = args.no_max_rejections;
    // Note: not const because modified if < 1
    int width_limit_samples = args.width_limit_samples;

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

    // if min/max k are not set, default is k
    if (min_k <= 0)
        min_k = k;
    if (max_k <= 0)
        max_k = k;

    if (no_max_k)
        max_k = n;

    double num_hyperedges;
    // Construct the arguments that will be passed to
    // my function one_sample via tbb:parallel_for_each
    // Fail if the maximum possible hyperedges for some (n,k)
    // Conditions:
    // nCk < m (not enough possible hyperedges)
    if (!config_model) {
        num_hyperedges = binom(n, k);
        if (num_hyperedges <= m) {
            std::cout << "Number of hyperedges binom(n = " << n << ", k = " << k << ") < m (" << num_hyperedges << " < " << m << "). Exiting." << std::endl;
            return 0;
        }
    }

    // Compute auto width limit if necessary
    if (auto_width_limit) {
        if (width_limit_samples < 1)
            width_limit_samples = num_samples*10; // Hard-coded sample size

        width_limit_exponent = compute_width_limit(n, m, gamma, k, min_k, max_k, width_limit_samples, config_model);
        std::cout << "Width Limit Exp: " << width_limit_exponent << std::endl;
    }

    if (!config_model)
        std::cout << "Sampling " << num_samples << " hypergraphs from uniform model with k=" << k << " n=" << n << " m=" << m << " min_k=" << min_k << " max_k=" << max_k << " and width_limit_exp=" << width_limit_exponent << std::endl;
    else
        std::cout << "Sampling " << num_samples << " hypergraphs from configuration model with k=" << k << " n=" << n << " gamma=" << gamma << " min_k=" << min_k << " max_k=" << max_k << " and width_limit_exp=" << width_limit_exponent << std::endl;

    // Construct output filename
    std::string filepath = args.filepath;
    std::string filename;
    if (!config_model)
        filename = filepath + construct_uniform_filename(n, m, k, num_samples, min_k, max_k);
    else
        filename = filepath + construct_config_filename(n, gamma, k, num_samples, min_k, max_k);

    if (!filter_isomorphic)
        filename += "_no-isomorphic-filter";

    if (width_limit_exponent > 0)
        filename += "_width-limit-exp-" + std::to_string(width_limit_exponent);

    filename += ".csv";

    std::vector<Params> loop_args; 
    // Construct arguments vector using Params helper struct
    for (int i = start_sample; i < num_samples; i++) {
        Params p(config_model, i, n, m, gamma, k, min_k, max_k, filter_isomorphic, sequential, width_limit_exponent);
        loop_args.push_back(p);
    }

    std::cout << "Constructed input arguments." << std::endl;

    // touch/empty the output file
    std::ofstream outfile;
    if (!append) {
        outfile.open(filename, std::ios::out | std::ios::trunc);
    } else {
        outfile.open(filename, std::ios::out | std::ios::app);
    }

    if (!outfile.is_open()) {
        std::cout << "Couldn't open file: " + filename << ". Exiting." << std::endl;
        return 0;
    } else {
        std::cout << "Output will be written to: " << filename << std::endl;
    }

    unsigned int num_rejections = 0;
    unsigned int max_rejections = num_samples*10;
    if (!sequential_samples) {
        tbb::spin_mutex feeder_mutex;
        tbb::parallel_for_each(loop_args.begin(), loop_args.end(),
                [&](Params p, tbb::feeder<Params>& feeder)
        {
            thread_local bool success;
            success = one_sample_write(p, outfile);
            if (!success) {
                {
                    tbb::spin_mutex::scoped_lock lock(feeder_mutex);
                    num_rejections += 1;
                    if (no_max_rejections || num_rejections < max_rejections) {
                        // Add another copy of p if this sample failed due to width limit
                        feeder.add(p);
                    }
                }
            }
        });
    } else {
        int num_successes = 0;
        bool success;
        while (num_successes < num_samples) {
            for(std::size_t i = 0; i < loop_args.size(); i++) {
                    success = one_sample_write(loop_args[i], outfile);
                if (success) {
                    num_successes += 1;
                    if (num_successes >= num_samples)
                        break;
                }
                else {
                    num_rejections += 1;
                    if (!no_max_rejections && num_rejections >= max_rejections)
                        break;
                }
            }
            if (!no_max_rejections && num_rejections >= max_rejections)
                break;
        }
    }
   
    if (width_limit_exponent > 0)
        std::cout << "Number of rejected samples: " << num_rejections << std::endl;

    outfile.close();
}

