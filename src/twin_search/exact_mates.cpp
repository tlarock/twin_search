#include <iostream>
#include <fstream>
#include <oneapi/tbb.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/global_control.h>
#include <boost/algorithm/string.hpp>

#include "argparse/argparse.hpp"
#include "hypergraph.hpp"
#include "factor_graph.hpp"
#include "projected_graph.hpp"
#include "twin_search.hpp"
#include "utils.hpp"

// struct to read and store command line arguments with argparse
struct MyArgs : public argparse::Args {
    std::string &filepath = kwarg("output-path", "Output filepath. Default is results/increasing_density/").set_default("results/increasing_density/");
    int &max_threads = kwarg("max-threads, tbb-max-allowed-parallelism", "Value passed to TBB to limit number of worker threads. If <= 0, ignored.").set_default(0);
};


struct TwinSearchChecker {
    //Hypergraph H;
    //ProjectedGraph proj;
    ublas::matrix<int> proj_mat;
    ublas::matrix<int> lg_mat;
    //UndirectedGraph lg_graph;
};

bool mateq(ublas::matrix<int> A, ublas::matrix<int> B) {
    // If the dimensions do not match, they are not equal
    if (A.size1() != B.size1() || A.size2() != B.size2())
        return false;

    for (std::size_t r = 0; r < A.size1(); ++r) {
        for (std::size_t c = 0; c < A.size2(); ++c) {
            if (A(r,c) != B(r,c))
                return false;
        }
    }
    return true;
}

std::vector<std::vector<int> > run_exact_mates_tests_parallel(std::vector<TwinSearchChecker> &gm_vect) {
    tbb::concurrent_vector<std::vector<int> > mate_pairs;
    if (gm_vect.size() < 2)
        return std::vector<std::vector<int> >(0);

    for(std::size_t i = 0; i < gm_vect.size()-1; i++) {
        tbb::parallel_for(std::size_t(i+1), gm_vect.size(), [&](std::size_t j){
                // check for projection and line graph equality
                if (mateq(gm_vect[i].proj_mat, gm_vect[j].proj_mat) && 
                        mateq(gm_vect[i].lg_mat, gm_vect[j].lg_mat)) {
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

// Takes a line from the exhaustive projections output and parses all of the
// hypergraphs into a vector of TwinSearchChecker objects
std::vector<TwinSearchChecker> parse_hypergraphs_line(std::string line) {
    // output structure
    std::vector<TwinSearchChecker> gm_vect;

    // Get vector of string representing hypergraphs
    std::vector<std::string> line_str;
    // split the line on /, then the
    // hypergraphs string will be in the last position
    boost::split(line_str, line, boost::is_any_of("/"));

    // Parse each hypergraph into a vector of strings
    std::vector<std::string> hypergraphs_str;
    // Split on semicolon to get each hypergraph
    boost::split(hypergraphs_str, line_str.back(), boost::is_any_of(";"));

    // for each hypergraph
    for (std::size_t hg_idx = 0; hg_idx < hypergraphs_str.size(); ++hg_idx) {
        // Parse the hyperedges into a vector of strings
        std::vector<std::string> hyperedges_strs;
        boost::split(hyperedges_strs, hypergraphs_str[hg_idx], boost::is_any_of("|"));

        // Construct hypergraph as vector of vector of ints
        std::vector<std::vector<int> > hyperedges;
        for (std::size_t he_idx = 0; he_idx < hyperedges_strs.size(); ++he_idx) {
            std::vector<std::string> he_str;
            boost::split(he_str, hyperedges_strs[he_idx], boost::is_any_of(":"));
            std::vector<int> he(he_str.size());
            // parse hypergraphs
            for (size_t i = 0; i < he_str.size(); i++) {
                he[i] = std::stoi( he_str[i] );
            }
            // Add hyperedge to hyperedges
            hyperedges.push_back(std::vector<int>(he));
        }

        // Construct the hypergraph object, from which we can compute
        // projections and line graphs
        Hypergraph H(hyperedges);
        ProjectedGraph proj(H, false);
        // Create the actual TwinSearchChecker object
        struct TwinSearchChecker gmc {.proj_mat = proj.proj_mat, .lg_mat = H.get_lg_mat()};

        // Push to the list
        gm_vect.push_back(gmc);
    }

    return gm_vect;
}

void write_exact_mates(std::vector<std::vector<int> > &mates_pairs, std::ofstream &outfile) {
    for (std::size_t i = 0; i < mates_pairs.size(); ++i)
    {
        outfile << mates_pairs[i][0] << "," << mates_pairs[i][1];
        if (i < mates_pairs.size()-1)
             outfile << ";";
        else
            outfile << std::endl;
    }
}

int main(int argc, char *argv[]) {
	// Set cout to be unbuffered
	std::cout.setf(std::ios::unitbuf);
    // Read input file name
    auto args = argparse::parse<MyArgs>(argc, argv);
    const std::string filepath = args.filepath;
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


    std::ifstream infile;
    infile.open(filepath, std::ios::in);
    std::ofstream outfile;
    outfile.open(filepath + ".exact_mates_pairs", std::ios::out | std::ios::trunc);
    std::string line;
    std::vector<TwinSearchChecker> gm_vect;
    unsigned int total_pairs = 0;
    // For each projection (line-by-line file read)
    while (std::getline(infile, line)) {
        // Get TwinSearchChecker objects corresponding to hypergraphs in this
        // projection
        gm_vect = parse_hypergraphs_line(line);
        // compute isomorphic and exact mates
        std::vector<std::vector<int> > mates_pairs = run_exact_mates_tests_parallel(gm_vect);
        total_pairs += mates_pairs.size();
        // Write to .mates_pairs file
        write_exact_mates(mates_pairs, outfile);
    }
    std::cout << "Total pairs of mates: " << total_pairs << std::endl;
    infile.close();
    outfile.close();
}
