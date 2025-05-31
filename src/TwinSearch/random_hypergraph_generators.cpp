#include "random_hypergraph_generators.hpp"
#include <cmath>
#include <algorithm>
#include <random>
#include <map>
#include <ctime>
#include <thread>
#include "utils.hpp"

// Generates a k-uniform hypergraph by sampling m hyperedges of size k
// uniformly at random. Guarantees that every node has non-zero degree.
Hypergraph sample_uniform_random(int n, int m, int k)
{
    int maximum_non_singletons = m*k;
    if (maximum_non_singletons < n) {
        return Hypergraph(std::vector<std::vector<int> > (0));
    }
    // Validate input: there must be at least m possible hyperedges
    double max_edges = binom(binom(n, k), m);
    if (m > max_edges)
    {
        return Hypergraph(std::vector<std::vector<int> > (0));
    }

    // Construct a thread_local RNG
    static thread_local std::mt19937 generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));

    std::uniform_int_distribution<int> distribution(0, n-1);
    std::set<int> nodes;
    std::set<std::vector<int> > hyperedges;
    std::vector<int> hyperedge;
    std::size_t nn = n;
    std::size_t mm = m;
    while (nodes.size() < nn) {
        nodes.clear();
        hyperedges.clear();
        while (hyperedges.size() < mm) {
            hyperedge = sample_hyperedge(k, generator, distribution);
            hyperedges.insert(hyperedge);
            for (int u : hyperedge) {
                nodes.insert(u);
            }
        }
    }

    return Hypergraph(hyperedges);
}

// Sample an individual hyperedge of size k from n nodes
// TODO: May be more efficient to add an output parameter rather than returning
// a new std::vector every time, though I will have to copy it after anyway.
std::vector<int> sample_hyperedge(int k, std::mt19937 &gen, std::uniform_int_distribution<int> &dist) {
    std::size_t kk = k;
    std::set<int> he_set;
    while (he_set.size() < kk) {
        he_set.insert(dist(gen));
    }
    std::vector<int> he_vect;
    for (int i : he_set) {
        he_vect.push_back(i);
    }
    sort(he_vect.begin(), he_vect.end());

    return he_vect;
}

Hypergraph uniform_hypergraph_configuration_model(int n, double gamma, int k, int max_degree) {
    std::map<int, int> node_degree_map = get_powerlaw_degrees(n, gamma, max_degree);
    return uniform_hypergraph_configuration_model(node_degree_map, k);
}

Hypergraph uniform_hypergraph_configuration_model(int n, double gamma, int k) {
    std::map<int, int> node_degree_map = get_powerlaw_degrees(n, gamma, n-1);
    return uniform_hypergraph_configuration_model(node_degree_map, k);
}

// k-uniform configuration model with degree distribution node_degrees
//
// NOTE: While this function is intended to implement the microcanonical model,
// the distribution in node_degrees is not gauranteed to be returned
// exactly. Firt, it may be modified if it is not hypergraphical, meaning
// the total degree is not divisible by k. Second, it may be modified
// if we reach the last hyperedge and have a repeated node. These modifications
// should be negligible for sparse hypergraphs.
Hypergraph uniform_hypergraph_configuration_model(std::map<int, int> node_degrees, int k) {
    // Construct a thread_local RNG
    static thread_local std::mt19937 gen(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));

    std::size_t kk = k;
    // Ensure that node_degrees is hypergraphical for a k-uniform hypergraph by
    // checking if the sum of the degrees is divisible by k
    int n = node_degrees.size();
    int total_degree = mapsum(node_degrees);
    int remainder = total_degree % k;
    if (remainder !=  0) {
        // If the sum is not divisible by k, we need to change some degrees to
        // make the distribution hypergraphical. For now, I just add 1 to the
        // first k-remainder nodes.
        int num_to_add = k - remainder;
        for(auto & [u, val] : node_degrees) {
            val += 1;
            num_to_add -= 1;
            if (num_to_add == 0)
                break;
        }
    }

    // Number of hyperedges is the total degree divided by the hyperedge size,
    // which should now be an integer
    total_degree = mapsum(node_degrees);
    int m  = total_degree / k;

    // Construct vector of stubs, containing node_degree[u] copies of u
    std::vector<int> stubs;
    for(const auto & [u, val] : node_degrees) {
        for(int i = 0; i < val; i++) {
            stubs.push_back(u);
        }
    }

    // Initialize int distribution between 0 and size of stubs
    std::uniform_int_distribution<int> dist(0, stubs.size()-1);

    // The hypergraph and hyperedge containers
    std::set<std::vector<int> > hypergraph;
    std::set<int> he_set;
    std::vector<int> he_vect(k);

    int uidx;
    while(!stubs.empty()) {
        // Make sure int distribution has max value of stubs
        dist = std::uniform_int_distribution<int>(0, stubs.size()-1);

        // Construct a hypergraph by sampling stubs, ensuring no repeated nodes
        he_set = std::set<int>();
        if(stubs.size() > kk) {
            while(he_set.size() < kk) {
                uidx = dist(gen);
                he_set.insert(stubs[uidx]);
            }
        } else if(stubs.size() == kk) {
            for(int u : stubs)
                he_set.insert(u);
        }
            
        if (he_set.size() == kk) { 
            // Turn the set hyperedge into a sorted vector hyperedge
            he_vect = std::vector<int>(0);
            std::copy(he_set.begin(), he_set.end(), std::back_inserter(he_vect));
            sort(he_vect.begin(), he_vect.end());

            // Check if this hyperedge is already in the hypergraph
            if (!hypergraph.contains(he_vect)) {
                // Add to the hypergraph
                hypergraph.insert(he_vect);

                // Remove one copy of each node from stubs
                for(std::size_t i = 0; i < he_vect.size(); i++) {
                    auto it = std::find(stubs.begin(), stubs.end(), he_vect[i]);
                    if (it != stubs.end()) {
                        stubs.erase(it);
                    }
                }
            }
        }

        // Deal with two potential problems adding the final hyperedge
        if (stubs.size() == kk) {
            // First: repeated node in final k stubs. We will replace it,
            // hopefully keeping our degree distribution as true to
            // node_degrees as possible.
            dist = std::uniform_int_distribution<int>(0, n-1);
            int new_stub;
            int num_stubs;
            for(std::size_t i = 0; i < stubs.size(); i++) {
                num_stubs = std::count(stubs.begin(), stubs.end(), stubs[i]);
                if (num_stubs > 1) {
                    // Replace item i
                    new_stub = dist(gen);
                    while(find(stubs.begin(), stubs.end(), new_stub) != stubs.end())
                        new_stub = dist(gen);
                    stubs[i] = new_stub;
                }
            }

            // Second: check if the hyperedge of the final k stubs is already
            // in the hypergraph. In this unlucky case, we will just resample
            // a totally new hyperedge to finish the process.
            he_vect = std::vector<int>(0);
            std::copy(stubs.begin(), stubs.end(), std::back_inserter(he_vect));
            sort(he_vect.begin(), he_vect.end());
            while(hypergraph.contains(he_vect)) {
                // need to sample a new set of stubs
                dist = std::uniform_int_distribution<int>(0, n-1);
                he_set = std::set<int>();
                while(he_set.size() < kk) {
                    he_set.insert(dist(gen));
                }
                he_vect = std::vector<int>(0);
                std::copy(he_set.begin(), he_set.end(), std::back_inserter(he_vect));
                sort(he_vect.begin(), he_vect.end());
            }

            // Insert the final hyperedge
            hypergraph.insert(he_vect);
            
            // Empty stubs
            stubs.clear();
        }
    }

    if (hypergraph.size() != static_cast<std::size_t> (m)) {
        m = static_cast<int>(hypergraph.size());
    }
    return Hypergraph(hypergraph, n, m);
}

// Constructs a non-uniform chung_lu_hypergraph with expected degree
// distribution node_degrees and expected hyperedge size distribution
// hyperedge_sizes.
// NOTE: This model offers no gaurantees about the sizes of hyperedges,
// and in fact hyperedges can be as large as m.
Hypergraph chung_lu_hypergraph(std::map<int, int> node_degrees, std::map<int, int> hyperedge_sizes) {
    // Construct a thread_local RNG
    static thread_local std::mt19937 generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_real_distribution<double> distribution(0, 1);

    // The loops go in decreasing order of degree and hyperedge size 
    std::vector<int> sorted_node_ids = sorted_keys_by_value(node_degrees);
    std::vector<int> sorted_edge_ids = sorted_keys_by_value(hyperedge_sizes);

    int n = node_degrees.size();
    int m = hyperedge_sizes.size();
    double total_degree = static_cast <double> (mapsum(node_degrees));

    // Need a data structure to store the hyperedges before turning them into
    // a hypergraph object. Would be better to pre-allocate if possible, but
    // need to ensure there are no duplicates as well
    std::map<int, std::set<int> > edge_map;
    std::map<std::set<int>, int> rev_edge_map;
    int edge; 
    int edge_id;
    double p;
    double log1p;
    double q;
    double r;
    // For each node u
    for(int u : sorted_node_ids) {
        // Starting from the first edge id
        edge_id = 0;
        edge = sorted_edge_ids[edge_id];
        // Compute the probability of adding node u to edge
        p = (node_degrees[u] * hyperedge_sizes[edge]) / total_degree;
        while(edge_id < m) {
            if(p < 1) {
                r = distribution(generator);
                log1p = std::log(1.0-p);
                // NOTE: xgi checks for log1p being 0...
                edge_id = edge_id + static_cast <int> (std::floor(std::log(r) / log1p));
            }
            // Ensure edge_id is still less than m after above modification
            if(edge_id < m) {
                edge = sorted_edge_ids[edge_id];
                q = (node_degrees[u] * hyperedge_sizes[edge]) / total_degree;
                r = distribution(generator);
                if(r < q / p) {
                    // DELETE partial hyperedge from rev_edge_map;
                    rev_edge_map.erase(edge_map[edge]);
                    // UPDATE new hyperedge
                    edge_map[edge].insert(u);
                    rev_edge_map[edge_map[edge]] = edge;
                }
                p = q;
                edge_id += 1;
            }
        }
    }

    // Using a set here to remove any duplicate hyperedges
    std::set<std::vector<int> > hyperedges;
    std::vector<int> he;
    for(const auto & [idx, he_set] : edge_map) {
        he = std::vector<int> (0);
        for (int u : he_set)
            he.push_back(u);

        sort(he.begin(), he.end());
        hyperedges.insert(he);
    }

    return Hypergraph(hyperedges, n, m);
}

std::map<int, int> get_powerlaw_degrees(int n, double gamma, int max_k) {
    static thread_local std::mt19937 generator(clock() + std::hash<std::thread::id>()(std::this_thread::get_id()));
    std::uniform_real_distribution<double> distribution(0, 1);
    std::map<int, int> degrees;
    double sampled_val;
    for(int u = 0; u < n; u++) {
        sampled_val = std::floor(sample_power_law(gamma, max_k, generator, distribution));
        degrees[u] = static_cast <int> (sampled_val);
    }
    return degrees;
}

double sample_power_law(double gamma, int max_k, std::mt19937 &gen, std::uniform_real_distribution<double> &dist) {
    double min_k = 1.0;
    double y = dist(gen);
    return std::pow((std::pow(max_k, gamma+1.0) - std::pow(min_k, gamma+1.0)) * y + std::pow(min_k, (1.0/(gamma+1.0))), 1.0 / (gamma+1.0));
}
