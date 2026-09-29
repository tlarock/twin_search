#include "hypergraph.hpp"
#include "random_hypergraph_generators.hpp"
#include "utils.hpp"
#include <gtest/gtest.h>

TEST(UniformHypergraphTest, HyperedgeSizeDistribution) {
    // Test pairwise graph
    int n = 10;
    int m = 10;

    for(int k = 3; k < 6; k++) {
        for (int trial = 0; trial < 10; trial++) {
            Hypergraph h = sample_uniform_random(n, m, k);

            for(const auto& [he_idx, he] : h.hyperedges) {
                EXPECT_EQ(he.size(), k);
                EXPECT_EQ(h.hyperedge_sizes[he_idx], k);
            }
        }
    }
}

TEST(UniformConfigurationHypergraphTest, KRegular) {
    // Test 2-uniform 2-regular graph
    int n = 50;
    int k = 2;
    std::map<int, int> node_degree_map;
    for(int node = 0; node < n; node++)
        node_degree_map[node] = k;
    int m = mapsum(node_degree_map) / k;
    Hypergraph h = uniform_hypergraph_configuration_model(node_degree_map, k);
    EXPECT_EQ(h.n, n);
    EXPECT_EQ(h.m, m);

    for(const auto& [he_idx, he] : h.hyperedges) {
        EXPECT_EQ(he.size(), k);
        EXPECT_EQ(h.hyperedge_sizes[he_idx], k);
    }

    // Test 3-uniform 3-regular
    k = 3;
    for(int node = 0; node < n; node++)
        node_degree_map[node] = k;
    m = mapsum(node_degree_map) / k;
    h = uniform_hypergraph_configuration_model(node_degree_map, k);
    EXPECT_EQ(h.n, n);
    EXPECT_EQ(h.m, m);

    for(const auto& [he_idx, he] : h.hyperedges) {
        EXPECT_EQ(he.size(), k);
        EXPECT_EQ(h.hyperedge_sizes[he_idx], k);
    }

    // Test 4-uniform 4-regular
    k = 4;
    for(int node = 0; node < n; node++)
        node_degree_map[node] = k;

    m = mapsum(node_degree_map) / k;
    h = uniform_hypergraph_configuration_model(node_degree_map, k);
    EXPECT_EQ(h.n, n);
    EXPECT_EQ(h.m, m);

    for(const auto& [he_idx, he] : h.hyperedges) {
        EXPECT_EQ(he.size(), k);
        EXPECT_EQ(h.hyperedge_sizes[he_idx], k);
    }
 
}

TEST(UniformConfigurationHypergraphTest, PowerLawDegrees) {
    // Test 3-uniform powerlaw distributed
    int total_degree;
    int remainder;
    int m;
    std::map<int, int> node_degree_map; 
    int n = 100;
    int max_deg = n-1;
    std::vector<double> gammas = {2.5, 3.0, 3.5};
    for(double gamma : gammas) {
        for(int k = 3; k < 6; k++) {
            for(int trial = 0; trial < 10; trial++) {
                node_degree_map = get_powerlaw_degrees(n, gamma, max_deg);
                total_degree = mapsum(node_degree_map);
                remainder = total_degree % k;
                if (remainder != 0) {
                    total_degree += (k-remainder);
                }

                m = total_degree / k;
                Hypergraph h = uniform_hypergraph_configuration_model(node_degree_map, k);
                EXPECT_EQ(h.n, n);
                EXPECT_EQ(h.m, m);

                for(const auto& [he_idx, he] : h.hyperedges) {
                    EXPECT_EQ(he.size(), k);
                    EXPECT_EQ(h.hyperedge_sizes[he_idx], k);
                }

                for(const auto& [node, membs] : h.node_memberships) {
                    EXPECT_GE(membs.size(), 1);
                    // TODO FIXME: It is useful to have some sort of test here to make
                    // sure that the degrees make sense beyond being non-zero.
                    // However, this particular test is very ad-hoc, could lead to
                    // non-deterministic failures. It is passing for the time
                    // being so leaving it here.
                    EXPECT_GE(membs.size(), std::max(1, node_degree_map[node]-k));
                }
            }
        }
    }
}

TEST(UniformConfigurationHypergraphTest, NoPrecompute) {
    // Test 3-uniform powerlaw distributed
    int n = 100;
    std::vector<double> gammas = {2.5, 3.0, 3.5};
    for(double gamma : gammas) {
        for(int k = 3; k < 6; k++) {
            for(int trial = 0; trial < 10; trial++) {
                Hypergraph h = uniform_hypergraph_configuration_model(n, gamma, k);
                EXPECT_EQ(h.n, n);

                for(const auto& [he_idx, he] : h.hyperedges) {
                    EXPECT_EQ(he.size(), k);
                    EXPECT_EQ(h.hyperedge_sizes[he_idx], k);
                }

                for(const auto& [node, membs] : h.node_memberships) {
                    EXPECT_GE(membs.size(), 1);
                }
            }
        }
    }
}

TEST(ChungLuHypergraphTest, NumNodesAndHyperedges) {
    // Test pairwise graph
    int n = 50;
    int m = 50;
    int k = 2;
    std::map<int, int> node_degree_map;
    for(int node = 0; node < n; node++)
        node_degree_map[node] = k;

    std::map<int, int> hyperedge_size_map;
    for(int edge_id = 0; edge_id < m; edge_id++)
        hyperedge_size_map[edge_id] = k;

    Hypergraph h = chung_lu_hypergraph(node_degree_map, hyperedge_size_map);
    EXPECT_EQ(h.n, n);
    EXPECT_EQ(h.m, m);

    for(const auto& [he_idx, he] : h.hyperedges) {
        if (he.size() == 0) {
            EXPECT_EQ(h.hyperedge_sizes[he_idx], 0);
        }
    }
}

TEST(ChungLuHypergraphTest, FromDistribution) {
    // Test pairwise graph
    int n = 50;
    int max_deg = n-1;
    int k = 3;
    double gamma = 3.0;
    std::map<int, int> node_degrees = get_powerlaw_degrees(n, gamma, max_deg);
    std::map<int, int> hyperedge_sizes;
    int total_degree = mapsum(node_degrees);
    int m = total_degree / k;
    for(int eid = 0; eid < m; eid++)
      hyperedge_sizes[eid] = k;

    Hypergraph h = chung_lu_hypergraph(node_degrees, hyperedge_sizes);
    EXPECT_EQ(h.n, n);
    for(const auto& [he_idx, he] : h.hyperedges) {
        if (he.size() == 0) {
            EXPECT_EQ(h.hyperedge_sizes[he_idx], 0);
        }
    }
}

TEST(PowerLawSamplingTests, BasicProperties) {
    int n = 100;
    int max_k = 100;
    double gamma = 3.0;
    std::map<int, int> degrees = get_powerlaw_degrees(n, gamma, max_k);
    for(const auto& [node, deg] : degrees) {
        EXPECT_GE(deg, 1);
        EXPECT_LE(deg, max_k);
    }
}

// Explicit-generator overloads exist so a run can be reproduced. The property
// that matters is that the sampled hypergraph is a function of the generator
// state alone - not of the calling thread, the thread count, or the order in
// which samples happen to be scheduled. count_twins_random relies on this to
// make --seed meaningful under tbb::parallel_for_each.
TEST(SeedingTest, SameSeedGivesTheSameHypergraph) {
    const int n = 8, m = 8, k = 3;
    for (unsigned int seed : {1u, 42u, 7919u}) {
        std::mt19937 a(seed), b(seed);
        Hypergraph ha = sample_uniform_random(n, m, k, a);
        Hypergraph hb = sample_uniform_random(n, m, k, b);
        EXPECT_EQ(ha.hyperedges, hb.hyperedges) << "seed " << seed;
        EXPECT_EQ(ha.n, hb.n);
        EXPECT_EQ(ha.m, hb.m);
    }
}

TEST(SeedingTest, DifferentSeedsGiveDifferentHypergraphs) {
    const int n = 8, m = 8, k = 3;
    std::mt19937 a(1u), b(2u);
    Hypergraph ha = sample_uniform_random(n, m, k, a);
    Hypergraph hb = sample_uniform_random(n, m, k, b);
    // Not guaranteed in principle, but the space here is astronomically large.
    EXPECT_NE(ha.hyperedges, hb.hyperedges);
}

// A generator passed by reference must be ADVANCED, so that consecutive draws
// from one stream differ. If it were taken by value, every sample in a loop
// would come out identical.
TEST(SeedingTest, GeneratorIsAdvancedBetweenDraws) {
    std::mt19937 gen(12345u);
    Hypergraph first = sample_uniform_random(8, 8, 3, gen);
    Hypergraph second = sample_uniform_random(8, 8, 3, gen);
    EXPECT_NE(first.hyperedges, second.hyperedges);
}

// The configuration model draws twice - the degree sequence and the stub
// matching - and both must come from the caller's generator, or seeding only
// half determines the result.
TEST(SeedingTest, ConfigurationModelIsFullySeeded) {
    std::mt19937 a(2024u), b(2024u);
    Hypergraph ha = uniform_hypergraph_configuration_model(10, 3.0, 3, 5, a);
    Hypergraph hb = uniform_hypergraph_configuration_model(10, 3.0, 3, 5, b);
    EXPECT_EQ(ha.hyperedges, hb.hyperedges);

    std::mt19937 c(2025u);
    Hypergraph hc = uniform_hypergraph_configuration_model(10, 3.0, 3, 5, c);
    EXPECT_NE(ha.hyperedges, hc.hyperedges);
}

// The historical overloads must keep working for callers that do not care.
TEST(SeedingTest, UnseededOverloadsStillProduceValidHypergraphs) {
    Hypergraph h = sample_uniform_random(8, 8, 3);
    EXPECT_EQ(h.m, 8);
    for (const auto &[idx, he] : h.hyperedges)
        EXPECT_EQ(he.size(), 3u);
}
