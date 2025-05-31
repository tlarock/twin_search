#ifndef RANDOM_HYPERGRAPH_GENERATORS_H
#define RANDOM_HYPERGRAPH_GENERATORS_H

#include <iostream>
#include <vector>
#include <set>
#include <random>
#include "hypergraph.hpp"

Hypergraph sample_uniform_random(int, int, int);
std::vector<int> sample_hyperedge(int, std::mt19937 &, std::uniform_int_distribution<int> &);
Hypergraph uniform_hypergraph_configuration_model(int n, double gamma, int k, int max_degree);
Hypergraph uniform_hypergraph_configuration_model(int n, double gamma, int k);
Hypergraph uniform_hypergraph_configuration_model(std::map<int, int> node_degrees, int k);
Hypergraph chung_lu_hypergraph(std::map<int, int> k1, std::map<int, int> k2);
std::map<int, int> get_powerlaw_degrees(int n, double gamma, int max_k);
double sample_power_law(double, int, std::mt19937 &, std::uniform_real_distribution<double> &);
#endif
