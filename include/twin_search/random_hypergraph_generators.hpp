#ifndef RANDOM_HYPERGRAPH_GENERATORS_H
#define RANDOM_HYPERGRAPH_GENERATORS_H

#include <iostream>
#include <vector>
#include <set>
#include <random>
#include "hypergraph.hpp"

// Each generator comes in two forms. The one taking an explicit std::mt19937&
// draws from that generator and nothing else, which is what makes a run
// reproducible: the caller decides the stream, so the result does not depend on
// which worker thread happens to run the sample. The form without a generator
// keeps the historical behaviour - a thread_local mt19937 seeded from the clock
// and the thread id - and is therefore NOT reproducible. Prefer the explicit
// form in anything whose output you intend to compare or publish.
Hypergraph sample_uniform_random(int, int, int, std::mt19937 &);
Hypergraph sample_uniform_random(int, int, int);
std::vector<int> sample_hyperedge(int, std::mt19937 &, std::uniform_int_distribution<int> &);
Hypergraph uniform_hypergraph_configuration_model(int n, double gamma, int k, int max_degree, std::mt19937 &);
Hypergraph uniform_hypergraph_configuration_model(int n, double gamma, int k, int max_degree);
Hypergraph uniform_hypergraph_configuration_model(int n, double gamma, int k);
Hypergraph uniform_hypergraph_configuration_model(std::map<int, int> node_degrees, int k, std::mt19937 &);
Hypergraph uniform_hypergraph_configuration_model(std::map<int, int> node_degrees, int k);
Hypergraph chung_lu_hypergraph(std::map<int, int> k1, std::map<int, int> k2, std::mt19937 &);
Hypergraph chung_lu_hypergraph(std::map<int, int> k1, std::map<int, int> k2);
std::map<int, int> get_powerlaw_degrees(int n, double gamma, int max_k, std::mt19937 &);
std::map<int, int> get_powerlaw_degrees(int n, double gamma, int max_k);

// The thread_local generator behind the non-reproducible overloads.
std::mt19937 &default_generator();
double sample_power_law(double, int, std::mt19937 &, std::uniform_real_distribution<double> &);
#endif
