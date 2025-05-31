#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/adjacency_list.hpp>

// wrapper type for an UndirectedGraph
// NOTE: Using boost::undirectedS with the boost::vf2_graph_iso does not
// properly handle self-loops, so make sure not to include them!
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> UndirectedGraph;

// convenience type for returning a map from k->set of cliques of size k
typedef std::map<int, std::set<std::vector<int> > > CliqueMap;

unsigned int factorial(unsigned int n);
unsigned int binom_exact(unsigned int, unsigned int);

// binomial coefficient based on boost::math::beta
double binom(double, double);

// sum of a 2d boost::ublas matrix<int>
int matsum(const boost::numeric::ublas::matrix<int> &m);

// maximum element of a 2d boost::ublas matrix<int>
int matmax(const boost::numeric::ublas::matrix<int> &m);

// Sum of values in a map of ints
int mapsum(std::map<int, int> &m);

// Takes a map of ints and returns vector of keys sorted by value
std::vector<int> sorted_keys_by_value(std::map<int, int> &m);

#endif
