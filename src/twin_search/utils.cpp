#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <float.h>
#include "utils.hpp"
#include <boost/math/special_functions/beta.hpp>

namespace ublas=boost::numeric::ublas;
namespace bmat=boost::math;

// Convenience function to sum a ublas matrix
// TODO: Research optimal way to iterate over boost matrices
int matsum(const ublas::matrix<int> &m) {
    int sum = 0;
    for (std::size_t i = 0; i < m.size1(); i++) {
        for (std::size_t j = 0; j < m.size2(); j++) {
            sum += m(i, j);
        }
    }
    return sum;
}

// Convenience function to sum a ublas matrix
// TODO: Research optimal way to iterate over boost matrices
int matmax(const ublas::matrix<int> &m) {
    int max = 0;
    for (std::size_t i = 0; i < m.size1(); i++) {
        for (std::size_t j = 0; j < m.size2(); j++) {
            if ( m(i,j) > max )
                max = m(i, j);
        }
    }
    return max;
}

// TODO: There is probably a better way to do this with std::sum and lambdas
int mapsum(std::map<int, int> &m) {
    int sum = 0;
    for(const auto& [u, val] : m)
        sum += val;
    return sum;
}

unsigned int factorial(unsigned int n) {
    if (n <= 1)
        return 1;
    else
        return n * factorial(n-1);
}

unsigned int binom_exact(unsigned int n, unsigned int k) {
    return factorial(n) / (factorial(k) * factorial(n-k));
}

// NOTE: I round the output here because it gives correct integer values after
// casting if needed, otherwise casting tends to makve values off by 1.
double binom(double n, double k) {
    return std::round(1 / ((n + 1) * bmat::beta(n - k + 1, k + 1)) );
}

std::vector<int> sorted_keys_by_value(std::map<int, int> &m) {
    std::vector<std::pair<int, int>> pairs;
    for (auto itr = m.begin(); itr != m.end(); ++itr)
        pairs.push_back(*itr);

    sort(pairs.begin(), pairs.end(), [=](std::pair<int, int>& a, std::pair<int, int>& b)
    {
        // TODO: Reverse this to get reversed list?
        return a.second < b.second;
    }
    );

    std::vector<int> sorted_keys;
    for(std::pair<int,int> p : pairs)
        sorted_keys.push_back(p.second);

    return sorted_keys;
}
