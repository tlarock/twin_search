#include "utils.hpp"
#include <gtest/gtest.h>

TEST(BinomTest, Correctness) {
    unsigned int exact;
    double beta_result;
    unsigned int beta_result_int;
    unsigned int max_n = 13;
    for(unsigned int n = 2; n < max_n; n++) {
        for(unsigned int k = 1; k < n+1; k++) {
            exact = binom_exact(n, k);
            beta_result = binom(n, k);
            beta_result_int = ( unsigned int )(beta_result);
            EXPECT_EQ(exact, beta_result_int);
        }
    }
}
