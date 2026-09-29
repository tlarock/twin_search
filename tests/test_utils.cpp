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

// binom() computes C(n,k) via boost::math::beta, which throws std::domain_error
// unless both its arguments are strictly positive. That made binom(n,k) abort
// for k > n instead of returning 0, and TwinSearch::compute_log_width_product
// calls it for every edge - so one unsatisfiable edge took the process down.
TEST(BinomTest, DegenerateArgumentsReturnZeroRatherThanThrowing) {
    EXPECT_NO_THROW({
        for (int n = 0; n < 8; n++)
            for (int k = n + 1; k < n + 6; k++)
                EXPECT_DOUBLE_EQ(binom(n, k), 0.0) << "n=" << n << " k=" << k;
    });

    EXPECT_NO_THROW({
        EXPECT_DOUBLE_EQ(binom(5, -1), 0.0);
        EXPECT_DOUBLE_EQ(binom(5, -3), 0.0);
        EXPECT_DOUBLE_EQ(binom(-1, 2), 0.0);
    });
}

// The guard must not disturb the values the drivers actually use.
TEST(BinomTest, BoundaryValuesUnchanged) {
    for (int n = 0; n < 13; n++) {
        EXPECT_DOUBLE_EQ(binom(n, 0), 1.0) << "n=" << n;
        EXPECT_DOUBLE_EQ(binom(n, n), 1.0) << "n=" << n;
        if (n >= 1)
            EXPECT_DOUBLE_EQ(binom(n, 1), static_cast<double>(n)) << "n=" << n;
    }
    EXPECT_DOUBLE_EQ(binom(6, 3), 20.0);
    EXPECT_DOUBLE_EQ(binom(20, 10), 184756.0);
}
