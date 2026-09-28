#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <numeric>
#include <discreture.hpp>
#include "combinations.hpp"

// combinations_of() replaces discreture::combinations() on the TwinSearch
// search hot path, because discreture::binomial() memoizes into an
// unsynchronized function-local static and is therefore unsafe to call from
// several threads at once. These tests pin the replacement to discreture's
// behaviour - same combinations, in the same order - so the swap cannot change
// any search result.

// Reference: the exact sequence discreture produces for (items, k), collected
// the same way TwinSearch::get_combinations used to collect it.
static std::vector<std::vector<int> > discreture_reference(const std::vector<int> &items, int k) {
    std::vector<std::vector<int> > out;
    auto combs = discreture::combinations(items, k);
    for (auto&& comb : combs) {
        std::vector<int> c;
        for (int x : comb)
            c.push_back(x);
        out.push_back(c);
    }
    return out;
}

// The equivalence that matters: identical sequences, element for element.
TEST(CombinationsTest, MatchesDiscretureExactlyOverManyNAndK) {
    for (int n = 0; n <= 14; n++) {
        std::vector<int> items(n);
        // Deliberately not 0..n-1: the generator must index into `items`
        // rather than assume the values are the indices.
        for (int i = 0; i < n; i++)
            items[i] = (i + 1) * 7;

        for (int k = 0; k <= n; k++) {
            std::vector<std::vector<int> > mine = combinations_of(items, k);
            std::vector<std::vector<int> > theirs = discreture_reference(items, k);
            ASSERT_EQ(mine.size(), theirs.size())
                << "count mismatch at n=" << n << " k=" << k;
            for (std::size_t i = 0; i < mine.size(); i++) {
                ASSERT_EQ(mine[i], theirs[i])
                    << "combination " << i << " differs at n=" << n << " k=" << k;
            }
        }
    }
}

// TwinSearch::get_combinations is called with weight = the residual matrix
// entry, which is always >= 1, but it is not otherwise bounded by the number
// of remaining clique-neighbours. k > n must yield nothing, not UB.
TEST(CombinationsTest, KGreaterThanNYieldsNothing) {
    std::vector<int> items {3, 1, 4, 1, 5};
    for (int k = static_cast<int>(items.size()) + 1; k < 20; k++) {
        EXPECT_TRUE(combinations_of(items, k).empty()) << "k=" << k;
        EXPECT_EQ(combinations_of(items, k).size(), discreture_reference(items, k).size());
    }
}

TEST(CombinationsTest, NegativeKYieldsNothing) {
    std::vector<int> items {1, 2, 3};
    EXPECT_TRUE(combinations_of(items, -1).empty());
    EXPECT_TRUE(combinations_of(items, -7).empty());
}

// k == 0 yields exactly one empty combination, which is what discreture does
// (binomial(n, 0) == 1). TwinSearch guards against ever calling with weight 0,
// but the two must agree regardless.
TEST(CombinationsTest, ZeroKYieldsOneEmptyCombination) {
    std::vector<int> items {1, 2, 3};
    std::vector<std::vector<int> > mine = combinations_of(items, 0);
    ASSERT_EQ(mine.size(), 1u);
    EXPECT_TRUE(mine[0].empty());
    EXPECT_EQ(mine, discreture_reference(items, 0));

    std::vector<int> empty_items;
    EXPECT_EQ(combinations_of(empty_items, 0).size(), 1u);
    EXPECT_TRUE(combinations_of(empty_items, 1).empty());
}

// Ordering is COLEXICOGRAPHIC in the *indices* into the input (largest index
// varies slowest), not lexicographic, and not ordered by value. discreture's
// default Combinations class is colex; its lexicographic generator is the
// separate LexCombinations class. Getting this wrong would still enumerate the
// same SET of twins, but would change which twin TwinSearch::filtered_twins
// keeps as each isomorphism class representative.
TEST(CombinationsTest, ColexicographicByIndexNotByValue) {
    std::vector<int> descending {50, 40, 30, 20};
    std::vector<std::vector<int> > mine = combinations_of(descending, 2);
    // index pairs (0,1) (0,2) (1,2) (0,3) (1,3) (2,3)
    std::vector<std::vector<int> > expected {
        {50, 40}, {50, 30}, {40, 30}, {50, 20}, {40, 20}, {30, 20}
    };
    EXPECT_EQ(mine, expected);
    EXPECT_EQ(mine, discreture_reference(descending, 2));
}

// Spell the colex index sequence out once at k=3 as well, so the property is
// pinned independently of discreture being available.
TEST(CombinationsTest, ColexIndexSequenceAtK3) {
    std::vector<int> items {0, 1, 2, 3, 4};
    std::vector<std::vector<int> > expected {
        {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3},
        {0,1,4}, {0,2,4}, {1,2,4}, {0,3,4}, {1,3,4}, {2,3,4}
    };
    EXPECT_EQ(combinations_of(items, 3), expected);
    EXPECT_EQ(combinations_of(items, 3), discreture_reference(items, 3));
}

TEST(CombinationsTest, CountsMatchBinomial) {
    for (int n = 0; n <= 14; n++) {
        std::vector<int> items(n);
        std::iota(items.begin(), items.end(), 0);
        for (int k = 0; k <= n; k++) {
            // exact expected count via Pascal's rule, no floating point
            std::vector<std::vector<long long> > C(n + 1, std::vector<long long>(n + 1, 0));
            for (int a = 0; a <= n; a++) {
                C[a][0] = 1;
                for (int b = 1; b <= a; b++)
                    C[a][b] = C[a-1][b-1] + (b <= a-1 ? C[a-1][b] : 0);
            }
            EXPECT_EQ(static_cast<long long>(combinations_of(items, k).size()), C[n][k])
                << "n=" << n << " k=" << k;
        }
    }
}

// num_combinations() only sizes a reserve(), so it is allowed to saturate,
// but it must be exact below the cap and must never over-report.
TEST(CombinationsTest, NumCombinationsExactBelowCapAndSaturates) {
    for (int n = 0; n <= 14; n++) {
        std::vector<int> items(n);
        std::iota(items.begin(), items.end(), 0);
        for (int k = 0; k <= n; k++) {
            std::size_t actual = combinations_of(items, k).size();
            std::size_t reserved = num_combinations(n, k);
            if (actual < 4096)
                EXPECT_EQ(reserved, actual) << "n=" << n << " k=" << k;
            else
                EXPECT_EQ(reserved, 4096u) << "n=" << n << " k=" << k;
        }
    }
    EXPECT_EQ(num_combinations(5, -1), 0u);
    EXPECT_EQ(num_combinations(5, 6), 0u);
    EXPECT_EQ(num_combinations(100, 50), 4096u);  // saturates, does not overflow
}

// The whole point of the replacement: no shared mutable state. Hammer it from
// many threads and require identical results to the single-threaded answer.
// (Under ThreadSanitizer this also asserts the absence of a data race.)
#include <thread>
#include <atomic>
TEST(CombinationsTest, ThreadSafeUnderConcurrentUse) {
    std::vector<std::vector<int> > expected;
    {
        std::vector<int> items(30);
        std::iota(items.begin(), items.end(), 0);
        expected = combinations_of(items, 2);
    }

    std::atomic<int> mismatches{0};
    std::vector<std::thread> threads;
    for (int t = 0; t < 8; t++) {
        threads.emplace_back([&, t]() {
            for (int rep = 0; rep < 200; rep++) {
                // vary n per thread, which is what drives discreture's memo
                // table to resize concurrently
                int n = 7 + ((t * 13 + rep * 7) % 30);
                std::vector<int> items(n);
                std::iota(items.begin(), items.end(), 0);
                std::vector<std::vector<int> > got = combinations_of(items, 2);
                std::size_t want = static_cast<std::size_t>(n) * (n - 1) / 2;
                if (got.size() != want)
                    mismatches++;
            }
            std::vector<int> items(30);
            std::iota(items.begin(), items.end(), 0);
            if (combinations_of(items, 2) != expected)
                mismatches++;
        });
    }
    for (auto &th : threads)
        th.join();
    EXPECT_EQ(mismatches.load(), 0);
}
