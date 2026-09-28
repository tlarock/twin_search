#ifndef TWIN_SEARCH_COMBINATIONS_H
#define TWIN_SEARCH_COMBINATIONS_H
#include <vector>
#include <cstddef>

// Self-contained combination generation.
//
// This replaces discreture::combinations() on the search hot path. The
// discreture implementation is correct but not thread safe: its
// Combinations(n, k) constructor calls discreture::binomial(), which memoizes
// into an unsynchronized function-local
//
//     static std::vector<std::vector<BigIntType>> B
//
// and grows it with B.resize() (Sequences.hpp:135). C++ guarantees thread-safe
// *initialization* of a function-local static, not thread-safe *mutation*.
// TwinSearch::get_combinations() runs at every node of the search tree on every
// worker thread, so two threads can resize and read that table concurrently.
// Because resize() to a smaller size destroys inner vectors another thread is
// still reading, this shows up as heap-use-after-free under AddressSanitizer
// and as intermittent segfaults in ordinary -O3 builds.
//
// Ordering note: combinations are produced in COLEXICOGRAPHIC order of the
// *indices* into `items`, i.e. ordered by largest index first, then next
// largest, and so on: for n=4, k=2 the index pairs come out as
//     (0,1) (0,2) (1,2) (0,3) (1,3) (2,3)
// This is what discreture::combinations() produces - note that discreture's
// lexicographic generator is a different class (LexCombinations), so "the
// obvious" lex order would NOT match. Matching matters beyond tidiness:
// TwinSearch::filtered_twins keeps the FIRST representative of each isomorphism
// class encountered, so a different enumeration order would silently change
// which twin is reported as the class representative.
// tests/test_combinations.cpp asserts the equivalence against discreture
// directly.

// Number of k-subsets of an n-set, for reserving output space only. Saturates
// at RESERVE_CAP rather than overflowing, so it must not be used where an
// exact count is needed (utils.hpp::binom is the exact version).
inline std::size_t num_combinations(int n, int k)
{
    // Large enough that realistic factor-graph neighbourhoods reserve exactly,
    // small enough that a pathological (n, k) cannot ask for a huge block.
    const std::size_t RESERVE_CAP = 4096;

    if (k < 0 || k > n)
        return 0;

    // C(n, k) == C(n, n-k); taking the smaller keeps the loop short.
    if (k > n - k)
        k = n - k;

    std::size_t result = 1;
    for (int i = 1; i <= k; i++) {
        // Exact at every step: the running value is always C(n-k+i, i).
        result = result * static_cast<std::size_t>(n - k + i)
                        / static_cast<std::size_t>(i);
        if (result >= RESERVE_CAP)
            return RESERVE_CAP;
    }

    return result;
}

// Streaming form: invokes fn(const std::vector<T>&) once per k-subset of
// `items`, in colexicographic order of the underlying indices, without
// materialising the whole sequence. The vector passed to fn is reused between
// calls, so copy it if you need to keep it.
//
// Returns without calling fn at all when k < 0 or k > items.size(), matching
// discreture's binomial(n, k) == 0 for those cases. k == 0 calls fn exactly
// once with an empty combination, again matching discreture.
template <typename T, typename F>
void for_each_combination(const std::vector<T> &items, int k, F &&fn)
{
    const int n = static_cast<int>(items.size());

    if (k < 0 || k > n)
        return;

    // indices is the current combination, held as strictly increasing offsets
    // into items, starting at the colex-smallest {0, 1, ..., k-1}.
    std::vector<int> indices(static_cast<std::size_t>(k));
    for (int i = 0; i < k; i++)
        indices[i] = i;

    std::vector<T> comb(static_cast<std::size_t>(k));

    while (true) {
        for (int i = 0; i < k; i++)
            comb[i] = items[indices[i]];
        fn(const_cast<const std::vector<T> &>(comb));

        // Colex successor: find the LEFTMOST index that can be incremented
        // without colliding with its right-hand neighbour (or with n, for the
        // last index), bump it, and reset everything to its left back to
        // {0, 1, ...}. When no index can move we have emitted the last one.
        int j = 0;
        while (j < k) {
            const int limit = (j + 1 < k) ? indices[j + 1] : n;
            if (indices[j] + 1 < limit)
                break;
            j++;
        }

        if (j == k)
            break;

        indices[j]++;
        for (int t = 0; t < j; t++)
            indices[t] = t;
    }
}

// All k-subsets of `items`, in colexicographic order of the underlying
// indices, materialised into a vector. Same contract as
// for_each_combination above; prefer that one when the caller only needs to
// walk the sequence once, since it avoids holding every combination at once.
template <typename T>
std::vector<std::vector<T> > combinations_of(const std::vector<T> &items, int k)
{
    std::vector<std::vector<T> > out;
    out.reserve(num_combinations(static_cast<int>(items.size()), k));
    for_each_combination(items, k,
        [&out](const std::vector<T> &c) { out.push_back(c); });
    return out;
}

#endif
