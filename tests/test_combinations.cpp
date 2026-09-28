#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <numeric>
#include <set>
#include <thread>
#include <atomic>
#include "combinations.hpp"

// combinations_of() replaced discreture::combinations() on the TwinSearch
// search hot path, because discreture::binomial() memoizes into an
// unsynchronized function-local static and so is unsafe to call from several
// threads at once (see combinations.hpp).
//
// The reference data below was GENERATED FROM discreture::combinations while
// that dependency was still present, by scripts/gen_golden.cpp. It therefore
// records discreture's actual behaviour rather than a restatement of our own.
// Regenerate it only if you deliberately intend to change the enumeration
// order - doing so silently changes which twin TwinSearch::filtered_twins keeps
// as each isomorphism class representative.

// ---------------------------------------------------------------------------
// BEGIN generated reference data - see scripts/gen_golden.cpp
// ---------------------------------------------------------------------------
// Explicit expected sequences for n = 0..5, all k.
// Generated from discreture::combinations (colex order).
struct ExplicitCase { int n; int k; std::vector<std::vector<int> > expected; };
static const std::vector<ExplicitCase> EXPLICIT_CASES = {
    {0, 0, {{}}},
    {1, 0, {{}}},
    {1, 1, {{0}}},
    {2, 0, {{}}},
    {2, 1, {{0},{1}}},
    {2, 2, {{0,1}}},
    {3, 0, {{}}},
    {3, 1, {{0},{1},{2}}},
    {3, 2, {{0,1},{0,2},{1,2}}},
    {3, 3, {{0,1,2}}},
    {4, 0, {{}}},
    {4, 1, {{0},{1},{2},{3}}},
    {4, 2, {{0,1},{0,2},{1,2},{0,3},{1,3},{2,3}}},
    {4, 3, {{0,1,2},{0,1,3},{0,2,3},{1,2,3}}},
    {4, 4, {{0,1,2,3}}},
    {5, 0, {{}}},
    {5, 1, {{0},{1},{2},{3},{4}}},
    {5, 2, {{0,1},{0,2},{1,2},{0,3},{1,3},{2,3},{0,4},{1,4},{2,4},{3,4}}},
    {5, 3, {{0,1,2},{0,1,3},{0,2,3},{1,2,3},{0,1,4},{0,2,4},{1,2,4},{0,3,4},{1,3,4},{2,3,4}}},
    {5, 4, {{0,1,2,3},{0,1,2,4},{0,1,3,4},{0,2,3,4},{1,2,3,4}}},
    {5, 5, {{0,1,2,3,4}}},
};

// FNV-1a digests of the full emitted sequence for n = 0..16, all k.
// Generated from discreture::combinations (colex order).
struct DigestCase { int n; int k; unsigned long long digest; };
static const std::vector<DigestCase> DIGEST_CASES = {
    { 0,  0,  9169359956432111873ULL},
    { 1,  0,  9169359956432111873ULL},
    { 1,  1,   987891489224502077ULL},
    { 2,  0,  9169359956432111873ULL},
    { 2,  1, 16505733451279107208ULL},
    { 2,  2,  2213671816526180420ULL},
    { 3,  0,  9169359956432111873ULL},
    { 3,  1,  1359596773354915716ULL},
    { 3,  2,  9223824645211818689ULL},
    { 3,  3, 18094442962846891832ULL},
    { 4,  0,  9169359956432111873ULL},
    { 4,  1,  6507723924854231741ULL},
    { 4,  2, 17323765313859517427ULL},
    { 4,  3,  4599602333607208385ULL},
    { 4,  4, 18083629052557239173ULL},
    { 5,  0,  9169359956432111873ULL},
    { 5,  1,  8720315516504994073ULL},
    { 5,  2, 12059660578509516887ULL},
    { 5,  3,   153563334521093101ULL},
    { 5,  4, 17373350174805731057ULL},
    { 5,  5,  6454123808866224257ULL},
    { 6,  0,  9169359956432111873ULL},
    { 6,  1,  7299984716848407236ULL},
    { 6,  2, 10461511192233672602ULL},
    { 6,  3, 12064967771767090331ULL},
    { 6,  4,  7418037735932461523ULL},
    { 6,  5,  2523518178647791056ULL},
    { 6,  6, 16377205101906446416ULL},
    { 7,  0,  9169359956432111873ULL},
    { 7,  1, 15364931235088993104ULL},
    { 7,  2,   753279790070195735ULL},
    { 7,  3, 10921994934536545614ULL},
    { 7,  4, 13931911954482663007ULL},
    { 7,  5,  6110185683549297918ULL},
    { 7,  6,  4539645725294609681ULL},
    { 7,  7,  9379327045201546508ULL},
    { 8,  0,  9169359956432111873ULL},
    { 8,  1,  7915671357966582889ULL},
    { 8,  2,  3435396058210993677ULL},
    { 8,  3, 12450694560915448165ULL},
    { 8,  4, 16275441688761434383ULL},
    { 8,  5,  8001531311547349181ULL},
    { 8,  6, 11929156983177319805ULL},
    { 8,  7, 14702734129687087777ULL},
    { 8,  8, 16548206944067404929ULL},
    { 9,  0,  9169359956432111873ULL},
    { 9,  1, 13411933070233747885ULL},
    { 9,  2, 12340876350088674461ULL},
    { 9,  3, 10448310858156811689ULL},
    { 9,  4, 15849918363237860519ULL},
    { 9,  5,   170197188340956365ULL},
    { 9,  6, 14624838571172300917ULL},
    { 9,  7,  2420186539017163333ULL},
    { 9,  8, 15715501274552624297ULL},
    { 9,  9, 17467573325362277565ULL},
    {10,  0,  9169359956432111873ULL},
    {10,  1, 15286279338524020480ULL},
    {10,  2, 12562288882362503872ULL},
    {10,  3,  5933941667826069449ULL},
    {10,  4, 13204360259948437799ULL},
    {10,  5, 15863730096809246667ULL},
    {10,  6,  6914887587943740639ULL},
    {10,  7,  4240762638218811981ULL},
    {10,  8,  6013122043199273789ULL},
    {10,  9, 13382105517342409712ULL},
    {10, 10,  6120692046392655372ULL},
    {11,  0,  9169359956432111873ULL},
    {11,  1, 10152792163723329820ULL},
    {11,  2,  4975575651856523805ULL},
    {11,  3,  4368395087284372860ULL},
    {11,  4,  8579687829826246511ULL},
    {11,  5,  8152162449979196175ULL},
    {11,  6,   672611682916182511ULL},
    {11,  7, 16454013999717951211ULL},
    {11,  8,  2374939787143168489ULL},
    {11,  9, 16621431361717163284ULL},
    {11, 10, 13728828302174908337ULL},
    {11, 11, 12302010304456495744ULL},
    {12,  0,  9169359956432111873ULL},
    {12,  1,  7889614659999033773ULL},
    {12,  2,  1875750773997116039ULL},
    {12,  3,  9596388154275545989ULL},
    {12,  4,  2441468007048041783ULL},
    {12,  5,   107104491420021557ULL},
    {12,  6,  7517200393151554297ULL},
    {12,  7,  4995453745965606133ULL},
    {12,  8, 11444385057471276695ULL},
    {12,  9, 16434226186375382553ULL},
    {12, 10, 14984308503530681647ULL},
    {12, 11, 11125696606752448905ULL},
    {12, 12,  3682929635546929797ULL},
    {13,  0,  9169359956432111873ULL},
    {13,  1, 11735152116091605769ULL},
    {13,  2, 14860769959948567747ULL},
    {13,  3,  7665212011001329445ULL},
    {13,  4,  4665287997572672027ULL},
    {13,  5,   891150559559457509ULL},
    {13,  6, 14634134681554695337ULL},
    {13,  7,  8511018644215916057ULL},
    {13,  8, 16960289979101706807ULL},
    {13,  9,  7196603592390344073ULL},
    {13, 10, 10780146451842920907ULL},
    {13, 11, 17479516078798596081ULL},
    {13, 12, 11695661828226122897ULL},
    {13, 13, 14934805713873366017ULL},
    {14,  0,  9169359956432111873ULL},
    {14,  1,  8316941436773853372ULL},
    {14,  2, 18248523403486310054ULL},
    {14,  3, 10424560277427983387ULL},
    {14,  4, 13493888146709590925ULL},
    {14,  5, 12279678167073782294ULL},
    {14,  6,  3374792239440220518ULL},
    {14,  7, 13598889845953573285ULL},
    {14,  8,  6654422932789791463ULL},
    {14,  9,  6397344378337090878ULL},
    {14, 10, 10437387955679889100ULL},
    {14, 11,  3728918338645085967ULL},
    {14, 12,  2620877688597523327ULL},
    {14, 13,  6174181483276173848ULL},
    {14, 14,   846721207508948904ULL},
    {15,  0,  9169359956432111873ULL},
    {15,  1, 11760617405398793720ULL},
    {15,  2,  2075411606301454219ULL},
    {15,  3,   688353835681707850ULL},
    {15,  4,  1224238670067250521ULL},
    {15,  5,  1261088031826609832ULL},
    {15,  6,  8339130185531798027ULL},
    {15,  7,   895470786910690194ULL},
    {15,  8, 13122493436373797903ULL},
    {15,  9,  7351535897245762766ULL},
    {15, 10, 14891815317837322633ULL},
    {15, 11,  4381464907052115932ULL},
    {15, 12, 14920549483589499679ULL},
    {15, 13,  9178085768102616198ULL},
    {15, 14,  3542570751771156089ULL},
    {15, 15,  6210518941028659668ULL},
    {16,  0,  9169359956432111873ULL},
    {16,  1, 16117530858796869081ULL},
    {16,  2,  9404986275092696001ULL},
    {16,  3, 16956289441887975889ULL},
    {16,  4,  9908481119728724757ULL},
    {16,  5,  6846322658722253089ULL},
    {16,  6, 15077299157577185593ULL},
    {16,  7,  1731056717458738153ULL},
    {16,  8, 12687661351240154263ULL},
    {16,  9, 10574994482152169673ULL},
    {16, 10, 14951047698320048097ULL},
    {16, 11,  4147096316218456529ULL},
    {16, 12, 16224362267681659685ULL},
    {16, 13,  4286016400889993817ULL},
    {16, 14, 17748289900386386337ULL},
    {16, 15, 16962201185399775337ULL},
    {16, 16,  8140303738287806913ULL},
};
// ---------------------------------------------------------------------------
// END generated reference data
// ---------------------------------------------------------------------------

// Must match gen_golden.cpp exactly.
static unsigned long long digest(const std::vector<std::vector<int> > &seq) {
    unsigned long long h = 1469598103934665603ULL;
    auto mix = [&h](unsigned long long v) {
        h ^= v; h *= 1099511628211ULL;
    };
    mix(static_cast<unsigned long long>(seq.size()));
    for (const auto &c : seq) {
        mix(0xFFFF);
        mix(static_cast<unsigned long long>(c.size()));
        for (int x : c)
            mix(static_cast<unsigned long long>(x) + 1);
    }
    return h;
}

static std::vector<int> iota_items(int n) {
    std::vector<int> items(n);
    std::iota(items.begin(), items.end(), 0);
    return items;
}

// Element-for-element equality with discreture, small n, fully spelled out.
TEST(CombinationsTest, MatchesDiscretureExplicitSmallN) {
    for (const ExplicitCase &c : EXPLICIT_CASES) {
        EXPECT_EQ(combinations_of(iota_items(c.n), c.k), c.expected)
            << "n=" << c.n << " k=" << c.k;
    }
}

// Same equality, carried out to n=16 via a digest of the whole sequence.
TEST(CombinationsTest, MatchesDiscretureDigestUpToN16) {
    for (const DigestCase &c : DIGEST_CASES) {
        EXPECT_EQ(digest(combinations_of(iota_items(c.n), c.k)), c.digest)
            << "n=" << c.n << " k=" << c.k;
    }
}

// Self-validating properties, independent of any reference data: the output
// must be exactly the set of k-subsets, each strictly increasing, with no
// duplicates, in strictly ascending colex order.
TEST(CombinationsTest, IsExactlyTheKSubsetsInStrictColexOrder) {
    auto colex_less = [](const std::vector<int> &a, const std::vector<int> &b) {
        // compare from the largest element downwards
        for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
            if (a[i] != b[i]) return a[i] < b[i];
        }
        return false;
    };

    for (int n = 0; n <= 12; n++) {
        for (int k = 0; k <= n; k++) {
            std::vector<std::vector<int> > got = combinations_of(iota_items(n), k);

            std::set<std::vector<int> > unique;
            for (const auto &c : got) {
                ASSERT_EQ(static_cast<int>(c.size()), k) << "n=" << n << " k=" << k;
                for (std::size_t i = 1; i < c.size(); i++)
                    ASSERT_LT(c[i-1], c[i]) << "not increasing, n=" << n << " k=" << k;
                for (int x : c)
                    ASSERT_TRUE(x >= 0 && x < n) << "out of range, n=" << n << " k=" << k;
                unique.insert(c);
            }
            EXPECT_EQ(unique.size(), got.size()) << "duplicates, n=" << n << " k=" << k;

            for (std::size_t i = 1; i < got.size(); i++)
                ASSERT_TRUE(colex_less(got[i-1], got[i]))
                    << "colex order violated at " << i << ", n=" << n << " k=" << k;

            // exact count via Pascal's rule, no floating point
            std::vector<std::vector<long long> > C(n + 1, std::vector<long long>(n + 1, 0));
            for (int a = 0; a <= n; a++) {
                C[a][0] = 1;
                for (int b = 1; b <= a; b++)
                    C[a][b] = C[a-1][b-1] + (b <= a-1 ? C[a-1][b] : 0);
            }
            EXPECT_EQ(static_cast<long long>(got.size()), C[n][k])
                << "n=" << n << " k=" << k;
        }
    }
}

// Ordering is colex in the INDICES, not in the values: combinations_of must
// index into `items` rather than assume the values are the indices.
TEST(CombinationsTest, ColexByIndexNotByValue) {
    std::vector<int> descending {50, 40, 30, 20};
    // index pairs (0,1) (0,2) (1,2) (0,3) (1,3) (2,3)
    std::vector<std::vector<int> > expected {
        {50, 40}, {50, 30}, {40, 30}, {50, 20}, {40, 20}, {30, 20}
    };
    EXPECT_EQ(combinations_of(descending, 2), expected);

    std::vector<int> spread {7, 14, 21, 28, 35};
    std::vector<std::vector<int> > got = combinations_of(spread, 3);
    std::vector<std::vector<int> > want {
        {7,14,21}, {7,14,28}, {7,21,28}, {14,21,28},
        {7,14,35}, {7,21,35}, {14,21,35}, {7,28,35}, {14,28,35}, {21,28,35}
    };
    EXPECT_EQ(got, want);
}

// TwinSearch::get_combinations passes weight = a residual matrix entry, which
// is >= 1 but is not bounded by the number of remaining clique-neighbours.
TEST(CombinationsTest, KGreaterThanNYieldsNothing) {
    std::vector<int> items {3, 1, 4, 1, 5};
    for (int k = static_cast<int>(items.size()) + 1; k < 20; k++)
        EXPECT_TRUE(combinations_of(items, k).empty()) << "k=" << k;
}

TEST(CombinationsTest, NegativeKYieldsNothing) {
    std::vector<int> items {1, 2, 3};
    EXPECT_TRUE(combinations_of(items, -1).empty());
    EXPECT_TRUE(combinations_of(items, -7).empty());
}

// k == 0 yields exactly one (empty) combination, matching discreture's
// binomial(n, 0) == 1. TwinSearch never calls with weight 0, but the contract
// is pinned regardless.
TEST(CombinationsTest, ZeroKYieldsOneEmptyCombination) {
    std::vector<int> items {1, 2, 3};
    std::vector<std::vector<int> > mine = combinations_of(items, 0);
    ASSERT_EQ(mine.size(), 1u);
    EXPECT_TRUE(mine[0].empty());

    std::vector<int> empty_items;
    EXPECT_EQ(combinations_of(empty_items, 0).size(), 1u);
    EXPECT_TRUE(combinations_of(empty_items, 1).empty());
}

// num_combinations() only sizes a reserve(), so it may saturate, but it must
// be exact below the cap and must never over-report.
TEST(CombinationsTest, NumCombinationsExactBelowCapAndSaturates) {
    for (int n = 0; n <= 14; n++) {
        for (int k = 0; k <= n; k++) {
            std::size_t actual = combinations_of(iota_items(n), k).size();
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

// The streaming form must agree with the materialising form.
TEST(CombinationsTest, ForEachMatchesMaterialisedForm) {
    for (int n = 0; n <= 10; n++) {
        for (int k = 0; k <= n; k++) {
            std::vector<std::vector<int> > streamed;
            for_each_combination(iota_items(n), k,
                [&](const std::vector<int> &c) { streamed.push_back(c); });
            EXPECT_EQ(streamed, combinations_of(iota_items(n), k))
                << "n=" << n << " k=" << k;
        }
    }
}

// The point of the replacement: no shared mutable state. Hammer it from many
// threads and require identical results to the single-threaded answer. Under
// ThreadSanitizer this also asserts the absence of a data race.
TEST(CombinationsTest, ThreadSafeUnderConcurrentUse) {
    std::vector<std::vector<int> > expected = combinations_of(iota_items(30), 2);

    std::atomic<int> mismatches{0};
    std::vector<std::thread> threads;
    for (int t = 0; t < 8; t++) {
        threads.emplace_back([&, t]() {
            for (int rep = 0; rep < 200; rep++) {
                // vary n per thread: this is what drove discreture's memo
                // table to resize concurrently
                int n = 7 + ((t * 13 + rep * 7) % 30);
                std::vector<std::vector<int> > got = combinations_of(iota_items(n), 2);
                std::size_t want = static_cast<std::size_t>(n) * (n - 1) / 2;
                if (got.size() != want) mismatches++;
            }
            if (combinations_of(iota_items(30), 2) != expected) mismatches++;
        });
    }
    for (auto &th : threads) th.join();
    EXPECT_EQ(mismatches.load(), 0);
}
