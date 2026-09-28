// One-shot generator: emits the golden reference table for
// tests/test_combinations.cpp, computed from discreture so the baked values
// record DISCRETURE's behaviour, not our reimplementation of it.
#include <cstdio>
#include <vector>
#include <discreture.hpp>

// FNV-1a over the flattened emitted sequence, with separators so that
// regroupings (e.g. {0,1},{2,3} vs {0},{1,2,3}) cannot collide.
static unsigned long long digest(const std::vector<std::vector<int> > &seq) {
    unsigned long long h = 1469598103934665603ULL;
    auto mix = [&h](unsigned long long v) {
        h ^= v; h *= 1099511628211ULL;
    };
    mix(static_cast<unsigned long long>(seq.size()));
    for (const auto &c : seq) {
        mix(0xFFFF);                       // combination separator
        mix(static_cast<unsigned long long>(c.size()));
        for (int x : c)
            mix(static_cast<unsigned long long>(x) + 1);
    }
    return h;
}

static std::vector<std::vector<int> > from_discreture(int n, int k) {
    std::vector<int> items(n);
    for (int i = 0; i < n; i++) items[i] = i;
    std::vector<std::vector<int> > out;
    auto combs = discreture::combinations(items, k);
    for (auto &&comb : combs) {
        std::vector<int> c;
        for (int x : comb) c.push_back(x);
        out.push_back(c);
    }
    return out;
}

int main() {
    // --- explicit sequences, small n, human-checkable ---
    std::printf("// Explicit expected sequences for n = 0..5, all k.\n");
    std::printf("// Generated from discreture::combinations (colex order).\n");
    std::printf("struct ExplicitCase { int n; int k; std::vector<std::vector<int> > expected; };\n");
    std::printf("static const std::vector<ExplicitCase> EXPLICIT_CASES = {\n");
    for (int n = 0; n <= 5; n++) {
        for (int k = 0; k <= n; k++) {
            auto seq = from_discreture(n, k);
            std::printf("    {%d, %d, {", n, k);
            for (std::size_t i = 0; i < seq.size(); i++) {
                std::printf("{");
                for (std::size_t j = 0; j < seq[i].size(); j++)
                    std::printf("%d%s", seq[i][j], j + 1 < seq[i].size() ? "," : "");
                std::printf("}%s", i + 1 < seq.size() ? "," : "");
            }
            std::printf("}},\n");
        }
    }
    std::printf("};\n\n");

    // --- digests, wider n, compact ---
    std::printf("// FNV-1a digests of the full emitted sequence for n = 0..16, all k.\n");
    std::printf("// Generated from discreture::combinations (colex order).\n");
    std::printf("struct DigestCase { int n; int k; unsigned long long digest; };\n");
    std::printf("static const std::vector<DigestCase> DIGEST_CASES = {\n");
    for (int n = 0; n <= 16; n++)
        for (int k = 0; k <= n; k++)
            std::printf("    {%2d, %2d, %20lluULL},\n", n, k, digest(from_discreture(n, k)));
    std::printf("};\n");
    return 0;
}
