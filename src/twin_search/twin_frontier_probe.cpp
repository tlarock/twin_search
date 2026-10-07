// Step 0 for checkpointing: how big is the search frontier, and does it grow?
//
// The straggler problem is that a single sample cannot be SPLIT. resume_cell.sh
// resumes at sample granularity, so a sample needing 30 hours never finishes on
// a 24-hour partition however often it is resubmitted. Serialising the search
// would fix that - EITHER by resuming one process, OR, better, by handing the
// frontier's independent subtrees to separate jobs that run concurrently.
//
// Whether any of that is worth building depends on one number: how many nodes
// are outstanding when the search is interrupted. If the frontier is small and
// plateaus, a checkpoint is cheap and can be taken often. If it grows without
// bound, checkpointing dies the same way memoising the search did - on memory
// rather than on principle. See docs/reuse-experiments.md, section 0: check the
// space cost BEFORE building the time optimisation.
//
// For each sample this drains at several wall-clock points and reports the
// frontier size, so growth is visible rather than inferred from one reading.
//
// Usage:
//   twin_frontier_probe -n 9 -m 37 -k 3 --seed 20260930 --indices 423 \
//                       --drain-at 0.25,0.5,1,2 [--max-threads 8] [--full]
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <cstdio>
#include <sstream>
#include <vector>

#include <oneapi/tbb.h>
#include <oneapi/tbb/global_control.h>

#include "argparse/argparse.hpp"

#include "hypergraph.hpp"
#include "projected_graph.hpp"
#include "twin_search.hpp"
#include "search_checkpoint.hpp"
#include "random_hypergraph_generators.hpp"
#include "utils.hpp"

struct ProbeArgs : public argparse::Args {
    int &n = kwarg("n, num-nodes", "Number of nodes.").set_default(9);
    int &m = kwarg("m", "Number of hyperedges.").set_default(37);
    int &k = kwarg("k", "Hyperedge dimension.").set_default(3);
    int &min_k = kwarg("min-k", "Minimum hyperedge size. <=0 means k.").set_default(0);
    int &max_k = kwarg("max-k", "Maximum hyperedge size. <=0 means k.").set_default(0);
    unsigned int &seed = kwarg("seed", "Master RNG seed.").set_default(20260930u);
    std::string &indices = kwarg("indices", "Comma-separated sample indices.").set_default(std::string("0"));
    std::string &drain_at = kwarg("drain-at", "Comma-separated wall-clock seconds at which to drain. One run per value.").set_default(std::string("1"));
    int &max_threads = kwarg("max-threads", "TBB parallelism limit.").set_default(0);
    bool &verify_resume = flag("verify-resume", "End-to-end check: run the sample to completion, then run it again draining at each --drain-at point, write a checkpoint, read it back, resume from it, and compare. The gate for the whole idea - an interrupted search must give the SAME answer, not a similar one.");
    std::string &ckpt_dir = kwarg("checkpoint-dir", "Where --verify-resume writes its checkpoint files.").set_default(std::string("/tmp"));
    bool &full = flag("full", "Also run the sample to completion first, for a denominator. Can be very slow on exactly the samples this is aimed at.");
};

static std::mt19937 sample_generator(unsigned int seed, int i, int attempt) {
    std::seed_seq seq{static_cast<std::uint32_t>(seed),
                      static_cast<std::uint32_t>(i),
                      static_cast<std::uint32_t>(attempt)};
    return std::mt19937(seq);
}

static bool draw(int n, int m, int k, unsigned int seed, int i, Hypergraph &out) {
    for (int attempt = 0; attempt < 1000; attempt++) {
        std::mt19937 gen = sample_generator(seed, i, attempt);
        Hypergraph h = sample_uniform_random(n, m, k, gen);
        if (h.m == m && h.n == n) { out = h; return true; }
    }
    return false;
}

static std::vector<double> parse_doubles(const std::string &csv) {
    std::vector<double> out;
    std::stringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) if (!tok.empty()) out.push_back(std::stod(tok));
    return out;
}

int main(int argc, char *argv[]) {
    std::cout.setf(std::ios::unitbuf);
    auto args = argparse::parse<ProbeArgs>(argc, argv);
    const int n = args.n, m = args.m, k = args.k;
    const int min_k = args.min_k > 0 ? args.min_k : k;
    const int max_k = args.max_k > 0 ? args.max_k : k;

    std::unique_ptr<oneapi::tbb::global_control> thread_limit;
    if (args.max_threads > 0)
        thread_limit = std::make_unique<oneapi::tbb::global_control>(
            oneapi::tbb::global_control::max_allowed_parallelism,
            static_cast<std::size_t>(args.max_threads));

    std::vector<int> idxs;
    {
        std::stringstream ss(args.indices);
        std::string tok;
        while (std::getline(ss, tok, ',')) if (!tok.empty()) idxs.push_back(std::stoi(tok));
    }
    const std::vector<double> drains = parse_doubles(args.drain_at);

    std::cout << "# twin_frontier_probe  n=" << n << " m=" << m << " k=" << k
              << " seed=" << args.seed << "\n";
    if (args.verify_resume)
        std::cout << "index\tdrain_s\tfrontier\tbanked\tref_twins\tres_twins"
                     "\tref_mates\tres_mates\tref_s\tresume_s\tverdict\n";
    else
        std::cout << "index\tdrain_s\twall_s\tfrontier\ttwins_so_far\tcnode_entries"
                     "\tbytes_cnode\tbytes_hyperedge\n";

    for (int i : idxs) {
        Hypergraph h;
        if (!draw(n, m, k, args.seed, i, h)) { std::cout << i << "\tDRAW-FAILED\n"; continue; }
        ProjectedGraph proj(h);

        if (args.verify_resume) {
            TwinSearch ref(proj, min_k, max_k, true, true, false, false);
            if (!ref.feasible) { std::cout << i << "\tINFEASIBLE\n"; continue; }
            const auto r0 = std::chrono::steady_clock::now();
            ref.parallel_search(true);
            const double ref_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - r0).count();
            std::vector<std::vector<int> > ref_t = ref.twins;
            for (std::vector<int> &g : ref_t) std::sort(g.begin(), g.end());
            std::sort(ref_t.begin(), ref_t.end());

            for (double d : drains) {
                TwinSearch part(proj, min_k, max_k, true, true, false, false);
                part.arm_drain(d);
                part.parallel_search(true);
                if (!part.drained) {
                    std::cout << "# i=" << i << " drain " << d << "s: completed first\n";
                    continue;
                }
                SearchCheckpoint cp;
                cp.n = n; cp.m = m; cp.k = k; cp.min_k = min_k; cp.max_k = max_k;
                cp.seed = args.seed; cp.index = i;
                cp.proj = twin_projection_key(proj.proj_mat);
                for (const std::vector<int> &g : part.twins)
                    cp.twins.push_back(part.inflate_cnodes(g));
                for (const std::vector<int> &g : part.frontier)
                    cp.frontier.push_back(part.inflate_cnodes(g));

                const std::string path = args.ckpt_dir + "/twin-ckpt-" + std::to_string(i);
                std::string err;
                if (!write_checkpoint(path, cp, err)) { std::cout << "# write failed: " << err << "\n"; continue; }
                SearchCheckpoint back;
                if (!read_checkpoint(path, back, err)) { std::cout << "# read failed: " << err << "\n"; continue; }
                if (back.proj != cp.proj) { std::cout << "# projection mismatch on reload\n"; continue; }

                TwinSearch res(proj, min_k, max_k, true, true, false, false);
                std::vector<std::vector<int> > fr, tw;
                for (const std::vector<std::vector<int> > &g : back.frontier) {
                    std::vector<int> cn;
                    for (const std::vector<int> &he : g) cn.push_back(res.cnode_for(he));
                    fr.push_back(cn);
                }
                for (const std::vector<std::vector<int> > &g : back.twins) {
                    std::vector<int> cn;
                    for (const std::vector<int> &he : g) cn.push_back(res.cnode_for(he));
                    tw.push_back(cn);
                }
                const auto s0 = std::chrono::steady_clock::now();
                const bool ok = res.parallel_search_from(fr, tw, true);
                const double res_s = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - s0).count();
                std::vector<std::vector<int> > res_t = res.twins;
                for (std::vector<int> &g : res_t) std::sort(g.begin(), g.end());
                std::sort(res_t.begin(), res_t.end());

                const bool same = ok && res_t == ref_t
                                  && res.mates.size() == ref.mates.size()
                                  && res.filtered_twins.size() == ref.filtered_twins.size();
                std::cout << i << "\t";
                std::printf("%.2f\t", d);
                std::cout << part.frontier.size() << "\t" << part.twins.size() << "\t"
                          << ref.twins.size() << "\t" << res.twins.size() << "\t"
                          << ref.mates.size() << "\t" << res.mates.size() << "\t";
                std::printf("%.2f\t%.2f\t", ref_s, res_s);
                std::cout << (same ? "IDENTICAL" : "MISMATCH") << "\n";
                std::remove(path.c_str());
            }
            continue;
        }

        if (args.full) {
            TwinSearch ts(proj, min_k, max_k, true, true, false, false);
            if (ts.feasible) {
                const auto t0 = std::chrono::steady_clock::now();
                ts.parallel_search(true);
                const double secs = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t0).count();
                std::cout << i << "\tFULL\t";
                std::printf("%.2f\t", secs);
                std::cout << "0\t" << ts.twins.size() << "\t0\t0\t0\n";
            }
        }

        for (double d : drains) {
            TwinSearch ts(proj, min_k, max_k, true, true, false, false);
            if (!ts.feasible) { std::cout << i << "\tINFEASIBLE\n"; break; }
            ts.arm_drain(d);
            const auto t0 = std::chrono::steady_clock::now();
            ts.parallel_search(true);
            const double secs = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            if (!ts.drained) {
                // Finished before the deadline: the frontier is empty because
                // there was nothing left, not because nothing was outstanding.
                std::cout << i << "\t";
                std::printf("%.2f\t%.2f\t", d, secs);
                std::cout << "COMPLETED\t" << ts.twins.size() << "\t0\t0\t0\n";
                break;
            }
            std::size_t entries = 0;
            for (const std::vector<int> &g : ts.frontier) entries += g.size();
            // Two serialisation sizes. cnode ids are compact but depend on
            // FactorGraph's construction order; explicit hyperedges are k bytes
            // each at n<=256 and survive any change to how ids are assigned,
            // which is what a checkpoint format should actually store.
            const std::size_t b_cnode = entries * 2;
            const std::size_t b_he = entries * static_cast<std::size_t>(k);
            std::cout << i << "\t";
            std::printf("%.2f\t%.2f\t", d, secs);
            std::cout << ts.frontier.size() << "\t" << ts.twins.size() << "\t"
                      << entries << "\t" << b_cnode << "\t" << b_he << "\n";
        }
    }
    return 0;
}
