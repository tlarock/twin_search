#ifndef TWIN_SEARCH_CHECKPOINT_H
#define TWIN_SEARCH_CHECKPOINT_H
#include <string>
#include <vector>
#include "projected_graph.hpp"

// On-disk form of a drained search.
//
// Hyperedges, not clique-node ids. Ids come from FactorGraph's construction
// order; they are probably stable, but "probably stable" is not a thing a file
// format should rest on, and a wrong id silently produces a valid-looking
// search of the wrong tree. Hyperedges cost k bytes instead of 2 and are
// self-describing.
//
// proj_rem is NOT stored for the frontier items either: it is recomputable
// from the projection and the chosen cliques (TwinSearch::residual_for), and
// it is several times the size of the item it belongs to.
struct SearchCheckpoint {
    int format = 1;
    std::string commit;           // refuse a checkpoint from a different build
    int n = 0, m = 0, k = 0, min_k = 0, max_k = 0;
    unsigned int seed = 0;
    int index = 0;                // sample index within the cell
    int attempt = 0;              // retry counter; nonzero under a width limit
    int width_limit_exp = 0;
    // Upper triangle of the projection, '-' separated. The integrity check
    // that matters: on resume the sample is redrawn, and if the redrawn
    // projection differs - wrong seed, wrong standard library, changed
    // generator - the frontier describes a different tree entirely.
    std::string proj;
    std::vector<std::vector<std::vector<int> > > twins;     // found so far
    std::vector<std::vector<std::vector<int> > > frontier;  // still to explore
};

// Upper triangle, '-' separated. Must stay byte-identical to the key the
// --dry-run path emits, or checkpoints and fingerprints disagree about what
// identifies a sample.
std::string twin_projection_key(const ProjMatT &M);

// Atomic: writes to <path>.tmp, flushes, then renames. A reader therefore sees
// either the previous checkpoint or the new one, never a half-written file -
// which matters because the common reason for writing one is that the process
// is about to die.
bool write_checkpoint(const std::string &path, const SearchCheckpoint &cp,
                      std::string &err);

// Rejects anything it cannot fully account for: truncated body, counts that
// disagree with the lists, missing end sentinel.
bool read_checkpoint(const std::string &path, SearchCheckpoint &cp,
                     std::string &err);

#endif
