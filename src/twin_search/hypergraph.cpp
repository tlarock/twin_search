#include "./hypergraph.hpp"



// Default constructor
Hypergraph::Hypergraph() {
    Hypergraph::n = 0;
    Hypergraph::m = 0;
}

// Print the hypergraph to the console, separating nodes with commas
// and hyperedges with spaces.
void Hypergraph::pretty_print() {
    // Named rather than per-line, so the whole dump is emitted as one unit.
    SyncStream out(std::cout);
    for (const auto& [he_idx, he]: hyperedges)
    {
        for (std::size_t i = 0; i < he.size()-1; i++)
        {
            out << he[i] << ",";
        }
        out << he.back() << " ";
    }
    out << "\n";
}

// Returns the bipartite representation of a hypergraph
UndirectedGraph Hypergraph::get_bipartite() {
    UndirectedGraph bipartite;
    // Node ids will map to themselves
    // hyperedges will map to he_idx+h.n
    for (const auto& [he_idx, he] : hyperedges) {
        // get the constiuent edge of he_idx
        for(int u : he) {
            boost::add_edge(u, he_idx+n, bipartite);
        }
    }
    return bipartite;
}

ublas::matrix<int> Hypergraph::get_incidence_matrix() {
    // The third argument is required. ublas::matrix(size1, size2) allocates
    // without initialising, so the entries this function does not explicitly
    // set to 1 would otherwise hold whatever was in that memory. That matters a
    // long way downstream: get_lg_mat() multiplies this matrix by its own
    // transpose, and get_line_graph() then does
    //     for (int v = 0; v < lg_mat(r,c); ++v) boost::add_edge(r, c, ...);
    // so a single large garbage entry asks boost to add billions of edges and
    // the process is killed. Whether it happens at all depends on whether the
    // allocator hands back fresh (kernel-zeroed) pages or dirty ones, which is
    // why it can lie dormant and then appear after an unrelated change.
    ublas::matrix<int> incidence(hyperedges.size(), n, 0);
    for (std::size_t row = 0; row < incidence.size1(); ++row) {
        for (auto node : hyperedges[row]) {
            incidence(row, node) = 1;
        }
    }
    return incidence;
}

// Returns an UndirectedGraph corresponding to the line
// graph of the hypergraph
// NOTE: Ignores self-loops because boost's isomorphism
// test does not deal with them correctly
UndirectedGraph Hypergraph::get_line_graph() {
    ublas::matrix<int> lg_mat = get_lg_mat();
    UndirectedGraph line_graph;
    for (std::size_t r = 0; r < lg_mat.size1(); ++r) {
        for (std::size_t c = 0; c < lg_mat.size2(); ++c) {
            if (r == c)
                continue;

            if (lg_mat(r, c) > 0) {
                for (int v = 0; v < lg_mat(r,c); ++v)
                    boost::add_edge(r, c, line_graph);
            }
        }
    }
    return line_graph;
}

ublas::matrix<int> Hypergraph::get_lg_mat() {
    ublas::matrix<int> incidence = get_incidence_matrix();
    ublas::matrix<int> lg_mat = ublas::prod(incidence, ublas::trans(incidence));
    return lg_mat;
}

// Stores input hyperedge `he` at index `he_idx`, sorted and with any repeated
// nodes removed, dropping it entirely if an identical hyperedge was already
// stored. See the declaration in hypergraph.hpp for why both repairs are
// required rather than merely tidy.
//
// Returns true if the hyperedge was stored.
bool Hypergraph::store_hyperedge(int he_idx, const std::vector<int> &he, std::set<int> &nodes,
                                 std::set<std::vector<int> > &seen, InputRepairs &repairs) {
    std::vector<int> cleaned(he.begin(), he.end());
    sort(cleaned.begin(), cleaned.end());

    const std::size_t size_before = cleaned.size();
    cleaned.erase(std::unique(cleaned.begin(), cleaned.end()), cleaned.end());
    if (cleaned.size() != size_before)
        repairs.hyperedges_with_repeated_nodes += 1;

    // Checked after de-duplication: {0,1,1} and {0,1} are the same hyperedge
    // once repeats are gone, so the collision only becomes visible here.
    if (!seen.insert(cleaned).second) {
        repairs.duplicate_hyperedges += 1;
        return false;
    }

    Hypergraph::hyperedges[he_idx] = cleaned;

    // Size and memberships must come from the cleaned hyperedge, not the input,
    // or a repeated node would still be counted twice.
    Hypergraph::hyperedge_sizes[he_idx] = static_cast<int>(cleaned.size());
    for (int node_id : cleaned) {
        Hypergraph::node_memberships[node_id].push_back(he_idx);
        nodes.insert(node_id);
    }
    return true;
}

// Tells the user their input was altered. Uses SyncStream so that the
// message is not interleaved with other threads' output: these constructors
// run inside TBB loops in count_twins_random and exhaustive_search_projections.
void Hypergraph::report_input_repairs(const InputRepairs &repairs) {
    if (!repairs.any())
        return;

    SyncStream out(std::cerr);   // see diagnostic() in utils.hpp
    out << "Warning: input hypergraph was not simple and has been modified.\n";
    if (repairs.hyperedges_with_repeated_nodes > 0) {
        out << "  - " << repairs.hyperedges_with_repeated_nodes
            << " hyperedge(s) contained a repeated node; the repeats were removed.\n";
    }
    if (repairs.duplicate_hyperedges > 0) {
        out << "  - " << repairs.duplicate_hyperedges
            << " hyperedge(s) duplicated an earlier hyperedge and were dropped.\n";
    }
    out << "  Results below describe the modified hypergraph, not the input as given."
        << std::endl;
}

// Throws if any node id falls outside 0,...,n-1.
//
// The constructors that take n from the caller cannot renumber to fit without
// silently changing what the results mean, so an out-of-range id is an error
// rather than something to repair. The error names the offending id and points
// at the constructor that does remap.
void Hypergraph::require_ids_in_range(const std::set<int> &nodes, int n) {
    if (nodes.empty())
        return;

    const int lowest = *nodes.begin();
    const int highest = *nodes.rbegin();
    if (lowest >= 0 && highest < n)
        return;

    throw std::out_of_range(std::format(
        "Hypergraph: node id {} is outside 0..{} for the given n={}. This "
        "constructor does not renumber nodes, because doing so would mean the "
        "node ids in any output no longer matched the ids supplied. Either "
        "pass an n that covers the ids, or use the single-argument Hypergraph "
        "constructor, which remaps ids to 0..n-1 and says so.",
        lowest < 0 ? lowest : highest, n - 1, n));
}

// Emits a note that node ids were renumbered. Called by whoever knows the ids
// came from a human, not by the constructor: the exhaustive driver builds a
// Hypergraph per candidate sub-hypergraph and most of those legitimately do not
// span every node, so reporting from the constructor produced hundreds of
// kilobytes of warnings on an ordinary run.
//
// This matters for provenance rather than correctness: the hypergraph is
// faithfully relabelled, but every node id written to an output file afterwards
// is a rank, not the id that was read in.
void Hypergraph::report_remapping() const {
    if (!nodes_were_remapped)
        return;

    diagnostic() << "Note: input node ids were not 0.." << (n - 1)
                 << ", so they have been renumbered to 0.." << (n - 1)
                 << " in ascending order. Node ids in any output refer to the "
                    "renumbered nodes, not the input ids."
                 << std::endl;
}
