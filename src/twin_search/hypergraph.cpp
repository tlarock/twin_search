#include "./hypergraph.hpp"

#include <syncstream>


// Default constructor
Hypergraph::Hypergraph() {
    Hypergraph::n = 0;
    Hypergraph::m = 0;
}

// Print the hypergraph to the console, separating nodes with commas
// and hyperedges with spaces.
void Hypergraph::pretty_print() {
    for (const auto& [he_idx, he]: hyperedges)
    {
        for (std::size_t i = 0; i < he.size()-1; i++)
        {
            std::cout << he[i] << ",";
        }
        std::cout << he.back() << " ";
    }
    std::cout << std::endl;
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
    ublas::matrix<int> incidence(hyperedges.size(), n);
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

// Tells the user their input was altered. Uses std::osyncstream so that the
// message is not interleaved with other threads' output: these constructors
// run inside TBB loops in count_twins_random and exhaustive_search_projections.
void Hypergraph::report_input_repairs(const InputRepairs &repairs) {
    if (!repairs.any())
        return;

    std::osyncstream out(std::cerr);
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
