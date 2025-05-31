#include "./hypergraph.hpp"


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
