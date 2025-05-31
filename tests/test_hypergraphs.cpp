#ifndef TEST_HYPERGRAPHS
#define TEST_HYPERGRAPHS

#include "hypergraph.hpp"

/*
 * util file for storing some simple example hypergraphs that
 * may be useful for testing.
 */

Hypergraph h1() {
    // non-uniform hypergraph
    // n: 6
    // m = 3
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 2};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 2, 3};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {1, 2, 4, 5};
    input_hyperedges.insert(he3);

    Hypergraph H(input_hyperedges);

    return H;
}

Hypergraph h2() {
    // 3-uniform hypergraph
    // n: 10
    // m: 7
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 2, 6};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 3, 6};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {4, 5, 8};
    input_hyperedges.insert(he3);
    std::vector<int> he4 {1, 2, 4};
    input_hyperedges.insert(he4);
    std::vector<int> he5 {1, 7, 8};
    input_hyperedges.insert(he5);
    std::vector<int> he6 {2, 7, 9};
    input_hyperedges.insert(he6);
    std::vector<int> he7 {4, 8, 9};
    input_hyperedges.insert(he7);

    Hypergraph H(input_hyperedges);

    return H;
}

Hypergraph h3() {
    // 3-uniform hypergraph
    // n: 10
    // m: 8
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {1, 3, 4};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 2, 9};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {1, 7, 9};
    input_hyperedges.insert(he3);
    std::vector<int> he4 {2, 5, 8};
    input_hyperedges.insert(he4);
    std::vector<int> he5 {0, 5, 6};
    input_hyperedges.insert(he5);
    std::vector<int> he6 {0, 7, 9};
    input_hyperedges.insert(he6);
    std::vector<int> he7 {2, 5, 9};
    input_hyperedges.insert(he7);
    std::vector<int> he8 {2, 6, 7};
    input_hyperedges.insert(he8);

    Hypergraph H(input_hyperedges);

    return H;
}

Hypergraph h4() {
    // 3-uniform hypergraph
    // n: 10
    // m: 6
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 4};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 3, 7};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {1, 4, 7};
    input_hyperedges.insert(he3);
    std::vector<int> he4 {2, 5, 7};
    input_hyperedges.insert(he4);
    std::vector<int> he5 {4, 5, 7};
    input_hyperedges.insert(he5);
    std::vector<int> he6 {6, 8, 9};
    input_hyperedges.insert(he6);

    Hypergraph H(input_hyperedges);

    return H;
}

Hypergraph h5() {
    // 3-uniform hypergraph
    // n: 10
    // m: 10
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 5};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 4, 7};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {2, 8, 9};
    input_hyperedges.insert(he3);
    std::vector<int> he4 {0, 1, 2};
    input_hyperedges.insert(he4);
    std::vector<int> he5 {1, 2, 9};
    input_hyperedges.insert(he5);
    std::vector<int> he6 {2, 4, 9};
    input_hyperedges.insert(he6);
    std::vector<int> he7 {0, 3, 9};
    input_hyperedges.insert(he7);
    std::vector<int> he8 {0, 6, 9};
    input_hyperedges.insert(he8);
    std::vector<int> he9 {1, 3, 4};
    input_hyperedges.insert(he9);
    std::vector<int> he10 {1, 4, 6};
    input_hyperedges.insert(he10);

    Hypergraph H(input_hyperedges);

    return H;
}

Hypergraph h6() {
    // CASE where node id 8 is missing
    // from the input so things need to be remapped.
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 2};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {1, 2, 3};
    input_hyperedges.insert(he3);
    std::vector<int> he4 {3, 4};
    input_hyperedges.insert(he4);
    std::vector<int> he5 {6, 7};
    input_hyperedges.insert(he5);
    std::vector<int> he6 {9, 10, 11};
    input_hyperedges.insert(he6);
    std::vector<int> he7 {11, 12};
    input_hyperedges.insert(he7);
    std::vector<int> he8 {5, 11, 12, 13, 14};
    input_hyperedges.insert(he8);
    std::vector<int> he9 {14, 15, 16};
    input_hyperedges.insert(he9);

    Hypergraph H(input_hyperedges);

    return H;

}

Hypergraph h7() {
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 6};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {1, 2, 8};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {1, 2, 4, 10};
    input_hyperedges.insert(he3);

    Hypergraph h(input_hyperedges);

    return h;
}

Hypergraph h8() {
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 2};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {0, 1, 3};
    input_hyperedges.insert(he2);

    Hypergraph h(input_hyperedges);

    return h;
}

Hypergraph GM() {
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1};
    input_hyperedges.insert(he1);
    std::vector<int> he2 {0, 1, 2, 3};
    input_hyperedges.insert(he2);
    std::vector<int> he3 {0, 2, 3, 4};
    input_hyperedges.insert(he3);
    std::vector<int> he4 {1, 2, 3, 4};
    input_hyperedges.insert(he4);
    std::vector<int> he5 {2, 4, 5};
    input_hyperedges.insert(he5);
    std::vector<int> he6 {3, 4, 5};
    input_hyperedges.insert(he6);

    Hypergraph h(input_hyperedges);
    return h;
}

Hypergraph h9() {
    std::set<std::vector<int> > input_hyperedges;
    std::vector<int> he1 {0, 1, 2};
    input_hyperedges.insert(he1);

    Hypergraph h(input_hyperedges);

    return h;
}

std::vector<Hypergraph> all_hypergraphs() {
    return std::vector<Hypergraph> {h1(), h2(), h3(), h4(), h5(), h6(), h7(), h8(), GM()};
}

#endif
