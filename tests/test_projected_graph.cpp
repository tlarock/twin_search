#include "hypergraph.hpp"
#include "projected_graph.hpp"
#include "test_hypergraphs.cpp"
#include <gtest/gtest.h>


TEST(ProjectedGraphTest, FindCliquesh1) {
    Hypergraph h = h1();
    ProjectedGraph proj(h);
    std::map<int, std::vector<std::vector<int > > >  gt;
    gt[2] = std::vector<std::vector<int>> {
        {0, 1},
        {0, 2},
        {1, 2},
        {1, 3},
        {2, 3},
        {1, 4},
        {1, 5},
        {2, 4},
        {2, 5},
        {4, 5}
    };

    gt[3] = std::vector<std::vector<int> > {
        {0, 1, 2},
        {1, 2, 3},
        {1, 2, 4},
        {1, 2, 5},
        {1, 4, 5},
        {2, 4, 5}
    };

    gt[4] = std::vector<std::vector<int> > {
        {1, 2, 4, 5}
    };

    CliqueMap cliques = proj.compute_cliques(2, h.n);
    EXPECT_EQ(cliques.size(), gt.size());
    for (std::size_t k = 2; k < 5; ++k) {
        EXPECT_EQ(cliques[k].size(), gt[k].size());
        for (std::size_t i = 0; i < gt[k].size(); ++i) {
            EXPECT_TRUE(cliques[k].contains(gt[k][i]));
        }
    }
}

TEST(ProjectedGraphTest, FindCliquesh2) {
    Hypergraph h = h2();
    ProjectedGraph proj(h);
    std::map<int, std::vector<std::vector<int > > >  gt;
    gt[2] = std::vector<std::vector<int>> { 
        { 0, 2 },
        { 0, 6 },
        { 1, 2 },
        { 1, 3 },
        { 1, 4 },
        { 1, 6 },
        { 1, 7 },
        { 1, 8 },
        { 2, 4 },
        { 2, 6 },
        { 2, 7 },
        { 2, 9 },
        { 3, 6 },
        { 4, 5 },
        { 4, 8 },
        { 4, 9 },
        { 5, 8 },
        { 7, 8 },
        { 7, 9 },
        { 8, 9 },
    };
    gt[3] = std::vector<std::vector<int>> { 
        { 0, 2, 6 },
        { 1, 2, 4 },
        { 1, 2, 6 },
        { 1, 2, 7 },
        { 1, 3, 6 },
        { 1, 4, 8 },
        { 1, 7, 8 },
        { 2, 4, 9 },
        { 2, 7, 9 },
        { 4, 5, 8 },
        { 4, 8, 9 },
        { 7, 8, 9 },
    };
    CliqueMap cliques = proj.compute_cliques(2, h.n);
    EXPECT_EQ(cliques.size(), gt.size());
    for (std::size_t k = 2; k < 5; ++k) {
        EXPECT_EQ(cliques[k].size(), gt[k].size());
        for (std::size_t i = 0; i < gt[k].size(); ++i) {
            EXPECT_TRUE(cliques[k].contains(gt[k][i]));
        }
    }
}

TEST(ProjectedGraphTest, FindCliquesh3) {
    Hypergraph h = GM();
    ProjectedGraph proj(h);
    std::map<int, std::vector<std::vector<int > > >  gt;
    gt[2] = std::vector<std::vector<int>> { 
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 0, 4 },
        { 1, 2 },
        { 1, 3 },
        { 1, 4 },
        { 2, 3 },
        { 2, 4 },
        { 2, 5 },
        { 3, 4 },
        { 3, 5 },
        { 4, 5 },
    };
    gt[3] = std::vector<std::vector<int>> { 
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 1, 4 },
        { 0, 2, 3 },
        { 0, 2, 4 },
        { 0, 3, 4 },
        { 1, 2, 3 },
        { 1, 2, 4 },
        { 1, 3, 4 },
        { 2, 3, 4 },
        { 2, 3, 5 },
        { 2, 4, 5 },
        { 3, 4, 5 },
    };
    gt[4] = std::vector<std::vector<int>> { 
        { 0, 1, 2, 3 },
        { 0, 1, 2, 4 },
        { 0, 1, 3, 4 },
        { 0, 2, 3, 4 },
        { 1, 2, 3, 4 },
        { 2, 3, 4, 5 },
    };
    gt[5] = std::vector<std::vector<int>> { 
        { 0, 1, 2, 3, 4 },
    };
    CliqueMap cliques = proj.compute_cliques(2, h.n);
    EXPECT_EQ(cliques.size(), gt.size());
    for (std::size_t k = 2; k < 5; ++k) {
        EXPECT_EQ(cliques[k].size(), gt[k].size());
        for (std::size_t i = 0; i < gt[k].size(); ++i) {
            EXPECT_TRUE(cliques[k].contains(gt[k][i]));
        }
    }
}


TEST(ProjectedGraphTest, FindCliquesh10) {
    Hypergraph h = h10();
    ProjectedGraph proj(h);
    std::map<int, std::vector<std::vector<int > > >  gt;
    gt[2] = std::vector<std::vector<int>> { 
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 0, 5 },
        { 1, 3 },
        { 1, 4 },
        { 1, 5 },
        { 2, 3 },
        { 2, 4 },
        { 2, 5 },
        { 3, 4 },
    };
    gt[3] = std::vector<std::vector<int>> { 
        { 0, 1, 3 },
        { 0, 1, 5 },
        { 0, 2, 3 },
        { 0, 2, 5 },
        { 1, 3, 4 },
        { 2, 3, 4 },
    };
    CliqueMap cliques = proj.compute_cliques(2, h.n);
    EXPECT_EQ(cliques.size(), gt.size());
    for (std::size_t k = 2; k < 5; ++k) {
        EXPECT_EQ(cliques[k].size(), gt[k].size());
        for (std::size_t i = 0; i < gt[k].size(); ++i) {
            EXPECT_TRUE(cliques[k].contains(gt[k][i]));
        }
    }
}

TEST(ProjectedGraphTest, FindCliquesh11) {
    Hypergraph h = h11();
    ProjectedGraph proj(h);
    std::map<int, std::vector<std::vector<int > > >  gt;
    gt[2] = std::vector<std::vector<int>> { 
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 0, 4 },
        { 0, 5 },
        { 1, 2 },
        { 1, 3 },
        { 1, 4 },
        { 1, 5 },
        { 2, 3 },
        { 2, 4 },
        { 3, 4 },
    };
    gt[3] = std::vector<std::vector<int>> { 
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 1, 4 },
        { 0, 1, 5 },
        { 0, 2, 3 },
        { 0, 2, 4 },
        { 0, 3, 4 },
        { 1, 2, 3 },
        { 1, 2, 4 },
        { 1, 3, 4 },
        { 2, 3, 4 },
    };
    gt[4] = std::vector<std::vector<int>> { 
        { 0, 1, 2, 3 },
        { 0, 1, 2, 4 },
        { 0, 1, 3, 4 },
        { 0, 2, 3, 4 },
        { 1, 2, 3, 4 },
    };
    gt[5] = std::vector<std::vector<int>> { 
        { 0, 1, 2, 3, 4 },
    };
    CliqueMap cliques = proj.compute_cliques(2, h.n);
    EXPECT_EQ(cliques.size(), gt.size());
    for (std::size_t k = 2; k < 5; ++k) {
        EXPECT_EQ(cliques[k].size(), gt[k].size());
        for (std::size_t i = 0; i < gt[k].size(); ++i) {
            EXPECT_TRUE(cliques[k].contains(gt[k][i]));
        }
    }
}
