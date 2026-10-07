#include "search_checkpoint.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>

std::string twin_projection_key(const ProjMatT &M) {
    std::string s;
    s.reserve(M.size1() * M.size1() * 2);
    for (std::size_t i = 0; i < M.size1(); i++)
        for (std::size_t j = i + 1; j < M.size2(); j++) {
            if (!s.empty()) s += '-';
            s += std::to_string(M(i, j));
        }
    return s;
}

namespace {

void write_graphs(std::ostream &os,
                  const std::vector<std::vector<std::vector<int> > > &gs) {
    for (const std::vector<std::vector<int> > &g : gs) {
        for (std::size_t i = 0; i < g.size(); i++) {
            if (i) os << '|';
            for (std::size_t j = 0; j < g[i].size(); j++) {
                if (j) os << ':';
                os << g[i][j];
            }
        }
        os << '\n';
    }
}

bool parse_graph(const std::string &line, std::vector<std::vector<int> > &out) {
    out.clear();
    std::stringstream ls(line);
    std::string he;
    while (std::getline(ls, he, '|')) {
        if (he.empty()) return false;
        std::vector<int> nodes;
        std::stringstream hs(he);
        std::string tok;
        while (std::getline(hs, tok, ':')) {
            if (tok.empty()) return false;
            try { nodes.push_back(std::stoi(tok)); }
            catch (const std::exception &) { return false; }
        }
        if (nodes.empty()) return false;
        out.push_back(nodes);
    }
    // An empty partial hypergraph is legitimate - it is the root of the search,
    // which is exactly what a drain at the very first node parks.
    return true;
}

}  // namespace

bool write_checkpoint(const std::string &path, const SearchCheckpoint &cp,
                      std::string &err) {
    const std::string tmp = path + ".tmp";
    {
        std::ofstream os(tmp, std::ios::out | std::ios::trunc);
        if (!os) { err = "cannot open " + tmp; return false; }
        os << "twin_search-checkpoint\t" << cp.format << '\n'
           << "commit\t" << cp.commit << '\n'
           << "n\t" << cp.n << '\n' << "m\t" << cp.m << '\n' << "k\t" << cp.k << '\n'
           << "min_k\t" << cp.min_k << '\n' << "max_k\t" << cp.max_k << '\n'
           << "seed\t" << cp.seed << '\n' << "index\t" << cp.index << '\n'
           << "attempt\t" << cp.attempt << '\n'
           << "width_limit_exp\t" << cp.width_limit_exp << '\n'
           << "proj\t" << cp.proj << '\n'
           << "twins\t" << cp.twins.size() << '\n'
           << "frontier\t" << cp.frontier.size() << '\n'
           << "#twins\n";
        write_graphs(os, cp.twins);
        os << "#frontier\n";
        write_graphs(os, cp.frontier);
        os << "#end\n";
        os.flush();
        if (!os) { err = "write failed for " + tmp; return false; }
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        err = "rename failed: " + tmp + " -> " + path;
        return false;
    }
    return true;
}

bool read_checkpoint(const std::string &path, SearchCheckpoint &cp,
                     std::string &err) {
    std::ifstream is(path);
    if (!is) { err = "cannot open " + path; return false; }
    cp = SearchCheckpoint();

    std::string line;
    std::size_t want_twins = 0, want_frontier = 0;
    bool seen_counts = false;
    while (std::getline(is, line)) {
        if (line == "#twins") { seen_counts = true; break; }
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) { err = "malformed header line: " + line; return false; }
        const std::string key = line.substr(0, tab), val = line.substr(tab + 1);
        try {
            if (key == "twin_search-checkpoint") cp.format = std::stoi(val);
            else if (key == "commit") cp.commit = val;
            else if (key == "n") cp.n = std::stoi(val);
            else if (key == "m") cp.m = std::stoi(val);
            else if (key == "k") cp.k = std::stoi(val);
            else if (key == "min_k") cp.min_k = std::stoi(val);
            else if (key == "max_k") cp.max_k = std::stoi(val);
            else if (key == "seed") cp.seed = static_cast<unsigned int>(std::stoul(val));
            else if (key == "index") cp.index = std::stoi(val);
            else if (key == "attempt") cp.attempt = std::stoi(val);
            else if (key == "width_limit_exp") cp.width_limit_exp = std::stoi(val);
            else if (key == "proj") cp.proj = val;
            else if (key == "twins") { want_twins = std::stoul(val); }
            else if (key == "frontier") { want_frontier = std::stoul(val); }
        } catch (const std::exception &) {
            err = "bad value for " + key + ": " + val;
            return false;
        }
    }
    if (!seen_counts) { err = "truncated before #twins"; return false; }
    if (cp.format != 1) {
        err = "unsupported checkpoint format " + std::to_string(cp.format);
        return false;
    }

    bool seen_frontier = false, seen_end = false;
    while (std::getline(is, line)) {
        if (line == "#frontier") { seen_frontier = true; continue; }
        if (line == "#end") { seen_end = true; break; }
        std::vector<std::vector<int> > g;
        if (!parse_graph(line, g)) { err = "malformed body line: " + line; return false; }
        (seen_frontier ? cp.frontier : cp.twins).push_back(g);
    }
    // The sentinel is the whole point of the atomic write: without it a file
    // truncated by a kill mid-write would parse as a shorter but valid
    // checkpoint, and the resume would silently drop part of the search tree.
    if (!seen_end) { err = "no #end sentinel - checkpoint is truncated"; return false; }
    if (!seen_frontier) { err = "no #frontier section"; return false; }
    if (cp.twins.size() != want_twins) {
        err = "twin count mismatch: header says " + std::to_string(want_twins)
            + ", body has " + std::to_string(cp.twins.size());
        return false;
    }
    if (cp.frontier.size() != want_frontier) {
        err = "frontier count mismatch: header says " + std::to_string(want_frontier)
            + ", body has " + std::to_string(cp.frontier.size());
        return false;
    }
    return true;
}
