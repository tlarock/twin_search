#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <map>
#include <set>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/adjacency_list.hpp>

// wrapper type for an UndirectedGraph
// NOTE: Using boost::undirectedS with the boost::vf2_graph_iso does not
// properly handle self-loops, so make sure not to include them!
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> UndirectedGraph;

// convenience type for returning a map from k->set of cliques of size k
typedef std::map<int, std::set<std::vector<int> > > CliqueMap;

// A synchronized output stream for warnings and diagnostics.
//
// Everything written through a SyncStream is buffered and emitted in one piece
// when it is destroyed, so messages cannot interleave even when several threads
// report at once. Unlike the per-file mutexes this replaces, that guarantee
// holds across translation units, which is what those mutexes could not do -
// each only excluded threads writing from its own file.
//
// This is std::osyncstream in all but name. It is written out by hand because
// Apple libc++ ships no <syncstream>, and this code has to build with the
// system compiler on macOS as well as with libstdc++ on Linux.
class SyncStream {
  public:
    explicit SyncStream(std::ostream &os) : out_(&os) {}
    SyncStream(SyncStream &&o) noexcept : out_(o.out_), buf_(std::move(o.buf_)) { o.out_ = nullptr; }
    SyncStream(const SyncStream &) = delete;
    SyncStream &operator=(const SyncStream &) = delete;
    ~SyncStream() { emit(); }

    template <typename T> SyncStream &operator<<(const T &v) { buf_ << v; return *this; }
    // Manipulators such as std::endl are function pointers, not values.
    SyncStream &operator<<(std::ostream &(*manip)(std::ostream &)) { buf_ << manip; return *this; }

    // Emits whatever has been buffered so far and resets the buffer. Called by
    // the destructor; public so a long-lived stream can flush early.
    void emit() {
        if (!out_)
            return;
        const std::string s = buf_.str();
        if (!s.empty()) {
            // One lock for every SyncStream in the process, which is what makes
            // the guarantee hold across translation units.
            static std::mutex m;
            const std::lock_guard<std::mutex> lock(m);
            *out_ << s;
            out_->flush();
        }
        buf_.str(std::string());
    }

  private:
    std::ostream *out_;
    std::ostringstream buf_;
};

// Returns a synchronized stream for warnings and diagnostics.
//
// Intended as a full statement, so the temporary is destroyed (and the message
// emitted) at the end of it:
//     diagnostic() << "something happened" << std::endl;
// For a multi-line report, name the stream so the whole report is emitted as a
// unit:
//     SyncStream out(std::cerr);
//
// Diagnostics go to stderr so that they cannot contaminate results written to
// stdout by the drivers.
inline SyncStream diagnostic() { return SyncStream(std::cerr); }

unsigned int factorial(unsigned int n);
unsigned int binom_exact(unsigned int, unsigned int);

// binomial coefficient based on boost::math::beta
double binom(double, double);

// sum of a 2d boost::ublas matrix<int>
int matsum(const boost::numeric::ublas::matrix<int> &m);

// maximum element of a 2d boost::ublas matrix<int>
int matmax(const boost::numeric::ublas::matrix<int> &m);

// Sum of values in a map of ints
int mapsum(std::map<int, int> &m);

// Takes a map of ints and returns vector of keys sorted by value
std::vector<int> sorted_keys_by_value(std::map<int, int> &m);

#endif
