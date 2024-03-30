/* A hashable edge. */

#include<stdexcept>
#include<unordered_set>
#include<utility>

#ifndef SCHENO__EDGE_H
#define SCHENO__EDGE_H

// Note: if 2^32 < n^2, increase this size and adjust SCHENO__MAX_EDGE_LABEL
//  (that's when n > 65536)
typedef uint32_t SCHENO__edge_int_type;
#define SCHENO__MAX_EDGE_LABEL 0xFFFFFFFF

#define EDGE(s, t, d) std::pair<int, int>((d || s < t ? s : t), (d || s < t ? t : s))

typedef std::pair<int, int> Edge;

// Note, if 2 * sizeof(int) > sizeof(long) and node ID's exceed
//      2^(sizeof(int) * 4), then this will have issues.
struct EdgeHash {
    size_t operator()(std::pair<int, int> const& p) const noexcept {
        if (2 * sizeof(int) > sizeof(long)) {
            return std::hash<size_t>{}(p.first << (sizeof(int) * 4) | p.second);
        }
        // Put the two integers in different halves of the long's bits.
        long combo = long(p.first) << (sizeof(int) * 8) | p.second;
        return std::hash<long>{}(combo);
    }
};

// Virtual base class for a set of edges.
class EdgeSetPair {
public:
    EdgeSetPair() {};

    virtual std::pair<std::unordered_set<Edge, EdgeHash>,
                      std::unordered_set<Edge, EdgeHash>>
                                edges_and_non_edges() =0;

    virtual long double heuristic_score() const =0;
};

#endif
