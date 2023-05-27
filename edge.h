/* A hashable edge. */

#include<stdexcept>
#include<unordered_set>
#include<utility>

#ifndef SYM__EDGE_H
#define SYM__EDGE_H

// Note: if 2^32 < n^2, increase this size and adjust SYM__MAX_EDGE_LABEL
//  (that's when n > 65536)
typedef uint32_t SYM__edge_int_type;
#define SYM__MAX_EDGE_LABEL 0xFFFFFFFF

#define EDGE(s, t, d) std::pair<int, int>((d || s < t ? s : t), (d || s < t ? t : s))

typedef std::pair<int, int> Edge;

struct EdgeHash {
    size_t operator()(std::pair<int, int> const& p) const noexcept {
        if (2 * sizeof(int) >= sizeof(long)) {
            return std::hash<size_t>{}(std::hash<int>{}(p.first) ^ p.second);
        }
        // Put the two integers in different halves of the long's bits.
        long combo = long(p.first) << sizeof(int) | p.second;
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
};

#endif
