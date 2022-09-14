/* A hashable edge. */

#include<stdexcept>
#include<utility>

#ifndef SYM__EDGE_H
#define SYM__EDGE_H

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

#endif
