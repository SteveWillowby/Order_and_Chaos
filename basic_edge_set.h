#include "nt_wrappers/nauty_traces.h"

#include<unordered_set>

#ifndef SCHENO__BASIC_EDGE_SET_H
#define SCHENO__BASIC_EDGE_SET_H

// This EdgeSet subclass is just a wrapper for unordered_set
//
// Note that it COPIES the unorderd_set
class BasicEdgeSetPair : public EdgeSetPair {
public:
    BasicEdgeSetPair(const std::unordered_set<Edge, EdgeHash>& e,
                     const std::unordered_set<Edge, EdgeHash>& ne) :
            sets(std::pair<std::unordered_set<Edge, EdgeHash>,
                           std::unordered_set<Edge, EdgeHash>>(e, ne)) {}

    std::pair<std::unordered_set<Edge, EdgeHash>,
              std::unordered_set<Edge, EdgeHash>>
                                edges_and_non_edges() { return sets; }

    long double heuristic_score() const { return 0.0; }

protected:
    std::pair<std::unordered_set<Edge, EdgeHash>,
              std::unordered_set<Edge, EdgeHash>> sets;
};

#endif
