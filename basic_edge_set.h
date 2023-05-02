#include "edge.h"

#include<unordered_set>

#ifndef SYM__BASIC_EDGE_SET_H
#define SYM__BASIC_EDGE_SET_H

// This EdgeSet subclass is just a wrapper for unordered_set
//
// Note that it COPIES the unorderd_set
class BasicEdgeSetPair : public EdgeSetPair {
public:
    BasicEdgeSetPair(const std::unordered_set<Edge, EdgeHash>& e,
                     const std::unordered_set<Edge, EdgeHash>& ne) :
            sets(std::pair<std::unordered_set<Edge, EdgeHash>,
                           std::unordered_set<Edge, EdgeHash>>(e, ne)) {}

    virtual std::pair<std::unordered_set<Edge, EdgeHash>,
                      std::unordered_set<Edge, EdgeHash>>
                                edges_and_non_edges() const { return sets; }

protected:
    std::pair<std::unordered_set<Edge, EdgeHash>,
              std::unordered_set<Edge, EdgeHash>> sets;
};

#endif
