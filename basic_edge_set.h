#include "edge.h"

#include<unordered_set>

#ifndef SYM__BASIC_EDGE_SET_H
#define SYM__BASIC_EDGE_SET_H

// This EdgeSet subclass is just a wrapper for unordered_set
//
// Note that it COPIES the unorderd_set
class BasicEdgeSet : public EdgeSet {
public:
    BasicEdgeSet(const std::unordered_set<Edge, EdgeHash>& e) : e_set(e) {}

    virtual const std::unordered_set<Edge, EdgeHash>& edges() { return e_set; }

protected:
    const std::unordered_set<Edge, EdgeHash> e_set;
};

#endif
