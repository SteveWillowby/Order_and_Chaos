#include "graph.h"
#include "sparse_graph.h"

#include<string>

#ifndef SYM__FILE_UTILS_H
#define SYM__FILE_UTILS_H

// This function assumes that the nodes are numbered 0 through the largest
//  node ID found in the edge list.
SparseGraph read_graph(const bool directed,
                       const std::string& edgelist_filename);

// Note that if the nodes in the node list are not labeled 0 through n-1, they
//  will be relabeled in sorted order as they are loaded.
SparseGraph read_graph(const bool directed,
                       const std::string& nodelist_filename,
                       const std::string& edgelist_filename);

// Note: if nodelist is empty then no nodelist is written
void write_graph(const Graph& g, const std::string& nodelist_filename,
                                 const std::string& edgelist_filename);

#endif
