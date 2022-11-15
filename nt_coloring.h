#include "coloring.h"

#include<vector>

#ifndef SYM__NT_COLORING_H
#define SYM__NT_COLORING_H

class NTColoring {
public:
    NTColoring(const std::vector<int>& color_list);
    NTColoring(const Coloring<int>& coloring);

    int* get_node_ids();
    int* get_partition_info();

protected:
    std::vector<int> node_ids;
    std::vector<int> partition_info;
};

#endif