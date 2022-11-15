#include "coloring.h"

#include<vector>

#ifndef SYM__NT_COLORING_H
#define SYM__NT_COLORING_H

class NTColoring {
public:
    // Takes colors for nodes 0 through n - 1.
    //  The colors list can have any values.
    NTColoring(const std::vector<int>& color_list);

    int* get_node_ids();
    int* get_partition_ints();

    // Call after Nauty or Traces has run to update the internal variables to
    //  match node_ids and partition_ints.
    //
    // Nauty and Traces may change those when running with this coloring.
    //
    // O(n)
    virtual void repair();

    virtual int num_cells() const;

    // Moves a node up one cell. If the node is already by itself in the
    //  highest cell, nothing happens.
    //
    // O(1)
    virtual void increment_cell(int node);

    // Moves a node down one cell. If the node is already by itself in the
    //  lowest cell, nothing happens.
    //
    // Runtime:
    //  O(n) if creates a new cell
    //  O(1) otherwise
    virtual void decrement_cell(int node);

    // Effectively calls increment or decrement cell repeatedly until the
    //  node has the proper cell value.
    //
    // Requires that -1 <= `cell` <= num cells.
    virtual void set_cell(int node, int cell);

    // O(1)
    virtual int get_cell(int node);

protected:
    std::vector<int> node_ids;
    std::vector<int> partition_ints;
    // A partition int at location i should be 0 iff i is the last listed node
    //  in the cell.

    // Example partition for {0, 2, 5}, {1}, {3, 4}
    // node_ids:       2, 5, 0, 1, 3, 4
    // partition_ints: 1, 1, 0, 0, 1, 0

    // node_locations[node_ids[i]] = i
    std::vector<int> node_locations;
    std::vector<int> node_cells;
    std::vector<int> cell_starts;

    // Cell ends is exclusive.
    std::vector<int> cell_ends;
};

#endif
