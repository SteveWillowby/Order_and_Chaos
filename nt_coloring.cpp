#include "nt_coloring.h"

#include<ordered_map>
#include<vector>

NTColoring::NTColoring(const std::vector<int>& color_list) {
    node_cells = std::vector<int>(color_list);
    colors = std::ordered_map<int>();

    // Count occurrences of each color, while sorting the colors.
    for (auto color = node_cells.begin(); color != node_cells.end(); color++) {
        auto itr = colors.find(*color);
        if (itr == colors.end()) {
            colors[*color] = 0;
        } else {
            colors[*color]++;
        }
    }

    partition_ints = std::vector<int>(node_cells.size(), 1);

    cell_starts = std::vector<int>(colors.size(), 0);
    cell_ends = std::vector<int>(colors.size(), 0);

    // Set up cell_starts.
    //  cell_ends will be used below to place nodes and thus will not have the
    //  correct values yet.
    int cumulative_sizes = 0;
    int color_count;
    int color;
    for (auto itr = colors.begin(); itr != colors.end(); itr++) {
        color = itr->first;
        color_count = itr->second;
        cell_starts[color] = cumulative_sizes;
        cell_ends[color] = cumulative_sizes;
        cumulative_sizes += color_count;
        partition_ints[cumulative_sizes - 1] = 0; // Mark the end of a partition
    }

    node_locations = std::vector<int>(node_cells.size(), 0);
    node_ids = std::vector<int>(node_cells.size(), 0);

    for (int node = 0; node < int(node_cells.size()); node++) {
        color = node_cells[node];
        node_locations[node] = cell_ends[color];
        node_ids[cell_ends[color]] = node;
        cell_ends[color]++;
    }
}

int* NTColoring::get_node_ids() {
    return &(node_ids[0]);
}

int* NTColoring::get_partition_ints() {
    return &(partition_ints[0]);
}

void NTColoring::repair() {
    cell_starts.clear()
    cell_starts.push_back(0);
    cell_ends.clear();

    int cell = 0;
    int node;
    for (int i = 0; i < int(node_cells.size()); i++) {
        node = node_ids[i];
        if (partition_ints[i] == 0) {
            cell_ends.push_back(i + 1);
            cell_starts.push_back(i + 1);
            cell++;
        }

        node_locations[node] = i;
        node_cells[node] = cell;
    }

    cell_starts.pop_back();
}

int NTColoring::num_cells() const {
    return cell_starts.size();
}

void NTColoring::increment_cell(int node) {

}

void NTColoring::decrement_cell(int node) {

}

void NTColoring::set_cell(int node, int cell) {

}

int NTColoring::get_cell(int node) {

}
