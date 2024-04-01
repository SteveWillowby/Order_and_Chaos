#include "nt_partition.h"

#include<cstddef>
#include<map>
#include<stdexcept>
#include<string>
#include<unordered_map>
#include<vector>

NTPartition::NTPartition(const std::vector<int>& cell_list)
                                        : _size(cell_list.size()) {
    std::map<int, size_t> colors = std::map<int, size_t>();

    // Count occurrences of each color, while sorting the colors.
    for (auto color = cell_list.begin(); color != cell_list.end(); color++) {
        auto itr = colors.find(*color);
        if (itr == colors.end()) {
            colors[*color] = 0;
        } else {
            colors[*color]++;
        }
    }

    partition_ints = new int[_size];
    node_ids = new int[_size];

    std::vector<int> cell_ends = std::vector<int>(colors.size(), 0);
    std::unordered_map<int, int> color_to_idx = std::unordered_map<int, int>();

    // cell_ends will be used below to place nodes.
    int cumulative_sizes = 0;
    int color_count;
    int color = 0;
    for (auto itr = colors.begin(); itr != colors.end(); itr++) {
        color_count = itr->second;
        cell_ends[color] = cumulative_sizes;
        cumulative_sizes += color_count;
        partition_ints[cumulative_sizes - 1] = 0; // Mark the end of a partition
        color_to_idx[itr->first] = color;
        color++;
    }

    for (int node = 0; node < int(_size); node++) {
        color = color_to_idx[cell_list[node]];
        node_ids[cell_ends[color]] = node;
        cell_ends[color]++;
    }
}

NTPartition::NTPartition(const NTPartition& ntp) : _size(ntp.size()) {
    node_ids = new int[_size];
    partition_ints = new int[_size];
    int* ntp_ni = ntp.get_node_ids();
    int* ntp_pi = ntp.get_partition_ints();

    for (size_t i = 0; i < _size; i++) {
        node_ids[i] = ntp_ni[i];
        partition_ints[i] = ntp_pi[i];
    }
}

NTPartition::NTPartition(const size_t s) : _size(s) {
    node_ids = new int[_size];
    partition_ints = new int[_size];

    for (int i = 0; i < int(_size - 1); i++) {
        node_ids[i] = i;
        partition_ints[i] = 1;
    }
    node_ids[_size - 1] = _size - 1;
    partition_ints[_size - 1] = 0;
}

NTPartition::~NTPartition() {
    delete node_ids;
    delete partition_ints;
}

NTPartition& NTPartition::operator=(const NTPartition& ntp) {
    if (this == &ntp) {
        return *this;
    }

    if (_size != ntp.size()) {
        _size = ntp.size();
        delete node_ids;
        delete partition_ints;
        node_ids = new int[_size];
        partition_ints = new int[_size];
    }

    for (size_t i = 0; i < _size; i++) {
        node_ids[i] = ntp.get_node_id(i);
        partition_ints[i] = ntp.get_partition_int(i);
    }

    return *this;
}

int* NTPartition::get_node_ids() const {
    return node_ids;
}

int* NTPartition::get_partition_ints() const {
    return partition_ints;
}

size_t NTPartition::size() const {
    return _size;
}

int NTPartition::get_node_id(size_t i) const {
    if (i >= _size) {
        throw std::range_error(std::string("Error! Index ")
                               + std::to_string(i)
                               + " outside the range of partition of size "
                               + std::to_string(_size));
    }
    return node_ids[i];
}

int NTPartition::get_partition_int(size_t i) const {
    if (i >= _size) {
        throw std::range_error(std::string("Error! Index ")
                               + std::to_string(i)
                               + " outside the range of partition of size "
                               + std::to_string(_size));
    }
    return partition_ints[i];
}

std::vector<int> NTPartition::get_cell_list() const {
    std::vector<int> cell_list = std::vector<int>(_size, 0);
    int cell = 0;
    int node;
    for (size_t i = 0; i < _size; i++) {
        node = node_ids[i];
        cell_list[node] = cell;

        if (partition_ints[i] == 0) {
            cell++;
        }
    }
    return cell_list;
}
