#include<cstddef>
#include<vector>

#ifndef SCHENO__NT_PARTITION_H
#define SCHENO__NT_PARTITION_H

class NTPartition {
public:
    // Takes colors (i.e. cells) for nodes 0 through n - 1.
    //  The colors list can have any values.
    NTPartition(const std::vector<int>& cell_list);

    // Puts all nodes in a single partition. O(s).
    NTPartition(const size_t s);

    NTPartition(const NTPartition& ntp);

    // Unless the compiler performs some sort of fancy trick, this copies all
    //  the data from `ntp`, which may be inefficient.
    //
    // TODO: Define move assignment operator or verify that it's implicitly
    //  defined already. Could add a bool `moved_or_deleted` (??)
    NTPartition& operator=(const NTPartition& ntp);

    ~NTPartition();

    // This is perhaps a deceptive use of `const`, in the sense that someone
    //  with access to this pointer can change the data in the class.
    //
    // However, this function itself does not change the class.
    int* get_node_ids() const;
    // Same as above.
    int* get_partition_ints() const;

    size_t size() const;

    int get_node_id(size_t i) const;
    int get_partition_int(size_t i) const;

    // Returns a list where the i'th entry is the color of node i.
    // O(n)
    virtual std::vector<int> get_cell_list() const;

protected:
    int* node_ids;
    int* partition_ints;
    size_t _size;
    // A partition int at location i should be 0 iff i is the last listed node
    //  in the cell.

    // Example partition for {0, 2, 5}, {1}, {3, 4}
    // node_ids:       2, 5, 0, 1, 3, 4
    // partition_ints: 1, 1, 0, 0, 1, 0
};

#endif
