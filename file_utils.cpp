#include "file_utils.h"
#include "graph.h"
#include "sparse_graph.h"

#include<fstream>
#include<set>
#include<stdexcept>
#include<string>
#include<unordered_map>
#include<vector>

std::vector<std::pair<int, int>> _read_edgelist(const std::string& filename);
SparseGraph _construct_graph(const bool directed,
                             const std::vector<std::pair<int, int>>& edgelist,
                             const std::set<int>& node_labels);

SparseGraph read_graph(const bool directed,
                       const std::string& edgelist_filename) {

    std::vector<std::pair<int, int>> el = _read_edgelist(edgelist_filename);
    if (el.size() == 0) {
        throw std::invalid_argument(
                std::string("Error! Cannot load an edgelist with no edges that")
                + " is not accompanied by a nodelist.");
    }

    int max_node = 0;
    for (auto pair_itr = el.begin(); pair_itr != el.end(); pair_itr++) {
        if (pair_itr->first > max_node) {
            max_node = pair_itr->first;
        }
        if (pair_itr->second > max_node) {
            max_node = pair_itr->second;
        }
    }
    std::vector<int> node_labels_list = std::vector<int>(max_node + 1, 0);
    for (int i = 0; i <= max_node; i++) {
         node_labels_list[i] = i;
    }
    // Hoping that this construction method is faster.
    std::set<int> node_labels = std::set<int>(node_labels_list.begin(),
                                              node_labels_list.end());

    return _construct_graph(directed, el, node_labels);
}

SparseGraph read_graph(const bool directed,
                       const std::string& nodelist_filename,
                       const std::string& edgelist_filename) {

    std::ifstream file = std::ifstream(nodelist_filename);
    if (!file) {
        throw std::invalid_argument("Error! Could not open file "
                                    + nodelist_filename);
    }

    std::string line;
    int node;
    std::set<int> nodes = std::set<int>();

    while (!file.bad() && !file.eof()) {
        std::getline(file, line);
        if (line.empty()) {
            break;
        }
        node = std::stoi(line);
        if (node < 0) {
            throw std::invalid_argument("Error! Negative node ID in edgelist "
                                        + nodelist_filename);
        }
        nodes.insert(node);
    }
    if (nodes.size() == 0) {
        throw std::invalid_argument("Error! No nodes in nodelist "
                                    + nodelist_filename);
    }

    file.close();

    return _construct_graph(directed, _read_edgelist(edgelist_filename), nodes);
}

void write_graph(const Graph& g, const std::string& nodelist_filename,
                                 const std::string& edgelist_filename) {
    if (!nodelist_filename.empty()) {
        // We have a nodelist filename
        std::ofstream nodelist(nodelist_filename);
        if (!nodelist) {
            throw std::invalid_argument("Error! Could not open file "
                                        + nodelist_filename);
        }
        for (int i = 0; i < int(g.num_nodes()); i++) {
            nodelist << i;
            if (i < int(g.num_nodes() - 1)) {
                nodelist << std::endl;
            }
        }
        nodelist.close();
    }

    std::ofstream edgelist(edgelist_filename);
    if (!edgelist) {
        throw std::invalid_argument("Error! Could not open file "
                                    + edgelist_filename);
    }
    if (g.directed) {
        size_t m = 0;
        for (int a = 0; a < int(g.num_nodes()); a++) {
            for (auto b_itr = g.out_neighbors(a).begin();
                      b_itr != g.out_neighbors(a).end(); b_itr++) {
                edgelist << a << " " << *b_itr;
                m++;
                if (m < g.num_edges()) {
                    edgelist << std::endl;
                }
            }
        }
    } else {
        // Undirected
        size_t m = 0;
        for (int a = 0; a < int(g.num_nodes()); a++) {
            for (auto b_itr = g.neighbors(a).begin();
                      b_itr != g.neighbors(a).end(); b_itr++) {
                if (a > *b_itr) {
                    continue;
                }
                edgelist << a << " " << *b_itr;
                m++;
                if (m < g.num_edges()) {
                    edgelist << std::endl;
                }
            }
        }
    }

    edgelist.close();
}

SparseGraph _construct_graph(const bool directed,
                             const std::vector<std::pair<int, int>>& edgelist,
                             const std::set<int>& node_labels) {

    SparseGraph g = SparseGraph(directed, node_labels.size());

    int last_node = node_labels.size() - 1;
    bool relabel = *(node_labels.begin()) != 0 ||
                   *(node_labels.rbegin()) != last_node;

    if (relabel) {
        std::unordered_map<int, int> relabeling = std::unordered_map<int, int>();
        int i = 0;
        for (auto itr = node_labels.begin();
                  itr != node_labels.end(); itr++) {
            relabeling[*itr] = i;
            i++;
        }

        for (auto edge_itr = edgelist.begin();
                  edge_itr != edgelist.end(); edge_itr++) {
            g.add_edge(relabeling[edge_itr->first],
                       relabeling[edge_itr->second]);
        }
    } else {
        for (auto edge_itr = edgelist.begin();
                  edge_itr != edgelist.end(); edge_itr++) {
            g.add_edge(edge_itr->first, edge_itr->second);
        }
    }

    return g;
}

std::vector<std::pair<int, int>> _read_edgelist(const std::string& filename) {
    std::ifstream file = std::ifstream(filename);
    if (!file) {
        throw std::invalid_argument("Error! Could not open file " + filename);
    }

    std::string line, s_a, s_b;
    int a, b;
    std::vector<std::pair<int, int>> result = std::vector<std::pair<int,int>>();

    while (!file.bad() && !file.eof()) {
        std::getline(file, line);
        if (line.size() < 3) {
            // A non-edge line.
            continue;
        }

        s_a = std::string(line.begin(), line.begin() + line.find(' '));
        s_b = std::string(line.begin() + line.find(' ') + 1, line.end());
        a = std::stoi(s_a);
        b = std::stoi(s_b);

        if (a < 0 || b < 0) {
            throw std::invalid_argument("Error! Negative node ID in edgelist "
                                        + filename);
        }
        result.push_back(std::pair<int, int>(a, b));
    }

    file.close();

    return result;
}

void make_nodelist(const std::string& edgelist_filename,
                   const std::string& nodelist_filename,
                   bool full_range) {

    std::vector<std::pair<int, int>> edges = _read_edgelist(edgelist_filename);
    std::set<int> nodes = std::set<int>();
    int a, b;
    for (auto itr = edges.begin(); itr != edges.end(); itr++) {
        a = itr->first;
        b = itr->second;
        nodes.insert(a);
        nodes.insert(b);
    }

    std::ofstream file = std::ofstream(nodelist_filename);
    if (full_range) {
        int max_node = 0;
        for (auto itr = nodes.begin(); itr != nodes.end(); itr++) {
            if (*itr > max_node) {
                max_node = *itr;
            }
        }

        for (int i = 0; i <= max_node; i++) {
            file << i;
            if (i < max_node) {
                file << std::endl;
            }
        }
    } else {
        for (auto itr = nodes.begin(); itr != nodes.end(); ) {
            file<<*itr;
            itr++;
            if (itr != nodes.end()) {
                file<<std::endl;
            }
        }
    }
    file.close();
}
