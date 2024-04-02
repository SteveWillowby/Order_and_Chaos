#include "nauty27r4_modified/nauty.h"
#include "nauty27r4_modified/nausparse.h"
#include "nauty27r4_modified/traces.h"

#include "edge.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"

#include<algorithm>
#include<cmath>
#include<set>
#include<unordered_map>
#include<unordered_set>
#include<vector>

TracesOptions default_traces_options() {
    // This line declares a TracesOptions object `to` with default settings.
    //  Apparently it sets something not mentioned in the manual.
    DEFAULTOPTIONS_TRACES(to);

    // NOTE: These comments are copied from the Nauty/Traces user guide v27,
    //  pages 21 and 22. The full pdf can be found in this repository at
    //  nauty27r4_modified/nug27.pdf
    //
    // Any modifications to the comments are in square brackets [].

    // If this is TRUE, the canonically labelled graph is produced as well as
    //  the automorphism group. Otherwise, only the automorphism group is
    //  determined.
    to.getcanon = false;

    // If this is TRUE, generators of the automorphism group will be written to
    //  the file outfile (see below). The format will depend on the settings of
    //  options cartesian and linelength (see below, again).
    // to.writeautoms = false;

    // If writeautoms = TRUE, the value of this option effects the format
    //  in which automorphisms are written. If cartesian = FALSE, the output is
    //  the usual cyclic representation of y, for example "(2 5 6)(3 4)". If
    //  cartesian = TRUE, the output for an automorphism L is the sequence of
    //  numbers "1y 2y . . . (n−1)y", for example "1 5 4 3 6 2".
    //to.cartesian = false;

    // Unused, must be FALSE. This release of Traces cannot handle di-graphs.
    //
    // [NOTE: In *this* repository, Traces *CAN* handle digraphs
    //  due to the use of nodes hidden from the user that are given colors.
    //  Regardless, this flag MUST always be false.]
    // to.digraph = false;

    // If this is TRUE, it is assumed that all vertices of the graph have the
    //  same colour (so the initial values of the parameters lab and ptn are
    //  ignored). If it is FALSE, the initial colouring of the vertices is
    //  determined by lab and ptn as described above [in the manual].
    to.defaultptn = false;

    // The value of this variable specifies the maximum number of characters
    //  per line (excluding end-of-line characters) which may be written to the
    //  file outfile (see below). Default 0.
    // to.linelength = 0;

    // This is the file to which the output selected by the options writeautoms
    //  and verbosity is sent. It must be already open and writable. The null
    //  pointer NULL is equivalent to stdout (the standard output).
    //  Default NULL.
    // to.outfile = NULL;

    // Must be 0 in this version.
    // to.strategy = 0;

    // A level of verbosity of messages while Traces is running. A value of
    //  0 means that no output will be written (except that automorphisms are
    //  written if the writeautoms option requests them). Larger values produce
    //  greater information about the execution, though its interpretation
    //  requires some knowledge of the algorithm. Default 0.
    to.verbosity = 0;

    // This can be used to provide known automorphisms to Traces and receive the
    //  automorphisms from Traces when it is finished. If it is NULL when Traces
    //  is called, Traces does not change it. If it is non-NULL, it is expected
    //  to point to a (perhaps empty) circular list of known automorphisms. (It
    //  is an error to give a permutation that is not an automorphism of the
    //  input coloured graph.) In this case, Traces will add automorphisms to
    //  the list so that the whole automorphism group is generated. See Section
    //  18 for detailed instructions and examples. Default NULL.
    // to.generators = NULL;

    // This is a pointer to a user-defined procedure which is to be called for
    //  each generator. Section 9 has details. No calls will be made if the
    //  value is NULL. Default NULL.
    // to.userautomproc = NULL;

    return to;
}

optionblk default_nauty_options() {
    DEFAULTOPTIONS_SPARSEGRAPH(no);

    no.getcanon = false;
    no.defaultptn = false;
    no.digraph = false;

    return no;
}

// 0 --> traces
// 1 --> nauty
// 2 --> fake_iso
NautyTracesResults __iso_program(int program, NTSparseGraph& g,
                                 const NautyTracesOptions& o,
                                 NTPartition& p);

NautyTracesResults nauty(NTSparseGraph& g, const NautyTracesOptions& o) {
    NTPartition partition = g.nauty_traces_coloring();
    return __iso_program(1, g, o, partition);
}

NautyTracesResults nauty(NTSparseGraph& g, const NautyTracesOptions& o,
                         NTPartition& p) {
    return __iso_program(1, g, o, p);
}

NautyTracesResults traces(NTSparseGraph& g, const NautyTracesOptions& o) {
    NTPartition partition = g.nauty_traces_coloring();
    return __iso_program(0, g, o, partition);
}

NautyTracesResults traces(NTSparseGraph& g, const NautyTracesOptions& o,
                          NTPartition& p) {
    return __iso_program(0, g, o, p);
}

NautyTracesResults fake_iso(NTSparseGraph& g, const NautyTracesOptions& o) {
    NTPartition partition = g.nauty_traces_coloring();
    return __iso_program(2, g, o, partition);
}

NautyTracesResults fake_iso(NTSparseGraph& g, const NautyTracesOptions& o,
                            NTPartition& p) {
    return __iso_program(2, g, o, p);
}

void __fake_iso(sparsegraph* g, int* partition_node_ids, int* partition_ints,
                int *orbits, long double* aut_base, long double* aut_exp);

// 0 --> traces
// 1 --> nauty
// 2 --> fake_iso
NautyTracesResults __iso_program(int program, NTSparseGraph& g,
                                 const NautyTracesOptions& o,
                                 NTPartition& p) {

    NautyTracesResults results;

    sparsegraph g_nt = g.as_nauty_traces_graph();

    sparsegraph* canon_rep = NULL;

    __nt_run_space.set_size(g_nt.nv, g_nt.nde);

    int* orbits = __nt_run_space.orbits;

    if (program == 0) {  // Traces
        TracesOptions to = default_traces_options();
        if (o.get_canonical_node_order) {
            to.getcanon = true;
            canon_rep = &__nt_run_space.g;
        }

        TracesStats ts;

        Traces(&g_nt, p.get_node_ids(), p.get_partition_ints(),
               orbits, &to, &ts, canon_rep);

        results.error_status = ts.errstatus;
        results.num_aut_base = ts.grpsize1;
        results.num_aut_exponent = ts.grpsize2;
    } else if (program == 1) {  // Nauty
        optionblk no = default_nauty_options();
        if (o.get_canonical_node_order) {
            no.getcanon = true;
            canon_rep = &__nt_run_space.g;
        }

        statsblk ns;

        sparsenauty(&g_nt, p.get_node_ids(), p.get_partition_ints(),
                    orbits, &no, &ns, canon_rep);

        results.error_status = ns.errstatus;
        results.num_aut_base = ns.grpsize1;
        results.num_aut_exponent = ns.grpsize2;
    } else {  // Fake ISO
        long double aut_base;
        long double aut_exp;

        __fake_iso(&g_nt, p.get_node_ids(), p.get_partition_ints(),
                   orbits, &aut_base, &aut_exp);

        results.error_status = 0;
        results.num_aut_base = aut_base;
        results.num_aut_exponent = aut_exp;
    }

    results.num_node_orbits = 0;
    results.num_edge_orbits = 0;

    int largest_orbit_id = 0;

    if (o.get_node_orbits) {
        results.node_orbits = Coloring<int>();
        for (size_t i = 0; i < g.num_nodes(); i++) {
            results.node_orbits.set(i, orbits[i]);
        }
        results.num_node_orbits = results.node_orbits.colors().size();
        if (o.get_edge_orbits) {
            largest_orbit_id = *(results.node_orbits.colors().rbegin());
        }
    }
    if (o.get_edge_orbits) {
        results.edge_orbits = Coloring<Edge, EdgeHash>();

        if (g.directed) {
            // Some edge nodes do not actually correspond to edges. So we look
            //  at the edges themselves and then look up the edge node ids.
            int edge_node;
            for (size_t a = 0; a < g.num_nodes(); a++) {
                for (auto b_itr = g.out_neighbors(a).begin();
                            b_itr != g.out_neighbors(a).end(); b_itr++) {
                    if (int(a) == *b_itr) {
                        continue;
                    }
                    edge_node = g.edge_node(a, *b_itr);
                    results.edge_orbits.set(EDGE(a, *b_itr, true),
                                            orbits[edge_node]);
                }
            }
        } else {
            // Undirected. All edge nodes correspond to actual edges.
            for (size_t i = g.num_nodes(); i < size_t(g_nt.nv); i++) {
                // Use details of the NTSparseGraph representation to speed
                //  things up.
                size_t loc = g_nt.v[i];
                results.edge_orbits.set(EDGE(g_nt.e[loc], g_nt.e[loc+1], false),
                                        orbits[i]);
            }
        }

        if (results.edge_orbits.size() > 0) {
            int c = *(results.edge_orbits.colors().rbegin());
            largest_orbit_id = (c > largest_orbit_id ? c : largest_orbit_id);
        }
        largest_orbit_id++; // Now strictly larger than any non-self-loop ids.

        for (size_t i = 0; i < g.num_nodes(); i++) {
            if (g.has_edge(i, i)) {
                results.edge_orbits.set(EDGE(i, i, g.directed),
                                        largest_orbit_id + orbits[i]);
            }
        }
        results.num_edge_orbits = results.edge_orbits.colors().size();

    }

    if (o.get_canonical_node_order) {
        // TODO: Verify whether the regular nodes are always listed first.
        //  (as opposed to the augmented edge nodes)
        results.canonical_node_order = std::vector<int>(g.num_nodes(), 0);
        int *canon = p.get_node_ids();
        for (size_t i = 0; i < g.num_nodes(); i++) {
            results.canonical_node_order[i] = canon[i];
        }
    }

    return results;
}

__NTRunSpace::__NTRunSpace() {
    SG_INIT(g);
    orbits = NULL;
    actual_n = 0;
    actual_elen = 0;
    g.wlen = 0;
    g.w = NULL;
}

__NTRunSpace::~__NTRunSpace() {
    if (actual_n > 0) {
        delete g.d;
        delete g.v;
        delete g.e;
        delete orbits;
    }
}

void __NTRunSpace::set_size(size_t n, size_t nde) {
    if (n > actual_n) {
        actual_n = (n * 3) / 2;
        if (actual_n > 0) {
            delete orbits;
            delete g.d;
            delete g.v;
        }
        orbits = new int[actual_n];
        g.d = new int[actual_n];
        g.v = new size_t[actual_n];

        g.dlen = actual_n;
        g.vlen = actual_n;
    }
    if (nde > actual_elen || nde == 0) {
        actual_elen = (nde * 3) / 2;
        if (actual_elen == 0) {
            actual_elen = 1;
        } else {
            delete g.e;
        }
        g.e = new int[actual_elen];
        g.elen = actual_elen;
    }
    g.nde = nde;
}

void __fake_iso(sparsegraph* g, int* partition_node_ids, int* partition_ints,
                int *orbits, long double* aut_base, long double* aut_exp) {

    long double log2_aut = 0.0;
    long naive_aut = 1;

    // Initialize the Data Structures

    int n = g->nv;
    int* degrees = g->d;
    size_t* nbr_start = g->v;
    int* neighbors = g->e;

    std::vector<int> node_to_label = std::vector<int>(n);
    std::unordered_map<int, std::unordered_set<int>> label_to_nodes;

    int node, neighbor;
    int label = 0;
    label_to_nodes.insert(
        std::pair<int, std::unordered_set<int>>(0, std::unordered_set<int>()));

    for (int idx = 0; idx < n; idx++) {
        node = partition_node_ids[idx];
        node_to_label[node] = label;
        label_to_nodes[label].insert(node);

        if (partition_ints[idx] == 0) {
            if (idx < n - 1) {
                label++;
                label_to_nodes.insert(std::pair<int, std::unordered_set<int>>(
                                        label, std::unordered_set<int>()));
            }
        }
    }

    std::set<int> affected_labels = std::set<int>();
    for (auto itr = label_to_nodes.begin(); itr != label_to_nodes.end(); itr++){
        if (itr->second.size() > 1) {
            affected_labels.insert(itr->first);
        }
    }

    std::vector<std::pair<std::vector<int>, int>> cell =
                    std::vector<std::pair<std::vector<int>, int>>();

    // Begin the WL Refinement Process
    int round = 0;
    int next_label = label_to_nodes.size();
    while (true) {  // Continue until all nodes have their own cell/label

        while (affected_labels.size() > 0) {  // Run WL

            // Pick the cell with the smallest ID and see if it splits.
            label = *(affected_labels.begin());

            // Fill Cell
            cell.clear();
            for (auto cell_itr = label_to_nodes[label].begin();
                      cell_itr != label_to_nodes[label].end(); cell_itr++) {
                node = *cell_itr;


                // Collect node's neighbor labels
                cell.push_back({std::vector<int>(), node});
                for (int k = nbr_start[node];
                         k < nbr_start[node] + degrees[node]; k++) {
                    neighbor = neighbors[k];
                    cell[cell.size() - 1].first.push_back(node_to_label[neighbor]);
                }

                // Sort node's neighbor labels
                std::sort(cell[cell.size() - 1].first.begin(),
                          cell[cell.size() - 1].first.end());

                // Append node's old label
                cell[cell.size() - 1].first.push_back(node_to_label[node]);
            }


            // Sort Cell
            std::sort(cell.begin(), cell.end());


            bool cell_split = false;
            // Relabel Cell Nodes
            for (size_t j = 0; j < cell.size(); j++) {
                node = cell[j].second;
                label_to_nodes[node_to_label[node]].erase(node);
                if (label_to_nodes[node_to_label[node]].size() == 0) {
                    label_to_nodes.erase(node_to_label[node]);
                }
                if (label_to_nodes.find(next_label) == label_to_nodes.end()) {
                    label_to_nodes.insert({next_label,
                                           std::unordered_set<int>()});
                }
                node_to_label[node] = next_label;
                label_to_nodes[next_label].insert(node);

                if (j < cell.size() - 1 && cell[j].first <
                                           cell[j + 1].first) {
                    next_label++;
                    cell_split = true;
                }
            }
            next_label++;


            affected_labels.erase(label);

            if (cell_split) {
                // Add affected cells (i.e. labels).
                for (size_t j = 0; j < cell.size(); j++) {
                    node = cell[j].second;
                    for (int k = nbr_start[node];
                             k < nbr_start[node] + degrees[node]; k++) {
                        neighbor = neighbors[k];
                        affected_labels.insert(node_to_label[neighbor]);
                    }
                }
            }
        }


        round++;
        if (round == 1) {
            // TODO: Consider zero-indexing these outputs.
            for (int node = 0; node < n; node++) {
                orbits[node] = node_to_label[node];
            }
        }

        if (label_to_nodes.size() == (size_t) n) {
            // Every node has its own cell/label
            break;
        }

        // Select a single node to have its own color.
        std::unordered_set<int> *cell_ptr;
        for (auto cell_itr = label_to_nodes.begin();
                  cell_itr != label_to_nodes.end(); cell_itr++) {
            if (cell_itr->second.size() == 1) {
                continue;
            }
            cell_ptr = &(cell_itr->second);
            break;
        }
        node = *(cell_ptr->begin());

        // Create new spot for node, but don't erase it from the old spot yet
        label_to_nodes.insert({next_label, std::unordered_set<int>()});
        label_to_nodes[next_label].insert(node);
        node_to_label[node] = next_label;
        next_label++;

        // Update automorphism estimate.
        log2_aut += std::log2l(cell_ptr->size());
        naive_aut *= cell_ptr->size();

        // Calculate Affected Cells
        for (auto node_itr = cell_ptr->begin();
                  node_itr != cell_ptr->end(); node_itr++) {
            for (int k = nbr_start[*node_itr];
                     k < nbr_start[*node_itr] + degrees[*node_itr]; k++) {
                neighbor = neighbors[k];
                affected_labels.insert(node_to_label[neighbor]);
            }
        }

        // Ok, now erase it.
        cell_ptr->erase(node);
    }

    // At the very end, convert the labeling into an NT partition
    std::vector<int> labels = std::vector<int>();
    for (auto lbl_itr = label_to_nodes.begin();
              lbl_itr != label_to_nodes.end(); lbl_itr++) {
        labels.push_back(lbl_itr->first);
    }
    std::sort(labels.begin(), labels.end());

    int idx = 0;
    for (auto lbl_itr = labels.begin(); lbl_itr != labels.end(); lbl_itr++) {
        label = *lbl_itr;
        for (auto node_itr = label_to_nodes[label].begin();
                  node_itr != label_to_nodes[label].end(); node_itr++) {
            node = *node_itr;
            partition_node_ids[idx] = node;
            partition_ints[idx] = 1;
            idx++;
        }
        partition_ints[idx - 1] = 0;
    }

    long double log10_aut = log2_aut / std::log2l(10);
    *aut_exp = std::floor(log10_aut);
    long double what_remains = log10_aut - *aut_exp;
    *aut_base = std::exp2l(what_remains * std::log2l(10));
}
