#include "edge.h"
#include "nauty27r4/nauty.h"
#include "nauty27r4/nausparse.h"
#include "nt_sparse_graph.h"
#include "traces.h"
#include "nauty27r4/traces.h"

#include<iostream> // TODO: Remove this and the cout statements.
#include<set>
#include<unordered_map>
#include<vector>

TracesOptions default_traces_options() {
    TracesOptions to;
    // NOTE: These comments are copied from the Nauty/Traces user guide v27,
    //  pages 21 and 22. The full pdf can be found in this repository at
    //  nauty27r4/nug27.pdf
    //
    // Any modifications to the comments are in square brackets [].

    // If this is TRUE, the canonically labelled graph is produced as well as
    //  the automorphism group. Otherwise, only the automorphism group is
    //  determined.
    to.getcanon = false;
    // If this is TRUE, generators of the automorphism group will be written to
    //  the file outfile (see below). The format will depend on the settings of
    //  options cartesian and linelength (see below, again).
    to.writeautoms = false;
    // If writeautoms = TRUE, the value of this option effects the format
    //  in which automorphisms are written. If cartesian = FALSE, the output is
    //  the usual cyclic representation of y, for example "(2 5 6)(3 4)". If
    //  cartesian = TRUE, the output for an automorphism L is the sequence of
    //  numbers "1y 2y . . . (n−1)y", for example "1 5 4 3 6 2".
    to.cartesian = false;
    // Unused, must be FALSE. This release of Traces cannot handle di-graphs.
    //
    // [NOTE: In *this* repository, Traces *CAN* handle digraphs
    //  due to the use of nodes hidden from the user that are given colors.
    //  Regardless, this flag must always be false.]
    to.digraph = false;
    // If this is TRUE, it is assumed that all vertices of the graph have the
    //  same colour (so the initial values of the parameters lab and ptn are
    //  ignored). If it is FALSE, the initial colouring of the vertices is
    //  determined by lab and ptn as described above [in the manual].
    to.defaultptn = false;
    // The value of this variable specifies the maximum number of characters
    //  per line (excluding end-of-line characters) which may be written to the
    //  file outfile (see below). Default 0.
    to.linelength = 0;
    // This is the file to which the output selected by the options writeautoms
    //  and verbosity is sent. It must be already open and writable. The null
    //  pointer NULL is equivalent to stdout (the standard output).
    //  Default NULL.
    to.outfile = NULL;
    // Must be 0 in this version.
    to.strategy = 0;
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
    to.generators = NULL;
    // This is a pointer to a user-defined procedure which is to be called for
    //  each generator. Section 9 has details. No calls will be made if the
    //  value is NULL. Default NULL.
    to.userautomproc = NULL;

    return to;
}

SYMTracesResults traces(NTSparseGraph& g, const SYMTracesOptions& o) {
    std::cout<<"Y"<<std::endl;

    SYMTracesResults results;

    std::cout<<"Z"<<std::endl;

    sparsegraph g_traces = g.as_nauty_traces_graph();

    std::cout<<"A"<<std::endl;

    int* orbits = new int[g_traces.nv];
    NTPartition partition = g.nauty_traces_coloring();

    std::cout<<"B"<<std::endl;

    TracesOptions to = default_traces_options();
    to.getcanon = o.get_canonical_node_order;

    TracesStats ts;

    std::cout<<"C"<<std::endl;

    Traces(&g_traces, partition.get_node_ids(), partition.get_partition_ints(),
           orbits, &to, &ts, NULL);

    std::cout<<"D"<<std::endl;

    results.error_status = ts.errstatus;
    results.num_aut_base = ts.grpsize1;
    results.num_aut_exponent = ts.grpsize2;
    results.num_node_orbits = 0;
    results.num_edge_orbits = 0;

    int largest_orbit_id = 0;

    if (o.get_node_orbits) {
        results.node_orbits = std::vector<int>(orbits, orbits + g.num_nodes());
        std::unordered_set<int> orbit_ids = std::unordered_set<int>();
        for (size_t i = 0; i < g.num_nodes(); i++) {
            orbit_ids.insert(orbits[i]);
        }
        results.num_node_orbits = orbit_ids.size();

        // We will need the largest orbit id to label the orbits of self-loops.
        if (o.get_edge_orbits) {
            for (auto id = orbit_ids.begin(); id != orbit_ids.end(); id++) {
                if (*id > largest_orbit_id) {
                    largest_orbit_id = *id;
                }
            }
        }
    }
    if (o.get_edge_orbits) {
        results.edge_orbits = std::unordered_map<Edge, int, EdgeHash>();
        std::unordered_set<int> orbit_ids = std::unordered_set<int>();
        if (g.directed) {
            // Some edge nodes do not actually correspond to edges. So we look
            //  at the edges themselves and then look up the edge node ids.
            int edge_node;
            for (size_t a = 0; a < g.num_nodes(); a++) {
                for (auto b_itr = g.out_neighbors(a).begin();
                            b_itr != g.out_neighbors(a).end(); b_itr++) {
                    edge_node = g.edge_node(a, *b_itr);
                    results.edge_orbits[EDGE(a, *b_itr, true)] =
                                                            orbits[edge_node];
                }
            }
        } else {
            // Undirected. All edge nodes correspond to actual edges.
            for (size_t i = g.num_nodes(); i < size_t(g_traces.nv); i++) {
                // Use details of the NTSparseGraph representation to speed
                //  things up.
                size_t loc = g_traces.v[i];
                results.edge_orbits[EDGE(g_traces.e[loc],
                                         g_traces.e[loc+1], false)] = orbits[i];
                orbit_ids.insert(orbits[i]);
            }
        }

        // Add self-loops.

        // We will need the largest orbit id to label the orbits of self-loops.
        for (auto id = orbit_ids.begin(); id != orbit_ids.end(); id++) {
            if (*id > largest_orbit_id) {
                largest_orbit_id = *id;
            }
        }
        largest_orbit_id++; // Now strictly larger than any non-self-loop ids.

        for (size_t i = 0; i < g.num_nodes(); i++) {
            if (g.has_edge(i, i)) {
                results.edge_orbits[EDGE(i, i, g.directed)] =
                                                  largest_orbit_id + orbits[i];
                orbit_ids.insert(largest_orbit_id + orbits[i]);
            }
        }
        results.num_edge_orbits = orbit_ids.size();
    }
    delete orbits;

    if (o.get_canonical_node_order) {
        // TODO: Verify whether the regular nodes are always listed first.
        results.canonical_node_order =
            std::vector<int>(partition.get_node_ids(),
                             partition.get_node_ids() + g.num_nodes());
    }

    return results;
}
