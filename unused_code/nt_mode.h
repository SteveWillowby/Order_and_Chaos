// Used to determine format of the Nauty/Traces input.

#ifndef NT_MODE_H
#define NT_MODE_H

class NTMode {

public:
    // default options: pt = true, nd = true
    NTMode();
    // custom options
    NTMode(const bool prefer_traces, const bool nauty_directed);

    const bool prefer_traces;
    const bool nauty_directed;

    /* These options are meant to be used in the following way:
     *
     * n = num nodes
     * m = num edges
     * m_undir = num edges in undirected version of graph
     *
     *
     * If the graph is undirected:
     *     If `prefer_traces`:
     *         Use Traces. Num nodes = n + int(has_edge_types) * m
     *     Else:
     *         Use Nauty. Num nodes = n + int(has_edge_types) * m
     *
     * Else:
     *     If `prefer_traces`:
     *         Use Traces. Num nodes = n + 2 * m_undir
     *     Else:
     *         If `nauty_directed`:
     *             Use Nauty. Num nodes = n + int(has_edge_types) * m
     *         Else:
     *             Use Nauty. Num nodes = n + 2 * m_undir
     */
};

#endif
