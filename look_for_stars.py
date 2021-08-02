from graph_sequencer_utility import GraphSequence

def find_stars(temporal_directed_edges):

    in_neighbors_by_timestamp = {}
    out_neighbors_by_timestamp = {}
    for (a, b, t) in temporal_edges:
        if t not in in_neighbors_by_timestamp:
            in_neighbors_by_timestamp[t] = {}
            out_neighbors_by_timestamp[t] = {}
        i_n_t = in_neighbors_by_timestamp[t]
        o_n_t = out_neighbors_by_timestamp[t]
        if c not in o_n_t:
            o_n_t[c] = set()
        if d not in o_n_t:
            o_n_t[d] = set()
        o_n_t[c].add(d)
        if d not in i_n_t:
            i_n_t[d] = set()
        if c not in i_n_t:
            i_n_t[c] = set()
        i_n_t[d].add(c)

    stars = []
    partial_stars = []
    for t, i_n_t in in_neighbors_by_timestamp.items():
        o_n_t = out_neighbors_by_timestamp[t]

        for node, in_neighbors in i_n_t.items():
            num_single_neighbors = 0
            for neighbor in in_neighbors:
                if len(o_n_t[neighbor]) + len(i_n_t[neighbor]) == 1:
                    num_single_neighbors += 1
            if num_single_neighbors == len(in_neighbors):
                stars.append((num_single_neighbors, node))
            elif num_single_neighbors >= 2:
                partial_stars.append((num_single_neighbors, node))

        for node, out_neighbors in o_n_t.items():
            num_single_neighbors = 0
            for neighbor in out_neighbors:
                if len(o_n_t[neighbor]) + len(i_n_t[neighbor]) == 1:
                    num_single_neighbors += 1
            if num_single_neighbors == len(in_neighbors):
                stars.append((num_single_neighbors, node))
            elif num_single_neighbors >= 2:
                partial_stars.append((num_single_neighbors, node))

    stars.sort(reverse=True)
    partial_stars.sort(reverse=True)
    return (stars, partial_stars)
