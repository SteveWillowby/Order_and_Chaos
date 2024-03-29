from ram_friendly_NT_session import RAMFriendlyNTSession

import math

# Assumes no self-loops
def idx_to_edge(i, num_nodes, directed):
    if directed:
        if i >= (num_nodes) * (num_nodes - 1):
            raise Exception(\
"Idx %d too large for idx_to_edge(%d, %d, %s)" % (i, i, num_nodes, directed))
        a = int(i / (num_nodes - 1))
        b = int(i % (num_nodes - 1))
        if b >= a:
            b += 1
        return (a, b)

    a = 0
    x = (num_nodes - 1)
    while x <= i:
        a += 1
        if (a + 1 == num_nodes):
            raise Exception(\
"Idx %d too large for idx_to_edge(%d, %d, %s)" % (i, i, num_nodes, directed))
        x += (num_nodes - (a + 1))
    x -= (num_nodes - (a + 1))
    return (a, a + (i - x) + 1)

def number_to_edges(num, num_nodes, directed):
    idx = 0
    edge_ints = []
    while (num >= 1):
        if (num & 0x1):
            edge_ints.append(idx)
        idx += 1
        num >>= 1

    edges = []
    for i in edge_ints:
        edges.append(idx_to_edge(i, num_nodes, directed))

    return edges

def edges_to_neighbors_collections(edges, num_nodes, directed):
    nc = [set() for _ in range(0, num_nodes)]

    for (a, b) in edges:
        nc[a].add(b)
        if (not directed):
            nc[b].add(a)

    return nc

def node_order_to_map(node_order):
    the_map = [0 for _ in range(0, len(node_order))]
    for i in range(0, len(node_order)):
        the_map[node_order[i]] = i
    return the_map

def canonical_form(edges, node_to_new_node, directed):
    if (directed):
        return \
            tuple(sorted([(node_to_new_node[a], \
                           node_to_new_node[b]) for (a, b) in edges]))

    new_edges = []
    for (a, b) in edges:
        new_a = node_to_new_node[a]
        new_b = node_to_new_node[b]
        if (new_a <= new_b):
            new_edges.append((new_a, new_b))
        else:
            new_edges.append((new_b, new_a))
    return tuple(sorted(new_edges))

if __name__ == "__main__":

    # N                       1  2  3   4    5     6       7     8      9
    #
    # Undirected Num Graphs:  1, 2, 4,  11,  34,   156,    1044, 12346, 274668
    #   Directed Num Graphs:  1, 3, 16, 218, 9608, 1540944

    DIRECTED = 1
    MAX_N = [8, 5][DIRECTED]

    for N in range(2, MAX_N + 1):
        MAX_E = int((N * (N - 1)) / (2 - DIRECTED))

        graphs = set()
        num_auts = []

        print("N = %d" % N)

        for i in range(0, 1 << MAX_E):
            edges = number_to_edges(i, N, DIRECTED)
            nc = edges_to_neighbors_collections(edges, N, DIRECTED)
            session = RAMFriendlyNTSession(directed=bool(DIRECTED), \
                                           has_edge_types=False, \
                                           neighbors_collections=nc)
            session.blank_coloring()
            node_order = session.get_canonical_order()
            num_aut = session.get_num_automorphisms()
            session.run()
            node_order = node_order.get()
            num_aut = num_aut.get()

            node_map = node_order_to_map(node_order)

            cf = canonical_form(edges, node_map, bool(DIRECTED))

            if (cf not in graphs):
                graphs.add(cf)
                num_auts.append(int(num_aut))

        n_fact = 1
        for i in range(1, N + 1):
            n_fact *= i

        rigid_matrices = (1 << MAX_E) / n_fact

        num_graphs = len(graphs)
        num_aut = sum(num_auts)

        print(rigid_matrices)
        print(num_graphs)
        print(num_aut)

        fancy_factor = num_graphs / rigid_matrices
        needed_factor = num_aut / rigid_matrices
        print("%f vs %f" % (fancy_factor, needed_factor))
        print("Needed power: %f" % (math.log2(needed_factor) / \
                                    math.log2(fancy_factor)))

        print("")
