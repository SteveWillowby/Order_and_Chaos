import math
from nauty_session import NautyTracesSession
import networkx as nx

def __get_graph_info__(nodes, given_edges, temporal, directed, \
                       use_color_direction=True, get_canon_order=False):
    if type(nodes) is not set:
        nodes = set(nodes)

    if directed and not use_color_direction:
        graph = nx.DiGraph()
    else:
        graph = nx.Graph()

    for node in nodes:
        graph.add_node(node)

    orbits = [list(nodes)]  # at first, everything in one orbit

    # Use only to store highlights for nodes not in `nodes`
    extra_highlights = []
    if temporal:
        edge_sets_by_timestamp = {}
        for (a, b, t) in given_edges:
            if t not in edge_sets_by_timestamp:
                edge_sets_by_timestamp[t] = set()
            edge_sets_by_timestamp[t].add((a, b))

        timestamps = sorted([t for t, _ in edge_sets_by_timestamp.items()])
        alt_nodes = []

        if directed and use_color_direction:
            dir_nodes_1 = []
            dir_nodes_2 = []
            extra_highlights = [dir_nodes_1, dir_nodes_2]

        for i in range(0, len(timestamps)):
            t = timestamps[i]
            for node in nodes:
                alt_nodes.append((node, t))
                graph.add_node((node, t))
                if i == 0:
                    graph.add_edge(node, (node, t))
                else:
                    graph.add_edge((node, timestamps[i - 1]), (node, t))
            for (a, b) in edge_sets_by_timestamp[t]:

                if directed and use_color_direction:
                    graph.add_node(((a, t), (b, t), 1))
                    graph.add_node(((a, t), (b, t), 2))
                    dir_nodes_1.append(((a, t), (b, t), 1))
                    dir_nodes_2.append(((a, t), (b, t), 2))
                    graph.add_edge((a, t), ((a, t), (b, t), 1))
                    graph.add_edge(((a, t), (b, t), 1), ((a, t), (b, t), 2))
                    graph.add_edge(((a, t), (b, t), 2), (b, t))

                    # only add edge once
                    if b not in graph.neighbors(a):
                        graph.add_edge(a, b)
                else:
                    graph.add_edge((a, t), (b, t))
    else:  # static
        if directed and use_color_direction:
            dir_nodes_1 = []
            dir_nodes_2 = []
            extra_highlights = [dir_nodes_1, dir_nodes_2]

        for (a, b) in given_edges:

            if directed and use_color_direction:
                graph.add_node((a, b, 1))
                graph.add_node((a, b, 2))
                dir_nodes_1.append((a, b, 1))
                dir_nodes_2.append((a, b, 2))
                graph.add_edge(a, (a, b, 1))
                graph.add_edge((a, b, 1), (a, b, 2))
                graph.add_edge((a, b, 2), b)

                # only add edge once
                if b not in graph.neighbors(a):
                    graph.add_edge(a, b)
            else:
                graph.add_edge(a, b)

    NTSession = NautyTracesSession(graph, mode="Nauty", sparse=True)
    NTSession.set_colors_by_highlights([list(nodes)] + extra_highlights)
    na = NTSession.get_num_automorphisms()
    if get_canon_order:
        co = NTSession.get_canonical_order()
    NTSession.complete()
    na = na.get()
    if get_canon_order:
        co = co.get()
        return (na, co)
    return (na, None)

# information_content()
#
# `only_consider_used_nodes` restricts the possible adjacency matrices to
#   only include nodes which had at least one edge.
#
# `extreme_only_consider_used_nodes` is only relevant for temporal graphs.
#   In addition to excluding nodes that never had an edge, this setting
#   excludes nodes at each timestamp that never had an edge before (or in) that
#   timestamp.
#
# `proportional_p` determines whether the Erdos-Renyi calculation uses p=0.5 or
#   p = |edges| / |possible edges|. The latter setting is used when
#   proportional_p=True.
#
# `extreme_proportional_p` is only relevant for temporal graphs. Rather than
#   using p = |edges| / |possible edges|, it does a separate calculation for
#   each timestamp.
#
# `num_timestamps` is only relevant if the graph is temporal. If used, it can
#   either have an integer or be set to "auto". This enables having temporal
#   snapshots with no edges in them.
#
# NOTE: Assumes that `edges` and `nodes` only contain each edge and node once.
#   This may cause problems if you pass a directed edge list but set
#   `directed` to False.
#
# NOTE: Do not pass self-loops.
def information_content(nodes, edges, proportional_p=False, \
                        directed=False, temporal=False, \
                        num_timestamps="auto", \
                        node_colors=None, only_consider_used_nodes=False, \
                        extreme_only_consider_used_nodes=False, \
                        extreme_proportional_p=False):

    if extreme_proportional_p:
        proportional_p = True

    if extreme_only_consider_used_nodes:
        only_consider_used_nodes = True

    if not temporal:
        extreme_only_consider_used_nodes = False
        extreme_proportional_p = False

    if only_consider_used_nodes:
        if len(edges) == 0:
            raise ValueError("Error! Cannot only consider nodes with edges " + \
                             "when there are no edges. Set " + \
                             "`only_consider_used_nodes` to False.")

    # Collect meta-data.
    if temporal:
        edges = [(t, a, b) for (a, b, t) in edges]
        edges.sort()
        edges = [(a, b, t) for (t, a, b) in edges]

        if only_consider_used_nodes:
            nodes = set()
        if extreme_only_consider_used_nodes:
            num_nodes_by_timestamp = []
        if extreme_proportional_p:
            num_edges_by_timestamp = []
            edges_in_timestamp = 0

        if len(edges) == 0:
            observed_timestamps = 0
        else:
            observed_timestamps = 1

        current_time = None
        for edge in edges:
            source = edge[0]
            target = edge[1]

            time = edge[2]

            if current_time is not None and time != current_time:
                observed_timestamps += 1
                if extreme_only_consider_used_nodes:
                    num_nodes_by_timestamp.append(len(nodes))
                if extreme_proportional_p:
                    num_edges_by_timestamp.append(edges_in_timestamp)
                    edges_in_timestamp = 0
            current_time = time

            if only_consider_used_nodes:
                nodes.add(source)
                nodes.add(target)

            if extreme_proportional_p:
                edges_in_timestamp += 1

        if extreme_only_consider_used_nodes:
            num_nodes_by_timestamp.append(len(nodes))
        if extreme_proportional_p:
            num_edges_by_timestamp.append(edges_in_timestamp)

        if num_timestamps == "auto":
            num_timestamps = observed_timestamps
        else:
            if num_timestamps < observed_timestamps:
                raise ValueError("Error! num_timestamps set to " + \
                                 ("%d, but " % num_timestamps) + \
                                 ("%d, were observed." % observed_timestamps))
            elif num_timestamps > observed_timestamps:
                for i in range(0, num_timestamps - observed_timestamps):
                    if extreme_only_consider_used_nodes:
                        num_nodes_by_timestamp.append(len(nodes))
                    if extreme_proportional_p:
                        num_edges_by_timestamp.append(0)

    else:  # static (not temporal)
        if only_consider_used_nodes:
            nodes = set()
            for edge in edges:
                nodes.add(edge[0])
                nodes.add(edge[1])

    # WE WORK IN LOG_2 SPACE
    #   otherwise, the numbers would get too large.

    # Compute the log-number of possible adjacency matrices.
    if temporal:
        if extreme_only_consider_used_nodes:
            log2_num_matrices = 0
            for n in num_nodes_by_timestamp:
                if directed:
                    log2_num_matrices += n * (n - 1)
                else:
                    log2_num_matrices += (n * (n - 1)) / 2
        else:
            num_nodes = len(nodes)
            if directed:
                log2_num_matrices = \
                    num_nodes * (num_nodes - 1) * num_timestamps
            else:
                log2_num_matrices = \
                    ((num_nodes * (num_nodes - 1)) / 2) * num_timestamps
    else:  # static
        num_nodes = len(nodes)
        if directed:
            log2_num_matrices = num_nodes * (num_nodes - 1)
        else:
            log2_num_matrices = (num_nodes * (num_nodes - 1)) / 2

    # Compute the log-probability of this adjacency matrix.
    if temporal and extreme_only_consider_used_nodes and \
            extreme_proportional_p:
        log2_prob_matrix = 0
        if extreme_proportional_p:
            for i in range(0, len(num_nodes_by_timestamp)):
                n = num_nodes_by_timestamp[i]
                e = num_edges_by_timestamp[i]
                if directed:
                    num_possible_edges = (n * (n - 1))
                else:
                    num_possible_edges = ((n * (n - 1)) / 2)

                if e == num_possible_edges or e == 0:
                    log2_prob_matrix = 0.0
                else:
                    prop = float(e) / num_possible_edges
                    log2_prop = math.log(prop, 2.0)
                    log2_non_prop = math.log(1.0 - prop, 2.0)
                    log2_prob_matrix += log2_prop * e + \
                                        log2_non_prop * (num_possible_edges - e)
    else:
        if proportional_p:
            # the number of matrices is 2^(edge slots)
            num_possible_edges = log2_num_matrices
            e = len(edges)
            if e == num_possible_edges or e == 0:
                log2_prob_matrix = 0.0
            else:
                prop = float(e) / num_possible_edges
                log2_prop = math.log(prop, 2.0)
                log2_non_prop = math.log(1.0 - prop, 2.0)
                log2_prob_matrix = log2_prop * e + \
                                   log2_non_prop * (num_possible_edges - e)
        else:
            # the number of matrices is 2^(edge slots)
            # the prob of a single matrix is 2^-(edge slots)
            log2_prob_matrix = -log2_num_matrices

    # Compute the log-number of automorphisms.
    (num_automorphisms, _) = __get_graph_info__(\
            nodes, edges, temporal, directed, use_color_direction=True, \
            get_canon_order=False)

    log2_num_automorphisms = math.log(num_automorphisms, 2.0)

    # Compute the log-number of nodes factorial:
    log2_n_fact = 0
    for i in range(1, len(nodes) + 1):
        log2_num = math.log(i, 2.0)
        log2_n_fact += log2_num

    # Information content is -log2(prob of graph | Erdos-Renyi)
    #   which equals:
    #       -log2(prob_of_this_matrix * num_matrices_that_make_this_graph)
    #   In turn, that is:
    #       -log2(prob_of_this_matrix * |nodes|! / num_automorphisms)
    return -1.0 * (log2_prob_matrix + log2_n_fact - log2_num_automorphisms)

# Used only for testing.
def __canonical_string__(nodes, edges, directed=False, temporal=False):
    (_, co) = __get_graph_info__(nodes, edges, temporal, directed, \
                       use_color_direction=True, get_canon_order=True)
    node_to_idx_map = {}
    for i in range(0, len(co)):
        node_to_idx_map[co[i]] = i

    nodes = [i for i in range(0, len(nodes))]
    if temporal:
        edges = [(node_to_idx_map[a], node_to_idx_map[b], t) for \
                    (a, b, t) in edges]
        if not directed:
            edges = [(min(a, b), max(a, b), t) for (a, b, t) in edges]
    else:
        edges = [(node_to_idx_map[a], node_to_idx_map[b]) for \
                    (a, b) in edges]
        if not directed:
            edges = [(min(a, b), max(a, b)) for (a, b) in edges]
    edges.sort()

    return str((nodes, edges))

# Used only for testing.
def __temporal_graph_mashup__(graphs, directed=False):
    from algorithmic_utils import get_all_k_permutations
    num_nodes = len(graphs[0][0])
    perms = get_all_k_permutations(num_nodes, num_nodes)
    canonical_strings = set()

    temporal_graphs = []

    for (nodes, first_edges) in graphs:
        for (_, second_edges) in graphs:
            for perm in perms:
                edges = []
                for (a, b) in first_edges:
                    edges.append((a, b, 1))
                for (a, b) in second_edges:
                    edges.append((nodes[perm[a-1]], nodes[perm[b-1]], 2))
                cs = __canonical_string__(nodes, edges, directed=directed,temporal=True)
                if cs not in canonical_strings:
                    canonical_strings.add(cs)
                    temporal_graphs.append((nodes, edges))
    return temporal_graphs

if __name__ == "__main__":
    prob_total = 0.0
    max_name_length = 23
    print("For the FOUR-node, UNDIRECTED graphs, we get:")
    for (nodes, edges, name) in [\
            ([1, 2, 3, 4], [], "Empty"), \
            ([1, 2, 3, 4], [(1, 2)], "Single Edge"), \
            ([1, 2, 3, 4], [(1, 2), (2, 3)], "Wedge"), \
            ([1, 2, 3, 4], [(1, 2), (3, 4)], "Two Disjoint Edges"), \
            ([1, 2, 3, 4], [(1, 2), (2, 3), (3, 1)], "Triangle"), \
            ([1, 2, 3, 4], [(1, 2), (2, 3), (3, 4)], "Three Chain"), \
            ([1, 2, 3, 4], [(1, 2), (1, 3), (1, 4)], "Three Star"), \
            ([1, 2, 3, 4], [(1, 2), (2, 3), (3, 4), (4, 1)], "Rectangle"), \
            ([1, 2, 3, 4], [(1, 2), (2, 3), (3, 4), (1, 3)], "Triangle with Dangler"), \
            ([1, 2, 3, 4], [(1, 2), (2, 3), (3, 4), (1, 3), (4, 1)], "Rectangle with Diagonal"), \
            ([1, 2, 3, 4], [(1, 2), (2, 3), (3, 4), (1, 3), (1, 4), (2, 4)], "Full")]:
        ic = information_content(nodes, edges, proportional_p=False, \
                        directed=False, temporal=False, \
                        num_timestamps="auto", \
                        node_colors=None, only_consider_used_nodes=False, \
                        extreme_only_consider_used_nodes=False, \
                        extreme_proportional_p=False)
        padding_string = ""
        for _ in range(0, max_name_length - len(name)):
            padding_string += " "
        print("\n  Information Content of %s graph: %s%f." % \
            (name, padding_string, ic))
        prob_total += math.pow(2.0, -ic)
    print("The total probabilities sum to ~%f." % prob_total)

    prob_total = 0.0
    max_name_length = 24
    print("\n\n\nFor the THREE-node, DIRECTED graphs, we get:")
    for (nodes, edges, name) in [\
            ([1, 2, 3], [], "Empty"), \
            ([1, 2, 3], [(1, 2)], "Single Edge"), \
            ([1, 2, 3], [(1, 2), (2, 3)], "Transitive Wedge"), \
            ([1, 2, 3], [(1, 2), (2, 1)], "Single Bi-Edge"), \
            ([1, 2, 3], [(1, 3), (2, 3)], "Inward Star"), \
            ([1, 2, 3], [(1, 2), (1, 3)], "Outward Star"), \
            ([1, 2, 3], [(1, 2), (2, 3), (3, 1)], "Triangle Cycle"), \
            ([1, 2, 3], [(1, 2), (2, 3), (1, 3)], "Triangle W/ Back-Edge"), \
            ([1, 2, 3], [(1, 2), (3, 2), (2, 1)], "Bi/In Star"), \
            ([1, 2, 3], [(1, 2), (2, 3), (2, 1)], "Bi/Out Star"), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 1)], "Bi-Wedge"), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (2, 3)], "In-Wedge with Bi-Base"), \
            ([1, 2, 3], [(1, 2), (2, 1), (3, 1), (3, 2)], "Out-Wedge with Bi-Base"), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 2)], "Trans-Wedge with Bi-Base"), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 1), (2, 3)], "All But One"), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)], "Full")]:
        ic = information_content(nodes, edges, proportional_p=False, \
                        directed=True, temporal=False, \
                        num_timestamps="auto", \
                        node_colors=None, only_consider_used_nodes=False, \
                        extreme_only_consider_used_nodes=False, \
                        extreme_proportional_p=False)
        padding_string = ""
        for _ in range(0, max_name_length - len(name)):
            padding_string += " "
        print("\n  Information Content of %s graph: %s%f." % \
            (name, padding_string, ic))
        prob_total += math.pow(2.0, -ic)
    print("The total probabilities sum to ~%f." % prob_total)

    prob_total = 0.0
    max_name_length = 19
    print("\n\n\nFor the THREE-node, UNDIRECTED, TEMPORAL graphs, we get:")
    for (nodes, edges, name) in [\
            ([1,2,3], [], "Empty"), \
            ([1,2,3], [(1,2,1)], "Single Edge | Empty"), \
            ([1,2,3], [(1,2,2)], "Empty | Single Edge"), \
            ([1,2,3], [(1,2,1),(1,3,1)], "Wedge | Empty"), \
            ([1,2,3], [(1,2,2),(1,3,2)], "Empty | Wedge"), \
            ([1,2,3], [(1,2,1),(1,2,2)], "Same Edge Twice"), \
            ([1,2,3], [(1,2,1),(1,3,2)], "Edge | Diff Edge"), \
            ([1,2,3], [(1,2,1),(1,3,1),(2,3,1)], "Triangle | Empty"), \
            ([1,2,3], [(1,2,1),(1,3,1),(2,3,1)], "Empty | Triangle"), \
            ([1,2,3], [(1,2,1),(2,3,1),(1,2,2)], "Wedge | Shared Edge"), \
            ([1,2,3], [(1,2,2),(2,3,2),(1,2,1)], "Shared Edge | Wedge"), \
            ([1,2,3], [(1,2,1),(2,3,1),(1,3,2)], "Wedge | Diff Edge"), \
            ([1,2,3], [(1,2,2),(2,3,2),(1,3,1)], "Diff Edge | Wedge"), \
            ([1,2,3], [(1,2,1),(2,3,1),(1,3,1),(1,2,2)], "Triangle | Edge"), \
            ([1,2,3], [(1,2,2),(2,3,2),(1,3,2),(1,2,1)], "Edge | Triangle"), \
            ([1,2,3], [(1,2,1),(2,3,1),(1,2,2),(2,3,2)], "Wedge | Same Wedge"), \
            ([1,2,3], [(1,2,1),(2,3,1),(1,2,2),(1,3,2)], "Wedge | Diff Wedge"), \
            ([1,2,3], [(1,2,1),(2,3,1),(3,1,1),(1,2,2),(1,3,2)], "Triangle | Wedge"), \
            ([1,2,3], [(1,2,2),(2,3,2),(3,1,2),(1,2,1),(1,3,1)], "Wedge | Triangle"), \
            ([1,2,3], [(1,2,1),(2,3,1),(3,1,1),(1,2,2),(2,3,2),(3,1,2)], "Full") \
            ]:
        ic = information_content(nodes, edges, proportional_p=False, \
                        directed=False, temporal=True, \
                        num_timestamps=2, \
                        node_colors=None, only_consider_used_nodes=False, \
                        extreme_only_consider_used_nodes=False, \
                        extreme_proportional_p=False)
        padding_string = ""
        for _ in range(0, max_name_length - len(name)):
            padding_string += " "
        print("\n  Information Content of %s graph: %s%f." % \
            (name, padding_string, ic))
        prob_total += math.pow(2.0, -ic)
    print("The total probabilities sum to ~%f." % prob_total)

    prob_total = 0.0
    print("\n\n\nFor the THREE-node, DIRECTED, TEMPORAL graphs, we get:")
    static_graphs = [\
            ([1, 2, 3], []), \
            ([1, 2, 3], [(1, 2)]), \
            ([1, 2, 3], [(1, 2), (2, 3)]), \
            ([1, 2, 3], [(1, 2), (2, 1)]), \
            ([1, 2, 3], [(1, 3), (2, 3)]), \
            ([1, 2, 3], [(1, 2), (1, 3)]), \
            ([1, 2, 3], [(1, 2), (2, 3), (3, 1)]), \
            ([1, 2, 3], [(1, 2), (2, 3), (1, 3)]), \
            ([1, 2, 3], [(1, 2), (3, 2), (2, 1)]), \
            ([1, 2, 3], [(1, 2), (2, 3), (2, 1)]), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 1)]), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (2, 3)]), \
            ([1, 2, 3], [(1, 2), (2, 1), (3, 1), (3, 2)]), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 2)]), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 1), (2, 3)]), \
            ([1, 2, 3], [(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)])]
    temporal_graphs = __temporal_graph_mashup__(static_graphs, directed=True)
    for (nodes, edges) in temporal_graphs:
        ic = information_content(nodes, edges, proportional_p=False, \
                        directed=True, temporal=True, \
                        num_timestamps=2, \
                        node_colors=None, only_consider_used_nodes=False, \
                        extreme_only_consider_used_nodes=False, \
                        extreme_proportional_p=False)
        prob_total += math.pow(2.0, -ic)
    print("The total probabilities sum to ~%f." % prob_total)
