from algorithmic_utils import get_reverse_dict
import bigfloat
from graph_types import GraphTypes
import math
from nauty_session import NautyTracesSession
import networkx as nx
import random

"""
log2_num_automorphisms(nodes, edges,
                       directed=False, temporal=False,
                       only_consider_used_nodes=False, \
                       node_colors=None)

  Returns the log-base-2 of the number of automorphisms.



log2_independent_edges_prob(nodes, edges, graph_type,
                            proportional_p=True,
                            all_timestamps="auto")

  Returns the log-base-2 of the probability of generating this graph if all
    edges are placed with equal, independent probability. If `proportional_p`
    is set, then the probability used is the ratio of edges to possible edges;
    otherwise the ratio used is 0.5.

    `graph_type` is from the GraphTypes class.

    `all_timestamps` can either be "auto" or a list/set of timestamps. This
        enables including timestamps in which there are no edges.



log2_ER_prob(nodes, edges, graph_type,
             all_timestamps="auto")

  Returns the log-base-2 of the probability of generating this graph with the
    Erdos-Renyi model.

    `graph_type` is from the GraphTypes class.

    `all_timestamps` can either be "auto" or a list/set of timestamps. This
        enables including timestamps in which there are no edges.
"""

# `node_times_if_relevant` must be a dict.
#   It can map nodes to single time values, which are assumed to be the time
#   at which the node is introduced, or it can map a node to a set of time
#   values, which are the times at which the node is considered to be part of
#   the graph. If using a single value, a node is considered to be part of all
#   snapshots during and after the time at which it is introduced.
#
# `timestamps_if_relevant` should be a list. Only needed if there are
#   timestamps not included in `node_times_if_relevant`
def __possible_num_edges__(nodes, graph_type, node_times_if_relevant=None, \
                           timestamps_if_relevant=None):
    if graph_type == GraphTypes.STATIC_UNDIRECTED or \
            graph_type == GraphTypes.STATIC_COLORED_UNDIRECTED:
        return (len(nodes) * (len(nodes) - 1)) / 2

    elif graph_type == GraphTypes.STATIC_DIRECTED or \
            graph_type == GraphTypes.STATIC_COLORED_DIRECTED:
        return len(nodes) * (len(nodes) - 1)

    elif graph_type == GraphTypes.BIPARTITE_UNDIRECTED:
        return len(nodes[0]) * len(nodes[1])

    elif graph_type == GraphTypes.BIPARTITE_DIRECTED:
        return 2 * len(nodes[0]) * len(nodes[1])

    elif graph_type == GraphTypes.NODE_JOINING_UNDIRECTED or \
            graph_type == GraphTypes.NODE_JOINING_DIRECTED or \
            graph_type == GraphTypes.NODE_JOINING_DIRECTED_ONE_WAY:

        assert node_times_if_relevant is not None

        intra_directions = 1 + int(\
            graph_type == GraphTypes.NODE_JOINING_DIRECTED or \
            graph_type == GraphTypes.NODE_JOINING_DIRECTED_ONE_WAY)
        extra_directions = 1 + int(\
            graph_type == GraphTypes.NODE_JOINING_DIRECTED)

        time_partitions = get_reverse_dict(node_times_if_relevant)
        if timestamps_if_relevant is None:
            times = [t for t, _ in time_partitions.items()]
        else:
            times = list(timestamps_if_relevant)
            for t in times:
                if t not in time_partitions:
                    time_partitions[t] = set()
        times.sort()
        num_possible_edges = 0
        num_previous_nodes = 0
        for t in times:
            num_nodes = len(time_partitions[t])

            intra_edges = ((num_nodes * (num_nodes - 1)) / 2) * \
                            intra_directions
            extra_edges = num_nodes * num_previous_nodes * extra_directions

            num_possible_edges += intra_edges + extra_edges
            num_previous_nodes += num_nodes
        return num_possible_edges
            
    elif graph_type == GraphTypes.TEMPORAL_UNDIRECTED or \
            graph_type == GraphTypes.TEMPORAL_COLORED_UNDIRECTED or \
            graph_type == GraphTypes.TEMPORAL_DIRECTED or \
            graph_type == GraphTypes.TEMPORAL_COLORED_DIRECTED:

        assert node_times_if_relevant is not None

        introduction_times = True
        for node, time_s in node_times_if_relevant.items():
            if type(time_s) is not set:
                break
            if type(time_s) is set and len(time_s) > 1:
                introduction_times = False
                break

        num_directions = 1 + int(\
            graph_type == GraphTypes.TEMPORAL_DIRECTED or \
            graph_type == GraphTypes.TEMPORAL_COLORED_DIRECTED)

        time_partitions = get_reverse_dict(node_times_if_relevant)
        if timestamps_if_relevant is None:
            times = [t for t, _ in time_partitions.items()]
        else:
            times = list(timestamps_if_relevant)
            for t in times:
                if t not in time_partitions:
                    time_partitions[t] = set()
        times.sort()
        num_possible_edges = 0
        cumulative_num_nodes = 0
        for t in times:
            num_nodes = len(time_partitions[t])
            cumulative_num_nodes += num_nodes
            if introduction_times:
                num_possible_edges += num_directions * \
                    ((cumulative_num_nodes * (cumulative_num_nodes - 1)) / 2)
            else:
                num_possible_edges += num_directions * \
                    ((num_nodes * (num_nodes - 1)) / 2)
        return num_possible_edges

# Options for `mode` include:
#   "cumulative"
#   "permanent"
#   "transitive"
#   "source_transitive" -- only give a node a timestamp t if it is a source
#       at time t
def __node_timestamps__(nodes, sorted_edges, \
                        all_timestamps="auto", mode="transitive"):
    if len(sorted_edges) == 0 and all_timestamps=="auto":
        all_timestamps = set([1])
    elif all_timestamps == "auto":
        all_timestamps = set()
        for edge in sorted_edges:
            all_timestamps.add(edge[2])
    sorted_timestamps = sorted(list(all_timestamps))

    assert mode == "permanent" or mode == "transitive" or \
           mode == "cumulative" or mode == "source_transitive"

    if mode == "permanent":
        # print("You are considering all timesteps for each node. Using " + \
        #       "introduction timestamps rather than sets of times.")

        return ({n: sorted_timestamps[0] for n in nodes}, sorted_timestamps)

    cumulative = 1
    transitive = 2
    source_transitive = 3
    if mode == "cumulative":
        mode = cumulative
    elif mode == "transitive":
        mode = transitive
    else:
        mode = source_transitive

    node_timestamps = {n: set() for n in nodes}
    prev_t = None
    if mode == cumulative:
        t_plus_sets = {sorted_timestamps[i]: \
                        set([t for t in sorted_timestamps[i:]]) \
                            for i in range(0, len(sorted_timestamps))}
    for edge in sorted_edges:
        source = edge[0]
        target = edge[1]
        t = edge[2]

        if prev_t is not None and t < prev_t:
            raise ValueError("Error! Edges were not sorted by timestamp!")
        if prev_t is None or t > prev_t:
            pass
        prev_t = t

        if mode == source_transitive:
            node_timestamps[source].add(t)
        elif mode == transitive:
            node_timestamps[source].add(t)
            node_timestamps[target].add(t)
        else:  # cumulative
            if len(node_timestamps[source]) == 0:
                node_timestamps[source] = t_plus_sets[t]
            if len(node_timestamps[target]) == 0:
                node_timestamps[target] = t_plus_sets[t]

    return (node_timestamps, sorted_timestamps)

def __possible_num_edges_with_assumptions__(nodes, sorted_edges, graph_type, \
                                            all_timestamps="auto"):

    if      graph_type == GraphTypes.STATIC_UNDIRECTED or \
            graph_type == GraphTypes.STATIC_COLORED_UNDIRECTED or \
            graph_type == GraphTypes.STATIC_DIRECTED or \
            graph_type == GraphTypes.STATIC_COLORED_DIRECTED or \
            graph_type == GraphTypes.BIPARTITE_UNDIRECTED or \
            graph_type == GraphTypes.BIPARTITE_DIRECTED:
        # The graph is "static"
        return __possible_num_edges__(nodes, graph_type, \
                    node_times_if_relevant=None, \
                    timestamps_if_relevant=None)

    elif    graph_type == GraphTypes.NODE_JOINING_UNDIRECTED or \
            graph_type == GraphTypes.NODE_JOINING_DIRECTED or \
            graph_type == GraphTypes.NODE_JOINING_DIRECTED_ONE_WAY or \
            graph_type == GraphTypes.TEMPORAL_UNDIRECTED or \
            graph_type == GraphTypes.TEMPORAL_COLORED_UNDIRECTED or \
            graph_type == GraphTypes.TEMPORAL_DIRECTED or \
            graph_type == GraphTypes.TEMPORAL_COLORED_DIRECTED:
        # The graph is "temporal"
        if      graph_type == GraphTypes.NODE_JOINING_UNDIRECTED or \
                graph_type == GraphTypes.NODE_JOINING_DIRECTED:
            (nt, at) = __node_timestamps__(nodes, sorted_edges, \
                                           all_timestamps=all_timestamps, \
                                           mode="transitive")
        elif    graph_type == GraphTypes.NODE_JOINING_DIRECTED_ONE_WAY:
            (nt, at) = __node_timestamps__(nodes, sorted_edges, \
                                           all_timestamps=all_timestamps, \
                                           mode="source_transitive")
        else:
            (nt, at) = __node_timestamps__(nodes, sorted_edges, \
                                           all_timestamps=all_timestamps, \
                                           mode="permanent")

        return __possible_num_edges__(nodes, graph_type, \
                    node_times_if_relevant=nt, \
                    timestamps_if_relevant=at)
    else:
        raise ValueError("Error! log2_..._prob() not coded for graph " + \
                         "type: " + GraphTypes.TYPE_NAME_STRINGS[graph_type])

def __sorted_edges__(edges):
    return [(a, b, t) for (t, a, b) in \
            sorted([(t, a, b) for (a, b, t) in edges])]
def log2_independent_edges_prob(nodes, edges, graph_type, \
                                proportional_p=True, \
                                all_timestamps="auto", \
                                node_colors=None):
    if node_colors is None:
        basic_highlights = [list(nodes)]
    else:
        assert len(node_colors) == len(nodes)

    if GraphTypes.IS_TEMPORAL[graph_type]:
        edges = __sorted_edges__(edges)
    pne = __possible_num_edges_with_assumptions__(nodes, edges, graph_type, \
                                                  all_timestamps=all_timestamps)

    if proportional_p:
        p = len(edges) / pne
    else:
        p = 0.5

    if p == 0.0:
        return 1.0
    elif p == 1.0:
        return 1.0
    log2_p = math.log(p, 2.0)
    log2_1_minus_p = math.log(1.0 - p, 2.0)
    log2_matrix_prob = len(edges) * log2_p + (pne - len(edges)) * log2_1_minus_p

    directed = GraphTypes.IS_DIRECTED[graph_type]
    temporal = GraphTypes.IS_TEMPORAL[graph_type]
    log2_na = log2_num_automorphisms(nodes, edges, \
                                     directed=directed, temporal=temporal, \
                                     only_consider_used_nodes=False, \
                                     node_colors=node_colors)

    if node_colors is None:
        color_partitions_sizes = [len(nodes)]
    else:
        color_partitions_sizes = {}
        for n in nodes:
            c = node_colors[n]
            if c not in color_partitions_sizes:
                color_partitions_sizes[c] = 0
            color_partitions_sizes[c] += 1
        color_partitions_sizes = \
            [count for c, count in color_partitions_sizes.items()]

    log2_max_permutations = 0.0
    for color_partition_size in color_partitions_sizes:
        for i in range(2, color_partition_size + 1):
            log2_max_permutations += math.log(i, 2.0)

    return (log2_max_permutations - log2_na) + log2_matrix_prob

def log2_ER_prob(nodes, edges, graph_type, \
                 all_timestamps="auto", node_colors=None):

    if len(edges) == 0:
        return 1.0

    if GraphTypes.IS_TEMPORAL[graph_type]:
        edges = __sorted_edges__(edges)
    pne = __possible_num_edges_with_assumptions__(nodes, edges, graph_type, \
                                                  all_timestamps=all_timestamps)

    if len(edges) == pne:
        return 1.0

    log2_matrix_prob = 0.0
    for i in range(0, len(edges)):
        n = len(edges) - i
        m = pne - i
        log2_matrix_prob += math.log(n / m, 2.0)

    directed = GraphTypes.IS_DIRECTED[graph_type]
    temporal = GraphTypes.IS_TEMPORAL[graph_type]
    log2_na = log2_num_automorphisms(nodes, edges, \
                                     directed=directed, temporal=temporal, \
                                     only_consider_used_nodes=False, \
                                     node_colors=node_colors)

    if node_colors is None:
        color_partitions_sizes = [len(nodes)]
    else:
        color_partitions_sizes = {}
        for n in nodes:
            c = node_colors[n]
            if c not in color_partitions_sizes:
                color_partitions_sizes[c] = 0
            color_partitions_sizes[c] += 1
        color_partitions_sizes = \
            [count for c, count in color_partitions_sizes.items()]
        
    log2_max_permutations = 0.0
    for color_partition_size in color_partitions_sizes:
        for i in range(2, color_partition_size + 1):
            log2_max_permutations += math.log(i, 2.0)

    return (log2_max_permutations - log2_na) + log2_matrix_prob

def log2_num_automorphisms(nodes, edges, \
                           directed=False, temporal=False, \
                           only_consider_used_nodes=False, \
                           node_colors=None):
    
    if only_consider_used_nodes:
        if len(edges) == 0:
            raise ValueError("Error! Cannot only consider nodes with edges " + \
                             "when there are no edges. Set " + \
                             "`only_consider_used_nodes` to False.")

    nauty_run_nodes = set()
    for edge in edges:
        nauty_run_nodes.add(edge[0])
        nauty_run_nodes.add(edge[1])

    if not only_consider_used_nodes:
        # Remove unused nodes to prevent factorial scaling of Nauty computation,
        #   then account for them after the run.
        if type(nodes) is set:
            unused_nodes = nodes - nauty_run_nodes
        else:
            unused_nodes = set(nodes) - nauty_run_nodes

        # Add a single unused node to the graph, just in case the graph has
        #   zero nodes otherwise.
        if len(unused_nodes) > 0:
            an_unused_node = unused_nodes.pop()
            unused_nodes.add(an_unused_node)
            nauty_run_nodes.add(an_unused_node)

    if node_colors is None:
        nauty_run_node_colors = None
    else:
        nauty_run_node_colors = {n: node_colors[n] for n in nauty_run_nodes}

    old_context = bigfloat.getcontext()
    i = 1
    done = False
    while not done:
        bf_context = bigfloat.Context(precision=2000*i, \
                                      emax=1000000*i, emin=-1000000*i)
        bigfloat.setcontext(bf_context)
        # Compute the log-number of automorphisms.
        (num_automorphisms, _) = __get_graph_info__(\
                nauty_run_nodes, edges, temporal, directed, \
                use_color_direction=True, \
                get_canon_order=False, \
                node_colors=nauty_run_node_colors)

        if num_automorphisms != bigfloat.BigFloat("inf"):
            done = True
        else:
            print("Not enough bits to store num of automorhisms.")
            print("Retrying.")
            
    bigfloat.setcontext(old_context)

    log2_num_automorphisms = bigfloat.log2(num_automorphisms)

    if not only_consider_used_nodes:
        if node_colors is None:
            partition_sizes = [len(unused_nodes)]
        else:
            partition_sizes = {}
            for node in unused_nodes:
                color = node_colors[node]
                if color not in partition_sizes:
                    partition_sizes[color] = 0
                partition_sizes[color] += 1
            partition_sizes = [count for col, count in partition_sizes.items()]
        # Add automorphisms for unused nodes:
        for partition_size in partition_sizes:
            for i in range(2, partition_size + 1):
                log2_num_automorphisms += math.log(i, 2.0)

    return log2_num_automorphisms

# Gives number of automorphisms and canonical node order.
def __get_graph_info__(nodes, given_edges, temporal, directed, \
                       use_color_direction=True, get_canon_order=False, \
                       node_colors=None):
    if type(nodes) is not set:
        nodes = set(nodes)

    if directed and not use_color_direction:
        graph = nx.DiGraph()
    else:
        graph = nx.Graph()

    for node in nodes:
        graph.add_node(node)

    orbits = [list(nodes)]  # at first, everything in one orbit

    if node_colors is None:
        basic_highlights = [list(nodes)]
    else:
        assert len(node_colors) == len(nodes)
        basic_highlights = []
        color_collection = [node_colors[n] for n in nodes]
        color_collection.sort()
        next_highlight_idx = 0
        color_to_highlight = {}
        for c in color_collection:
            if c not in color_to_highlight:
                color_to_highlight[c] = next_highlight_idx
                next_highlight_idx += 1
                basic_highlights.append([])
        for node in nodes:
            basic_highlights[color_to_highlight[node_colors[node]]].append(node)
        

    # Use only to store highlights for nodes not in `nodes`
    extra_highlights = []
    if temporal:
        if directed and use_color_direction:
            # Since we're already adding nodes to indicate direction,
            #   we can color those same nodes to set timesamps.
            timestamp_sets_by_edge = {}
            for (a, b, t) in given_edges:
                if (a, b) not in timestamp_sets_by_edge:
                    timestamp_sets_by_edge[(a, b)] = set()
                timestamp_sets_by_edge[(a, b)].add(t)
            timestamp_tuples_by_edge = {}
            timestamp_tuples = set()
            for (a, b), s in timestamp_sets_by_edge.items():
                timestamp_tuples_by_edge[(a, b)] = tuple(sorted(list(s)))
                timestamp_tuples.add(timestamp_tuples_by_edge[(a, b)])
            del timestamp_sets_by_edge

            # Need to sort the tuples to get consistent highlight order.
            #   Otherwise, the canonical node ordering may be changed.
            extra_highlights = [[]]
            timestamp_tuples = sorted(list(timestamp_tuples))
            next_highlight_idx = 1
            timestamp_highlights = {}
            for tup in timestamp_tuples:
                timestamp_highlights[tup] = next_highlight_idx
                next_highlight_idx += 1
                extra_highlights.append([])

            for (a, b), tup in timestamp_tuples_by_edge.items():
                color_idx = timestamp_highlights[tup]
                graph.add_node((a, b, 1))
                graph.add_node((a, b, 2))
                extra_highlights[0].append((a, b, 1))
                extra_highlights[color_idx].append((a, b, 2))

                graph.add_edge(a, (a, b, 1))
                graph.add_edge((a, b, 1), (a, b, 2))
                graph.add_edge((a, b, 2), b)
                graph.add_edge(a, b)
            # print(extra_highlights)
        else:
            edge_sets_by_timestamp = {}
            for (a, b, t) in given_edges:
                if t not in edge_sets_by_timestamp:
                    edge_sets_by_timestamp[t] = set()
                edge_sets_by_timestamp[t].add((a, b))

            timestamps = \
                sorted([t for t, _ in edge_sets_by_timestamp.items()])

            for i in range(0, len(timestamps)):
                t = timestamps[i]
                for node in nodes:
                    graph.add_node((node, t))
                    if i == 0:
                        graph.add_edge(node, (node, t))
                    else:
                        graph.add_edge((node, timestamps[i - 1]), (node, t))
                for (a, b) in edge_sets_by_timestamp[t]:

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

    if use_color_direction:
        NTSession = NautyTracesSession(graph, mode="Traces")
    else:
        NTSession = NautyTracesSession(graph, mode="Nauty", sparse=True)

    NTSession.set_colors_by_highlights(basic_highlights + extra_highlights)
    na = NTSession.get_num_automorphisms()
    if get_canon_order:
        co = NTSession.get_canonical_order()
    NTSession.complete()
    na = na.get()
    if get_canon_order:
        co = co.get()
        return (na, co)
    return (na, None)

# average_ER gives log of average random NA
def measure_of_structure(nodes_lists, edges, graph_type, \
                         all_timestamps="auto", \
                         ER_try_count=10, \
                         node_colors=None, \
                         average_ER=False):
    directed = GraphTypes.IS_DIRECTED[graph_type]
    temporal = GraphTypes.IS_TEMPORAL[graph_type]
    all_nodes = []
    for l in nodes_lists:
        all_nodes += l
    assert len(set(all_nodes)) == len(all_nodes)

    min_ER_log2_na = None
    ER_log2_na_values = []

    if GraphTypes.IS_TEMPORAL[graph_type]:
        if all_timestamps == "auto":
            all_timestamps = set()
            for edge in edges:
                all_timestamps.add(edge[2])
            all_timestamps = sorted(list(all_timestamps))

        if graph_type == GraphTypes.NODE_JOINING_UNDIRECTED or \
                graph_type == GraphTypes.NODE_JOINING_DIRECTED or \
                graph_type == GraphTypes.NODE_JOINING_DIRECTED_ONE_WAY:
            assert "This Code" == "Not Implemented"
            nt = {n: set() for n in all_nodes}
            for edge in edges:
                nt[edge[0]].add(edge[2])
                nt[edge[1]].add(edge[2])
        else:
            nt = None

        done = False
        while not done:
            done = True
            for i in range(0, ER_try_count):
                ER_graph = __create_ER_graph__(\
                        nodes_lists, len(edges), graph_type, \
                        num_timestamps_if_relevant=len(all_timestamps), \
                        node_timestamps_if_relevant=nt)
                (ER_nodes, ER_edges) = (ER_graph[0][0], ER_graph[1])
                na_value = log2_num_automorphisms(ER_nodes, ER_edges, \
                                         directed=directed, temporal=temporal, \
                                         only_consider_used_nodes=False, \
                                         node_colors=node_colors)
                ER_log2_na_values.append(na_value)
                if float(na_value) == 0.0:
                    min_ER_log2_na = na_value
                    if not average_ER:
                        break
                elif min_ER_log2_na is None or na_value < min_ER_log2_na:
                    min_ER_log2_na = na_value
                    if not average_ER:
                        done = False
                        break

    elif graph_type == GraphTypes.BIPARTITE_UNDIRECTED or \
            graph_type == GraphTypes.BIPARTITE_DIRECTED:
        assert "Needs graph colors for num-automorphisms" == "True"
        done = False
        while not done:
            done = True
            for i in range(0, ER_try_count):
                ER_graph = __create_ER_graph__(\
                            nodes_lists, len(edges), graph_type, \
                            num_timestamps_if_relevant=None, \
                            node_timestamps_if_relevant=None)
                ER_nodes = ER_graph[0][0] | ER_graph[0][1]
                ER_edges = ER_graph[1]
                colors = "NOWHERE TO BE FOUND RIGHT NOW"
                assert type(colors) is dict
                na_value = log2_num_automorphisms(ER_nodes, ER_edges, \
                                     directed=directed, temporal=temporal, \
                                     only_consider_used_nodes=False, \
                                     node_colors=node_colors)
                ER_log2_na_values.append(na_value)

                if float(na_value) == 0.0:
                    min_ER_log2_na = na_value
                    if not average_ER:
                        break
                elif min_ER_log2_na is None or na_value < min_ER_log2_na:
                    min_ER_log2_na = na_value
                    if not average_ER:
                        done = False
                        break
    else:
        done = False
        while not done:
            done = True
            for i in range(0, ER_try_count):
                ER_graph = __create_ER_graph__(\
                            nodes_lists, len(edges), graph_type, \
                            num_timestamps_if_relevant=None, \
                            node_timestamps_if_relevant=None)
                (ER_nodes, ER_edges) = (ER_graph[0][0], ER_graph[1])
                na_value = log2_num_automorphisms(ER_nodes, ER_edges, \
                                     directed=directed, temporal=temporal, \
                                     only_consider_used_nodes=False, \
                                     node_colors=node_colors)
                ER_log2_na_values.append(na_value)

                if float(na_value) == 0.0:
                    min_ER_log2_na = na_value
                    if not average_ER:
                        break
                elif min_ER_log2_na is None or na_value < min_ER_log2_na:
                    min_ER_log2_na = na_value
                    if not average_ER:
                        done = False
                        break

    graph_log2_num_automorphisms = log2_num_automorphisms(all_nodes, edges, \
                                        directed=directed, temporal=temporal, \
                                        only_consider_used_nodes=False, \
                                        node_colors=node_colors)

    if not average_ER:
        return (graph_log2_num_automorphisms, min_ER_log2_na)

    # Take log of average or average log?
    assert len(ER_log2_na_values) == ER_try_count
    avg = 0.0
    overflow_tracker = 0.0
    for log_val in ER_log2_na_values:
        avg += bigfloat.pow(2.0, log_val)
        if avg < overflow_tracker:
            raise ValueError("Error! Bigfloat environment not large enough " + \
                             "for averaging number of automorphisms.")
        overflow_tracker = avg
    avg /= len(ER_log2_na_values)
    return (graph_log2_num_automorphisms, float(bigfloat.log2(avg)))

# In a node_joining network, the node_timestamps are the join times of the
#   nodes.
#
# Returns a Prepped_Fixed_Graph
#
# Otherwise, nodes are assumed to be present throughout.
def __create_ER_graph__(nodes_lists, num_edges, graph_type, \
                        num_timestamps_if_relevant=None, \
                        node_timestamps_if_relevant=None):

    # highlights = list(nodes_lists)
    edges = set()
    if graph_type == GraphTypes.NODE_JOINING_UNDIRECTED:
        assert node_timestamps_if_relevant is not None
        assert num_timestamps_if_relevant is None
        assert "this code" == "not implemented"

    elif graph_type == GraphTypes.NODE_JOINING_DIRECTED:
        assert node_timestamps_if_relevant is not None
        assert num_timestamps_if_relevant is None
        assert "this code" == "not implemented"

    elif graph_type == GraphTypes.NODE_JOINING_DIRECTED_ONE_WAY:
        assert node_timestamps_if_relevant is not None
        assert num_timestamps_if_relevant is None
        assert "this code" == "not implemented"

    elif GraphTypes.IS_TEMPORAL[graph_type]:
        assert node_timestamps_if_relevant is None
        assert num_timestamps_if_relevant is not None
        assert len(nodes_lists) == 1

        # TODO: Change return type
        """
        if GraphTypes.IS_DIRECTED[graph_type]:
            next_A = max([max(nl) for nl in nodes_lists]) + 1
            next_B = next_A + num_edges
            highlights += [[], []]
        else:
            next_A = max([max(nl) for nl in nodes_lists]) + 1
            highlights += [[]]

        flattened_edges = set()
        """

        for i in range(0, num_edges):
            valid = False
            while not valid:
                source = random.randint(0, len(nodes_lists[0]) - 1)
                # This way of generating target ensures no self-loops.
                target = random.randint(0, len(nodes_lists[0]) - 2)
                if target >= source:
                    target += 1
                time = random.randint(1, num_timestamps_if_relevant)

                source = nodes_lists[0][source]  # Move from indices to labels
                target = nodes_lists[0][target]

                if GraphTypes.IS_DIRECTED[graph_type]:
                    if (source, target, time) not in edges:
                        edges.add((source, target, time))
                        valid = True
                else:
                    s = min(source, target)
                    t = max(source, target)
                    if (s, t, time) not in edges:
                        edges.add((s, t, time))
                        valid = True

    else:
        assert node_timestamps_if_relevant is None
        assert num_timestamps_if_relevant is None

        if graph_type == GraphTypes.BIPARTITE_UNDIRECTED or \
                graph_type == GraphTypes.BIPARTITE_DIRECTED:
            # Bipartite
            assert len(nodes_lists) == 2

            for i in range(0, num_edges):
                valid = False
                while not valid:
                    source = random.randint(0, len(nodes_lists[0]) - 1)
                    target = random.randint(0, len(nodes_lists[1]) - 1)
                    source = nodes_lists[0][source]
                    target = nodes_lists[1][target]

                    if GraphTypes.IS_DIRECTED[graph_type]:
                        swap = random.randint(0, 1)
                        if swap == 1:
                            swap = source
                            source = target
                            target = swap
                    if (source, target) not in edges:
                        edges.add((source, target))
                        valid = True

        else:
            assert len(nodes_lists) == 1

            for i in range(0, num_edges):
                valid = False
                while not valid:
                    source = random.randint(0, len(nodes_lists[0]) - 1)
                    # This way of generating target ensures no self-loops.
                    target = random.randint(0, len(nodes_lists[0]) - 2)
                    if target >= source:
                        target += 1
                    source = nodes_lists[0][source]
                    target = nodes_lists[0][target]

                    if GraphTypes.IS_DIRECTED[graph_type]:
                        if (source, target) not in edges:
                            edges.add((source, target))
                            valid = True
                    else:
                        s = min(source, target)
                        t = max(source, target)
                        if (s, t) not in edges:
                            edges.add((s, t))
                            valid = True

    return (nodes_lists, edges)


# Used only for testing.
def __canonical_string__(nodes, edges, directed=False, temporal=False,\
                         use_color_direction=True, \
                         node_coloring=None):
    (_, co) = __get_graph_info__(nodes, edges, temporal, directed, \
                       use_color_direction=use_color_direction, \
                       get_canon_order=True, \
                       node_colors=node_coloring)
    node_to_idx_map = {}
    for i in range(0, len(co)):
        node_to_idx_map[co[i]] = i

    nodes = [i for i in range(0, len(nodes))]

    if node_coloring is None:
        colors = []
    else:
        colors = [node_coloring[co[i]] for i in range(0, len(nodes))]

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
    return str((nodes, colors, edges))

# Used only for testing.
def __temporal_graph_mashup__(graphs, directed=False, node_coloring=None):
    from algorithmic_utils import get_all_k_permutations
    num_nodes = len(graphs[0][0])
    perms = get_all_k_permutations(num_nodes, num_nodes)
    canonical_strings = set()

    temporal_graphs = []

    for (nodes, first_edges) in graphs:
        for (_, second_edges) in graphs:
            for perm_A in perms:
                for perm_B in perms:
                    edges = []
                    for (a, b) in first_edges:
                        edges.append((nodes[perm_A[a-1]], nodes[perm_A[b-1]],1))
                    for (a, b) in second_edges:
                        edges.append((nodes[perm_B[a-1]], nodes[perm_B[b-1]],2))
                    cs = __canonical_string__(nodes, edges, directed=directed,\
                                              temporal=True, \
                                              use_color_direction=True, \
                                              node_coloring=node_coloring)
                    if cs not in canonical_strings:
                        canonical_strings.add(cs)
                        temporal_graphs.append((nodes, edges))
    return temporal_graphs

class Prepped_Fixed_Graph:

    # Pass [in-neighbors, out-neighbors] if directed
    # Pass [neighbors] is undirected
    def __init__(self, nodes, neighbors_sets, directed, highlights=None):
        self.directed = directed
        if directed:
            self.predecessors = neighbors_sets[0]
            self.successors = neighbors_sets[1]
        else:
            self.neighbors = neighbors_sets[0]
        if highlights is None:
            self.highlights = [n for n in nodes]
        else:
            self.highlights = highlights

        self.nodes = nodes

    def nodes(self):
        return self.nodes

    def is_directed(self):
        return self.directed

    def neighbors(self, node):
        if self.directed:
            raise ValueError("Can only call neighbors() on undirected graph.")
        return self.neighbors[node]

    def successors(self, node):
        if not self.directed:
            raise ValueError("Can only call successors() on directed graph.")
        return self.successors[node]

    def predecessors(self, node):
        if not self.directed:
            raise ValueError("Can only call predecessors() on directed graph.")
        return self.predecessors[node]

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
        ic = log2_independent_edges_prob(nodes, edges, GraphTypes.STATIC_UNDIRECTED, \
                                         proportional_p=False, \
                                         all_timestamps="auto")
        padding_string = ""
        for _ in range(0, max_name_length - len(name)):
            padding_string += " "
        print("\n  Information Content of %s graph: %s%f." % \
            (name, padding_string, ic))
        prob_total += math.pow(2.0, ic)
    print("The total probabilities sum to ~%f." % prob_total)

    prob_total = 0.0
    max_name_length = 26
    print("For the THREE-node, UNDIRECTED, COLORED graphs, we get:")
    for (nodes, edges, coloring, name) in [\
            ([1,2,3], [], {1:1,2:1,3:1}, "Empty, All Same"), \
            ([1,2,3], [], {1:1,2:1,3:2}, "Empty, Two Same"), \
            ([1,2,3], [], {1:1,2:2,3:3}, "Empty, All Diff"), \
            ([1,2,3], [(1, 2)], {1:1,2:1,3:1}, "Single Edge, All Same"), \
            ([1,2,3], [(1, 2)], {1:1,2:1,3:2}, "Single Edge, Two Same Same"), \
            ([1,2,3], [(1, 2)], {1:1,2:2,3:1}, "Single Edge, Two Diff Same"), \
            ([1,2,3], [(1, 2)], {1:1,2:2,3:3}, "Single Edge, All Diff 3"), \
            ([1,2,3], [(1, 2)], {1:1,2:3,3:2}, "Single Edge, All Diff 2"), \
            ([1,2,3], [(1, 2)], {1:2,2:3,3:1}, "Single Edge, All Diff 1"), \
            ([1,2,3], [(1, 3), (2, 3)], {1:1,2:1,3:1}, "Wedge, All Same"), \
            ([1,2,3], [(1, 3), (2, 3)], {1:1,2:1,3:2}, "Wedge, Two Same Same"), \
            ([1,2,3], [(1, 3), (2, 3)], {1:1,2:2,3:1}, "Wedge, Two Diff Same"), \
            ([1,2,3], [(1, 2), (2, 3)], {1:1,2:2,3:3}, "Wedge, All Diff 3"), \
            ([1,2,3], [(1, 2), (2, 3)], {1:1,2:3,3:2}, "Wedge, All Diff 2"), \
            ([1,2,3], [(1, 2), (2, 3)], {1:2,2:3,3:1}, "Wedge, All Diff 1"), \
            ([1,2,3], [(1, 2), (2, 3), (3, 1)], {1:1,2:1,3:1}, "Triangle, All Same"), \
            ([1,2,3], [(1, 2), (2, 3), (3, 1)], {1:1,2:1,3:2}, "Triangle, Two Same"), \
            ([1,2,3], [(1, 2), (2, 3), (3, 1)], {1:1,2:2,3:3}, "Triangle, All Diff")]:
        ic = log2_independent_edges_prob(nodes, edges, GraphTypes.STATIC_UNDIRECTED, \
                                         proportional_p=False, \
                                         all_timestamps="auto", \
                                         node_colors=coloring)
        padding_string = ""
        for _ in range(0, max_name_length - len(name)):
            padding_string += " "
        print("\n  Information Content of %s graph: %s%f." % \
            (name, padding_string, ic))
        prob_total += math.pow(2.0, ic)
    print("We had 3 different colorings: all same, one diff, all diff.")
    print("  Thus, the total probabilities sum to ~%f." % prob_total)

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
        ic = log2_independent_edges_prob(nodes, edges, GraphTypes.STATIC_DIRECTED, \
                                         proportional_p=False, \
                                         all_timestamps="auto")
        padding_string = ""
        for _ in range(0, max_name_length - len(name)):
            padding_string += " "
        print("\n  Information Content of %s graph: %s%f." % \
            (name, padding_string, ic))
        prob_total += math.pow(2.0, ic)
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
        ic = log2_independent_edges_prob(nodes, edges, GraphTypes.TEMPORAL_UNDIRECTED, \
                                         proportional_p=False, \
                                         all_timestamps=[1, 2])
        padding_string = ""
        for _ in range(0, max_name_length - len(name)):
            padding_string += " "
        print("\n  Information Content of %s graph: %s%f." % \
            (name, padding_string, ic))
        prob_total += math.pow(2.0, ic)
    print("The total probabilities sum to ~%f." % prob_total)

    print("\n\n\nFor the THREE-node, COLORED, DIRECTED, TEMPORAL graphs, we get:")
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
    for coloring in [{1: 1, 2: 1, 3: 1}, \
                     {1: 1, 2: 1, 3: 2}, \
                     {1: 1, 2: 2, 3: 3}]:
        prob_total = 0.0
        temporal_graphs = __temporal_graph_mashup__(static_graphs, directed=True, \
                                                    node_coloring=coloring)
        print("  With coloring %s, we get %d graphs total." % \
                (coloring, len(temporal_graphs)))
        for (nodes, edges) in temporal_graphs:
            ic = log2_independent_edges_prob(nodes, edges, GraphTypes.TEMPORAL_DIRECTED, \
                                             proportional_p=False, \
                                             all_timestamps=[1, 2], \
                                             node_colors=coloring)
            prob_total += math.pow(2.0, ic)
        print("    The total probabilities sum to ~%f." % prob_total)
