from algorithmic_utils import get_reverse_dict
import bigfloat
from graph_types import GraphTypes
import math

"""
log2_num_automorphisms(nodes, edges,
                       directed=False, temporal=False,
                       only_consider_used_nodes=False)
"""

# node_times_if_relevant must be a dict.
#   It can map nodes to single time values, which are assumed to be the time
#   at which the node is introduced, or it can map nodes to a set of time
#   values, which are the times at which the node is considered to be part of
#   the graph. If using a single value, a node is considered to be part of all
#   snapshots during and after the time at which it is introduced.
def __possible_num_edges__(nodes, graph_type, node_times_if_relevant=None):
    if graph_type == GraphType.STATIC_UNDIRECTED or \
            graph_type == GraphType.STATIC_COLORED_UNDIRECTED:
        return (len(nodes) * (len(nodes) - 1)) / 2

    elif graph_type == GraphType.STATIC_DIRECTED or \
            graph_type == GraphType.STATIC_COLORED_DIRECTED:
        return len(nodes) * (len(nodes) - 1)

    elif graph_type == GraphType.BIPARTITE_UNDIRECTED:
        return len(nodes[0]) * len(nodes[1])

    elif graph_type == GraphType.BIPARTITE_DIRECTED:
        return 2 * len(nodes[0]) * len(nodes[1])

    elif graph_type == GraphType.NODE_JOINING_UNDIRECTED or \
            graph_type == GraphType.NODE_JOINING_DIRECTED or \
            graph_type == GraphType.NODE_JOINING_DIRECTED_ONE_WAY:

        assert node_times_if_relevant is not None

        intra_directions = 1 + int(\
            graph_type == GraphType.NODE_JOINING_DIRECTED or \
            graph_type == GraphType.NODE_JOINING_DIRECTED_ONE_WAY)
        extra_directions = 1 + int(\
            graph_type == GraphType.NODE_JOINING_DIRECTED)

        time_partitions = get_reverse_dict(node_times_if_relevant)
        times = [t for t, _ in time_partitions.items()]
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
            
    elif graph_type == GraphType.TEMPORAL_UNDIRECTED or \
            graph_type == GraphType.TEMPORAL_COLORED_UNDIRECTED or \
            graph_type == GraphType.TEMPORAL_DIRECTED or \
            graph_type == GraphType.TEMPORAL_COLORED_DIRECTED:

        assert node_times_if_relevant is not None

        introduction_times = True
        for node, time_s in node_times_if_relevant:
            if type(time_s) is not set:
                break
            if type(time_s) is set and len(time_s) > 1:
                introduction_times = False
                break

        num_directions = 1 + int(\
            graph_type == GraphType.TEMPORAL_DIRECTED or \
            graph_type == GraphType.TEMPORAL_COLORED_DIRECTED)

        time_partitions = get_reverse_dict(node_times_if_relevant)
        times = [t for t, _ in time_partitions.items()]
        times.sort()
        num_possible_edges = 0
        cumulative_num_nodes = 0
        for t in times:
            num_nodes = len(time_partitions[t])
            cumulative_num_nodes += num_nodes
            if introduction_times:
                num_possible_edges += num_directions * \
                    ((cumulative_nodes * (cumulative_nodes - 1)) / 2)
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
    if all_timestamps == "auto":
        all_timestamps = set()
        for edge in sorted_edges:
            all_timestamps.add(edge[2])
    sorted_timestamps = sorted(list(all_timestamps))

    assert mode == "permanent" or mode == "transitive" or \
           mode == "cumulative" or mode == "source_transitive"

    if mode == "permanent":
        print("You are considering all timesteps for each node. Perhaps your" +\
              " code could be more efficient if you have a special case.")

        # Note that this does not use excessive space, because all_timestamps is
        #   a pointer.
        return {n: all_timestamps for n in nodes}

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

    return node_timestamps

def log2_ER_prob(nodes, edges, graph_type, \
                 proportional_p=True, \
                 node_timestamps=None):
    pass

# TODO: Finish/re-write. Instead of all these flags, have graph types.
def __matrix_meta_info__(num_nodes, edges, num_timestamps="auto", \
                         directed=False, temporal=False, \
                         only_consider_used_nodes=False, \
                         extreme_only_consider_used_nodes=False, \
                         proportional_p=True, extreme_proportional_p=True):

    if extreme_proportional_p:
        proportional_p = True
    if extreme_only_consider_used_nodes:
        only_consider_used_nodes = True

    if temporal:
        if (not extreme_proportional_p) and \
                not extreme_only_consider_used_nodes:
            if num_timestamps == "auto" or only_consider_used_nodes:
                timestamps = set()
                nodes = set()
                for edge in edges:
                    if num_timestamps == "auto":
                        timestamps.add(edge[2])
                    if only_consider_used_nodes:
                        nodes.add(edge[0])
                        nodes.add(edge[1])
                if num_timestamps == "auto":
                    num_timestamps = len(timestamps)
                    del timestamps
                if only_consider_used_nodes:
                    num_nodes = len(nodes)
                    del nodes

            num_possible_edges = (num_nodes * (num_nodes-1) * num_timestamps) /\
                                    (2 - int(directed))

            return {"temporal": True, "directed": directed, \
                    "num_timestamps": num_timestamps, \
                    "num_possible_edges": num_possible_edges, \
                    "num_edges": len(edges), \
                    "num_nodes": num_nodes}

        edges = [(t, a, b) for (a, b, t) in edges]
        edges.sort()
        edges = [(a, b, t) for (t, a, b) in edges]

        if only_consider_used_nodes:
            nodes = set()
        if extreme_only_consider_used_nodes:
            num_nodes_by_timestamp = []
        if extreme_proportional_p:
            edges_by_timestamp = []
            edges_in_timestamp = set()

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
                    edges_by_timestamp.append(edges_in_timestamp)
                    edges_in_timestamp = set()
            current_time = time

            if only_consider_used_nodes:
                nodes.add(source)
                nodes.add(target)

            if extreme_proportional_p:
                edges_in_timestamp.add(edge)

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
                        edges_by_timestamp.append(set())

    else:  # static (not temporal)
        if only_consider_used_nodes:
            nodes = set()
            for edge in edges:
                nodes.add(edge[0])
                nodes.add(edge[1])

def log2_ER_prob_of_matrix(nodes, edges, num_timestamps="auto", \
                           directed=False, temporal=False, \
                           only_consider_used_nodes=False, \
                           proportional_p=True):

    if temporal and num_timestamps == "auto":
        timestamps = set()
        for edge in edges:
            timetamps.add(edge[2])
        num_timestamps = len(timestamps)


def log2_num_automorphisms(nodes, edges, \
                           directed=False, temporal=False, \
                           only_consider_used_nodes=False):
    
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

    # Compute the log-number of automorphisms.
    (num_automorphisms, _) = __get_graph_info__(\
            nauty_run_nodes, edges, temporal, directed, \
            use_color_direction=True, \
            get_canon_order=False)

    log2_num_automorphisms = float(bigfloat.log2(num_automorphisms))

    if not only_consider_used_nodes:
        # Add automorphisms for unused nodes:
        for i in range(2, len(unused_nodes) + 1):
            log2_num_automorphisms += math.log(i, 2.0)

    return log2_num_automorphisms
