from graph_info_content import information_content, min_information_content_limit, max_information_content_limit
from graph_types import GraphTypes
import matplotlib.pyplot as plt
import random
from structure_measures import measure_of_structure


# graph_sequence should be a list of (nodes, edges) tuples.
def plot_graph_IC_sequence(graph_sequence, \
                           directed=False, temporal=False, \
                           sequence_name="A Graph Sequence", \
                           set_num_timestamps=None, \
                           normalize=True, \
                           proportional_p=False):
    graph_indices = [i for i in range(0, len(graph_sequence))]
    ic = []
    min_ic_limit = []
    max_ic_limit = []
    for (nodes, edges) in graph_sequence:

        num_timestamps = 1
        if temporal:
            if set_num_timestamps is not None:
                num_timestamps = set_num_timestamps
            else:
                timestamps = set()
                for (a, b, t) in edges:
                    timestamps.add(t)
                num_timestamps = len(timestamps)

        ic.append(information_content(nodes, edges, \
                        proportional_p=proportional_p, \
                        directed=directed, temporal=temporal, \
                        num_timestamps=num_timestamps, \
                        node_colors=None, only_consider_used_nodes=False, \
                        extreme_only_consider_used_nodes=False, \
                        extreme_proportional_p=False))

        min_ic_limit.append(\
            min_information_content_limit(len(nodes), \
                                          directed=directed, \
                                          num_timestamps=num_timestamps))
        max_ic_limit.append(\
            max_information_content_limit(len(nodes), \
                                          directed=directed, \
                                          num_timestamps=num_timestamps))
        # print("  Min: %f\n  IC: %f\n  Max: %f" % \
        #         (min_ic_limit[-1], ic[-1], max_ic_limit[-1]))

    if normalize:
        for i in range(0, len(ic)):
            if max_ic_limit[i] == min_ic_limit[i]:
                min_ic_limit[i] = 1
                max_ic_limit[i] = 1
                ic[i] = 1
            else:
                ic[i] = (ic[i] - min_ic_limit[i]) / \
                            (max_ic_limit[i] - min_ic_limit[i])
            min_ic_limit[i] = 0
            max_ic_limit[i] = 1

    plt.plot(graph_indices, min_ic_limit, color="red")
    plt.plot(graph_indices, max_ic_limit, color="red")
    plt.plot(graph_indices, ic, color="blue")
    plt.title("Info Content of %s" % sequence_name)
    plt.xlabel("Graph Sequence Index")
    plt.ylabel("Information Content")
    plt.savefig("figs/IC_of_%s.png" % sequence_name)
    plt.close()

def plot_graph_SM_sequence(graph_sequence, \
                           directed=False, temporal=False, \
                           sequence_name="A Graph Sequence"):
    if directed and temporal:
        graph_type = GraphTypes.TEMPORAL_DIRECTED
    elif temporal:
        graph_type = GraphTypes.TEMPORAL_UNDIRECTED
    elif directed:
        graph_type = GraphTypes.STATIC_DIRECTED
    else:
        graph_type = GraphTypes.STATIC_UNDIRECTED

    graph_indices = [i for i in range(0, len(graph_sequence))]
    sm = []
    for (nodes, edges) in graph_sequence:
        sm.append(measure_of_structure([nodes], edges, graph_type, \
                                       all_timestamps="auto", \
                                       ER_try_count=10))

    plt.plot(graph_indices, sm, color="blue")
    plt.title("Structure Measure of %s" % sequence_name)
    plt.xlabel("Graph Sequence Index")
    plt.ylabel("Structure Measure")
    plt.savefig("figs/SM_of_%s.png" % sequence_name)
    plt.close()

# Reads the edge list with potentially duplicated edges and returns two lists:
#   nodes and edges
def __read_edge_list__(filename, directed, temporal):
    f = open(filename, "r")
    edges = set()
    nodes = set()
    self_loops = 0
    for line in f:
        line = line.strip()
        line = line.split(" ")
        assert len(line) == 2 + int(temporal)
        source = int(line[0])
        target = int(line[1])

        nodes.add(source)
        nodes.add(target)

        if source == target:
            self_loops += 1
            continue

        if directed:
            edge = [source, target]
        else:
            edge = [min(source, target), max(source, target)]

        if temporal:
            timestamp = line[2]
            if timestamp.isnumeric():
                timestamp = int(timestamp)
            edge.append(timestamp)
        edges.add(tuple(edge))

    f.close()
    print("Input file had %d self-loops, all of which (if any) were removed." %\
        self_loops)
    return (list(nodes), list(edges))

def __triangles_sequence__(sequence_length):
    sequence = []
    nodes_list = []
    edges_list = []
    for i in range(0, sequence_length):
        nodes_list.append(i*3)
        nodes_list.append(i*3 + 1)
        nodes_list.append(i*3 + 2)
        edges_list.append((i*3, i*3 + 1))
        edges_list.append((i*3 + 1, i*3 + 2))
        edges_list.append((i*3 + 2, i*3))
        sequence.append((list(nodes_list), list(edges_list)))
    return sequence

def __triangles_sequence_with_all_nodes_always__(sequence_length):
    sequence = []
    nodes_list = []
    edges_list = []
    for i in range(0, sequence_length):
        nodes_list.append(i*3)
        nodes_list.append(i*3 + 1)
        nodes_list.append(i*3 + 2)
        edges_list.append((i*3, i*3 + 1))
        edges_list.append((i*3 + 1, i*3 + 2))
        edges_list.append((i*3 + 2, i*3))
        sequence.append((nodes_list, list(edges_list)))
    return sequence

def __random_edge_addition__(num_nodes, edges_per_iter, sequence_length, \
                             directed=False):
    sequence = []
    nodes = [i for i in range(0, num_nodes)]
    edges_set = set()
    for _ in range(0, sequence_length):
        num_added = 0
        while num_added < edges_per_iter:
            a = random.randint(0, num_nodes - 1)
            b = random.randint(0, num_nodes - 1)
            if a == b:
                continue
            if (a, b) not in edges_set and \
                    (directed or (b, a) not in edges_set):
                edges_set.add((a, b))
                num_added += 1
        sequence.append((nodes, list(edges_set)))
    return sequence

def __binary_tree_sequence__(num_trees):
    sequence = []
    nodes_list = [1]
    edges_list = []
    num_nodes = 1
    sequence.append(([1], []))
    for _ in range(0, num_trees - 1):
        new_num_nodes = num_nodes * 2 + 1
        for n in range(num_nodes + 1, new_num_nodes + 1):
            nodes_list.append(n)
            edges_list.append((int(n / 2), n))
        sequence.append((list(nodes_list), list(edges_list)))
        num_nodes = new_num_nodes
    return sequence

# Note, can change the order of elements within `edges`.
#
# O(|edges| log |edges| + |timestamps| * target_num_buckets)
def __bucket_temporal_edges__(edges, target_num_buckets):
    assert type(edges) is list
    for i in range(0, len(edges)):
        edges[i] = (edges[i][2], edges[i][0], edges[i][1])
    edges.sort()

    possible_bucket_start_positions = []

    current_time = None
    for i in range(0, len(edges)):
        (t, a, b) = edges[i]
        if current_time is None or t != current_time:
            current_time = t
            possible_bucket_start_positions.append(i)
        edges[i] = (a, b, t)

    if target_num_buckets > len(possible_bucket_start_positions):
        print(("Note! It is not possible to have %d " % target_num_buckets) + \
              ("with this edge list - only %d timestamps." % \
                len(possible_bucket_start_positions)))
        buckets = [(possible_bucket_start_positions[i], \
                    possible_bucket_start_positions[i + 1]) for \
                      i in range(0, len(possible_bucket_start_positions)-1)] + \
                  [(possible_bucket_start_positions[-1], len(edges))]
        edge_buckets = []
        for (start, end) in buckets:
            edge_buckets.append(edges[start:end])
        return edge_buckets

    possible_bucket_start_positions.append(len(edges))
    target_edges_per_bucket = int(len(edges) / target_num_buckets)
    pbsp = possible_bucket_start_positions
    tbs = target_edges_per_bucket
    # The first target-minus-1 buckets are as small as possible.
    candidate_bucketing = [(i, i + 1) for i in range(0, target_num_buckets - 1)]
    # The last bucket uses the rest of the edges.
    candidate_bucketing.append((target_num_buckets - 1, len(pbsp) - 1))
    # Sum up the differences between bucket size and target bucket size.
    #   Lower is better.
    bucket_score = sum([abs(tbs - (pbsp[end] - pbsp[start])) for \
                            (start, end) in candidate_bucketing])
    updated_bucket_score = bucket_score  # used for testing
    done = False
    while not done:
        done = True
        # Consider moving any bucket start position but the first one.
        for i in range(1, target_num_buckets):
            idx = target_num_buckets - i
            (start, end) = candidate_bucketing[idx]
            if start == end - 1:
                continue  # can't move start
            # See if moving bucket[idx] start pos further to the right helps.
            (prev_start, prev_end) = candidate_bucketing[idx - 1]
            assert prev_end == start
            current_idx_score_contribution = \
                abs(tbs - (pbsp[end] - pbsp[start])) + \
                abs(tbs - (pbsp[prev_end] - pbsp[prev_start]))
            new_idx_score_contribution = \
                abs(tbs - (pbsp[end] - pbsp[start + 1])) + \
                abs(tbs - (pbsp[prev_end + 1] - pbsp[prev_start]))
            if new_idx_score_contribution <= current_idx_score_contribution:
                candidate_bucketing[idx] = (start + 1, end)
                candidate_bucketing[idx - 1] = (prev_start, prev_end + 1)
                updated_bucket_score += new_idx_score_contribution - \
                                        current_idx_score_contribution
                done = False
                break

    bucket_score = sum([abs(tbs - (pbsp[end] - pbsp[start])) for \
                            (start, end) in candidate_bucketing])
    print("Average deviation from target bucket size: %f" % \
            (bucket_score / target_num_buckets))

    assert updated_bucket_score == bucket_score
    
    edge_buckets = []
    for (start, end) in candidate_bucketing:
        edge_buckets.append(edges[pbsp[start]:pbsp[end]])
    return edge_buckets


def __plot_temporal_sequence__(filename, directed=True, num_buckets=None, \
                               use_all_nodes=True):
    (nodes, edges) = __read_edge_list__(filename, directed, True)
    edges_by_timestamp = {}
    timestamps = []
    for (a, b, t) in edges:
        if a == b:
            continue
        if t not in edges_by_timestamp:
            edges_by_timestamp[t] = set()
            timestamps.append(t)
        if directed:
            edges_by_timestamp[t].add((a, b))
        else:
            edges_by_timestamp[t].add((min(a, b), max(a, b)))
    timestamps.sort()
    if num_buckets is not None:
        # Get buckets, then replace timestamps.
        bucketed_edges = __bucket_temporal_edges__(edges, num_buckets)
        edges = []
        for i in range(0, len(bucketed_edges)):
            bucket_set = set()
            for (a, b, _) in bucketed_edges[i]:
                if directed:
                    bucket_set.add((a, b, i))
                else:
                    bucket_set.add((min(a, b), max(a, b), i))
            edges += list(bucket_set)

    graph_sequence = []
    edge_sub_list = []
    current_t = None

    for (a, b, t) in edges:
        if current_t is not None and current_t != t:
            graph_sequence.append((nodes, list(edge_sub_list)))
        edge_sub_list.append((a, b, t))
        current_t = t
    graph_sequence.append((nodes, edge_sub_list))

    if num_buckets is not None:
        num_timestamps_to_use = len(graph_sequence)
    else:
        num_timestamps_to_use = len(timestamps)

    plot_graph_SM_sequence(graph_sequence, directed=directed, temporal=True, \
                           sequence_name=(filename.split("/")[-1] + \
                                          (" - %d buckets" % num_buckets)))
                           # set_num_timestamps=None)  # this just for IC

if __name__ == "__main__":

    plot_graph_SM_sequence(__triangles_sequence__(9), \
                        directed=False, temporal=False, \
                        sequence_name="Triangles Sequence")
    plot_graph_SM_sequence(__triangles_sequence_with_all_nodes_always__(9), \
                        directed=False, temporal=False, \
                        sequence_name="Triangles Sequence with All Nodes")
    plot_graph_SM_sequence(__random_edge_addition__(num_nodes=9*3, \
                                                 edges_per_iter=3, \
                                                 sequence_length=9, \
                                                 directed=False), \
                        directed=False, temporal=False, \
                        sequence_name="Randomly Adding 3 Edges Each Time")

    plot_graph_SM_sequence(__binary_tree_sequence__(num_trees=7), \
                        directed=False, temporal=False, \
                        sequence_name="Binary Trees")

    __plot_temporal_sequence__("datasets/college-temporal.g", \
                               directed=True, num_buckets=10)
    __plot_temporal_sequence__("datasets/eucore-temporal.g", \
                               directed=True, num_buckets=10)

    # __plot_temporal_sequence__("datasets/college-temporal.g", \
    #                            directed=True, num_buckets=100)
    # __plot_temporal_sequence__("datasets/eucore-temporal.g", \
    #                            directed=True, num_buckets=100)

    # __plot_temporal_sequence__("datasets/wiki-en-additions.g", \
    #                            directed=True, num_buckets=10)
    # __plot_temporal_sequence__("datasets/wiki-en-additions.g", \
    #                            directed=True, num_buckets=100)
