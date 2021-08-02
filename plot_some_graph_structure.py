from graph_types import GraphTypes
import matplotlib.pyplot as plt
import random
from structure_measures import measure_of_structure, log2_independent_edges_prob, log2_ER_prob

class GraphSequence:

    def __init__(self):
        self.__has_next__ = False
        self.sequence_type = None
        self.name = "A Graph Sequence"

    def set_with_list(self, l):
        self.sequence_type = "list"
        self.next_idx = 0
        self.sequence_list = l
        self.__has_next__ = self.next_idx < len(self.sequence_list)

    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def set_with_temporal_graph_file(self, filename, directed=True, \
                                     num_buckets=None):
        self.sequence_type = "temporal_file"

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

        edges = [(a, b, t) for (t, a, b) in \
                    sorted([(t, a, b) for (a, b, t) in edges])]

        self.full_sorted_edges = edges
        self.cumulative_edges = []
        self.nodes = nodes
        self.next_edge_idx = 0
        self.__has_next__ = self.next_edge_idx < len(self.full_sorted_edges)

        if num_buckets is not None:
            self.name = filename.split("/")[-1] + \
                                   (" - %d buckets" % num_buckets)
        else:
            self.name = filename.split("/")[-1] + " - no buckets"

    def set_window_sequence_with_temporal_file(self, filename, \
                                                     time_numbers_per_unit, \
                                                     unit_name, \
                                                     units_per_window, \
                                                     directed=True, \
                                                     start_offset_units=0, \
                                                     flatten_window=False, \
                                                     weight_repeats=False, \
                                                     windows_overlap=True):
        self.sequence_type = "windows_of_temporal_file"
        self.windows_overlap = windows_overlap
        self.flatten_window = flatten_window
        self.weight_repeats = weight_repeats

        (nodes, edges) = __read_edge_list__(filename, directed, True)
        # Sort timestamps.
        edges = [(a, b, t) for (t, a, b) in \
                    sorted([(t, a, b) for (a, b, t) in edges])]
        # Relabel timestamps with the basic unit.
        start_time = edges[0][2] + start_offset_units
        interval_end_time = start_time + time_numbers_per_unit
        replacement_timestamp = 1
        edges_lists = []
        edge_dict = {}  # Keep track of how often an edge appears in a timestamp
        nodes_sets = []
        nodes_set = set()
        for i in range(0, len(edges)):
            (a, b, t) = edges[i]
            while t >= interval_end_time:
                replacement_timestamp += 1
                interval_end_time += time_numbers_per_unit
                edges_lists.append([(a, b, t, w) for \
                                        (a, b, t), w in edge_dict.items()])
                edge_dict = {}
                if len(nodes_set) == 0:
                    print("Note: The %dth %s was empty." % \
                            (replacement_timestamp - 1, unit_name))
                nodes_sets.append(nodes_set)
                nodes_set = set()
            new_edge = (a, b, replacement_timestamp)
            if new_edge not in edge_dict:
                edge_dict[new_edge] = 0
            edge_dict[new_edge] += 1
            nodes_set.add(a)
            nodes_set.add(b)
        edges_lists.append([(a, b, t, w) for \
                               (a, b, t), w in edge_dict.items()])
        nodes_sets.append(nodes_set)

        self.full_edges_lists = edges_lists
        self.full_nodes_sets = nodes_sets
        self.units_per_window = units_per_window
        self.current_window_idx = 0
        self.last_window_idx = len(self.full_edges_lists) - units_per_window
        self.__has_next__ = self.current_window_idx <= self.last_window_idx

        if flatten_window:
            flattened_str = " flattened"
        else:
            flattened_str = ""
        if windows_overlap:
            overlapping_str = " overlapping"
        else:
            overlapping_str = ""
        if weight_repeats:
            weighted_str = " weighted"
        else:
            weighted_str = ""
            
        self.name = filename.split("/")[-1] + \
                        (" -%s%s%s %d %s windows" % \
                           (weighted_str, flattened_str, overlapping_str, \
                            units_per_window, unit_name))

    def has_next(self):
        return self.__has_next__

    def next(self):
        assert self.has_next()

        if self.sequence_type == "list":
            self.next_idx += 1
            self.__has_next__ = self.next_idx < len(self.sequence_list)
            return self.sequence_list[self.next_idx - 1]

        elif self.sequence_type == "temporal_file":
            t = self.full_sorted_edges[self.next_edge_idx][2]
            while self.next_edge_idx < len(self.full_sorted_edges) and \
                    self.full_sorted_edges[self.next_edge_idx][2] == t:
                self.cumulative_edges.append(\
                    self.full_sorted_edges[self.next_edge_idx])
                self.next_edge_idx += 1
            self.__has_next__ = self.next_edge_idx < len(self.full_sorted_edges)
            return (self.nodes, self.cumulative_edges)

        elif self.sequence_type == "windows_of_temporal_file":
            nodes = set()
            if self.flatten_window:
                edges = set()
            else:
                edges = []
            for i in range(self.current_window_idx, \
                           self.current_window_idx + self.units_per_window):
                nodes |= self.full_nodes_sets[i]
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                # TODO: ALLOW EDGE WEIGHTS IN CALCULATIONS!!!!!!!!!!!!
                if self.flatten_window:
                    edges |= set([(a, b, 1) for \
                                  (a, b, t, w) in self.full_edges_lists[i]])
                else:
                    edges += [(a, b, t) for \
                              (a, b, t, w) in self.full_edges_lists[i]]
            if self.windows_overlap:
                self.current_window_idx += 1
            else:
                self.current_window_idx += self.units_per_window
            self.__has_next__ = self.current_window_idx <= self.last_window_idx
            return (list(nodes), edges)

        

# graph_sequence should be a list of (nodes, edges) tuples.
def plot_graph_IC_sequence(graph_sequence, \
                           directed=False, temporal=False):
    if directed and temporal:
        graph_type = GraphTypes.TEMPORAL_DIRECTED
    elif temporal:
        graph_type = GraphTypes.TEMPORAL_UNDIRECTED
    elif directed:
        graph_type = GraphTypes.STATIC_DIRECTED
    else:
        graph_type = GraphTypes.STATIC_UNDIRECTED

    sequence_name = graph_sequence.get_name()

    gi = 0
    graph_indices = []
    ic = []
    while graph_sequence.has_next():
        (nodes, edges) = graph_sequence.next()
        if len(nodes) == 0:
            nodes = [1]
        gi += 1
        graph_indices.append(gi)

        num_timestamps = 1
        if temporal:
            if set_num_timestamps is not None:
                num_timestamps = set_num_timestamps
            else:
                timestamps = set()
                for (a, b, t) in edges:
                    timestamps.add(t)
                num_timestamps = len(timestamps)

        ic.append(-log2_ER_prob(nodes, edges, graph_type))

    plt.plot(graph_indices, ic, color="blue")
    plt.title("Info Content of %s" % sequence_name)
    plt.xlabel("Graph Sequence Index")
    plt.ylabel("Information Content")
    plt.savefig("figs/IC_of_%s.png" % sequence_name.replace(" ", "_"))
    plt.close()

def plot_graph_SM_sequence(graph_sequence, \
                           directed=False, temporal=False, \
                           add_nodes_edges_plot=False):
    if directed and temporal:
        graph_type = GraphTypes.TEMPORAL_DIRECTED
    elif temporal:
        graph_type = GraphTypes.TEMPORAL_UNDIRECTED
    elif directed:
        graph_type = GraphTypes.STATIC_DIRECTED
    else:
        graph_type = GraphTypes.STATIC_UNDIRECTED

    sequence_name = graph_sequence.get_name()

    gi = 0
    graph_indices = []
    real_graph = []
    rand_graph = []
    n = []
    m = []
    # sm = []
    while graph_sequence.has_next():
        (nodes, edges) = graph_sequence.next()
        if add_nodes_edges_plot:
            n.append(len(nodes))
            m.append(len(edges))
        if len(nodes) == 0:
            nodes = [1]
        gi += 1
        graph_indices.append(gi)
        (real, rand) = measure_of_structure([nodes], edges, graph_type, \
                                       all_timestamps="auto", \
                                       ER_try_count=20)
        real_graph.append(real)
        rand_graph.append(rand)
        plt.plot([gi, gi], [rand, real], color="gray")

    # plt.plot(graph_indices, sm, color="blue")
    plt.scatter(graph_indices, real_graph, color="blue")
    plt.scatter(graph_indices, rand_graph, color="red")
    plt.title("Structure Measure of %s" % sequence_name)
    plt.xlabel("Graph Sequence Index")
    plt.ylabel("Structure Measure")
    plt.savefig("figs/SM_of_%s.png" % sequence_name.replace(" ", "_"))
    plt.close()

    if add_nodes_edges_plot:
        plt.scatter(graph_indices, n, color="gray")
        plt.scatter(graph_indices, m, color="purple")
        plt.title("Nodes and Edges of %s" % sequence_name)
        plt.xlabel("Graph Sequence Index")
        plt.ylabel("|V| (gray), |E|, (purple)")
        plt.savefig("figs/Nodes_Edges_of_%s.png" % \
                        sequence_name.replace(" ", "_"))
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


if __name__ == "__main__":

    """
    list_GS = GraphSequence()
    list_GS.set_with_list(__triangles_sequence__(9))
    list_GS.set_name("Triangles Sequence")
    plot_graph_SM_sequence(list_GS, directed=False, temporal=False)
    list_GS.set_with_list(__triangles_sequence__(9))
    list_GS.set_name("Triangles Sequence")
    plot_graph_IC_sequence(list_GS, directed=False, temporal=False)

    list_GS.set_with_list(__triangles_sequence_with_all_nodes_always__(9))
    list_GS.set_name("Triangles Sequence with All Nodes")
    plot_graph_SM_sequence(list_GS, directed=False, temporal=False)
    list_GS.set_with_list(__triangles_sequence_with_all_nodes_always__(9))
    list_GS.set_name("Triangles Sequence with All Nodes")
    plot_graph_IC_sequence(list_GS, directed=False, temporal=False)

    rand_edge_addition = __random_edge_addition__(num_nodes=9*3, \
                                                   edges_per_iter=3, \
                                                   sequence_length=9, \
                                                   directed=False)
    list_GS.set_with_list(rand_edge_addition)
    list_GS.set_name("Randomly Adding 3 Edges Each Time")
    plot_graph_SM_sequence(list_GS, directed=False, temporal=False)
    list_GS.set_with_list(rand_edge_addition)
    list_GS.set_name("Randomly Adding 3 Edges Each Time")
    plot_graph_IC_sequence(list_GS, directed=False, temporal=False)

    list_GS.set_with_list(__binary_tree_sequence__(num_trees=7))
    list_GS.set_name("Binary Trees")
    plot_graph_SM_sequence(list_GS, directed=False, temporal=False)
    list_GS.set_with_list(__binary_tree_sequence__(num_trees=7))
    list_GS.set_name("Binary Trees")
    plot_graph_IC_sequence(list_GS, directed=False, temporal=False)
    """

    """
    file_GS = GraphSequence()
    file_GS.set_with_temporal_graph_file("datasets/college-temporal.g", \
                                         directed=True, \
                                         num_buckets=10)
    plot_graph_SM_sequence(file_GS, directed=True, temporal=True)

    file_GS.set_with_temporal_graph_file("datasets/eucore-temporal.g", \
                                         directed=True, \
                                         num_buckets=10)
    plot_graph_SM_sequence(file_GS, directed=True, temporal=True)
    """

    # TODO: Add support for edge weights in calcs

    # Day, Week resolution
    #
    # The -28561 makes the start time 12 AM on April 15th, 2004
    #   (That's Pacific Daylight Time)
    #   Choosing -10561 would be 5 AM instead.
    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60*60*24), \
        unit_name="days",
        units_per_window=7, \
        start_offset_units=-28561, \
        windows_overlap=False, \
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True, \
                           add_nodes_edges_plot=True)

    # Hour, Day resolution
    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60*60), \
        unit_name="hours",
        units_per_window=24, \
        start_offset_units=-28561, \
        windows_overlap=False, \
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True, \
                           add_nodes_edges_plot=True)

    # Hour, Hour resolution
    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60*60), \
        unit_name="hour",
        units_per_window=1, \
        start_offset_units=-28561, \
        windows_overlap=False, \
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True, \
                           add_nodes_edges_plot=True)

    """
    # Minute, Hour resolution
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60), \
        unit_name="minutes",
        units_per_window=60, \
        directed=True)
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True)
    """
