# init()
#
# has_next()
# next()
#
# set_name(name)
# get_name()
#
# set_with_list(l)
#   l should contain (nodes, edges) tuples
#
# set_with_temporal_graph_file(filename, directed=True, num_buckets=None)
#
# set_window_sequence_with_temporal_file(filename,
#                                        time_numbers_per_unit,
#                                        unit_name,
#                                        units_per_window,
#                                        directed=True,
#                                        start_offset_number=0,
#                                        flatten_window=False,
#                                        weight_repeats=False,
#                                        windows_overlap=True)
#
#   time_numbers_per_unit -- number of raw timestamp increments per unit of
#       time to be analyzed. E.g., set to 60 when analyzing minutes with raw
#       timestamps of seconds
#
#   unit_name -- e.g. "minutes"
#
#   units_per_window -- any integer >= 1
#
#   directed -- whether the edges in the graph are directed
#
#   start_offset_number -- any integer <= 0 -- If 0, the first timestamp is
#       considered the start of the first unit. Otherwise, the unit begins
#       at the first timestamp plus this number.
#
#       For example, if the first timestamp is 5 minutes and 23 seconds into a
#       game, (value: 323) but you want the minutes to start at actual minutes
#       set start_offset_number=-23 to start at the 5th minute or
#       start_offset_number=-323 to start at the beginning of the game.
#
#   flatten_window -- within a given window, squash all timestamps to a single
#       timestamp
#
#   weight_repeats -- when re-labeling timestamps, either via units_per_window,
#       flatten_window, or both, setting this to True makes the edges have
#       weights according to how often they were repeated in the unit or window.
#
#   windows_overlap -- If True, they each k-unit window shares k-1 units with
#       the previous window. Otherwise, each k-unit window is disjoint, back-to-
#       back.


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
                                                     start_offset_number=0, \
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
        start_time = edges[0][2] + start_offset_number
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
