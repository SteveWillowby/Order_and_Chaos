import bigfloat
from graph_sequencer_utility import GraphSequence
from graph_types import GraphTypes
from look_for_stars import find_stars
import math
import matplotlib.pyplot as plt
import random
from structure_measures import measure_of_structure, log2_num_automorphisms, __create_ER_graph__



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


def plot_graph_percentile_sequence(graph_sequence, \
                                   directed=True, temporal=True, \
                                   num_ER_graphs=1999, \
                                   num_timestamps="auto", \
                                   agg_NM="median", \
                                   show_plot=True):
    if directed and temporal:
        graph_type = GraphTypes.TEMPORAL_DIRECTED
    elif temporal:
        graph_type = GraphTypes.TEMPORAL_UNDIRECTED
    elif directed:
        graph_type = GraphTypes.STATIC_DIRECTED
    else:
        graph_type = GraphTypes.STATIC_UNDIRECTED

    old_context = bigfloat.getcontext()
    bf_context = bigfloat.Context(precision=2000, \
                                  emax=1000000, emin=-1000000)
    bigfloat.setcontext(bf_context)

    log2_na = []
    n = []
    m = []
    if num_timestamps == "auto":
        first = True
    elif temporal:
        T = num_timestamps
    else:
        T = None
    while graph_sequence.has_next():
        (nodes, edges) = graph_sequence.next()
        if len(edges) == 0 and len(nodes) == 0:
            nodes = [1]
        elif len(nodes) == 0:
            raise ValueError("Error! Cannot have edges without nodes!")
        n.append(len(nodes))
        m.append(len(edges))

        if temporal and num_timestamps == "auto":
            if first:
                T = len(set([edge[2] for edge in edges]))
                first = False
            else:
                T = max(T, len(set([edge[2] for edge in edges])))

        log2_na.append(log2_num_automorphisms(nodes, edges, \
                       directed=directed, temporal=temporal,
                       only_consider_used_nodes=False, \
                       node_colors=None))

    if agg_NM == "median":
        N = sorted(n)[int(len(n) / 2)]
        M = sorted(m)[int(len(m) / 2)]
    elif agg_NM == "mean":
        N = int(float(sum(n)) / len(n) + 0.5)
        M = int(float(sum(m)) / len(m) + 0.5)
    else:
        raise ValueError(("Error! Unknown aggregation method '%s'" % agg_NM) + \
                         " for |V|, |E|. Enter 'median' or 'mean'.")

    ER_log2_na = []
    for _ in range(0, num_ER_graphs):
        ER_graph = __create_ER_graph__(\
                        [[j for j in range(0, N)]], M, graph_type, \
                        num_timestamps_if_relevant=T, \
                        node_timestamps_if_relevant=None)
        (ER_nodes, ER_edges) = (ER_graph[0][0], ER_graph[1])

        ER_log2_na.append(log2_num_automorphisms(ER_nodes, ER_edges, \
                          directed=directed, temporal=temporal,
                          only_consider_used_nodes=False, \
                          node_colors=None))

    bigfloat.setcontext(old_context)

    as_dict = {}
    for val in ER_log2_na:
        if val not in as_dict:
            as_dict[val] = 0
        as_dict[val] += 1
    ER_log2_na = sorted([(val, count) for val, count in as_dict.items()])
    log2_na = sorted([(log2_na[i], i) for i in range(0, len(log2_na))])
    real_idx = 0
    ER_idx = 0
    cumulative_count = 0
    percentiles = [None for i in range(0, len(log2_na))]
    prev_ER_val = 0.0
    ER_percentiles = []
    while real_idx < len(log2_na):
        (ER_val, count) = ER_log2_na[ER_idx]
        (real_val, sequence_idx) = log2_na[real_idx]
        if real_val == ER_val:
            # Set the percentile as the "middle" of the set of matching values..
            percentile = 100 * (cumulative_count + (count + 1) / 2.0) / \
                                    float(num_ER_graphs + 1)
            percentiles[sequence_idx] = percentile
            real_idx += 1
        elif prev_ER_val < real_val and real_val < ER_val:
            # Sandwich between two ER values
            percentile = 100 * (cumulative_count + 1) / \
                                    float(num_ER_graphs + 1)
            percentiles[sequence_idx] = percentile
            real_idx += 1
        elif real_val == prev_ER_val:
            assert float(real_val) == 0.0
            percentiles[sequence_idx] = 0.0
            real_idx += 1
        else:
            assert real_val > ER_val
            if ER_idx == len(ER_log2_na) - 1:
                percentiles[sequence_idx] = 100.0
                real_idx += 1
            else:
                ER_percentile = 100 * (cumulative_count + (count / 2.0)) / \
                                            num_ER_graphs
                ER_percentiles.append(ER_percentile)
                                    
                ER_idx += 1
                prev_ER_val = ER_val
                cumulative_count += count

    # Data collected. Now plot.
    sequence_name = graph_sequence.get_name()
    plt.scatter([0 for _ in ER_percentiles], ER_percentiles, color="red")
    plt.plot([i for i in range(1, len(percentiles) + 1)], percentiles, \
                color="blue")
    plt.title("Randomness Percentile of %s" % sequence_name)
    plt.xlabel("Graph Sequence Index")
    plt.ylabel("Randomness Percentile for (N,M,T): (%d, %d, %d)" % (N, M, T))
    plt.savefig("figs/RP_of_%s.png" % sequence_name.replace(" ", "_"))
    if show_plot:
        plt.show()
    plt.close()

# `window_size` is only relevant if using windows from the
#   graph_sequencer_utility. It ensures that ER graphs use all the timestamps
#   available in the window. DO NOT use in combination with flattened graphs.
def plot_graph_SM_sequence(graph_sequence, \
                           directed=False, temporal=False, \
                           add_nodes_edges_plot=False, \
                           graph_flattened=False, \
                           windows_overlap=False, \
                           show_plot=False, \
                           include_star_amt=False, \
                           include_pair_amt=False, \
                           relative=True, \
                           window_size=None, \
                           specified_NE=None):

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
    if not relative:
        real_graph = []
        rand_graph = []
    n = []
    m = []
    if relative:
        sm = []
    if include_star_amt or include_pair_amt:
        star_sm = []
        pair_sm = []
    while graph_sequence.has_next():
        (nodes, edges) = graph_sequence.next()
        if add_nodes_edges_plot:
            n.append(len(nodes))
            m.append(len(edges))
        if len(nodes) == 0:
            nodes = [1]
        gi += 1
        graph_indices.append(gi)
        if window_size is None:
            timestamps = "auto"
        else:
            # the graph sequencer 1-indexes timestamps
            if graph_flattened:
                timestamps = [1]
            elif windows_overlap:
                timestamps = [gi + i for i in range(0, window_size)]
            else:
                timestamps = [(gi - 1) * window_size + i for \
                                    i in range(1, window_size + 1)]

            for edge in edges:
                assert edge[-1] in timestamps

        (real, rand) = measure_of_structure([nodes], edges, graph_type, \
                                        all_timestamps=timestamps, \
                                        ER_try_count=1000, \
                                        average_ER=False)
        if float(rand) != 0.0:
            print(("Did not obtain rigid graph for %d edges" % len(edges)) + \
                  (" spread over %d nodes" % len(nodes)))
        if relative:
            sm.append(real - rand)
            plt.plot([gi, gi], [real - rand, 0.0], color="lightgray")
        else:
            real_graph.append(real)
            rand_graph.append(rand)
            plt.plot([gi, gi], [rand, real], color="lightgray")
        if include_star_amt or include_pair_amt:
            (stars, partial_stars, pairs) = \
                find_stars(edges, directed=directed, temporal=temporal)
            s_au = 0.0
            for (size, _) in stars + partial_stars:
                for value in range(2, size + 1):
                    s_au += math.log(value, 2.0)
            pair_au = 0.0
            for value in range(2, len(pairs) + 1):
                pair_au += math.log(value, 2.0)

            if relative and real > 0:
                # The % of the log2-difference between real and rand as
                #   contributed by stars, pairs
                star_sm.append((s_au / real) * (real - rand))
                pair_sm.append((pair_au / real) * (real - rand))
            else:
                star_sm.append(s_au)
                pair_sm.append(pair_au)

    if relative:
        plt.plot(graph_indices, sm, color="blue")
    else:
        plt.plot(graph_indices, real_graph, color="blue")
        plt.plot(graph_indices, rand_graph, color="brown")
    if include_star_amt:
        plt.plot(graph_indices, star_sm, color="green")
    if include_pair_amt:
        plt.plot(graph_indices, pair_sm, color="orange")

    plt.title("%sStructure Measure of %s" % \
        (["", "Rel. "][int(relative)], sequence_name))
    plt.xlabel("Graph Sequence Index")
    plt.ylabel("Structure Measure")
    plt.savefig("figs/%sSM_of_%s.png" % (\
        ["", "relative_"][int(relative)], sequence_name.replace(" ", "_")))
    if show_plot:
        plt.show()
    plt.close()

    if add_nodes_edges_plot:
        plt.plot(graph_indices, n, color="gray")
        plt.plot(graph_indices, m, color="purple")
        plt.title("Nodes and Edges of %s" % sequence_name)
        plt.xlabel("Graph Sequence Index")
        plt.ylabel("|V| (gray), |E|, (purple)")
        plt.savefig("figs/Nodes_Edges_of_%s.png" % \
                        sequence_name.replace(" ", "_"))
        if show_plot:
            plt.show()
        plt.close()


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

def squashed_email_graph_filename(email_graph_filename, timestamp_range=None):
    f = open(email_graph_filename, "r")
    lines = f.readlines()
    f.close()
    tmp_filename = "/tmp/" + email_graph_filename.split("/")[-1]
    f = open(tmp_filename, "w")
    for i in range(0, len(lines)):
        line = lines[i].strip().split(" ")
        sender = line[0]
        timestamp = line[-1]
        if timestamp_range is not None and \
                (int(timestamp) < timestamp_range[0] or \
                 int(timestamp) > timestamp_range[1]):
            continue
        for j in range(1, len(line) - 1):
            recipient = line[j]
            if recipient != sender:
                f.write("%s %s %s" % (sender, recipient, timestamp))
                if j < len(line) - 2 or i < len(lines) - 1:
                    f.write("\n")
    f.close()
    return tmp_filename

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

    Jan_1_1980 = 62483932800 + 8 * 60 * 60
    Jan_1_2000 = int(365.25 * 20) * 24 * 60 * 60 + Jan_1_1980
    Jan_1_2001 = int(365.25 * 20 + 366) * 24 * 60 * 60 + Jan_1_1980
    Jan_1_2003 = (365 * 2) * 24 * 60 * 60 + Jan_1_2001
    squashed = squashed_email_graph_filename(\
                    "datasets/enron_dataset/exclusive_emails.txt", \
                    timestamp_range=[Jan_1_2001, Jan_1_2003])
    # Hour, Day resolution
    flatten = False
    overlap = True
    units_per_window = 24
    if flatten:
        window_size = None
    else:
        window_size = units_per_window

    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename=squashed, \
        time_numbers_per_unit=(60*60), \
        unit_name="hour", \
        units_per_window=units_per_window, \
        start_offset_number=0, \
        windows_overlap=overlap, \
        flatten_window=flatten, \
        weight_repeats=False, \
        directed=True)

    """
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True, \
                           add_nodes_edges_plot=True, \
                           graph_flattened=flatten, \
                           windows_overlap=overlap, \
                           show_plot=True, \
                           include_star_amt=True, \
                           include_pair_amt=True, \
                           relative=True, \
                           window_size=window_size)
    """
    plot_graph_percentile_sequence(graph_sequence=window_GS, \
                                   directed=True, temporal=True, \
                                   num_ER_graphs=3999, \
                                   num_timestamps=7, \
                                   agg_NM="median", \
                                   show_plot=True)

    # TODO: Add support for edge weights in calcs


    """
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
        start_offset_number=-28561, \
        windows_overlap=False, \
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True, \
                           add_nodes_edges_plot=True, \
                           graph_flattened=True, \
                           show_plot=False, \
                           include_star_amt=True, \
                           include_pair_amt=True, \
                           relative=True, \
                           window_size=None)

    # Hour, Day resolution
    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60*60), \
        unit_name="hours",
        units_per_window=24, \
        start_offset_number=-28561, \
        windows_overlap=False, \
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True, \
                           add_nodes_edges_plot=True, \
                           graph_flattened=True, \
                           show_plot=False, \
                           include_star_amt=True, \
                           include_pair_amt=True, \
                           relative=True, \
                           window_size=None)

    # Hour, Hour resolution
    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60*60), \
        unit_name="hour",
        units_per_window=1, \
        start_offset_number=-28561, \
        windows_overlap=False, \
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)
    plot_graph_SM_sequence(window_GS, directed=True, temporal=True, \
                           add_nodes_edges_plot=True, \
                           graph_flattened=True, \
                           show_plot=False, \
                           include_star_amt=True, \
                           include_pair_amt=True, \
                           relative=True, \
                           window_size=None)
    """

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
