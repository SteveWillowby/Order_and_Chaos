import bigfloat
import math

"""
log2_num_automorphisms(nodes, edges,
                       directed=False, temporal=False,
                       only_consider_used_nodes=False)
"""

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
