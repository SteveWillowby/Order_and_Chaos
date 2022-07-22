import bigfloat
from graph_types import GraphTypes
from nauty_session import NautyTracesSession
import networkx as nx
from structure_measures import measure_of_structure
import sys

def get_node_contributions(edges, directed, ER_try_count=100):
    if directed:
        graph_type = GraphTypes.STATIC_DIRECTED
    else:
        graph_type = GraphTypes.STATIC_UNDIRECTED

    if type(edges) is not set:
        edges = set(edges)
    neighbors = {}
    nodes = set()
    for (a, b) in edges:
        if a not in nodes:
            neighbors[a] = set()
            nodes.add(a)
        if b not in nodes:
            neighbors[b] = set()
            nodes.add(b)
        neighbors[a].add(b)
        neighbors[b].add(a)

    print("Formatted graph")
    sys.stdout.flush()

    (real, rand) = measure_of_structure([list(nodes)], edges, \
                                graph_type, \
                                all_timestamps="auto", \
                                ER_try_count=ER_try_count*2, \
                                node_colors=None, \
                                average_ER=False)
    basic_structure_amount = real - rand
    print("Basic Structure: %f" % basic_structure_amount)
    sys.stdout.flush()

    node_contributions = []
    for node in nodes:
        relevant_edges = edges - set([(node, n) for n in neighbors[node]] + \
                                     [(n, node) for n in neighbors[node]])
        relevant_nodes = set([a for (a, b) in relevant_edges] + \
                             [b for (a, b) in relevant_edges])
        (real, rand) = measure_of_structure([list(relevant_nodes)], \
                                            relevant_edges, \
                                            graph_type, \
                                            all_timestamps="auto", \
                                            ER_try_count=ER_try_count, \
                                            node_colors=None, \
                                            average_ER=False)
        node_contributions.append(((real - rand) - basic_structure_amount,node))
        print("%d %f" % (node, node_contributions[-1][0]))
        sys.stdout.flush()
    node_contributions.sort(reverse=True)
    return node_contributions

def get_centered_aut_logs(nodes, edges, directed):
    if type(nodes) is not list:
        nodes = list(nodes)

    extra_nodes = []
    extra_edges = []
    extra_highlights = []
    if directed:
        extra_highlights = [[], []]
        for (s, t) in edges:
            extra_nodes.append((s, t, 1))
            extra_nodes.append((s, t, 2))
            extra_highlights[0].append((s, t, 1))
            extra_highlights[1].append((s, t, 2))
            extra_edges.append((s, (s, t, 1)))
            extra_edges.append(((s, t, 1), (s, t, 2)))
            extra_edges.append(((s, t, 2), t))

    edges += extra_edges
    full_nodes = nodes + extra_nodes

    start_graph = nx.Graph()
    for node in full_nodes:
        start_graph.add_node(node)
    for (a, b) in edges:
        start_graph.add_edge(a, b)

    
    old_context = bigfloat.getcontext()
    bf_context = bigfloat.Context(precision=2000, \
                                  emax=1000000, emin=-1000000)
    bigfloat.setcontext(bf_context)

    NTSession = NautyTracesSession(start_graph, mode="Traces", \
                                   sparse=True, \
                                   allow_edits=False)
    na = NTSession.get_num_automorphisms()
    ao = NTSession.get_automorphism_orbits()
    NTSession.complete()
    na = na.get()
    ao = ao.get()
    print(na)
    ao_excluding_edge_nodes = []
    for o in ao:
        if o[0] in nodes:
            ao_excluding_edge_nodes.append(o)
    print(len(ao_excluding_edge_nodes))
    exit()

    log_auts = []
    batch_size = 150
    num_batches = int((len(nodes) + batch_size - 1) / batch_size)
    print("For %d nodes we have %d batches of size %d (or smaller)" % \
            (len(nodes), num_batches, batch_size))
    for i in range(0, num_batches):
        print("Beginning batch %d" % (i + 1))
        sys.stdout.flush()
        start = batch_size * i
        end = min(len(nodes), (i + 1) * batch_size)
        NTSession = NautyTracesSession(start_graph, mode="Traces", \
                                       sparse=True, \
                                       allow_edits=False)
        for j in range(start, end):
            node = nodes[j]
            if j == 0:
                highlights = [[node], nodes[1:]]
            elif j == len(nodes) - 1:
                highlights = [[node], nodes[:-1]]
            else:
                highlights = [[node], nodes[:j] + nodes[j+1:]]

            highlights += extra_highlights
            NTSession.set_colors_by_highlights(highlights)
            log_auts.append(NTSession.get_num_automorphisms())

        NTSession.complete()
    log_auts = [(float(bigfloat.log2(log_auts[i].get())), nodes[i]) for \
                    i in range(0, len(nodes))]

    bigfloat.setcontext(old_context)
    log_auts.sort()
    values_of_specials = []
    for (log_na, node) in log_auts:
        print("%d %s" % (node, log_na))
        if node <= 64:
            values_of_specials.append((node, log_na))
    print("-------------------")
    print("Specials:")
    print(values_of_specials)

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

if __name__ == "__main__":
    (nodes, edges) = __read_edge_list__("datasets/twitter_subsample_network.g",\
                                        directed=True, temporal=False)
    get_centered_aut_logs(nodes, edges, directed=True)
    exit(0)

    print("Loaded graph.")
    nc = get_node_contributions(edges, directed=True, ER_try_count=30)
    print(nc)

    print("Start: %s" % str(nc[0]))
    print("End:   %s" % str(nc[-1]))
    for i in range(0, len(nc)):
        (contribution, node) = nc[i]
        if 1 <= node and node <= 64:
            print("At index %d we have %s" % (i, str((contribution, node))))
