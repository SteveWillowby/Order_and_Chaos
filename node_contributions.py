from graph_types import GraphTypes
from structure_measures import measure_of_structure

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

    basic_structure_amount = measure_of_structure([list(nodes)], edges, \
                                graph_type, \
                                all_timestamps="auto", \
                                ER_try_count=ER_try_count*10, \
                                node_colors=None, \
                                average_ER=False)

    node_contributions = []
    for node in nodes:
        print(node)
        relevant_edges = edges - set([(node, n) for n in neighbors[node]] + \
                                     [(n, node) for n in neighbors[node]])
        relevant_nodes = set([a for (a, b) in relevant_edges] + \
                             [b for (a, b) in relevant_edges])
        contribution = measure_of_structure([list(relevant_nodes)], \
                                            relevant_edges, \
                                            graph_type, \
                                            all_timestamps="auto", \
                                            ER_try_count=ER_try_count, \
                                            node_colors=None, \
                                            average_ER=False)
        node_contributions.append((contribution - basic_structure_amount, node))
    node_contributions.sort(reverse=True)
    return node_contributions

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
    nc = get_node_contributions(edges, directed=True, ER_try_count=100)
    print(nc)

    print("Start: %s" % str(nc[0]))
    print("End:   %s" % str(nc[-1]))
    for i in range(0, len(nc)):
        (contribution, node) = nc[i]
        if 1 <= node and node <= 64:
            print("At index %d we have %s" % (i, str((contribution, node))))
