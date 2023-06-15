import sys

if __name__ == "__main__":
    node_file = sys.argv[1]
    graph_file = sys.argv[2]
    result_file = sys.argv[3]

    print("THIS ASSUMES UNDIRECTED GRAPH!!!")

    f = open(graph_file, "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip().split(" ") for l in lines if len(l.strip()) > 0]
    edges = [(int(x[0]), int(x[1])) for x in lines]
    edges = [(min(a, b), max(a, b)) for (a, b) in edges]

    f = open(node_file, "r")
    lines = f.readlines()
    f.close()

    nodes = set([int(l.strip()) for l in lines if len(l.strip()) > 0])

    edge_nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])
    assert nodes == edge_nodes
    node_map = sorted(list(nodes))
    node_map = {i: node_map[i] for i in range(0, len(node_map))}

    f = open(result_file, "r")
    lines = f.readlines()
    f.close()

    i = 0
    while lines[i].strip() != "With edges:":
        i += 1

    l = lines[i + 1]
    l = l.split("), (")
    l = [x.replace(")", "").replace("(", "") for x in l]
    l = [x.split(", ") for x in l]
    flipped_edges = [(node_map[int(x[0])], node_map[int(x[1])]) for x in l]
    flipped_edges = [(min(a, b), max(a, b)) for (a, b) in flipped_edges]

    flipped_nodes = set([a for (a, b) in flipped_edges] + \
                        [b for (a, b) in flipped_edges])

    assert flipped_nodes - nodes == set()

    edges = set(edges)
    meta_node = max(nodes) + 1
    flipped_edges = set(flipped_edges)

    noise = flipped_edges
    origin_graph = (edges - flipped_edges) | (flipped_edges - edges)

    origin_nodes = set([a for (a, b) in origin_graph] + \
                       [b for (a, b) in origin_graph])

    f = open("merged_graph.csv", "w")
    f.write("Source,Target,Weight\n")
    for (a, b) in sorted(list(origin_graph)):
        f.write("%d,%d,%d\n" % (a, b, 1))

    for node in sorted(list(nodes - origin_nodes)):
        f.write("%d,%d,%d\n" % (node, meta_node, 1))
    if len(nodes - edge_nodes) > 0:
        for node in sorted(list(nodes & origin_nodes)):
            f.write("%d,%d,%d\n" % (node, meta_node, 1))
            break
    f.close()

    print("There were %d additions and %d removals." % \
            (len(flipped_edges - edges), len(flipped_edges & edges)))
    print("There are %d flipped edges and %d disconnected nodes." \
            % (len(flipped_edges), len(nodes - origin_nodes)))
    print("%d of flipped edges involve a disconnected node." % \
            sum([1 for (a, b) in flipped_edges if a in (nodes - origin_nodes) or \
                                                  b in (nodes - origin_nodes)]))

    f = open("merged_graph_with_noise.csv", "w")
    f.write("Source,Target,Weight\n")
    for (a, b) in sorted(list(origin_graph - noise)):
        f.write("%d,%d,%d\n" % (a, b, 1))
    for (a, b) in sorted(list(noise)):
        if (a, b) in origin_graph:
            f.write("%d,%d,%d\n" % (a, b, 1))
        else:
            f.write("%d,%d,%d\n" % (a, b, 1))

    for node in sorted(list(nodes - (origin_nodes | flipped_nodes))):
        f.write("%d,%d,%d\n" % (node, meta_node, 1))
    if len(nodes - origin_nodes) > 0:
        for node in sorted(list(nodes & (origin_nodes | flipped_nodes))):
            f.write("%d,%d,%d\n" % (node, meta_node, 1))
            break
    f.close()
