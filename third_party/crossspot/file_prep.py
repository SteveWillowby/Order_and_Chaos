import sys

if __name__ == "__main__":
    filename = sys.argv[1]
    undirected = len(sys.argv) > 2

    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    lines = [s.strip().replace(" ", ",") + ",1\n" for s in lines]

    edges = []
    node_to_node = {}

    f = open("prepared.csv", "w")
    for line in lines:

        s = line.split(",")
        n1 = int(s[0])
        n2 = int(s[1])
        if n1 not in node_to_node:
            node_to_node[n1] = []
        node_to_node[n1].append(n2)

        edges.append((n1, n2))
        if undirected:
            if n2 not in node_to_node:
                node_to_node[n2] = []
            node_to_node[n2].append(n1)

        f.write(line)
    f.close()

    edge_to_node = sorted(edges)
    edge_to_node = {edge_to_node[i]: i for i in range(0, len(edge_to_node))}

    n_to_n_reverse = {}
    for a, s in node_to_node.items():
        for b in s:
            if b not in n_to_n_reverse:
                n_to_n_reverse[b] = []
            n_to_n_reverse[b].append(a)

    new_edges = []

    edges = set(edges)

    for (a, b) in edges:
        for c in node_to_node[a]:
            if (a, c) in edges and c != b:
                (x, y) = (edge_to_node[(a, b)], edge_to_node[(a, c)])
                if undirected:
                    new_edges.append((min(x, y), max(x, y)))
                else:
                    new_edges.append((x, y))
        if a in n_to_n_reverse:
            for c in n_to_n_reverse[a]:  # Anything pointing to a
                if (c, a) in edges and c != b:  # Might not be present
                    (x, y) = (edge_to_node[(a, b)], edge_to_node[(c, a)])
                    if undirected:
                        new_edges.append((min(x, y), max(x, y)))
                    else:
                        new_edges.append((x, y))

        if b in node_to_node:
            for c in node_to_node[b]:
                if (b, c) in edges and a != c:
                    (x, y) = (edge_to_node[(a, b)], edge_to_node[(b, c)])
                    if undirected:
                        new_edges.append((min(x, y), max(x, y)))
                    else:
                        new_edges.append((x, y))

        if b in n_to_n_reverse:
            for c in n_to_n_reverse[b]:
                if (c, b) in edges and c != a:
                    (x, y) = (edge_to_node[(a, b)], edge_to_node[(c, b)])
                    if undirected:
                        new_edges.append((min(x, y), max(x, y)))
                    else:
                        new_edges.append((x, y))

    new_edges = set(new_edges)

    f = open("line_graph.csv", "w")
    for (a, b) in new_edges:
        f.write("%d,%d,1\n" % (a, b))
    f.close()
