import sys

if __name__ == "__main__":
    graph  = sys.argv[1]
    output = sys.argv[2]
    directed = sys.argv[3]
    if directed == "u":
        directed = False
    else:
        directed = True

    f = open(output, "r")
    lines = f.readlines()
    f.close()

    s1 = lines[3].strip()[1:-1].replace(" ","").split(",")
    s2 = lines[4].strip()[1:-1].replace(" ","").split(",")

    s1 = set([int(x) for x in s1])
    s2 = set([int(x) for x in s2])

    f = open(graph, "r")
    lines = f.readlines()
    f.close()

    lines = [x.strip().split(" ") for x in lines]
    graph_edges = set()
    nodes = set()
    for line in lines:
        if len(line) > 0:
            a = int(line[0])
            b = int(line[1])
            nodes.add(a)
            nodes.add(b)
            if directed:
                graph_edges.add((a, b))
            else:
                graph_edges.add((min(a, b), max(a, b)))

    un_chosen_1 = nodes - s1
    un_chosen_2 = nodes - s2
    if len(un_chosen_1) < len(s1):
        print("Set 1 bigger than complement")
        # print("Using complement for set 1")
        # s1 = un_chosen_1
    if len(un_chosen_2) < len(s2):
        print("Set 2 bigger than complement")
        # print("Using complement for set 2")
        # s2 = un_chosen_2

    empty = False
    if len(s1) == 0:
        print("Set 1 has no elements")
        empty = True
    if len(s2) == 0:
        print("Set 2 has no elements")
        empty = True

    if empty:
        print("Quitting")
        f = open("chosen_edges.txt", "w")
        f.write("\n")
        # f.write("1 2\n")
        f.close()
        exit(0)

    edges = set()
    for n in s1:
        for n2 in s2:
            if directed:
                edges.add((n, n2))
            else:
                edges.add((min(n, n2), max(n, n2)))

    # Treat the edges as the structure
    noise_edges = graph_edges - edges
    # Consider the following
    noise_edges = (graph_edges - edges) | (edges - graph_edges)

    f = open("chosen_edges.txt", "w")
    for (a, b) in noise_edges:
        f.write("%d %d\n" % (a, b))
    f.close()

    f = open("nodes.txt", "w")
    for n in nodes:
        f.write("%d\n" % n)
    f.close()

    # print(s1)
    # print(s2)
    # print(nodes)
