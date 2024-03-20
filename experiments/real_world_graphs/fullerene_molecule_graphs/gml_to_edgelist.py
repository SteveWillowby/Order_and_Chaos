import sys

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Error! Need gml filename as input.")
        exit(1)

    filename = sys.argv[1]
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip().split(" ") for l in lines]
    nodes = set()
    edges = set()
    s = None
    for line in lines:
        if line[0] == "id":
            nodes.add(int(line[1]))
        elif line[0] == "source":
            s = int(line[1])
        elif line[0] == "target":
            t = int(line[1])
            edges.add((s, t))

    filename = filename.replace(" ", "_")

    edge_nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])
    if len(edge_nodes) < len(nodes):
        print("Extra Nodes -- Making Nodelist")
        out_filename = filename[:-4] + "_nodelist.txt"
        f = open(out_filename, "w")
        for n in nodes:
            f.write("%d\n" % n)
        f.close()

    out_filename = filename[:-4] + ".txt"
    f = open(out_filename, "w")
    for (a, b) in sorted(list(edges)):
        f.write("%d %d\n" % (a, b))
    f.close()
