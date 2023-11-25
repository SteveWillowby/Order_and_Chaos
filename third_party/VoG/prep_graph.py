import sys


if __name__ == "__main__":

    filename = sys.argv[1]

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip().split(" ") for l in lines]
    lines = [(int(x[0]), int(x[1])) for x in lines]
    lines = [(a + 1, b + 1) for (a, b) in lines]
    # The graph is undirected, so remove duplicate edges
    lines = [(min(a, b), max(a, b)) for (a, b) in lines]
    lines = list(set(lines))

    lines = [str(a) + "," + str(b) + ",1\n" for (a, b) in lines]
    f = open("DATA/input_graph.txt", "w")
    for l in lines:
        f.write(l)
    f.close()
