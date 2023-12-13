import sys

if __name__ == "__main__":

    filename = sys.argv[1]   # edgelist
    dir_str = sys.argv[2]    # should be true or false

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip().split(" ") for l in lines]
    edges = [(int(x[0]), int(x[1])) for x in lines]
    edges.sort()
    nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])

    f = open("graph_file.json", "w")
    f.write("[\n")
    for n in nodes:
        f.write("\t{ \"vertex\": {\n")
        f.write("\t\t\"id\": \"%d\"\n" % n)
        f.write("\t  }},\n")
    for i in range(0, len(edges)):
        (a, b) = edges[i]

        f.write("\t{ \"edge\": {\n")
        f.write("\t\t\"id\": \"%d\",\n" % i)
        f.write("\t\t\"source\": \"%d\",\n" % a)
        f.write("\t\t\"target\": \"%d\",\n" % b)
        f.write("\t\t\"directed\": \"%s\"\n" % dir_str)
        f.write("\t  }}")

        if i < len(edges) - 1:
            f.write(",\n")
        else:
            f.write("\n")
    f.write("]")
