import sys

# Usage: python3 make_node_list.py <edgelist_filename/path>
#
# Expects edges to be pairs of integers separated by a space or comma
#   (but all spaces or all commas, not a mix in the same file)
if __name__ == "__main__":
    edgelist_name = sys.argv[1]

    if "edges" in edgelist_name:
        nodelist_name = edgelist_name.replace("edges", "nodes")
    else:
        nodelist_name = edgelist_name.split(".")
        if len(nodelist_name) == 1:
            nodelist_name = nodelist_name[0] + "_nodes"
        else:
            nodelist_name[-2] += "_nodes"
            nln = nodelist_name[0]
            for i in range(1, len(nodelist_name)):
                nln += "." + nodelist_name[i]
            nodelist_name = nln

    edge_file = open(edgelist_name, "r")
    edges = edge_file.read()
    edge_file.close()
    
    assert ("," in edges) or (" " in edges)

    spaces = " " in edges

    edges = edges.strip().split("\n")

    if spaces:
        edges = [e.strip().split(" ") for e in edges]
    else:
        edges = [e.strip().split(",") for e in edges]

    edges = [(int(a), int(b)) for (a, b) in edges]
    nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])
    nodes = sorted(list(nodes))

    node_file = open(nodelist_name, "w")
    for n in nodes:
        node_file.write("%d\n" % n)
    node_file.close()
