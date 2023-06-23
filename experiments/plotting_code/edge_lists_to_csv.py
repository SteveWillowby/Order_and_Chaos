import sys

# Pass in noise file, hypothesis graph, output filename
#
# Combines edge lists and appends file idx as a label
#
# Specifically, labels the edge with the FIRST file it sees it in

if __name__ == "__main__":

    if (len(sys.argv) != 4):
        print("Error! Requires exactly 3 arguments!")
        exit(1)

    A_file = sys.argv[1]
    B_file = sys.argv[2]
    output_file = sys.argv[-1]

    A_edges = set()
    B_edges = set()

    with open(A_file, "r") as f:
        A_edges = set([l.strip().replace(" ", ",") for l in f.readlines()])
    with open(B_file, "r") as f:
        B_edges = set([l.strip().replace(" ", ",") for l in f.readlines()])

    f = open(output_file, "w")
    f.write("Source,Target,Weight\n")
    both_edges = A_edges & B_edges
    A_edges = A_edges - both_edges
    B_edges = B_edges - both_edges
    # f.write("100000000,100000001,0.000001\n")  # One dummy edge
    for l in A_edges:  # A only -- deleted edges
        f.write(l + (",%f\n" % 0.01))
    for l in both_edges:  # both -- added edges
        f.write(l + (",%f\n" % 1))
    for l in B_edges:  # B only -- original edges
        f.write(l + (",%f\n" % 2))
    f.close()

    # types = ["A", "B", "C", "D", "E", "F", "G", "H"]
    """
    input_files = sys.argv[1:-1]

    f2 = open(output_file, "w")
    f2.write("Source,Target,Weight\n")
    for i in range(0, len(input_files)):
        input_file = input_files[i]
        f1 = open(input_file, "r")
        edge_list = f1.readlines()
        f1.close()

        edge_list = [l.strip().replace(" ", ",") for l in edge_list]
        new_edge_list = set(edge_list) - seen_edges
        old_edge_list = set(edge_list) & seen_edges
        seen_edges |= new_edge_list
        new_edge_list = sorted(list(new_edge_list))
        old_edge_list = sorted(list(old_edge_list))
        old_edge_list = [l + (",%d\n" % (i*2 + 5)) for l in old_edge_list]
        new_edge_list = [l + (",%d\n" % (i*2 + 6)) for l in new_edge_list]

        for line in old_edge_list:
            f2.write(line)
        for line in new_edge_list:
            f2.write(line)
    f2.close()
    """
