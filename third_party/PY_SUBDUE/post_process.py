import os
import sys

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Need four arguments:\n" + \
              "\tsubdue-file graph-file temp-file\n" + \
              "\tbest-noise and directed/undirected")
        exit()

    sub_file = sys.argv[1]
    graph_file = sys.argv[2]
    tmp_file = sys.argv[3]
    best_noise_file = sys.argv[4]
    directed = sys.argv[5]
    if directed.lower() in ["true", "1", "directed"]:
        directed = True
        dir_str = "-d"
    else:
        directed = False
        dir_str = "-u"

    f = open(graph_file, "r")
    lines = f.readlines()
    f.close()

    graph_nodes = set()
    graph_edges = set()
    lines = [l.strip().split(" ") for l in lines]
    for l in lines:
        a = int(l[0])
        b = int(l[1])
        graph_nodes.add(a)
        graph_nodes.add(b)
        if directed:
            graph_edges.add((a, b))
        else:
            graph_edges.add((min(a, b), max(a, b)))

    f = open("nodes.txt", "w")
    graph_nodes_list = list(graph_nodes)
    for i in range(0, len(graph_nodes_list)):
        f.write("%d" % graph_nodes_list[i])
        if i < len(graph_nodes_list) - 1:
            f.write("\n")
    f.close()

    f = open(sub_file, "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip() for l in lines]

    # Get edges in the order the file lists them
    edge_list = []

    for l in lines:
        if len(l) > 4 and l[:4] == "edge":
            l = l.split("--")
            a = l[0].split("(")[1]
            b = l[1].split(")")[0]
            a = int(a)
            b = int(b)

            if not directed:
                (a, b) = (min(a, b), max(a, b))

            edge_list.append((a, b))

    best_score = None
    best_edge_set = set()

    edge_set = set(edge_list)
    # for (a, b) in edge_list:
    #     old_len = len(edge_set)
    #     edge_set.add((a, b))
    #     if len(edge_set) == old_len:
    #         # No new edges
    #         continue

    noise_set = (graph_edges - edge_set) | (edge_set - graph_edges)
    f = open(tmp_file, "w")
    noise_set_list = list(noise_set)
    for i in range(0, len(noise_set_list)):
        (c, d) = noise_set_list[i]
        f.write("%d %d" % (c, d))
        if i < len(noise_set_list) - 1:
            f.write("\n")
    f.close()

    os.system(("../executables/score_only -graph %s " % graph_file) + \
              ("-edges %s -nodes nodes.txt %s > testing/score.txt" % \
                    (tmp_file, dir_str)))

    f = open("testing/score.txt", "r")
    lines = f.readlines()
    f.close()

    l = lines[-1].strip()
    score = float(l)

    if best_score is None or best_score < score:
        best_score = score
        best_noise_set = set(noise_set)

    print("Best Score: %f" % score)
    print("Num Edges in Best Noise Set: %d" % len(best_noise_set))
    # print("Best Noise Set: %s" % str(best_noise_set))

    f = open(best_noise_file, "w")
    best_noise_set = list(best_noise_set)
    for i in range(0, len(best_noise_set)):
        (a, b) = best_noise_set[i]
        f.write("%d %d" % (a, b))
        if i < len(best_noise_set) - 1:
            f.write("\n")
    f.close()
