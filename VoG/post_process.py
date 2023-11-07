import os
import random
import sys

def rand_edge_set(graph_edges, nodes, num_added, num_removed):
    ge_list = list(graph_edges)
    removed = set()
    while len(removed) < num_removed:
        removed.add(ge_list[random.randint(0, len(ge_list) - 1)])

    added = set()
    nodes_list = list(nodes)
    while len(added) < num_added:
        a = nodes_list[random.randint(0, len(nodes) - 1)]
        b = nodes_list[random.randint(0, len(nodes) - 1)]
        if a == b:
            continue
        e = (min(a, b), max(a, b))
        if e in graph_edges:
            continue
        added.add(e)

    return (added, removed)

def run_scorer(graph_edges, nodes, noise_edges):
    f = open("nodes.txt", "w")
    for n in nodes:
        f.write("%d\n" % n)
    f.close()

    f = open("graph.txt", "w")
    for (a, b) in graph_edges:
        f.write("%d %d\n" % (a, b))
    f.close()

    f = open("noise_edges.txt", "w")
    for (a, b) in noise_edges:
        f.write("%d %d\n" % (a, b))
    f.close()
    os.system("../executables/score_only -graph graph.txt " + \
              "-edges noise_edges.txt -nodes nodes.txt > result.txt")

    f = open("result.txt", "r")
    lines = f.readlines()
    f.close()

    l = lines[-1].strip()
    score = float(l)
    return score


if __name__ == "__main__":

    model_file = sys.argv[1]
    graph_file = sys.argv[2]
    chosen_model_lines = sys.argv[3]
    output_base = sys.argv[4]
    output_name = output_base + ".txt"
    output_graph = output_base + "_og_graph.txt"
    output_noise = output_base + "_og_noise.txt"
    output_nodes = output_base + "_og_nodes.txt"

    f = open(model_file, "r")
    structures = f.readlines()
    structures = [l.strip() for l in structures]
    f.close()

    f = open(chosen_model_lines, "r")
    relevant_lines = f.readlines()
    f.close()
    relevant_lines = [s.strip() for s in relevant_lines]
    relevant_line_ints = []
    for s in relevant_lines:
        if len(s) > 0:
            relevant_line_ints.append(int(s) - 1)

    # print("Num Structures: %d" % len(structures))
    # print("Lines: %s" % str(relevant_line_ints))
    structures = [structures[x] for x in relevant_line_ints]

    f = open(graph_file, "r")
    graph_edges = f.readlines()
    graph_edges = [l.strip().split(",") for l in graph_edges]
    graph_edges = [(int(x[0]), int(x[1])) for x in graph_edges]
    graph_edges = [(min(a, b), max(a, b)) for (a, b) in graph_edges]
    graph_edges = set(graph_edges)
    f.close()

    nodes = set([a for (a, b) in graph_edges] + \
                [b for (a, b) in graph_edges])

    struct_edges = set()
    noise_edges  = set()

    start_score = None
    best_score = None
    best_noise_edges  = set()
    best_struct_edges = set()

    print("##########################")
    print(output_name)

    for structure in ["bb    ", "cc     "] + structures:
        structure = structure.strip()

        if len(structure) == 0:
            continue

        s = structure.split(" ")
        s_type = s[0]
        s = structure[2:]
        print("")
        if s_type == "bb":
            print("No Noise:")
            struct_edges = graph_edges
        elif s_type == "cc":
            print("All Noise:")
            struct_edges = set()
        elif s_type == "st":
            # Star
            print("Star")
            v = s.split(", ")
            hub = v[0]
            spokes = v[1].strip().split(" ")
            se = [(int(hub), int(sp)) for sp in spokes]
            se = [(min(a, b), max(a, b)) for (a, b) in se]
            struct_edges |= set(se)
        elif s_type == "nb" or s_type == "bc":
            # (Near-)Bipartite Core
            print("Bipartite Core")
            v = s.split(", ")
            g1 = v[0].strip().split(" ")
            g2 = v[1].strip().split(" ")
            se = []
            for a in g1:
                for b in g2:
                    se.append((int(a), int(b)))
            se = [(min(a, b), max(a, b)) for (a, b) in se]
            struct_edges |= set(se) & graph_edges
        elif s_type == "ch":
            # Chain -- NOT Cycle
            print("Chain")
            v = s.strip().split(" ")
            se = [(v[i], v[i + 1]) for i in range(0, len(v) - 1)]
            se = [(int(a), int(b)) for (a, b) in se]
            se = [(min(a, b), max(a, b)) for (a, b) in se]
            struct_edges |= set(se)
        elif s_type == "fc" or s_type == "nc":
            # Full/Near-Clique
            print("Clique")
            if s_type == "nc":
                s = s.split(", ")[1]
            v = s.strip().split(" ")
            v = [int(x) for x in v]
            se = []
            for i in range(0, len(v)):
                for j in range(i + 1, len(v)):
                    se.append((v[i], v[j]))
            se = [(min(a, b), max(a, b)) for (a, b) in se]
            struct_edges |= set(se) & graph_edges

        noise_edges  = graph_edges - struct_edges
        noise_edges |= struct_edges - graph_edges

        print("Num Noise Edges: \t%d" % len(noise_edges))
        print("Num Orig. Edges: \t%d" % len(graph_edges & struct_edges))
        print("Num Struct Edges: \t%d" % len(struct_edges))

        score = run_scorer(graph_edges, nodes, noise_edges)
        print(score)

        if s_type == "bb":
            start_score = score
        if s_type == "cc":
            all_noise_score = score

        if True: # or best_score is None or score > best_score:
            best_score = score
            best_noise_edges = set(noise_edges)
            best_struct_edges = set(struct_edges)

    f = open(output_graph, "w")
    ge_list = list(graph_edges)
    for i in range(0, len(ge_list)):
        f.write("%d %d" % ge_list[i])
        if i < len(ge_list) - 1:
            f.write("\n")
    f.close()

    f = open(output_noise, "w")
    n_list = list(noise_edges)
    for i in range(0, len(n_list)):
        f.write("%d %d" % n_list[i])
        if i < len(n_list) - 1:
            f.write("\n")
    f.close()

    f = open(output_nodes, "w")
    n_list = list(nodes)
    for i in range(0, len(n_list)):
        f.write("%d" % n_list[i])
        if i < len(n_list) - 1:
            f.write("\n")
    f.close()

    num_rand_scores = 10
    rand_score_avg = 0
    for k in range(0, num_rand_scores):
        (added, removed) = \
            rand_edge_set(graph_edges, nodes, \
                          len(noise_edges - graph_edges), \
                          len(noise_edges & graph_edges))
        rand_score = run_scorer(graph_edges, nodes, added | removed)
        rand_score_avg += rand_score
    rand_score_avg /= num_rand_scores


    f = open(output_name, "w")
    f.write("All-is-Structure Score: \t%f\n" % start_score)
    f.write("All-is-Noise Score: \t%f\n" % all_noise_score)
    f.write("Best Score: \t%f\n" % best_score)
    f.write("Avg. Rand Noise Score: \t%f\n" % rand_score_avg)
    f.write("# Edges in Best Noise: \t%f\n" % len(noise_edges))
    f.write("# Edges in Best Structure: \t%f\n" % len(struct_edges))
    f.write("# Edges in Orig. Graph: \t%f\n" % len(graph_edges))
    f.write("Noise  Edges: \t%s\n" % \
str(list(best_noise_edges)).replace("),", ")").replace("[","").replace("]",""))
    f.write("Struct Edges: \t%s\n" % \
str(list(best_struct_edges)).replace("),",")").replace("[","").replace("]",""))
    f.close()
