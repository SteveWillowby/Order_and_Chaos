import matplotlib.pyplot as plt
from utils import *
import sys

__graphs_base__ = "../experiments/real_world_graphs/"

__graphs_list__ = ["karate.txt", "season_4_undirected_edges.txt", \
                   "maayan-foodweb.txt",       "pol_blogs.txt", \
                   "eucore.txt",       "cora.txt"]

__graph_idxs__ = {"karate": 0, "season_4": 1, "foodweb": 2, \
                  "pol_blogs": 3, "eucore": 4, "cora": 5}

__nodes_list__ =  [None,         None, \
                   "maayan-foodweb_nodes.txt", "pol_blogs_nodes.txt", \
                   "eucore_nodes.txt", "cora_nodes.txt"]

__dir_list__ =    [False,        False, \
                   True,                       True, \
                   True,               True]


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Needs a graph argument:\n\tkarate\n\tseason_4\n\tfoodweb\n" + \
              "\tpol_blogs\n\teucore\n\tcora")

    graph_name = sys.argv[1].lower()
    if graph_name in __graph_idxs__:
        idx = __graph_idxs__[graph_name]
    else:
        print("Needs a graph argument:\n\tkarate\n\tseason_4\n\tfoodweb\n" + \
              "\tpol_blogs\n\teucore\n\tcora")
        exit(1)

    graph_file = __graphs_base__ + __graphs_list__[idx]
    directed = __dir_list__[idx]

    # Hard reset directed to False
    directed = False

    edges = get_edgeset(graph_file, directed)
    if __nodes_list__[idx] is not None:
        nodes_file = __graphs_base__ + __nodes_list__[idx]
        nodes = get_nodeset(nodes_file)
    else:
        nodes = edges_to_nodes(edges)

    noise_edges = set()
    struct_edges = set(edges)

    ga_itrs_per_cycle = 60
    num_cycles = 6

    scores = []
    for _ in range(0, num_cycles):
        vog_edges = struct_edges

        if len(vog_edges) == 0:
            print("Quitting Early -- No Edges Left")
            break

        (struct_edges, noise_edges) = \
            run_VoG(vog_edges, directed=directed)
        scores.append(run_scorer(edges, nodes, noise_edges, directed))

        seed_noise = (edges - struct_edges) | (struct_edges - edges)

        (struct_edges, noise_edges) = \
            run_GA(edges, directed=directed, nodes=nodes, \
                   n_itr=ga_itrs_per_cycle, seed_noise=seed_noise)

        scores.append(run_scorer(edges, nodes, noise_edges, directed))

    print("SCORES! %s" % str(scores))
    plt.plot([i for i in range(0, len(scores))], scores)
    plt.title("Evolution of Decomp. Score")
    plt.xlabel("VoG + GA Decomps")
    plt.ylabel("Score")
    plt.savefig("%s_fig.png" % graph_name)
    plt.close()
