import matplotlib.pyplot as plt
from utils import *


__graphs_base__ = "../experiments/real_world_graphs/"

__graphs_list__ = ["karate.txt", "season_4_undirected_edges.txt", \
                   "maayan-foodweb.txt",       "pol_blogs.txt", \
                   "eucore.txt",       "cora.txt"]

__nodes_list__ =  [None,         None, \
                   "maayan-foodweb_nodes.txt", "pol_blogs_nodes.txt", \
                   "eucore_nodes.txt", "cora_nodes.txt"]

__dir_list__ =    [False,        False, \
                   True,                       True, \
                   True,               True]

__name_list__ =   ["karate", "season_4", "foodweb", "pol_blogs", "eucore", \
                   "cora"]


if __name__ == "__main__":
    for i in range(0, len(__graphs_list__)):
        graph_file = __graphs_base__ + __graphs_list__[i]
        directed = __dir_list__[i]
        edges = get_edgeset(graph_file, directed)
        if __nodes_list__[i] is None:
            nodes = edges_to_nodes(edges)
        else:
            nodes_file = __graphs_base__ + __nodes_list__[i]
            nodes = get_nodeset(nodes_file)

        core_data = []
        truss_data = []
        rand_data = []

        k = 2
        while len(core_data) == 0 or core_data[-1][1] < len(edges):
            (_, noise_edges) = run_k_core(edges, directed=directed, k=k)
            score = run_scorer(edges, nodes, noise_edges, directed)
            core_data.append((score, len(noise_edges)))
            k = k + 1

        k = 3
        while len(truss_data) == 0 or truss_data[-1][1] < len(edges):
            (_, noise_edges) = run_k_truss(edges, directed=directed, k=k)
            score = run_scorer(edges, nodes, noise_edges, directed)
            truss_data.append((score, len(noise_edges)))
            k = k + 1

        noise_sizes = [size for (score, size) in core_data] + \
                      [size for (score, size) in truss_data]

        noise_sizes = sorted(list(set(noise_sizes)))

        for size in noise_sizes:
            num_rand_scores = 10
            avg_rand_score = 0
            for _ in range(0, num_rand_scores):
                rand_noise = rand_noise_set(edges, nodes, 0, size)
                avg_rand_score += run_scorer(edges, nodes, rand_noise, directed)
            avg_rand_score /= num_rand_scores

            rand_data.append((avg_rand_score, size))


        core_x  = [float(size) / len(edges) for (score, size) in core_data]
        truss_x = [float(size) / len(edges) for (score, size) in truss_data]
        rand_x  = [float(size) / len(edges) for (score, size) in rand_data]
        core_y =  [score for (score, size) in core_data]
        truss_y = [score for (score, size) in truss_data]
        rand_y =  [score for (score, size) in rand_data]
        plt.plot(core_x, core_y, \
                 truss_x, truss_y, \
                 rand_x, rand_y, marker="o")
        plt.legend(["Core Scores", "Truss Scores", "Rand Scores"])
        plt.title("Decomp. Score as Edges are Removed")
        plt.xlabel("Fraction of Edges Removed")
        plt.ylabel("Score")
        graph_name = __name_list__[i]
        plt.savefig("results/core_n_truss_%s_fig.png" % (graph_name))
        plt.close()
