import matplotlib.pyplot as plt
from utils import *
import sys

__graphs_base__ = "../experiments/real_world_graphs/"

__graphs_list__ = ["karate.txt", "season_4_undirected_edges.txt", \
                   "maayan-foodweb.txt",       "pol_blogs.txt", \
                   "eucore.txt",       "cora.txt", \
                   "test_graph.txt", "collins_yeast_interactome.txt"]

__graph_idxs__ = {"karate": 0, "season_4": 1, "foodweb": 2, \
                  "pol_blogs": 3, "eucore": 4, "cora": 5, \
                  "test": 6, "yeast": 7}

__nodes_list__ =  [None,         None, \
                   "maayan-foodweb_nodes.txt", "pol_blogs_nodes.txt", \
                   "eucore_nodes.txt", "cora_nodes.txt",
                   None,             None]

__dir_list__ =    [False,        False, \
                   True,                       True, \
                   True,               True, \
                   False,            False]

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Needs a graph argument:\n\tkarate\n\tseason_4\n\tfoodweb\n" + \
              "\tpol_blogs\n\teucore\n\tcora\n\ttest\n\tyeast")

        print("Also needs an alg name: vog, subdue, or kcore")

    graph_name = sys.argv[1].lower()
    if graph_name in __graph_idxs__:
        idx = __graph_idxs__[graph_name]
    else:
        print("Needs a graph argument:\n\tkarate\n\tseason_4\n\tfoodweb\n" + \
              "\tpol_blogs\n\teucore\n\tcora\n\ttest\n\tyeast")
        exit(1)

    alg_name = sys.argv[2].lower()
    if alg_name == "vog":
        third_party_runner = run_VoG
    elif alg_name == "subdue":
        third_party_runner = run_C_SUBDUE
    elif alg_name == "kcore":
        third_party_runner = (lambda edges, directed=False : \
                                run_k_core(edges, directed=directed, k=3))
    else:
        print("Needs an alg name: vog, subdue, or kcore")

    if graph_name == "test":
        graph_file = __graphs_list__[idx]
    else:
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

    all_noise_score  = run_scorer(edges, nodes, edges, directed)
    all_struct_score = run_scorer(edges, nodes, set(), directed)

    noise_edges = set()
    struct_edges = set(edges)

    ga_itrs_per_cycle = 60
    num_cycles = 6

    scores = []
    struct_size = []
    num_added = []
    num_removed = []

    for _ in range(0, num_cycles):
        third_party_edges = struct_edges

        if len(third_party_edges) == 0:
            print("Quitting Early -- No Edges Left")
            break

        (struct_edges, noise_edges) = \
            third_party_runner(third_party_edges, directed=directed)
        scores.append(run_scorer(edges, nodes, noise_edges, directed))
        struct_size.append(len(struct_edges))
        print("There are %d struct edges after 3rd party." % len(struct_edges))
        print(struct_edges)
        print(noise_edges)

        num_added.append(len(noise_edges - edges))
        num_removed.append(len(edges & noise_edges))

        seed_noise = (edges - struct_edges) | (struct_edges - edges)

        (struct_edges, noise_edges) = \
            run_GA(edges, directed=directed, nodes=nodes, \
                   n_itr=ga_itrs_per_cycle, seed_noise=seed_noise)

        scores.append(run_scorer(edges, nodes, noise_edges, directed))
        struct_size.append(len(struct_edges))
        print("There are %d struct edges after GA." % len(struct_edges))
        num_added.append(len(noise_edges - edges))
        num_removed.append(len(edges & noise_edges))

    num_runs = 10
    rand_scores = []
    for i in range(0, len(scores)):
        na = num_added[i]
        nr = num_removed[i]
        avg_score = 0
        for _ in range(0, num_runs):
            rand_edges = rand_noise_set(edges, nodes, na, nr)
            avg_score += run_scorer(edges, nodes, rand_edges, directed)
        avg_score /= num_runs
        rand_scores.append(avg_score)

    print("SCORES! %s" % str(scores))
    print("SIZES!  %s" % str(struct_size))
    print("Rand Scores! %s" % str(rand_scores))
    x_axis = [i + 1 for i in range(0, len(scores))]
    score_bot_ref = min(all_noise_score, all_struct_score)
    score_top_ref = max(all_noise_score, all_struct_score)
    plot_scores = [(s - score_bot_ref) / (score_top_ref - score_bot_ref) \
                        for s in scores]
    plot_rand_scores = [(s - score_bot_ref) / (score_top_ref - score_bot_ref) \
                        for s in rand_scores]
    plot_edges =  [ne / len(edges) for ne in struct_size]

    plt.plot(x_axis, plot_scores, \
             x_axis, plot_edges, \
             x_axis, plot_rand_scores)
    plt.legend(["Rel. Score", "S. E / O. E", "Rand. SS"])
    plt.title("Evolution of Decomp. Score")
    plt.xlabel("VoG + GA Decomps")
    plt.ylabel("Score")
    plt.savefig("%s_fig.png" % (graph_name + "_" + alg_name))
    plt.close()
