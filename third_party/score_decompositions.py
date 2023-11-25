from utils import *
import sys

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

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Error! Need to pass an algorithm as input: VoG or GA")
        exit(1)

    algorithm = sys.argv[1]

    always_undirected = False
    if algorithm.lower() == "vog":
        always_undirected = True
        decomp_fn = run_VoG
    elif algorithm.lower() == "ga":
        decomp_fn = run_GA
    else:
        print("Error! Need to pass an algorithm as input: VoG or GA")
        exit(1)


    for i in range(0, len(__graphs_list__)):
        graph_file = __graphs_base__ + __graphs_list__[i]
        directed = __dir_list__[i] and not always_undirected
        edges = get_edgeset(graph_file, directed)
        if __nodes_list__[i] is None:
            nodes = edges_to_nodes(edges)
        else:
            nodes_file = __graphs_base__ + __nodes_list__[i]
            nodes = get_nodeset(nodes_file)

        (struct_edges, noise_edges) = decomp_fn(edges, directed=directed)

        assert len(struct_edges) + len(noise_edges) >= len(edges)

        print("###")
        print("#Edges:  %d" % len(edges))
        print("#Struct: %d" % len(struct_edges))
        print("#Noise:  %d" % len(noise_edges))

        score = run_scorer(edges, nodes, noise_edges, directed)

        # Get random scores
        num_added =   len(noise_edges - edges)
        num_removed = len(edges & noise_edges)

        num_rand_scores = 10
        avg_rand_score = 0
        for _ in range(0, num_rand_scores):
            rand_noise = rand_noise_set(edges, nodes, num_added, num_removed)
            avg_rand_score += run_scorer(edges, nodes, rand_noise, directed)
        avg_rand_score /= num_rand_scores

        if i == 1:
            break
