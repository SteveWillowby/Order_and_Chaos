# Usage:
#
# python3 score_decompositions.py <alg_name> <core>
#
#   `alg_name` should be one of the following:
#       vog, kcore, ktruss, subdue, ga
#
#   `core` is an optional argument. If you set it to core_only then
#       the program will run the alg once, get all the non-singleton nodes,
#       then re-run on the portion of the graph with only those nodes.

from utils import *

import matplotlib.pyplot as plt
import random
import sys

__graphs_base__ = "../experiments/real_world_graphs/"

__graphs_list__ = ["karate.txt", "season_4_undirected_edges.txt", \
                   "maayan-foodweb.txt",       "pol_blogs.txt", \
                   "eucore.txt",       "cora.txt", \
                   "flickr.txt", "epinions.txt", "enron_email.txt", \
                   "Fullerene_C180.txt", "Fullerene_C720.txt", \
                   "Fullerene_C6000.txt"]

__nodes_list__ =  [None,         None, \
                   "maayan-foodweb_nodes.txt", "pol_blogs_nodes.txt", \
                   "eucore_nodes.txt", "cora_nodes.txt", \
                   None,         None,           None, \
                   None,                 None, \
                   None]

__dir_list__ =    [False,        False, \
                   True,                       True, \
                   True,               True, \
                   False,        False,           True, \
                   False,                False, \
                   False]

__name_list__ =   ["karate", "season_4", "foodweb", "pol_blogs", "eucore", \
                   "cora", "flickr", "epinions", "enron", \
                   "fullerene_c180", "fullerene_c720", "fullerene_c6000"]

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Error! Need to pass an algorithm as input:\n" + \
              "\tVoG, SUBDUE, GA, kcore, or ktruss")
        exit(1)

    algorithm = sys.argv[1]

    core_only = False
    rand_graphs = False
    for idx in range(2, len(sys.argv)):
        core_only   |= sys.argv[idx].lower() in ["core_only"]
        rand_graphs |= sys.argv[idx].lower() in ["rand", "rand_graphs"]

    rand_rounds = 100
    rand_min_p = 0.01
    rand_max_p = 0.5
    rand_n = 200
    rand_results = []

    preprocess = False # If True, runs the GA before running a third-party alg.

    always_undirected = False
    if algorithm.lower() == "vog":
        always_undirected = True
        decomp_fn = run_VoG
    elif algorithm.lower() == "ga":
        decomp_fn = run_GA
    elif algorithm.lower() == "subdue":
        #   Run SUBDUE but swap noise and structure
        # decomp_fn = (lambda edges, directed=False : \
        #              [(y, x) for (x, y) in \
        #                 [run_C_SUBDUE(edges, min_size=3, max_size=10, \
        #                               iterations=0, directed=directed)]][0])

        # Run SUBDUE like normal
        decomp_fn = (lambda edges, directed=False : \
                        run_C_SUBDUE(edges, min_size=3, max_size=10, \
                                      iterations=0, directed=directed))
    # elif algorithm.lower() == "grami":
    #     decomp_fn = (lambda edges, directed=False : \
    #                     run_GraMi(edges, directed=directed, \
    #                               min_support="auto"))
    elif algorithm.lower() == "kcore":
        k = 3
        print("Using k = %d" % k)
        decomp_fn = (lambda edges, directed=False : \
                        run_k_core(edges, directed=directed, k=k))
    elif algorithm.lower() == "ktruss":
        k = 3
        print("Using k = %d" % k)
        decomp_fn = (lambda edges, directed=False : \
                        run_k_truss(edges, directed=directed, k=k))
    else:
        print("Error! Need to pass an algorithm as input: VoG, SUBDUE, GA, or kcore")
        exit(1)

    num_rounds = len(__graphs_list__)
    if rand_graphs:
        num_rounds = rand_rounds
    for i in range(0, num_rounds):

        if rand_graphs:
            directed = False
            rand_p = random.random() * (rand_max_p - rand_min_p) + rand_min_p
            (nodes, edges) = get_ER_rand_graph(rand_n, rand_p, directed)

        else:
            if __name_list__[i] in ["flickr", "epinions"]:
                continue
            if __name_list__[i] in ["enron"] and \
                    algorithm.lower() not in ["vog"]:
                continue
            if __name_list__[i] in ["fullerene_c6000"] and \
                    algorithm.lower() not in ["subdue"]:
                continue

            graph_file = __graphs_base__ + __graphs_list__[i]
            directed = __dir_list__[i] and not always_undirected

            edges = get_edgeset(graph_file, directed)
            if __nodes_list__[i] is None:
                nodes = edges_to_nodes(edges)
            else:
                nodes_file = __graphs_base__ + __nodes_list__[i]
                nodes = get_nodeset(nodes_file)

        if preprocess and algorithm.lower() != "ga":
            (edges, _) = run_GA(edges, directed=directed)

        (struct_edges, noise_edges) = decomp_fn(edges, directed=directed)

        if core_only:
            nodes = edges_to_nodes(struct_edges)
            if len(nodes) == 0:
                continue
            edges = limit_edges_by_nodeset(edges, nodes)
            noise_edges = limit_edges_by_nodeset(noise_edges, nodes)
            # struct_edges is already limited to that set of nodes

            # (struct_edges, noise_edges) = decomp_fn(edges, directed=directed)


        assert len(struct_edges) + len(noise_edges) >= len(edges)

        print("###")
        if rand_graphs:
            print("#  Random Graph with n = %d, p = %f" % (rand_n, rand_p))
        else:
            print("#  %s" % __name_list__[i])
        print("#Edges:  %d" % len(edges))
        print("#Struct: %d" % len(struct_edges))
        print("#Noise:  %d" % len(noise_edges))

        score = run_scorer(edges, nodes, noise_edges, directed)
        all_noise = run_scorer(edges, nodes, set(edges), directed)
        no_noise = run_scorer(edges, nodes, set(), directed)

        # Get random scores
        num_added =   len(noise_edges - edges)
        num_removed = len(edges & noise_edges)

        print("#Added:  %d" % num_added)
        print("#Del:    %d" % num_removed)
        print("\t#All Noise:      %f" % all_noise[0])
        print("\t#All Structure:  %f" % no_noise[0])
        print("\t\t#From Aut (no singletons):    %f" % no_noise[1])
        print("\t\t#From Aut (singletons only):  %f" % no_noise[2])
        print("\t\t#From AO (no sing. swaps):    %f" % no_noise[3])
        print("\t\t#From AO (sing. swaps only):  %f" % no_noise[4])
        print("\t\t#From noise size probability: %f" % no_noise[5])
        print("\t#Score:          %f" % score[0])
        print("\t#W/O Singletons: %f" % (score[1] + score[3] + score[5]))
        print("\t\t#From Aut (no singletons):    %f" % score[1])
        print("\t\t#From Aut (singletons only):  %f" % score[2])
        print("\t\t#From AO (no sing. swaps):    %f" % score[3])
        print("\t\t#From AO (sing. swaps only):  %f" % score[4])
        print("\t\t#From noise size probability: %f" % score[5])

        num_rand_scores = 10
        avg_rand_score = [0, 0, 0, 0, 0, 0]
        for _ in range(0, num_rand_scores):
            rand_noise = rand_noise_set(edges, nodes, num_added, num_removed)
            rand_score = run_scorer(edges, nodes, rand_noise, directed)
            for j in range(0, 6):
                avg_rand_score[j] += rand_score[j]
        for j in range(0, 6):
            avg_rand_score[j] /= num_rand_scores

        print("\t#AR Score:       %f" % avg_rand_score[0])
        print("\t#W/O Singletons: %f" % \
                (avg_rand_score[1] + avg_rand_score[3] + avg_rand_score[5]))
        print("\t\t#From Aut (no singletons):    %f" % avg_rand_score[1])
        print("\t\t#From Aut (singletons only):  %f" % avg_rand_score[2])
        print("\t\t#From AO (no sing. swaps):    %f" % avg_rand_score[3])
        print("\t\t#From AO (sing. swaps only):  %f" % avg_rand_score[4])
        print("\t\t#From noise size probability: %f" % avg_rand_score[5])

        bottoms      = [no_noise[5], score[5], avg_rand_score[5]]
        aut_no_sing  = [no_noise[1], score[1], avg_rand_score[1]]
        aut_sing     = [no_noise[2], score[2], avg_rand_score[2]]
        ao_no_sing   = [no_noise[3], score[3], avg_rand_score[3]]
        ao_sing      = [no_noise[4], score[4], avg_rand_score[4]]

        if rand_graphs:
            rand_results.append((rand_p, all_noise[0], no_noise[0], \
                                 score[0], avg_rand_score[0]))

        else:

            bar_names = ["No Noise", "Alg's Choice", "Random"]
            graph_name = __name_list__[i]
            graph_disp_name = graph_name[0].upper() + graph_name[1:]
            dir_str = ["Und", "D"][int(directed)] + "irected"
            stacks = [aut_no_sing, aut_sing, ao_no_sing, ao_sing]
            labels = ["Aut - No Sing.","Aut - Sing.","AO - No Sing.","AO - Sing."]
            plt.figure(figsize=(9,6))
            for j in range(0, len(stacks)):
                s = stacks[j]
                plt.bar(bar_names, s, bottom=bottoms, label=labels[j])
                for i in range(0, len(bottoms)):
                    bottoms[i] += s[i]
            plt.legend()
            plt.suptitle("Score Decompositions - %s" % algorithm, size=20)
            plt.title("%s %s" % (graph_disp_name, dir_str), size=18)
            plt.xlabel("Models", size=18, labelpad=10)
            plt.ylabel("Gains Above Edit Cost", size=18)
            plt.savefig("results/%s%s_decomp_%s_%s%s.png" % \
                            (["", "core_only/"][int(core_only)], \
                             algorithm.lower(), graph_name, dir_str.lower(), \
                             ["", "_core"][int(core_only)]), \
                        bbox_inches='tight', pad_inches=0.5)
            plt.close()

    if rand_graphs:
        dir_str = ["Und", "D"][int(directed)] + "irected"

        # Sort by rand prob used
        rand_results.sort()
        x_axis = [p for (p, _, __, ___, ____) in rand_results]
        y_axis_1 = [alg - all_struct for \
                        (_, __, all_struct, alg, ___) in rand_results]
        y_axis_2 = [rand - all_struct for \
                        (_, __, all_struct, ___, rand) in rand_results]
        plt.figure(figsize=(9,6))
        plt.scatter(x_axis, y_axis_1, label=algorithm)
        plt.scatter(x_axis, y_axis_2, label="Random Changes")
        plt.legend()
        plt.suptitle("Score Decompositions of Noise - %s" % algorithm, size=20)
        plt.title("%d-Node ER Graphs" % rand_n, size=18)
        plt.xlabel("Edge Probability", size=18, labelpad=10)
        plt.ylabel("Gains Above All Structure", size=18)
        plt.savefig("results/%s%s_rand_decomp_%s%s.png" % \
                        (["", "core_only/"][int(core_only)], \
                         algorithm.lower(), dir_str.lower(), \
                         ["", "_core"][int(core_only)]), \
                    bbox_inches='tight', pad_inches=0.5)
