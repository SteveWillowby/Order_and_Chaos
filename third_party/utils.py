import random
import os


# Returns the score
def run_scorer(graph_edges, nodes, noise_edges, directed):
    nodes_file =  "/tmp/score_nodes.txt"
    graph_file =  "/tmp/score_graph.txt"
    noise_file =  "/tmp/score_noise.txt"
    result_file = "/tmp/score_result.txt"

    write_nodeset(nodes_file, nodes)
    write_edgeset(graph_file, graph_edges)
    write_edgeset(noise_file, noise_edges)

    dir_str = ["-u", "-d"][int(directed)]

    os.system(("../executables/score_only -graph %s " % graph_file) + \
              ("-edges %s -nodes %s %s > %s" % \
                (noise_file, nodes_file, dir_str, result_file)))

    f = open(result_file, "r")
    lines = f.readlines()
    f.close()

    l = lines[-1].strip()
    score = float(l)
    return score

# Returns:
#   (struct_edges, noise_edges)
def run_GA(edges, directed=False, nodes=None, n_itr=120, seed_noise=None):

    tmp_graph = "/tmp/GA_tmp_graph.txt"
    tmp_nodes = "/tmp/GA_tmp_nodes.txt"
    tmp_out   = "/tmp/GA_tmp_out"
    tmp_seed  = "/tmp/GA_tmp_seed.txt"
    dir_str = ["-u", "-d"][int(directed)]

    if nodes is None:
        nodes = edges_to_nodes(edges)

    # Zero-index everything
    nodes = sorted(list(nodes))
    new_nodes = set([i for i in range(0, len(nodes))])
    old_to_new = {nodes[i]: i for i in range(0, len(nodes))}
    new_to_old = {i: nodes[i] for i in range(0, len(nodes))}
    edges = set([(old_to_new[a], old_to_new[b]) for (a, b) in edges])

    write_nodeset(tmp_nodes, new_nodes)
    write_edgeset(tmp_graph, edges)

    print("GA --- THERE ARE %d NODES" % len(nodes))
    print("GA --- THERE ARE %d EDGES" % len(edges))

    if seed_noise is None:
        os.system("../executables/main -graph %s -nodes %s -o %s -n_itr %d %s" % \
                    (tmp_graph, tmp_nodes, tmp_out, n_itr, dir_str))
    else:
        seed_noise = set([(old_to_new[a], old_to_new[b]) \
                                for (a, b) in seed_noise])
        write_edgeset(tmp_seed, seed_noise)
        os.system("../executables/main -graph %s -nodes %s -seed %s -o %s -n_itr %d %s" % \
                    (tmp_graph, tmp_nodes, tmp_seed, tmp_out, n_itr, dir_str))

    noise = get_edgeset(tmp_out + "_noise.txt", directed)
    struct_edges = (edges - noise) | (noise - edges)

    # Restore the node labels.
    return (set([(new_to_old[a], new_to_old[b]) for (a, b) in struct_edges]), \
            set([(new_to_old[a], new_to_old[b]) for (a, b) in noise]))

# Returns:
#   (struct_edges, noise_edges)
def run_VoG(edges, directed=False):

    # Make sure the edges are undirected
    #   (and no self-loops)
    assert len(edges) * 2 == len(edges | set([(b, a) for (a, b) in edges]))

    file_base = "input_graph"
    tmp_graph = "VoG/DATA/%s.txt" % file_base

    nodes = edges_to_nodes(edges)
    print("VOG -- THERE ARE %d NODES" % len(nodes))
    print("VOG -- THERE ARE %d EDGES" % len(edges))

    # 1-index everything
    nodes = sorted(list(nodes))
    new_nodes = set([i for i in range(1, len(nodes) + 1)])
    old_to_new = {nodes[i]: (i + 1) for i in range(0, len(nodes))}
    new_to_old = {(i + 1): nodes[i] for i in range(0, len(nodes))}
    edges = set([(old_to_new[a], old_to_new[b]) for (a, b) in edges])

    write_VoG_edgeset(tmp_graph, edges)

    os.system("cd VoG; matlab -r run_structureDiscovery")
    os.system(("python2.7 VoG/MDL/greedySearch_nStop.py %s " % tmp_graph) + \
              ("VoG/DATA/%s_orderedALL.model >/dev/null 2>&1" % file_base))
    os.system("mv heuristic* VoG/DATA/")

    (s_e, n_e) = read_VoG_decomposition(edges, new_nodes, file_base)

    # Restore the node labels.
    return (set([(new_to_old[a], new_to_old[b]) for (a, b) in s_e]), \
            set([(new_to_old[a], new_to_old[b]) for (a, b) in n_e]))


# Returns:
#   (struct_edges, noise_edges)
def read_VoG_decomposition(graph_edges, nodes, file_base):
    model_file = ("VoG/DATA/%s_orderedALL.model" % file_base)
    lines_file = ("VoG/DATA/heuristicSelection_nStop_ALL_%s_orderedALL.model" \
                    % file_base)

    f = open(model_file, "r")
    structures = f.readlines()
    structures = [l.strip() for l in structures]
    f.close()

    f = open(lines_file, "r")
    relevant_lines = f.readlines()
    f.close()
    relevant_lines = [s.strip() for s in relevant_lines]
    relevant_line_ints = []
    for s in relevant_lines:
        if len(s) > 0:
            relevant_line_ints.append(int(s) - 1)

    structures = [structures[x] for x in relevant_line_ints]

    struct_edges = set()
    noise_edges  = set()

    for structure in structures:
        structure = structure.strip()

        if len(structure) == 0:
            continue

        s = structure.split(" ")
        s_type = s[0]
        s = structure[2:]
        if s_type == "st":
            # Star
            v = s.split(", ")
            hub = v[0]
            spokes = v[1].strip().split(" ")
            se = [(int(hub), int(sp)) for sp in spokes]
            se = [(min(a, b), max(a, b)) for (a, b) in se]
            struct_edges |= set(se)
        elif s_type == "nb" or s_type == "bc":
            # (Near-)Bipartite Core
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
            v = s.strip().split(" ")
            se = [(v[i], v[i + 1]) for i in range(0, len(v) - 1)]
            se = [(int(a), int(b)) for (a, b) in se]
            se = [(min(a, b), max(a, b)) for (a, b) in se]
            struct_edges |= set(se)
        elif s_type == "fc" or s_type == "nc":
            # Full/Near-Clique
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

    return (struct_edges, noise_edges)


def get_edgeset(edge_file, directed, remove_self_loops=True):
    f = open(edge_file, "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip() for l in lines]
    new_lines = []
    for l in lines:
        l = l.split(" ")
        if len(l) == 2:
            new_lines.append(l)
    lines = [(int(x[0]), int(x[1])) for x in new_lines]

    edges = set()
    for (a, b) in lines:
        if remove_self_loops and a == b:
            continue
        if directed:
            edges.add((a, b))
        else:
            edges.add((min(a, b), max(a, b)))
    return edges


def edges_to_nodes(edgeset):
    return set([a for (a, b) in edgeset] + [b for (a, b) in edgeset])


def get_nodeset(node_file):
    f = open(node_file, "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip() for l in lines]
    nodes = set()
    for l in lines:
        if len(l) > 0:
            nodes.add(int(l))
    return nodes


def write_edgeset(edge_file, edges):
    f = open(edge_file, "w")
    edges = list(edges)
    for i in range(0, len(edges)):
        f.write("%d %d" % edges[i])
        if i < len(edges) - 1:
            f.write("\n")
    f.close()


def write_nodeset(node_file, nodes):
    f = open(node_file, "w")
    nodes = list(nodes)
    for i in range(0, len(nodes)):
        f.write(str(nodes[i]))
        if i < len(nodes) - 1:
            f.write("\n")
    f.close()


def write_VoG_edgeset(edge_file, edges):
    lines = [str(a) + "," + str(b) + ",1\n" for (a, b) in edges]
    f = open(edge_file, "w")
    for l in lines:
        f.write(l)
    f.close()


# Returns a random set of noise edges
def rand_noise_set(graph_edges, nodes, num_added, num_removed, directed=False):
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

        if directed:
            e = (a, b)
        else:
            e = (min(a, b), max(a, b))

        if e in graph_edges:
            continue
        added.add(e)

    return added | removed
