import random
import os


# Returns a six-tuple containing:
#   The total score
#   The part of the score due to Aut of connected components
#   The part of the score due to Aut of singletons
#   The part of the score due to AO size excepting singleton swaps
#   The part of the score due to AO size from singleton swaps
#   The part of the score due to noise size probability
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

    score_info = (float(lines[-1].strip()),
                  float(lines[-11].strip()),
                  float(lines[-9].strip()),
                  float(lines[-7].strip()),
                  float(lines[-5].strip()),
                  float(lines[-3].strip()))
    return score_info

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


# Returns:
#   (struct_edges, noise_edges)
def run_C_SUBDUE(edges, directed=False, \
                 min_size=3, max_size=6, iterations=0, \
                 temp_in_filename="C_SUBDUE/testing/graph_file.txt", \
                 temp_out_filename="C_SUBDUE/testing/output.txt"):

    nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])

    nodes = sorted(list(nodes))
    node_to_id = {nodes[i]: i + 1 for i in range(0, len(nodes))}
    id_to_node = {new: old for (old, new) in node_to_id.items()}

    nodes = [node_to_id[n] for n in nodes]
    edges = set([(node_to_id[a], node_to_id[b]) for (a, b) in edges])

    dir_str = ["u", "d"][int(directed)]

    f = open(temp_in_filename, "w")
    for n in nodes:
        f.write("v %d node_type_1\n" % n)
    for (a, b) in edges:
        f.write("%s %d %d edge_type_1\n" % (dir_str, a, b))
    f.close()

    os.system("rm " + temp_out_filename)
    os.system(("C_SUBDUE/bin/subdue -compress -eval 1 -minsize %d " % min_size)\
              + ("-maxsize %d -iterations %d " % (max_size, iterations)) \
              + "-nsubs 1 -output 3 " + temp_in_filename + " >> " + \
              temp_out_filename)

    chosen_edges = \
        read_C_SUBDUE_output_maybe(temp_out_filename, edges, nodes, directed)

    noise_edges = (edges - chosen_edges) | (chosen_edges - edges)

    chosen_edges = set([(id_to_node[a], id_to_node[b]) \
                                for (a, b) in chosen_edges])
    noise_edges  = set([(id_to_node[a], id_to_node[b]) \
                                for (a, b) in noise_edges])

    return (chosen_edges, noise_edges)

def read_C_SUBDUE_output_maybe(filename, edges, nodes, directed=False):

    assert min(nodes) == 1
    assert max(nodes) == len(nodes)

    ids_to_nodes = {n: set([n]) for n in nodes}

    f = open(filename)
    lines = f.readlines()
    f.close()

    lines = [l.strip().split(" ") for l in lines]
    struct_info = []
    instance_on = False
    for l in lines:
        if l[0] == "Best":
            struct_info.append([])  # New set of structures
            instance_on = False
        elif l[0] == "Instance":
            # New structure instance: nodes and edges
            struct_info[-1].append([[], []])
            instance_on = True
        elif instance_on and l[0] == "v":
            node = int(l[1])
            struct_info[-1][-1][0].append(node)
        elif instance_on and (l[0] == "u" or l[0] == "d" or l[0] == "e"):
            node_a = int(l[1])
            node_b = int(l[2])
            struct_info[-1][-1][1].append((node_a, node_b))

    # Follow the heirarchical collapse of nodes into single nodes
    #   Every time an edge is removed, add it to chosen_edges
    chosen_edges = set()

    remaining_nodenames = set(nodes)
    for struct_round in struct_info:
        selected_nodes = set()
        for struct in struct_round:
            for n in struct[0]:
                selected_nodes.add(n)
        unselected_nodes = remaining_nodenames - selected_nodes

        old_ids_to_nodes = dict(ids_to_nodes)
        unselected_nodes = sorted(list(unselected_nodes))

        ids_to_nodes = {}

        idx = 1
        for struct in struct_round:
            struct_nodes = set(struct[0])
            struct_edges = struct[1]

            for (a, b) in struct_edges:
                real_nodes_a = old_ids_to_nodes[a]
                real_nodes_b = old_ids_to_nodes[b]
                for x in real_nodes_a:
                    for y in real_nodes_b:
                        if directed:
                            if (x, y) in edges:
                                chosen_edges.add((x, y))
                        else:
                            (x_, y_) = (min(x, y), max(x, y))
                            if (x_, y_) in edges:
                                chosen_edges.add((x_, y_))

            # Assume the new IDs begin with the structures listed
            #
            # This seems to be the case.
            real_graph_nodes = set()
            for n in struct_nodes:
                real_graph_nodes |= old_ids_to_nodes[n]

            ids_to_nodes[idx] = real_graph_nodes
            idx = idx + 1

        # Assume all the un-chosen nodes are given the largest IDs
        #
        # This seems to be the case.
        for i in range(0, len(unselected_nodes)):
            ids_to_nodes[i + 1 + len(struct_round)] = \
                                old_ids_to_nodes[unselected_nodes[i]]
        remaining_nodenames = set([n for (n, s) in ids_to_nodes.items()])

    return chosen_edges


# Returns:
#   (struct_edges, noise_edges)
def run_PY_SUBDUE(edges, directed=False, \
                  temp_in_filename="PY_SUBDUE/testing/graph_file.txt", \
                  temp_out_filename="PY_SUBDUE/testing/output.txt"):

    nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])

    dir_str = ["False", "True"][int(directed)]

    # Prepare graph
    f = open(temp_filename, "w")
    f.write("[\n")
    for n in nodes:
        f.write("\t{ \"vertex\": {\n")
        f.write("\t\t\"id\": \"%d\"\n" % n)
        f.write("\t  }},\n")
    for i in range(0, len(edges)):
        (a, b) = edges[i]

        f.write("\t{ \"edge\": {\n")
        f.write("\t\t\"id\": \"%d\",\n" % i)
        f.write("\t\t\"source\": \"%d\",\n" % a)
        f.write("\t\t\"target\": \"%d\",\n" % b)
        f.write("\t\t\"directed\": \"%s\"\n" % dir_str)
        f.write("\t  }}")

        if i < len(edges) - 1:
            f.write(",\n")
        else:
            f.write("\n")
    f.write("]")

    # Actually run SUBDUE
    os.system("python3 SUBDUE/src/Subdue.py --minsize 2 --maxsize 10 " + \
              "--numbest 1 --overlap none --iterations 0 testing/graph_file.json " + \
              ("-out %s" % temp_out_filename))

    # Get the result.


"""
def run_GraMi(edges, directed=False, min_support="auto"):
    nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])

    relabels = {n: n for n in nodes}
    if min(nodes) != 1 or max(nodes) != len(nodes):
        nodes = sorted(list(nodes))
        relabels = {nodes[i]: i for i in range(0, len(nodes))}
        nodes = set([i for i in range(0, len(nodes))])
        edges = [(relabels[a], relabels[b]) for (a, b) in edges]

    if type(min_support) is str and min_support == "auto":
        # TODO: Update this logic
        min_support = int(len(edges) / len(nodes)) + int(len(nodes) / 2)
    elif type(min_support) is int:
        pass
    else:
        raise Exception("Error! min_support cannot have a non-integer value other than \"auto\"")

    in_filename = "GraMi_input.lg"
    nodes = sorted(list(nodes))
    f = open("GraMi/Datasets/" + in_filename, "w")
    f.write("# t 1\n")
    for n in nodes:
        f.write("v %d 1\n" % n)
    for (a, b) in sorted(list(edges)):
        f.write("e %d %d 1\n" % (a, b))
    # f.write("t # -1\n")
    f.close()

    out_filename = "/tmp/GraMi_output.txt"

    os.system("cd GraMi; ./grami -f %s -s %d -t %d -p 1 > %s" % \
                (in_filename, min_support, int(directed), out_filename))

    return (edges, set())
"""

def run_k_core(edges, directed=False, k=3):
    # The code is actually the same whether directed or undirected.
    #   The `directed` flag is included just to match a template of runners.

    nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])
    nodes = list(nodes)

    out_neighbors = {n: set() for n in nodes}
    in_neighbors = {n: set() for n in nodes}

    for (a, b) in edges:
        out_neighbors[a].add(b)
        in_neighbors[b].add(a)

    done = False
    while not done:
        done = True
        for i in range(0, len(nodes)):
            node = nodes[i]
            degree = len(out_neighbors[node]) + len(in_neighbors[node])
            if degree > 0 and degree < k:
                done = False
                for nbr in out_neighbors[node]:
                    in_neighbors[nbr].remove(node)
                out_neighbors[node] = set()
                for nbr in in_neighbors[node]:
                    out_neighbors[nbr].remove(node)
                in_neighbors[node] = set()

    struct_edges = set()
    for node in nodes:
        for nbr in out_neighbors[node]:
            struct_edges.add((node, nbr))

    noise_edges = edges - struct_edges
    return (struct_edges, noise_edges)


# Keep maximal set of edges each part of at least k - 2 triangles.
def run_k_truss(edges, directed=False, k=3):
    # The code is actually the same whether directed or undirected.
    #   The `directed` flag is included just to match a template of runners.

    if k < 3:
        k = 3

    nodes = set([a for (a, b) in edges] + [b for (a, b) in edges])

    # Ensure undirected
    und_edges = set([(min(a, b), max(a, b)) for (a, b) in edges])

    # Store all neighbors with larger IDs
    greater_neighbors = {n: set() for n in nodes}

    for (a, b) in und_edges:
        assert a != b
        greater_neighbors[a].add(b)


    # Find all (undirected) triangles

    edge_to_triangles = {}
    nodes = sorted(list(nodes))

    for a in nodes:
        nbr_list = sorted(list(greater_neighbors[a]))
        for i in range(0, len(nbr_list) - 1):
            b = nbr_list[i]

            for j in range(i + 1, len(nbr_list)):
                c = nbr_list[j]
                if (b, c) in und_edges:
                    tri_edges = [(a, b), (a, c), (b, c)]
                    triangle = (a, b, c)

                    for edge in tri_edges:
                        if edge not in edge_to_triangles:
                            edge_to_triangles[edge] = set()

                        edge_to_triangles[edge].add(triangle)

    # Purge bad (undirected) edges

    done = False
    while not done:
        done = True
        for edge in und_edges:
            if edge not in edge_to_triangles:
                continue

            triangles = edge_to_triangles[edge]
            if len(triangles) == 0 or len(triangles) >= (k - 2):
                continue

            done = False

            triangles = list(triangles)
            for (a, b, c) in triangles:
                tri_edges = [(a, b), (a, c), (b, c)]
                for tri_edge in tri_edges:
                    edge_to_triangles[tri_edge].remove((a, b, c))


    # Convert back to graph
    struct_edges = set()
    for (a, b), triangles in edge_to_triangles.items():
        assert len(triangles) == 0 or len(triangles) >= (k - 2)
        if len(triangles) == 0:
            continue
        if (a, b) in edges:
            struct_edges.add((a, b))
        if (b, a) in edges:
            struct_edges.add((b, a))

    noise_edges = edges - struct_edges
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


def limit_edges_by_nodeset(edges, nodes):
    new_edges = set()
    for (a, b) in edges:
        if a in nodes and b in nodes:
            new_edges.add((a, b))
    return new_edges
