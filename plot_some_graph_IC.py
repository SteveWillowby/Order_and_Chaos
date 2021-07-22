from graph_info_content import information_content, min_information_content_limit, max_information_content_limit
import matplotlib.pyplot as plt
import random


# graph_sequence should be a list of (nodes, edges) tuples.
def plot_graph_sequence(graph_sequence, \
                        directed=False, temporal=False, \
                        sequence_name="A Graph Sequence"):
    graph_indices = [i for i in range(0, len(graph_sequence))]
    ic = []
    min_ic_limit = []
    max_ic_limit = []
    for (nodes, edges) in graph_sequence:
        ic.append(information_content(nodes, edges, proportional_p=False, \
                        directed=directed, temporal=temporal, \
                        num_timestamps="auto", \
                        node_colors=None, only_consider_used_nodes=False, \
                        extreme_only_consider_used_nodes=False, \
                        extreme_proportional_p=False))

        num_timestamps = 1
        if temporal:
            timestamps = set()
            for (a, b, t) in edges:
                timestamps.add(t)
            num_timestamps = len(timestamps)
        min_ic_limit.append(\
            min_information_content_limit(len(nodes), \
                                          directed=directed, \
                                          num_timestamps=num_timestamps))
        max_ic_limit.append(\
            max_information_content_limit(len(nodes), \
                                          directed=directed, \
                                          num_timestamps=num_timestamps))

    # plt.clear()
    plt.plot(graph_indices, min_ic_limit, color="red")
    plt.plot(graph_indices, max_ic_limit, color="red")
    plt.plot(graph_indices, ic, color="blue")
    plt.title("Information Content of %s" % sequence_name)
    plt.xlabel("Graph Sequence Index")
    plt.ylabel("Information Content")
    plt.show()

def __triangles_sequence__(sequence_length):
    sequence = []
    nodes_list = []
    edges_list = []
    for i in range(0, sequence_length):
        nodes_list.append(i*3)
        nodes_list.append(i*3 + 1)
        nodes_list.append(i*3 + 2)
        edges_list.append((i*3, i*3 + 1))
        edges_list.append((i*3 + 1, i*3 + 2))
        edges_list.append((i*3 + 2, i*3))
        sequence.append((list(nodes_list), list(edges_list)))
    return sequence

def __triangles_sequence_with_all_nodes_always__(sequence_length):
    sequence = []
    nodes_list = []
    edges_list = []
    for i in range(0, sequence_length):
        nodes_list.append(i*3)
        nodes_list.append(i*3 + 1)
        nodes_list.append(i*3 + 2)
        edges_list.append((i*3, i*3 + 1))
        edges_list.append((i*3 + 1, i*3 + 2))
        edges_list.append((i*3 + 2, i*3))
        sequence.append((nodes_list, list(edges_list)))
    return sequence

def __random_edge_addition__(num_nodes, edges_per_iter, sequence_length, \
                             directed=False):
    sequence = []
    nodes = [i for i in range(0, num_nodes)]
    edges_set = set()
    for _ in range(0, sequence_length):
        num_added = 0
        while num_added < edges_per_iter:
            a = random.randint(0, num_nodes - 1)
            b = random.randint(0, num_nodes - 1)
            if a == b:
                continue
            if (a, b) not in edges_set and \
                    (directed or (b, a) not in edges_set):
                edges_set.add((a, b))
                num_added += 1
        sequence.append((nodes, list(edges_set)))
    return sequence

def __binary_tree_sequence__(num_trees):
    sequence = []
    nodes_list = [1]
    edges_list = []
    num_nodes = 1
    sequence.append(([1], []))
    for _ in range(0, num_trees - 1):
        new_num_nodes = num_nodes * 2 + 1
        for n in range(num_nodes + 1, new_num_nodes + 1):
            nodes_list.append(n)
            edges_list.append((int(n / 2), n))
        sequence.append((list(nodes_list), list(edges_list)))
        num_nodes = new_num_nodes
    return sequence

if __name__ == "__main__":

    plot_graph_sequence(__triangles_sequence__(9), \
                        directed=False, temporal=False, \
                        sequence_name="Triangles Sequence")
    plot_graph_sequence(__triangles_sequence_with_all_nodes_always__(9), \
                        directed=False, temporal=False, \
                        sequence_name="Triangles Sequence with All Nodes")
    plot_graph_sequence(__random_edge_addition__(num_nodes=9*3, \
                                                 edges_per_iter=3, \
                                                 sequence_length=9, \
                                                 directed=False), \
                        directed=False, temporal=False, \
                        sequence_name="Randomly Adding 3 Edges Each Time")

    plot_graph_sequence(__binary_tree_sequence__(num_trees=7), \
                        directed=False, temporal=False, \
                        sequence_name="Binary Trees")
