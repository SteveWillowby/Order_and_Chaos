import math
import random

from prob_calc_utils import *

# Assumes there are no duplicate edges in edge_list
#
# Returns how many more bits this subgraph surprises us with than we'd expect
#   a random subgraph of the same size to surprise us with.
#
# Currently this suffers when n is large.
#
# This function assumes that the m edges forming the subgraph you are looking
#   at were randomly sampled from your overall set of edges, which in turn is
#   assumed to be an ER graph.
#
# Things might work better if I had a way to calculate P(occ at least once)
def subgraph_structure_score(n, directed, edge_list, auto_solver_class):
    m = len(edge_list)
    nc = edge_list_to_neighbors_collections(edge_list, directed)

    # fn = subgraph_prob
    fn = subgraph_surprisal
    # fn = expected_subgraph_count

    reverse_divide = (lambda x: x[1] / x[0])
    divide = (lambda x: x[0] / x[1])
    subtract = (lambda x: x[0] - x[1])
    fn_only = (lambda x: x[0])

    combiner = subtract

    fn_value = fn(n, m, nc, directed, auto_solver_class)
    expected_fn_value = expected_fn_estimate(n, m, directed, fn, auto_solver_class)
    # expected_fn_value = expected_number_of_subgraphs(n, m, directed)
    # print("%f vs %f" % (fn_value, expected_fn_value))

    return combiner((fn_value, expected_fn_value))
    # return 2.0 ** -subgraph_surprisal(n, m, nc, directed, auto_solver_class)

# The expected number of times you would find an m-edge pattern isomorphic
#   to this m-edge pattern in an 0.5 edge-prob ER graph on n nodes.
def expected_subgraph_count(n, m, nc, directed, auto_solver_class):
    nc = graph_without_singletons(nc)

    log2_full_factorial = log2_factorial(n)
    log2_edge_choice = float(-m)  # log2(0.5^m)

    num_singletons = n - len(nc)
    log2_aut = log2_automorphisms(directed, has_edge_types=False, \
                                  neighbors_collections=nc, \
                                  auto_solver_class=auto_solver_class)
    log2_aut += log2_factorial(num_singletons)

    # log2( P(a given m-edge set) * # of relevant m-edge sets )
    log2_expected_count = log2_edge_choice + (log2_full_factorial - log2_aut)

    expected_count = 2.0 ** log2_expected_count
    assert expected_count > 0
    return expected_count

def expected_number_of_subgraphs(n, m, directed):
    if directed:
        max_m = n * (n - 1)
    else:
        max_m = int((n * (n - 1)) / 2)

    total_expected = 0.0
    for curr_m in range(m, max_m + 1):
        # log2(Prob there are curr_m edges)
        log2_p_curr_m = (-max_m) + log2_a_choose_b(max_m, curr_m)
        # log2(# m-edge sets among curr_m edges)
        log2_m_in_curr_m = log2_a_choose_b(curr_m, m)

        total_expected += 2.0 ** (log2_p_curr_m + log2_m_in_curr_m)

    return total_expected

# Negative log-probability that if you randomly sample |edge_list| edges on
#   n nodes that you get an edge list isomorphic to this one.
def subgraph_surprisal(n, m, nc, directed, auto_solver_class):
    nc = graph_without_singletons(nc)

    if directed:
        max_m = n * (n - 1)
    else:
        max_m = int((n * (n - 1)) / 2)

    log2_full_factorial = log2_factorial(n)
    log2_edge_choice = log2_a_choose_b(max_m, m)
    # log2_edge_choice = float(-m)  # log2(0.5^m)
    # log2_edge_choice = 0.0

    num_singletons = n - len(nc)
    log2_aut = log2_automorphisms(directed, has_edge_types=False, \
                                  neighbors_collections=nc, \
                                  auto_solver_class=auto_solver_class)
    log2_aut += log2_factorial(num_singletons)

    return -((log2_full_factorial - log2_aut) - log2_edge_choice)

# Probability that if you randomly sample |edge_list| edges on
#   n nodes that you get an edge list isomorphic to this one.
def subgraph_prob(n, m, nc, directed, auto_solver_class):
    return 2.0 ** -subgraph_surprisal(n, m, nc, directed, auto_solver_class)

def expected_fn_estimate(n, m, directed, fn, auto_solver_class):
    if (n, m, directed, fn) in NMD_SUBGRAPH_EXPECTED_FN_VALUES:
        return NMD_SUBGRAPH_EXPECTED_FN_VALUES[(n, m, directed, fn)]

    values = []
    TARGET_NUM_SAMPLES = 1000000

    CHUNK_SIZE = int(math.sqrt(TARGET_NUM_SAMPLES))
    NUM_CHUNKS = CHUNK_SIZE

    for _ in range(0, NUM_CHUNKS):
        chunk_values = []
        for _ in range(0, CHUNK_SIZE):
            nc = ER_graph(n, m, directed)
            chunk_values.append(fn(n, m, nc, directed, auto_solver_class))
        values.append(sum(chunk_values) / len(chunk_values))

    expected_value = sum(values) / len(values)

    NMD_SUBGRAPH_EXPECTED_FN_VALUES[(n, m, directed, fn)] = expected_value
    return expected_value

# Used as a reserve bank.
NMD_SUBGRAPH_EXPECTED_FN_VALUES = {}

def edge_list_to_neighbors_collections(edge_list, directed):
    nodes = set([a for (a, b) in edge_list] + [b for (a, b) in edge_list])
    nodes = sorted(list(nodes))
    relabeled = {nodes[i]: i for i in range(0, len(nodes))}

    nc = [set() for _ in nodes]
    if directed:
        for (a, b) in edge_list:
            nc[relabeled[a]].add(relabeled[b])
    else:
        for (a, b) in edge_list:
            nc[relabeled[a]].add(relabeled[b])
            nc[relabeled[b]].add(relabeled[a])
    return nc

def graph_without_singletons(neighbors_collections):
    nodes = set()
    for nc in neighbors_collections:
        nodes |= nc

    nodes = sorted(list(nodes))
    relabeled = {nodes[i]: i for i in range(0, len(nodes))}

    return [set([relabeled[nbr] for nbr in neighbors_collections[node]]) \
                    for node in nodes]

def connected_components(neighbors_collections, directed):

    n = len(neighbors_collections)

    if directed:
        # Create an undirected version.
        orig_neighbors_collections = neighbors_collections
        neighbors_collections = [set(s) for s in neighbors_collections]
        for a in range(0, n):
            for b in neighbors_collections[a]:
                if b != a:
                    neighbors_collections[b].add(a)

    node_pools = []
    new_nodes = set(range(0, n))
    while len(new_nodes) > 0:
        new_node = new_nodes.pop()
        frontier = set([new_node])
        visited = set([new_node])
        while len(frontier) > 0:
            new_frontier = set()
            for node in frontier:
                new_frontier |= (neighbors_collections[node] - visited)
            visited |= new_frontier
            frontier = new_frontier
        node_pools.append(visited)
        new_nodes -= visited

    new_NCs = []
    if directed:
        del neighbors_collections
        neighbors_collections = orig_neighbors_collections

    for node_pool in node_pools:
        node_pool = sorted(list(node_pool))
        relabeled = {node_pool[i]: i for i in range(0, len(node_pool))}
        new_NCs.append([set([relabeled[nbr] for nbr in neighbors_collections[node]]) for node in node_pool])
    return new_NCs

def ER_graph(n, m, directed):
    if directed:
        max_m = n * (n - 1)
    else:
        max_m = int((n * (n - 1)) / 2)

    if m > int(max_m / 2):
        target_m = max_m - m
    else:
        target_m = m

    edges = set()
    if directed:
        while len(edges) < target_m:
            a = random.randint(0, n - 1)
            b = random.randint(0, n - 2)
            b += int(b >= a)
            edges.add((a, b))
    else:
        while len(edges) < target_m:
            a = random.randint(0, n - 1)
            b = random.randint(0, n - 2)
            b += int(b >= a)
            edges.add((min(a, b), max(a, b)))

    if target_m < m:
        neighbors_collections = [set(range(0, n)) for _ in range(0, n)]
        if directed:
            for (a, b) in edges:
                neighbors_collections[a].remove(b)
        else:
            for (a, b) in edges:
                neighbors_collections[a].remove(b)
                neighbors_collections[b].remove(a)
    else:
        neighbors_collections = [set() for _ in range(0, n)]
        if directed:
            for (a, b) in edges:
                neighbors_collections[a].add(b)
        else:
            for (a, b) in edges:
                neighbors_collections[a].add(b)
                neighbors_collections[b].add(a)

    del edges

    # if len(connected_components(graph_without_singletons(neighbors_collections), directed)) > 1:
    #     return ER_graph(n, m, directed)

    return neighbors_collections

if __name__ == "__main__":
    from py_NT_session import PyNTSession

    # print(expected_surprisal_estimate(n=6, m=6, directed=False, \
    #                                   auto_solver_class=PyNTSession))

    wedge_nc = [set([1, 2]), set([0]), set([0])]
    print(expected_subgraph_count(n=3, m=2, directed=False, nc=wedge_nc, auto_solver_class=PyNTSession))
    print(expected_number_of_subgraphs(n=3, m=2, directed=False))
    print(6.0 / 8.0)

    six_ring = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    six_ring_with_chord = six_ring + [(0, 3)]
    eight_ring = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)]
    print("six-ring")
    print(subgraph_structure_score(n=6, directed=False, \
                                   edge_list=six_ring, \
                                   auto_solver_class=PyNTSession))
    print("six-ring with chord")
    print(subgraph_structure_score(n=6, directed=False, \
                                   edge_list=six_ring_with_chord, \
                                   auto_solver_class=PyNTSession))

    print("Setting n = 8")

    print("eight-ring")
    print(subgraph_structure_score(n=8, directed=False, \
                                   edge_list=eight_ring, \
                                   auto_solver_class=PyNTSession))
