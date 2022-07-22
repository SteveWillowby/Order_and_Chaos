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
    return subgraph_surprisal(n, directed, edge_list, auto_solver_class) - \
           expected_surprisal_estimate(n, m, directed, auto_solver_class)

# Assumes there are no duplicate edges in edge_list
def subgraph_surprisal(n, directed, edge_list, auto_solver_class):
    m = len(edge_list)

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

    if directed:
        max_m = n * (n - 1)
    else:
        max_m = int((n * (n - 1)) / 2)

    log2_full_factorial = log2_factorial(n)
    log2_edge_choice = log2_a_choose_b(max_m, m)

    num_singletons = n - len(nc)
    log2_aut = log2_automorphisms(directed, has_edge_types=False, \
                                  neighbors_collections=nc, \
                                  auto_solver_class=auto_solver_class)
    log2_aut += log2_factorial(num_singletons)

    return -((log2_full_factorial - log2_aut) - log2_edge_choice)

# AKA Entropy Estimate
def expected_surprisal_estimate(n, m, directed, auto_solver_class):
    if (n, m, directed) in NMD_SUBGRAPH_EXPECTED_SURPRISAL:
        return NMD_SUBGRAPH_EXPECTED_SURPRISAL[(n, m, directed)]

    if directed:
        max_m = n * (n - 1)
    else:
        max_m = int((n * (n - 1)) / 2)

    log2_full_factorial = log2_factorial(n)
    log2_edge_choice = log2_a_choose_b(max_m, m)

    surprisals = []
    NUM_RAND_GRAPHS = 1000
    for _ in range(0, NUM_RAND_GRAPHS):
        nc = ER_graph(n, m, directed)
        nc = graph_without_singletons(nc)
        num_singletons = n - len(nc)

        log2_aut = log2_automorphisms(directed, has_edge_types=False, \
                                      neighbors_collections=nc, \
                                      auto_solver_class=auto_solver_class)
        log2_aut += log2_factorial(num_singletons)

        log2_ER_prob = (log2_full_factorial - log2_aut) - log2_edge_choice

        surprisals.append(-log2_ER_prob)

    expected_surprisal = sum(surprisals) / len(surprisals)

    NMD_SUBGRAPH_EXPECTED_SURPRISAL[(n, m, directed)] = expected_surprisal
    return expected_surprisal

# Used as a reserve bank.
NMD_SUBGRAPH_EXPECTED_SURPRISAL = {}

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
        new_nodes -= node_pools

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

    return neighbors_collections

if __name__ == "__main__":
    from py_NT_session import PyNTSession

    print(expected_surprisal_estimate(n=6, m=6, directed=False, \
                                      auto_solver_class=PyNTSession))

    six_ring = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    six_ring_with_chord = six_ring + [(0, 3)]
    eight_ring = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 0)]
    print(subgraph_structure_score(n=20, directed=False, \
                                   edge_list=six_ring, \
                                   auto_solver_class=PyNTSession))
    print(subgraph_structure_score(n=20, directed=False, \
                                   edge_list=six_ring_with_chord, \
                                   auto_solver_class=PyNTSession))
    print(subgraph_structure_score(n=20, directed=False, \
                                   edge_list=eight_ring, \
                                   auto_solver_class=PyNTSession))
