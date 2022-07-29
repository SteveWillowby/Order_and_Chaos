def edge_difference(nc_A, nc_B, directed):
    s = 0
    for i in range(0, len(nc_B)):
        for j in nc_B[i]:
            if (not directed) and j > i:
                continue
            if j not in nc_A[i]:
                s += 1
    return s

if __name__ == "__main__":
    from enumerate_graphs import AdjEnumerator
    from subgraph_scorer import subgraph_structure_score, edge_list_to_neighbors_collections

    from py_NT_session import PyNTSession

    graph_count_values = {False: {1: 1, 2: 2, 3: 8, 4: 64, 5: 1024, 6: 32768, 7: 2097152, 8: 268435456}, \
                           True: {1: 1, 2: 4, 3: 64, 4: 4096, 5: 1048576, 6: 1073741824}}

    directed = False

    test_names = []
    test_graphs = []

    test_names.append("5-chain on six nodes")
    test_graphs.append([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)])

    test_names.append("5-clique on five nodes")
    test_graphs.append([(0, 1), (0, 2), (0, 3), (0, 4), \
                        (1, 2), (1, 3), (1, 4), \
                        (2, 3), (2, 4), \
                        (3, 4)])

    test_names.append("grid on 9 nodes")
    test_graphs.append([(0, 1), (1, 2), \
                        (3, 4), (4, 5), \
                        (6, 7), (7, 8), \
                        (0, 3), (3, 6), \
                        (1, 4), (4, 7), \
                        (2, 5), (5, 8)]) 

    test_names.append("grid on 9 nodes with one corner edge missing")
    test_graphs.append(list(test_graphs[-1]))
    test_graphs[-1].pop()

    test_names.append("peterson graph")
    test_graphs.append([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), \
                        (5, 6), (6, 7), (7, 8), (8, 9), (9, 5), \
                        (0, 5), (1, 7), (2, 9), (3, 6), (4, 8)])

    for i in range(0, len(test_graphs)):
        edge_list = test_graphs[i]
        test_name = test_names[i]

        best_score = None
        best_NCs = [None, None, None, None, None]

        nc_A = edge_list_to_neighbors_collections(edge_list, directed=directed)

        graph_count_max = graph_count_values[directed][len(nc_A)]
        percent_done = 0

        print("\n\nRunning for %s" % test_name)

        first = True
        n = len(nc_A)
        AE = AdjEnumerator(n, directed)
        while True:
            nc_B = AE.next_graph()
            if nc_B is None:
                break

            if first:  # skip the empty graph
                first = False
                continue

            count = AE.graph_count()
            if int((count * 100.0) / graph_count_max) > percent_done:
                percent_done = int((count * 100.0) / graph_count_max)
                if percent_done % 10 == 0:
                    print("%d percent done -- %d" % (percent_done, count))

            # print(nc_B)

            score = subgraph_structure_score(n, directed, nc_B, auto_solver_class=PyNTSession)
            if best_score is None or best_score < score:
                score -= edge_difference(nc_A, nc_B, directed)
                if best_score is None or best_score < score:
                    best_score = score
                    best_NCs = [(score, nc_B)] + [best_NCs[i] for i in range(1, len(best_NCs))]

        print("Original graph:")
        print(nc_A)
        print("Top Graphs:")
        for (score, nc_B) in best_NCs:
            print("\nScore: %f -- Graph: %s" % (score, nc_B))
