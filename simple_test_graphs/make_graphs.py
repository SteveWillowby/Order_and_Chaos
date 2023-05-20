def write_graph(edges, base_name):
    nodes = sorted(list(set([a for (a, b) in edges] + \
                            [b for (a, b) in edges])))

    f = open("%s_%d_nodes.txt" % (base_name, len(nodes)), "w")
    for node in nodes:
        f.write("%d\n" % node)
    f.close()
    f = open("%s_%d_edges.txt" % (base_name, len(nodes)), "w")
    for (a, b) in edges:
        f.write("%d %d\n" % (a, b))
    f.close()

# n numbers between 0 and n
def johnson_indices_to_number(indices, n):
    num = 0
    factor = 1
    for idx in indices:
        num += idx * factor
        factor *= n
    return num

def less_all(a, l):
    for b in l:
        if a >= b:
            return False
    return True

def johnson_node_indices(n, k):
    if k == 1:
        return [(i, ) for i in range(0, n)]
    sub = johnson_node_indices(n, k - 1)
    indices = []
    for i in range(0, n):
        indices += [tuple([i] + list(x)) for x in sub if less_all(i, x)]
    return indices

if __name__ == "__main__":
    NUM_NODES = 128

    TEST_ONLY = False

    # Binary Tree
    edges = []
    total_n = 1
    layer = 1
    while total_n + layer * 2 <= NUM_NODES:
        for i in range(total_n + 1 - layer, total_n + 1):
            edges.append((i, i * 2))
            edges.append((i, (i * 2) + 1))
        layer *= 2
        total_n += layer

    if not TEST_ONLY:
        write_graph(edges, "binary_tree")

    # Ring
    edges = [(i, (i + 1) % NUM_NODES) for i in range(0, NUM_NODES)]
    if not TEST_ONLY:
        write_graph(edges, "ring")

    # Wreath
    edges = []
    degree_right = 4
    degree = -1 + 2 * degree_right
    for i in range(0, NUM_NODES):
        for d in range(1, degree_right + 1):
            edges.append((i, (i + d) % NUM_NODES))

    edges = [(i, (i + 1) % NUM_NODES) for i in range(0, NUM_NODES)]
    if not TEST_ONLY:
        write_graph(edges, "wreath_d%d" % degree)

    # Johnson Graph
    edges = []
    k = 3
    candidate_n = 1
    prev_n_choose_k = None
    while True:
        n_choose_k = 1
        for i in range(0, k):
            n_choose_k *= (candidate_n - i)
        for i in range(1, k + 1):
            n_choose_k /= i

        if n_choose_k > NUM_NODES:
            candidate_n -= 1
            break
        else:
            candidate_n += 1

        prev_n_choose_k = n_choose_k

    # print("%d choose %d is %d" % (candidate_n, k, prev_n_choose_k))

    jni = johnson_node_indices(candidate_n, k)
    for x in jni:
        for idx in range(0, k):
            y = x[idx]
            for value in range(y + 1, candidate_n):
                in_tuple = False
                for z in ((list(x))[idx + 1:]):
                    if value == z:
                        in_tuple = True
                        break
                if not in_tuple:
                    next_tuple = list(x)
                    next_tuple[idx] = value
                    next_tuple = tuple(sorted(next_tuple))
                    edges.append((johnson_indices_to_number(x, candidate_n), \
                                  johnson_indices_to_number(next_tuple, candidate_n)))

    if not TEST_ONLY:
        write_graph(edges, "johnson_%d_%d" % (candidate_n, k))
