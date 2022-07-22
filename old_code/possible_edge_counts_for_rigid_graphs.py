
# Begin with asymmetric trees.
def tree_counts():
    n_m_combos = {1: {0: 1}, 7: {6: 1}, 8: {7: 1}, 9: {8: 3}, 10: {9: 6}, \
                  11: {10: 15}, 12: {11: 29}, 13: {12: 67}, 14: {13: 139}, \
                  15: {14: 310}, 16: {15: 667}, 17: {16: 1480}, 18: {17: 3244},\
                  19: {18: 7241}, 20: {19: 16104}, 21: {20: 36192}, \
                  22: {21: 81435}, 23: {22: 184452}, 24: {23: 418870}, \
                  25: {24: 955860}, 26: {25: 2187664}, 27: {26: 5025990}, \
                  28: {27: 11580130}, 29: {28: 26765230}, 30: {29: 62027433}, \
                  31: {30: 144133676}, 32: {31: 335731381}, \
                  33: {32: 783859852}, 34: {33: 1834104934}, \
                  35: {34: 4300433063}, 36: {35: 10102854473}, \
                  37: {36: 23778351222} \
                  }
    return n_m_combos

def add_complements(combos):
    for n, m_dict in combos.items():
        max_m = int((n * (n - 1)) / 2)
        l = [(m, count) for m, count in m_dict.items()]
        for (m, count) in l:
            if max_m - m not in m_dict:
                m_dict[max_m - m] = count

def foundation():
    d = tree_counts()
    d[6] = {6: 1}
    add_complements(d)
    return d

def min_m_for_rigid(num_nodes):
    if num_nodes == 1:
        return 0
    elif num_nodes == 6:
        return 6
    assert num_nodes >= 7

    total_m = 0
    total_n = 0
    list_form = [(n, m_dict[n - 1]) for n, m_dict in tree_counts().items()]
    list_form.sort()
    for n, c in list_form:
        if n == 1:
            next_n = 7
        else:
            next_n = n + 1

        for i in range(0, c):
            if total_n + n == num_nodes:
                # We reach the total exactly.
                total_m += n - 1
                return total_m
            elif i < c - 1 and total_n + 2 * n > num_nodes:
                # There are at least two unused n-graphs.
                # However, adding both would get us over the target.
                total_m += (num_nodes - total_n) - 1
                return total_m
            elif i == c - 1 and total_n + n + next_n > num_nodes:
                # There is just one unused n-graph.
                # Adding it and the next graph would take us over the target.
                total_m += (num_nodes - total_n) - 1
                return total_m
            else:
                # Not near the total yet.
                total_n += n
                total_m += n - 1
                assert total_n < num_nodes
    raise ValueError(("Error! Cannot find min_m for a rigid %d" % num_nodes) + \
                     "-node graph. Add larger options to tree_counts()")

def min_m_for_rigid_range(max_num_nodes, multiples_of=1):
    if max_num_nodes == 1:
        if multiples_of == 1:
            return [0]
        else:
            return []
    if multiples_of == 1:
        list_of_edge_counts = [0, 6]
    elif multiples_of in [2, 3, 6]:
        list_of_edge_counts = [6]
    else:
        list_of_edge_counts = []
    if max_num_nodes == 6:
        return list_of_edge_counts

    assert max_num_nodes >= 7

    total_m = 0
    total_n = 0
    target_n = 7
    list_form = [(n, m_dict[n - 1]) for n, m_dict in tree_counts().items()]
    list_form.sort()
    for n, c in list_form:
        if n == 1:
            next_n = 7
        else:
            next_n = n + 1

        for i in range(0, c):
            done = False
            while not done:
                if total_n + n == target_n:
                    # We reach the total exactly.
                    if target_n % multiples_of == 0:
                        list_of_edge_counts.append(total_m + n - 1)
                    target_n += 1
                elif i < c - 1 and total_n + 2 * n > target_n:
                    # There are at least two unused n-graphs.
                    # However, adding both would get us over the target.
                    if target_n % multiples_of == 0:
                        list_of_edge_counts.append(total_m + \
                                                   (target_n - total_n) - 1)
                    target_n += 1
                elif i == c - 1 and total_n + n + next_n > target_n:
                    # There is just one unused n-graph.
                    # Adding it and the next graph would take us over the target.
                    if target_n % multiples_of == 0:
                        list_of_edge_counts.append(total_m + \
                                                   (target_n - total_n) - 1)
                    target_n += 1
                else:
                    # Not near the total yet.
                    total_n += n
                    total_m += n - 1
                    assert total_n < target_n
                    done = True
                if target_n > max_num_nodes:
                    return list_of_edge_counts
    raise ValueError(("Error! Cannot find min_m for a rigid %d" % target_n) + \
                     "-node graph. Add larger options to tree_counts()")

def paired_up(start_dict, node_max=None):
    new_dict = foundation()

    old_total_combos = sum([len(m_dict) for n, m_dict in start_dict.items()])
    new_total_combos = sum([len(m_dict) for n, m_dict in new_dict.items()])

    list_form = [(n, m_dict) for n, m_dict in start_dict.items()]
    list_form.sort(reverse=True)
    for i in range(0, len(list_form)):
        for j in range(i, len(list_form)):
            (n1, m_dict_1) = list_form[i]
            (n2, m_dict_2) = list_form[j]
            n3 = n1 + n2
            if node_max is not None and n3 > node_max:
                continue
            sub_list_form_1 = [(m, c) for m, c in m_dict_1.items()]
            sub_list_form_2 = [(m, c) for m, c in m_dict_2.items()]
            for k in range(0, len(sub_list_form_1)):
                for l in range(int(i == j) * k, len(sub_list_form_2)):
                    (m1, c1) = sub_list_form_1[k]
                    (m2, c2) = sub_list_form_2[l]
                    num_combos = c1 * c2
                    if i == j and k == l:
                        num_combos -= 1
                    m3 = m1 + m2
                    max_m3 = int((n3*(n3 - 1)) / 2)
                    if n3 not in new_dict:
                        new_dict[n3] = {}
                    if m3 not in new_dict[n3]:
                        new_total_combos += 1
                        new_dict[n3][m3] = 0
                    new_dict[n3][m3] += num_combos
                    complement_edges = max_m3 - m3
                    # TODO: Loosen this condition
                    if complement_edges != m3:
                        if complement_edges not in new_dict[n3]:
                            new_total_combos += 1
                            new_dict[n3][complement_edges] = 0
                        new_dict[n3][complement_edges] += num_combos

    assert new_total_combos == sum([len(m_dict) for n, m_dict in new_dict.items()])
    print("\nStarting Combos Before: %d" % old_total_combos)
    print("Starting Combos After:  %d" % new_total_combos)
    return new_dict

def max_hypothetically_possible_combos(max_n):
    assert max_n == 1 or max_n >= 6
    mins = min_m_for_rigid_range(max_n, multiples_of=1)
    total = 0
    for i in range(0, len(mins)):
        if i == 0:
            n = 1
        else:
            n = i + 5
        max_m = int((n * (n - 1)) / 2)
        total += (max_m - mins[i]) + 1 - mins[i]
    return total

if __name__ == "__main__":
    """
    skip_size = 10000
    mins = min_m_for_rigid_range(skip_size * 1000, multiples_of=skip_size)
    for i in range(0, len(mins)):
        n = (i + 1) * skip_size
        print("%d: \t%d, \t%f" % (n, mins[i], n / mins[i]))
    """

    combos = {}
    combos = paired_up(combos, node_max=100)
    combos = paired_up(combos, node_max=100)
    combos = paired_up(combos, node_max=100)
    combos = paired_up(combos, node_max=100)
    # TODO: Find the bug.
    #
    # I believe it is because paired_up() forgets that some combos already used
    #   up all the available graphs of a given node count. For example, in the
    #   first call, it may combine a 7-node tree graph with a 10-node graph;
    #   then, it the second call, it may combine that 17-node graph with the 7-
    #   node tree graph again, even though there's only one 7-node tree graph.
    print(max_hypothetically_possible_combos(100))
