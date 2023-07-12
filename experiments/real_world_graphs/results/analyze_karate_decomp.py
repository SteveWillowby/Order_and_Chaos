

if __name__ == "__main__":
    f = open("../backup_results/karate_graph.txt", "r")
    edges = f.readlines()
    f.close()

    edges = [tuple(l.strip().split(" ")) for l in edges]

    f = open("../karate_labels.txt", "r")
    membership = f.readlines()
    f.close()

    membership = [l.strip().split(" ") for l in membership]
    membership = {x[0]: x[1] for x in membership}

    degrees = {n: 0 for (n, mem) in membership.items()}
    for (a, b) in edges:
        degrees[a] += 1
        degrees[b] += 1

    degrees = sorted([(d, n) for (n, d) in degrees.items()], reverse=True)
    print(degrees[:4])
    top_4_nodes = [n for (d, n) in degrees[:4]]

    x = top_4_nodes[0]
    y = None
    for z in top_4_nodes[1:]:
        if (x, z) in edges or (z, x) in edges:
            y = z
            break
    assert y is not None

    group_1 = [x, y]
    group_2 = []
    for q in top_4_nodes:
        if q not in group_1:
            group_2.append(q)

    group_assignments = {group_1[0]: membership[group_1[0]], \
                         group_1[1]: membership[group_1[0]], \
                         group_2[0]: membership[group_2[0]], \
                         group_2[1]: membership[group_2[1]]}

    for (a, b) in edges:
        if (a in group_1):
            group_assignments[b] = membership[group_1[0]]
        elif (b in group_1):
            group_assignments[a] = membership[group_1[0]]
        elif (a in group_2):
            group_assignments[b] = membership[group_2[0]]
        elif (b in group_2):
            group_assignments[a] = membership[group_2[0]]
        else:
            raise Exception("Error! Neither %s nor %s are group hubs " % (a, b))

    agreed = 0
    disagreed = 0
    for (x, mem) in membership.items():
        if group_assignments[x] == mem:
            agreed += 1
        else:
            disagreed += 1
            print("Disagreed on %s" % x)

    print("Percent Agreed: %f" % (agreed / (agreed + disagreed)))
