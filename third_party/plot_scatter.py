import matplotlib.pyplot as plt
import math
import sys

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Error! Needs filename and algorithm name.")
        exit(1)

    filename = sys.argv[1]

    alg_name = sys.argv[2]

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    tok_seq = ["#Directed", "#Nodes", "#Edges", \
               "#Und. Edges", "#Struct", "#Gain", "#AR Gain"]
    directed  = []
    nodes     = []
    edges     = []
    und_edges = []
    struct    = []
    gain      = []
    ar_gain   = []
    types = [bool, int, int, int, int, float, float]
    lists = [directed, nodes, edges, und_edges, struct, gain, ar_gain]

    # Collect Information
    idx = 0
    for l in lines:
        l = l.strip()

        next_tok = tok_seq[idx]

        if l[:min(len(next_tok), len(l))] == next_tok:
            lists[idx].append(types[idx](l.split(" ")[-1].split("\t")[-1]))
            idx = (idx + 1) % len(tok_seq)

    # Plot
    fns = [(lambda d, n, e, ue, s: (e * 2) / n), \
           (lambda d, n, e, ue, s: ((e * (int(not d) + 1)) / (n * (n - 1)))), \
           (lambda d, n, e, ue, s: \
            math.log10(((e * (int(not d) + 1)) / (n * (n - 1))))), \
           (lambda d, n, e, ue, s: \
            math.log10(((ue * 2) / (n * (n - 1))))), \
           (lambda d, n, e, ue, s: ((s + 1) * 2) / n), \
           (lambda d, n, e, ue, s: \
            (((s + 1) * (int(not d) + 1)) / (n * (n - 1)))), \
           (lambda d, n, e, ue, s: \
            math.log10((((s + 1) * (int(not d) + 1)) / (n * (n - 1)))))]
    terms = ["Avg. Degree", "Density", "Log Density", "Und. Log Density", \
             "Rem. Avg. Degree", "Rem. Density", "Log Rem. Density"]

    indices = [1, 2, 3, 5, 6]

    gain_aug    = list(gain)
    ar_gain_aug = list(ar_gain)
    for i in range(0, len(gain_aug)):
        if gain_aug[i] == 0.0:
            gain_aug[i] = 0.0001
        if ar_gain_aug[i] == 0.0:
            ar_gain_aug[i] = 0.0001

    log10_gain    = [(g / abs(g)) * math.log10(abs(g)) for g in gain_aug]
    log10_ar_gain = [(g / abs(g)) * math.log10(abs(g)) for g in ar_gain_aug]

    for idx in indices:
        fn = fns[idx]
        plt.figure(figsize=(9,6))
        x_axis = [fn(directed[i], nodes[i], edges[i], und_edges[i], struct[i]) \
                        for i in range(0, len(nodes))]
        plt.scatter(x_axis, gain, label=alg_name, marker="x", c="#000000",  s=100)
        # plt.scatter(x_axis, gain, label=alg_name, marker="o", c="#BBBBBB",  s=80)
        # plt.scatter(x_axis, ar_gain, label="random", marker="x", c="#000000", s=100)
        # plt.legend()
        plt.title("Effectiveness of %s vs. %s" % (alg_name, terms[idx]), size=18)
        plt.xlabel(terms[idx], size=14)
        plt.ylabel("Gain above all Structure", size=14)
        plt.show()
