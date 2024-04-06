import matplotlib.pyplot as plt
import math
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error! Needs filename.")
        exit(1)

    filename = sys.argv[1]

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    tok_seq = ["#Directed", "#Nodes", "#Edges", \
               "#Struct", "#Gain", "#AR Gain"]
    directed = []
    nodes    = []
    edges    = []
    struct   = []
    gain     = []
    ar_gain  = []
    types = [bool, int, int, int, float, float]
    lists = [directed, nodes, edges, struct, gain, ar_gain]

    # Collect Information
    idx = 0
    for l in lines:
        l = l.strip()

        next_tok = tok_seq[idx]

        if l[:min(len(next_tok), len(l))] == next_tok:
            lists[idx].append(types[idx](l.split(" ")[-1].split("\t")[-1]))
            idx = (idx + 1) % len(tok_seq)

    # Plot
    fns = [(lambda d, n, e, s: (e * 2) / n), \
           (lambda d, n, e, s: ((e * (int(not d) + 1)) / (n * (n - 1)))), \
           (lambda d, n, e, s: \
            math.log10(((e * (int(not d) + 1)) / (n * (n - 1))))), \
           (lambda d, n, e, s: ((s + 1) * 2) / n), \
           (lambda d, n, e, s: \
            (((s + 1) * (int(not d) + 1)) / (n * (n - 1)))), \
           (lambda d, n, e, s: \
            math.log10((((s + 1) * (int(not d) + 1)) / (n * (n - 1)))))]
    terms = ["Avg. Degree", "Density", "Log Density", \
             "Rem. Avg. Degree", "Rem. Density", "Log Rem. Density"]

    for idx in range(0, len(fns)):
        fn = fns[idx]
        plt.figure(figsize=(9,6))
        x_axis = [fn(directed[i], nodes[i], edges[i], struct[i]) \
                        for i in range(0, len(nodes))]
        plt.scatter(x_axis, gain, label="algorithm", marker="o", c="#BBBBBB",  s=80)
        plt.scatter(x_axis, ar_gain, label="random", marker="x", c="#000000", s=100)
        plt.legend()
        plt.title("Effectiveness of Algorithm vs. Graph Property", size=18)
        plt.xlabel(terms[idx], size=14)
        plt.ylabel("Gain above all Structure", size=14)
        plt.show()
