import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error! Needs filename.")
        exit(1)

    filename = sys.argv[1]

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    tok_seq = ["#Nodes", "#Edges", "#Struct", "#Gain", "#AR Gain"]
    nodes   = []
    edges   = []
    struct  = []
    gain    = []
    ar_gain = []
    types = [int, int, int, float, float]
    lists = [nodes, edges, struct, gain, ar_gain]

    # Collect Information
    idx = 0
    for l in lines:
        l = l.strip()

        next_tok = tok_seq[idx]

        if l[:min(len(next_tok), len(l))] == next_tok:
            lists[idx].append(types[idx](l.split(" ")[-1].split("\t")[-1]))
            idx = (idx + 1) % len(tok_seq)

    # Plot
    fns = [(lambda n, e, s: (e * 2) / n), \
           (lambda n, e, s: ((e * 2) / (n * (n - 1)))), \
           (lambda n, e, s: (s * 2) / n), \
           (lambda n, e, s: ((s * 2) / (n * (n - 1))))]
    terms = ["Avg. Degree", "Density", "Rem. Avg. Degree", "Rem. Density"]

    for idx in range(0, len(fns)):
        fn = fns[idx]
        plt.figure(figsize=(9,6))
        x_axis = [fn(nodes[i], edges[i], struct[i]) for i in range(0, len(nodes))]
        plt.scatter(x_axis, gain, label="algorithm")
        plt.scatter(x_axis, ar_gain, label="random")
        plt.legend()
        plt.title("Effectiveness of Algorithm vs. Graph Property", size=18)
        plt.xlabel(terms[idx], size=14)
        plt.ylabel("Gain above all Structure", size=14)
        plt.show()
