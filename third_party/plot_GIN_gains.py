import matplotlib.pyplot as plt

if __name__ == "__main__":

    alg_file_name = "GIN"
    alg_name      = "GIN"
    graph_file_names = ["karate", "season_4", "foodweb", "pol_blogs", "eucore", "cora"]
    graph_names      = ["Karate", "Football", "Foodweb", "PolBlogs",  "EUCore", "Cora"]

    NUM_FIGS = 6
    FIGS_AVAILABLE = 6  # Number of completely processed data files available

    MAX_NUM_POINTS  = 400
    MIN_X_INCREMENT = 2.0 / MAX_NUM_POINTS

    font = "Georgia"
    plt.rcParams["font.family"] = font
    
    # fig, axs = plt.subplots(NUM_FIGS, sharex=True, figsize=(6, 12))
    fig = plt.figure()
    fig.set_size_inches((6, 12))
    gs = fig.add_gridspec(NUM_FIGS, hspace=0.35)
    axs = gs.subplots(sharex=True)

    fig.suptitle('%s\'s SCHENO Score Gains over "All is Structure" Hypothesis' % alg_name, y=0.94, x=0.46)
    # fig.xlabel("Schema Size / Original Edge Count")
    for i in range(0, NUM_FIGS):
        graph_idx = i % FIGS_AVAILABLE

        f = open("results/pan_threshold_%s_scores_%s.txt" % \
                    (alg_file_name, graph_file_names[graph_idx]))
        lines = f.readlines()
        f.close()
        values = [tuple(l.strip().split(" ")) for l in lines]

        schema_to_edges_ratios = [float(a) for (a, _, __) in values]
        alg_gains              = [float(b) for (_, b, __) in values]
        rand_gains             = [float(c) for (_, __, c) in values]

        indices = [0]
        for j in range(1, len(schema_to_edges_ratios)):
            if schema_to_edges_ratios[j] - schema_to_edges_ratios[indices[-1]] >= MIN_X_INCREMENT:
                indices.append(j)

        schema_to_edges_ratios = [schema_to_edges_ratios[idx] for idx in indices]
        alg_gains              = [alg_gains[idx]              for idx in indices]
        rand_gains             = [rand_gains[idx]             for idx in indices]

        axs[i].plot([-1, 3], [0, 0], color="gray", linewidth=0.5)
        axs[i].plot(schema_to_edges_ratios, alg_gains,  label=alg_name, \
                    color="black")
        axs[i].plot(schema_to_edges_ratios, rand_gains, label="Random", \
                    color="black", linestyle=":")
        axs[i].legend()
        axs[i].set_xlim([-0.1, 2.1])

        axs[i].set_title(graph_names[graph_idx])
        # axs[i].xaxis.grid()
        # axs[i].yaxis.grid()

    for ax in axs.flat:
        ax.set(xlabel="Schema Size / Original Edge Count", ylabel="log2(Score Gain)")
    for ax in fig.get_axes():
        ax.label_outer()

    plt.savefig("results/%s_all_graphs.svg" % alg_name, bbox_inches='tight', pad_inches=0.1)
