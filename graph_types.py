# Simply a namespace of graph classifications

class GraphTypes:
    STATIC_UNDIRECTED = 0
    STATIC_DIRECTED = 1

    BIPARTITE_UNDIRECTED = 2
    BIPARTITE_DIRECTED = 3

    STATIC_COLORED_UNDIRECTED = 4
    STATIC_COLORED_DIRECTED = 5

    # NODE_JOINING graphs are temporal graphs:
    #   New nodes join, and can only connect to nodes from
    #   the same timestep or previous timesteps.
    #
    #   Nodes only form connections when they appear.
    NODE_JOINING_UNDIRECTED = 6
    NODE_JOINING_DIRECTED = 7
    #   Nodes add only out-edges when they join.
    NODE_JOINING_DIRECTED_ONE_WAY = 8

    TEMPORAL_UNDIRECTED = 9
    TEMPORAL_DIRECTED = 10

    TEMPORAL_COLORED_UNDIRECTED = 11
    TEMPORAL_COLORED_DIRECTED = 12
