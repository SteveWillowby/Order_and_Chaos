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

    TYPE_NAME_STRINGS = {\
        STATIC_UNDIRECTED: "Static Undirected", \
        STATIC_DIRECTED: "Static Directed", \
        BIPARTITE_UNDIRECTED: "Bipartite Undirected", \
        BIPARTITE_DIRECTED: "Bipartite Directed", \
        STATIC_COLORED_UNDIRECTED: "Static Colored Undirected", \
        STATIC_COLORED_DIRECTED: "Static Colored Directed", \
        NODE_JOINING_UNDIRECTED: "Node-Joining Undirected", \
        NODE_JOINING_DIRECTED: "Node-Joining Directed", \
        NODE_JOINING_DIRECTED_ONE_WAY: "Node-Joining Directed One-Way", \
        TEMPORAL_UNDIRECTED: "Temporal Undirected", \
        TEMPORAL_DIRECTED: "Temporal Directed", \
        TEMPORAL_COLORED_UNDIRECTED: "Temporal Colored Undirected", \
        TEMPORAL_COLORED_DIRECTED: "Temporal Colored Directed"  \
        }

    IS_DIRECTED = {\
        STATIC_UNDIRECTED: False, \
        STATIC_DIRECTED: True, \
        BIPARTITE_UNDIRECTED: False, \
        BIPARTITE_DIRECTED: True, \
        STATIC_COLORED_UNDIRECTED: False, \
        STATIC_COLORED_DIRECTED: True, \
        NODE_JOINING_UNDIRECTED: False, \
        NODE_JOINING_DIRECTED: True, \
        NODE_JOINING_DIRECTED_ONE_WAY: True, \
        TEMPORAL_UNDIRECTED: False, \
        TEMPORAL_DIRECTED: True, \
        TEMPORAL_COLORED_UNDIRECTED: False, \
        TEMPORAL_COLORED_DIRECTED: True \
        }

    IS_TEMPORAL = {\
        STATIC_UNDIRECTED: False, \
        STATIC_DIRECTED: False, \
        BIPARTITE_UNDIRECTED: False, \
        BIPARTITE_DIRECTED: False, \
        STATIC_COLORED_UNDIRECTED: False, \
        STATIC_COLORED_DIRECTED: False, \
        NODE_JOINING_UNDIRECTED: True, \
        NODE_JOINING_DIRECTED: True, \
        NODE_JOINING_DIRECTED_ONE_WAY: True, \
        TEMPORAL_UNDIRECTED: True, \
        TEMPORAL_DIRECTED: True, \
        TEMPORAL_COLORED_UNDIRECTED: True, \
        TEMPORAL_COLORED_DIRECTED: True \
        }

