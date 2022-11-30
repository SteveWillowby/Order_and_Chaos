from ram_friendly_NT_session import RAMFriendlyNTSession 
from graph_loader import read_graph
import sys

if __name__ == "__main__":
    if (len(sys.argv) == 1):
        print("Requires a graph file name (e.g. cora.g)")
        exit(1)

    directed = True
    if (len(sys.argv) >= 3):
        if (sys.argv[2] == "false" or sys.argv[2] == "False"):
            directed = False

    graph = read_graph("../real_world_graphs/" + sys.argv[1], directed=directed)
    edges = graph[0][3]
    node_labels = graph[1]
    num_self_loops = sum(node_labels)

    num_nodes = len(edges)
    num_edges = sum([len(s) for s in edges]) + num_self_loops
    print("ACTUAL NUM NODES AND EDGES:")
    print("Num Nodes: %d    Num Edges: %d" % (num_nodes, num_edges))

    nt_session = RAMFriendlyNTSession(mode="Traces",
                                      directed=directed,
                                      has_edge_types=False,
                                      neighbors_collections=edges,
                                      dreadnaut_call="../nauty27r4/dreadnaut",
                                      kill_py_graph=True,
                                      only_one_call=True,
                                      print_notes=True)

    # Include self-loops
    nt_session.set_colors_by_coloring(node_labels)

    orbits = nt_session.get_automorphism_orbits()
    num_aut = nt_session.get_num_automorphisms()
    canon_order = nt_session.get_canonical_order()

    nt_session.run()

    orbits = orbits.get()
    num_aut = num_aut.get()
    canon_order = canon_order.get()

    print("Num Aut: %s" % num_aut)
    print("Num Orbits: %d" % len(orbits))
