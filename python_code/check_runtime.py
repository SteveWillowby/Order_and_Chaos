from ram_friendly_NT_session import RAMFriendlyNTSession 
from graph_loader import read_graph

if __name__ == "__main__":

    edges = read_graph("../real_world_graphs/cora.g", directed=True)[0][3]

    num_nodes = len(edges)
    num_edges = sum([len(s) for s in edges])
    print("Num Nodes: %d    Num Edges: %d" % (num_nodes, num_edges))

    nt_session = RAMFriendlyNTSession(directed=True,
                                      has_edge_types=False,
                                      neighbors_collections=edges,
                                      dreadnaut_call="../nauty27r4/dreadnaut",
                                      kill_py_graph=True,
                                      only_one_call=True)

    orbits = nt_session.get_automorphism_orbits()
    num_aut = nt_session.get_num_automorphisms()

    nt_session.run()

    orbits = orbits.get()
    num_aut = num_aut.get()

    print("Num Aut: %s" % num_aut)
    print("Num Orbits: %d" % len(orbits))
