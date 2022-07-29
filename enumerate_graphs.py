# Iterates through all distinct adjacency matrices on n nodes.
#   Returns them in the "neighbors collections" format.
#
# No self-loops.
#
# O(2^(n choose 2) ) if undirected
# O(2^(n * (n - 1)) ) if directed
class AdjEnumerator:

    def __init__(self, n, directed):
        assert n >= 1
        self.__n__ = n
        self.__directed__ = directed

        self.__counters__ = [0 for _ in range(0, n)]
        if self.__directed__:
            self.__maxs__ = [self.__n__ - 1 for _ in range(0, self.__n__)]
        else:
            self.__maxs__ = [self.__n__ - (i + 1) for i in range(0, self.__n__)]
        self.__maxs__ = [int(2**v) for v in self.__maxs__]

        self.__masks__ = [int(2**i) for i in range(0, self.__n__)]

        self.__done__ = False
        self.__graph_count__ = 0

    # Returns the number of graphs returned thus far.
    def graph_count(self):
        return self.__graph_count__

    # Returns a new graph, or None if there are no more graphs to find.
    def next_graph(self):
        if self.__done__:
            return None

        nc = []
        if self.__directed__:
            for i in range(0, self.__n__):
                s = set()
                c = self.__counters__[i]
                for j in range(0, self.__n__):
                    if j == i:
                        continue
                    j_idx = j - int(j > i)

                    if c & self.__masks__[j_idx]:
                        s.add(j)
                nc.append(s)
        else:
            for i in range(0, self.__n__):
                s = set()
                c = self.__counters__[i]
                for j_idx in range(0, self.__n__ - (i + 1)):
                    if c & self.__masks__[j_idx]:
                        j = i + j_idx + 1
                        s.add(j)
                nc.append(s)

        self.__done__ = True
        for i in range(0, self.__n__):
            if self.__counters__[i] < self.__maxs__[i] - 1:
                self.__counters__[i] += 1
                self.__done__ = False
                break
            else:
                self.__counters__[i] = 0

        self.__graph_count__ += 1
        return nc

# Iterates through all isomorphically distinct graphs on n nodes. 
#
# No self-loops.
#
# O(2^(n choose 2) * automorphism stuff) if undirected -- or better
# O(2^(n * (n - 1)) * automorphism stuff) if directed -- or better
class GraphEnumerator:

    def __init__(self, n, directed, auto_solver_class, hashing=False):
        assert n >= 1
        self.__n__ = n
        self.__directed__ = directed

        self.__curr_G__ = None
        self.__curr_orbit_pair__ = None
        self.__curr__ = []
        self.__next__ = {}
        self.__graph_count__ = 0
        self.__auto_solver_class__ = auto_solver_class
        self.__hashing__ = hashing
        self.__first__ = True

        empty_graph = [set() for _ in range(0, n)]
        (canon, node_to_orbit, orbits, orbitwise_orbits) = \
            self.__canonical_form_and_orbits__(empty_graph)
        self.__next__[canon] = \
            (empty_graph, node_to_orbit, orbits, orbitwise_orbits)

    # Returns the number of graphs returned thus far.
    def graph_count(self):
        return self.__graph_count__

    # Returns a new graph, or None if there are no more graphs to find.
    def next_graph(self):
        if self.__first__:
            self.__first__ = False
            self.__graph_count__ += 1
            return [set() for _ in range(0, self.__n__)]

        while True:
            if self.__curr_G__ is None:
                if len(self.__curr__) == 0:
                    if len(self.__next__) == 0:
                        return None
                    self.__curr__ = [(nc, n_to_o, o, o_o) for \
                                        ((key, (nc, n_to_o, o, o_o))) \
                                            in self.__next__.items()]
                    self.__next__ = {}
                self.__curr_G__ = self.__curr__.pop()
                self.__curr_orbit_pair__ = (0, 0)

            (nc, node_to_orbit, orbits, orbitwise_orbits) = self.__curr_G__
            num_orbits = len(orbits)

            (a, b) = self.__curr_orbit_pair__
            if a >= len(orbits):
                self.__curr_G__ = None
                continue

            a_orbit = orbits[a]
            a_orbit_size = len(a_orbit)
            if a_orbit_size == 1:
                if b >= len(orbits):
                    self.__curr_orbit_pair__ = (a + 1, 0)
                    continue
                elif b == a:
                    self.__curr_orbit_pair__ = (a, b + 1)
                    continue
                b_orbit = orbits[b]
            else:
                if b >= len(orbitwise_orbits[a]):
                    self.__curr_orbit_pair__ = (a + 1, 0)
                    continue
                b_orbit = orbitwise_orbits[a][b]

            node_a = a_orbit[0]
            node_b = b_orbit[0]

            self.__curr_orbit_pair__ = (a, b + 1)

            if node_b == node_a:
                assert len(b_orbit) == 1
                continue

            # TODO: Consider making this more efficient by somehow pre-
            #   determining which nodes are not already neighbors.
            if node_b in nc[node_a]:
                continue

            if not self.__directed__:
                if min(a_orbit) > min(orbits[node_to_orbit[node_b]]):
                    continue

            new_graph = [set(s) for s in nc]
            new_graph[node_a].add(node_b)
            if not self.__directed__:
                new_graph[node_b].add(node_a)

            (canon, n_to_o, o, o_o) = \
                self.__canonical_form_and_orbits__(new_graph)
            if canon not in self.__next__:
                self.__next__[canon] = (new_graph, n_to_o, o, o_o)
                self.__graph_count__ += 1
                return new_graph

    def __canonical_form_and_orbits__(self, nc):
        session = self.__auto_solver_class__(directed=self.__directed__, \
                                             has_edge_types=False, \
                                             neighbors_collections=nc)
        canon_order = session.get_canonical_order()
        orbits = session.get_automorphism_orbits()
        session.run()
        canon_order = canon_order.get()
        orbits = orbits.get()
        num_orbits = len(orbits)

        if num_orbits == self.__n__:
            orbitwise_orbits = None
        else:
            orbitwise_orbits = []
            for orbit_idx in range(0, num_orbits):
                orbit = orbits[orbit_idx]
                if len(orbit) == 1:
                    orbitwise_orbits.append(None)
                    continue

                a = orbit[0]
                orbit_minus_a = []
                for b in orbit:
                    if b != a:
                        orbit_minus_a.append(b)
                session.set_colors_by_partitions(\
                    [[a], orbit_minus_a] + \
                    [orbits[i] for i in range(0, orbit_idx)] + \
                    [orbits[i] for i in range(orbit_idx + 1, num_orbits)])

                sub_orbits = session.get_automorphism_orbits()
                session.run()
                sub_orbits = sub_orbits.get()
                orbitwise_orbits.append(sub_orbits)

        session.end_session()

        node_to_orbit = [None for _ in range(0, self.__n__)]
        for orbit_idx in range(0, len(orbits)):
            for node in orbits[orbit_idx]:
                node_to_orbit[node] = orbit_idx

        order_to_node = canon_order
        node_to_order = [(canon_order[i], i) for i in range(0, self.__n__)]
        node_to_order.sort()
        node_to_order = [i for (node, i) in node_to_order]

        if self.__hashing__ or self.__directed__:
            canon_edgelist = \
                tuple([tuple(sorted([node_to_order[b] for b in nc[a]])) \
                                        for a in order_to_node])
        else:
            canon_edgelist = []
            for i in range(0, self.__n__):
                greater_neighbors = []
                for n in nc[order_to_node[i]]:
                    j = node_to_order[n]
                    if j >= i:
                        greater_neighbors.append(j)
                canon_edgelist.append(tuple(sorted(greater_neighbors)))
            canon_edgelist = tuple(canon_edgelist)

        if self.__hashing__:
            canon_form = hash(canon_edgelist)
        else:
            canon_form = canon_edgelist

        return (canon_form, node_to_orbit, orbits, orbitwise_orbits)


if __name__ == "__main__":
    from py_NT_session import PyNTSession

    correct_values = {False: {1: 1, 2: 2, 3: 4, 4: 11, 5: 34, 6: 156, 7: 1044, 8: 12346, 9: 274668}, \
                       True: {1: 1, 2: 3, 3: 16, 4: 218, 5: 9608, 6: 1540944}}

    hashing = False
    print("Testing Graph Enumerator with%s Hashing" % {True: "", False: "out"}[hashing])
    for directed in [False, True]:
        print("Directed = %s" % directed)
        for n in range(1, {False: 9, True: 6}[directed]):
            GE = GraphEnumerator(n=n, directed=directed, auto_solver_class=PyNTSession, \
                                 hashing=hashing)
            while GE.next_graph() is not None:
                pass
            print(GE.graph_count())
            assert GE.graph_count() == correct_values[directed][n]

    correct_values = {False: {1: 1, 2: 2, 3: 8, 4: 64, 5: 1024, 6: 32768, 7: 2097152, 8: 268435456}, \
                       True: {1: 1, 2: 4, 3: 64, 4: 4096, 5: 1048576, 6: 1073741824}}

    print("Testing Adjacency Matrix Enumerator")
    for directed in [False, True]:
        print("Directed = %s" % directed)
        for n in range(1, {False: 8, True: 6}[directed]):
            AE = AdjEnumerator(n=n, directed=directed)
            while AE.next_graph() is not None:
                pass
            print(AE.graph_count())
            assert AE.graph_count() == correct_values[directed][n]
