from graph_sequencer_utility import GraphSequence

def find_stars(edges, directed=True, temporal=True):

    if not directed:
        if temporal:
            edges = set([(min(a, b), max(a, b), t) for (a, b, t) in edges])
        else:
            edges = set([(min(a, b), max(a, b)) for (a, b) in edges])

    edge_colors = {}
    if not temporal:
        edge_colors = {(a, b): 0 for (a, b) in edges}
    else:
        timestamps_per_edge = {}
        for (a, b, t) in edges:
            if (a, b) not in timestamps_per_edge:
                timestamps_per_edge[(a, b)] = set()
            timestamps_per_edge[(a, b)].add(t)

        timestamps_to_color = {}
        next_color = 0
        for edge, times in timestamps_per_edge.items():
            times = tuple(sorted(list(times)))
            if times not in timestamps_to_color:
                timestamps_to_color[times] = next_color
                next_color += 1
            edge_colors[edge] = timestamps_to_color[times]

    # edge_colors just contains the directed information
    #
    # connection_types factors in back-edges for colors
    connection_types = {}
    if not directed:
        connection_types = edge_colors
    elif not temporal:
        for (a, b), _ in edge_colors.items():
            if (b, a) in edge_colors:
                # Bidirected
                connection_types[(a, b)] = 2
                connection_types[(b, a)] = 2
            else:
                # Out
                connection_types[(a, b)] = 0
                # In
                connection_types[(b, a)] = 1
    else:
        # Directed and temporal
        color_pairs_to_connection_type = {}
        next_connection_type = 0
        for (a, b), c in edge_colors.items():
            if (b, a) in edge_colors:
                c_alt = edge_colors[(b, a)]
            else:
                c_alt = -1
            if (c, c_alt) not in color_pairs_to_connection_type:
                color_pairs_to_connection_type[(c,c_alt)] = next_connection_type
                next_connection_type += 1
            if (c_alt, c) not in color_pairs_to_connection_type:
                color_pairs_to_connection_type[(c_alt,c)] = next_connection_type
                next_connection_type += 1
            connection_types[(a, b)] = color_pairs_to_connection_type[(c,c_alt)]
            connection_types[(b, a)] = color_pairs_to_connection_type[(c_alt,c)]

    # At this point, connection_types is "bidirected"
    for (a, b), _ in connection_types.items():
        assert (b, a) in connection_types

    neighbors_by_type = {}
    num_neighbors = {}
    for (a, b), type_value in connection_types.items():
        if a not in num_neighbors:
            num_neighbors[a] = 0
            neighbors_by_type[a] = {}
        if type_value not in neighbors_by_type[a]:
            neighbors_by_type[a][type_value] = set()

        num_neighbors[a] += 1
        neighbors_by_type[a][type_value].add(b)

    stars = []
    partial_stars = []
    pairs = []
    for node, n_b_t in neighbors_by_type.items():
        for type_value, neighbors in n_b_t.items():
            num_singleton_neighbors = 0
            for neighbor in neighbors:
                if num_neighbors[neighbor] == 1:
                    num_singleton_neighbors += 1
            if num_singleton_neighbors == num_neighbors[node]:
                if num_singleton_neighbors == 1:
                    if connection_types[(node, neighbor)] == \
                            connection_types[(neighbor, node)] and \
                                node < neighbor:
                        pairs.append((node, neighbor))
                else:
                    stars.append((num_singleton_neighbors, node))
            elif num_singleton_neighbors > 1:
                partial_stars.append((num_singleton_neighbors, node))

    stars.sort(reverse=True)
    partial_stars.sort(reverse=True)
    return (stars, partial_stars, pairs)

if __name__ == "__main__":
    # Day, Week resolution
    #
    # The -28561 makes the start time 12 AM on April 15th, 2004
    #   (That's Pacific Daylight Time)
    #   Choosing -10561 would be 5 AM instead.
    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60*60*24), \
        unit_name="days",
        units_per_window=7, \
        start_offset_number=-28561, \
        windows_overlap=False, \
        flatten_window=False, \
        weight_repeats=False, \
        directed=True)

    print("Weekly Analysis")
    i = 1
    largest_star = 0
    largest_week = None
    while window_GS.has_next():
        print("Week %d" % i)
        (nodes, edges) = window_GS.next()
        (stars, partial_stars, _) = \
            find_stars(edges, temporal=True, directed=True)
        if len(stars) > 0:
            print("  Largest Star: %d" % stars[0][0])
            if stars[0][0] > largest_star:
                largest_star = stars[0][0]
                largest_week = i
        if len(partial_stars) > 0:
            print("  Largest Partial Star: %d" % partial_stars[0][0])
            if partial_stars[0][0] > largest_star:
                largest_star = partial_stars[0][0]
                largest_week = i
        if len(stars) == 0 and len(partial_stars) == 0:
            print("  No Stars At All!")
        i += 1
    print("Largest star of %d at week %d" % (largest_star, largest_week))
        

    # Hour, Day resolution
    window_GS = GraphSequence()
    window_GS.set_window_sequence_with_temporal_file(\
        filename="datasets/college-temporal.g", \
        time_numbers_per_unit=(60*60), \
        unit_name="hours",
        units_per_window=24, \
        start_offset_number=-28561, \
        windows_overlap=False, \
        flatten_window=False, \
        weight_repeats=False, \
        directed=True)
    print("Daily Analysis")
    i = 1
    largest_star = 0
    largest_day = None
    while window_GS.has_next():
        print("Day %d" % i)
        (nodes, edges) = window_GS.next()
        (stars, partial_stars, _) = \
            find_stars(edges, temporal=True, directed=True)
        if len(stars) > 0:
            print("  Largest Star: %d" % stars[0][0])
            if stars[0][0] > largest_star:
                largest_star = stars[0][0]
                largest_day = i
        if len(partial_stars) > 0:
            print("  Largest Partial Star: %d" % partial_stars[0][0])
            if partial_stars[0][0] > largest_star:
                largest_star = partial_stars[0][0]
                largest_day = i
        if len(stars) == 0 and len(partial_stars) == 0:
            print("  No Stars At All!")

        i += 1
    print("Largest star of %d at day %d" % (largest_star, largest_day))
