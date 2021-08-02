from graph_sequencer_utility import GraphSequence

def find_stars(temporal_directed_edges):

    in_neighbors_by_timestamp = {}
    out_neighbors_by_timestamp = {}
    for (a, b, t) in temporal_directed_edges:
        if t not in in_neighbors_by_timestamp:
            in_neighbors_by_timestamp[t] = {}
            out_neighbors_by_timestamp[t] = {}
        i_n_t = in_neighbors_by_timestamp[t]
        o_n_t = out_neighbors_by_timestamp[t]
        if a not in o_n_t:
            o_n_t[a] = set()
            i_n_t[a] = set()
        if b not in o_n_t:
            o_n_t[b] = set()
            i_n_t[b] = set()
        o_n_t[a].add(b)
        i_n_t[b].add(a)

    stars = []
    partial_stars = []
    bidir_pairs = []
    for t, i_n_t in in_neighbors_by_timestamp.items():
        o_n_t = out_neighbors_by_timestamp[t]

        # Inward star (e.g. followers)
        for node, in_neighbors in i_n_t.items():
            if len(in_neighbors) == 0:
                continue
            num_single_neighbors = 0
            for neighbor in in_neighbors:
                if len(o_n_t[neighbor]) + len(i_n_t[neighbor]) == 1:
                    num_single_neighbors += 1
            if num_single_neighbors == len(in_neighbors):
                stars.append((num_single_neighbors, node))
            elif num_single_neighbors >= 2:
                partial_stars.append((num_single_neighbors, node))

        # Outward star (e.g. text recipients)
        for node, out_neighbors in o_n_t.items():
            if len(out_neighbors) == 0:
                continue
            num_single_neighbors = 0
            for neighbor in out_neighbors:
                if len(o_n_t[neighbor]) + len(i_n_t[neighbor]) == 1:
                    num_single_neighbors += 1
            if num_single_neighbors == len(in_neighbors):
                stars.append((num_single_neighbors, node))
            elif num_single_neighbors >= 2:
                partial_stars.append((num_single_neighbors, node))

        # Bidirectional star
        for node, out_neighbors in o_n_t.items():
            if len(out_neighbors & i_n_t[node]) == 0:
                continue
            num_single_neighbors = 0
            for neighbor in out_neighbors & i_n_t[node]:
                if len(o_n_t[neighbor]) + len(i_n_t[neighbor]) == 2:
                    num_single_neighbors += 1
            if num_single_neighbors == len(out_neighbors & i_n_t[node]):
                if num_single_neighbors == 1 and len(out_neighbors) == 1:
                    bidir_pairs.append((2, node))
                else:
                    stars.append((num_single_neighbors, node))
            elif num_single_neighbors >= 2:
                partial_stars.append((num_single_neighbors, node))

    stars.sort(reverse=True)
    partial_stars.sort(reverse=True)
    return (stars, partial_stars, bidir_pairs)

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
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)

    print("Weekly Analysis")
    i = 1
    largest_star = 0
    largest_week = None
    while window_GS.has_next():
        print("Week %d" % i)
        (nodes, edges) = window_GS.next()
        (stars, partial_stars, _) = find_stars(temporal_directed_edges=edges)
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
        flatten_window=True, \
        weight_repeats=False, \
        directed=True)
    print("Daily Analysis")
    i = 1
    largest_star = 0
    largest_day = None
    while window_GS.has_next():
        print("Day %d" % i)
        (nodes, edges) = window_GS.next()
        (stars, partial_stars, _) = find_stars(temporal_directed_edges=edges)
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
