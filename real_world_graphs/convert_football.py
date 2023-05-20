def date_to_season(d):
    (year, month, day) = d
    if month < 8:
        return year - 2005
    return year - 2006

if __name__ == "__main__":

    f = open("cfb_team.csv", "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip().split(",") for l in lines if len(l.strip()) > 0]
    team_name_ids = list(set([(l[2], l[0]) for l in lines]))
    name_to_id = {a: b for (a, b) in team_name_ids}
    id_to_name = {b: a for (a, b) in team_name_ids}
    ids = set([b for (a, b) in team_name_ids])

    f = open("id_to_team.txt", "w")
    for (i, n) in id_to_name.items():
        f.write("%s %s\n" % (i, n))
    f.close()

    f = open("cfb_game.csv", "r")
    lines = f.readlines()
    f.close()

    lines = [l.strip().split(",") for l in lines if len(l.strip()) > 0]
    date_ta_tb_result = [(tuple([int(s) for s in l[1].split("-")]), \
                          l[2], l[3], l[5]) for l in lines]
    season_ta_tb_result = [(date_to_season(d), a, b, c) for (d, a, b, c) in \
                            date_ta_tb_result]
    # start_years = sorted(list(set([x[0][0] for x in date_ta_tb_result])))[:-1]

    seasons = sorted(list(set([x[0] for x in season_ta_tb_result])))
    season_lists = [[] for _ in seasons]
    for (s, a, b, c) in season_ta_tb_result:
        season_lists[s].append((a, b, c))

    for i in range(0, len(season_lists)):
        l = season_lists[i]
        edges = []
        for (ta, tb, result) in l:
            edges.append((ta, tb))
            if (result == "NEUTRAL"):
                edges.append((tb, ta))
        edges = sorted(list(set(edges)))
        nodes = sorted(list(set([a for (a, b) in edges] + \
                                [b for (a, b) in edges])))

        f = open("season_%d_directed_edges.txt" % i, "w")
        for (a, b) in edges:
            f.write("%s %s\n" % (a, b))
        f.close()
        f = open("season_%d_directed_nodes.txt" % i, "w")
        for x in nodes:
            f.write("%s\n" % x)
        f.close()

        undir_edges = sorted(list(set([(min(a, b), max(a, b)) for (a, b) in edges])))
        f = open("season_%d_undirected_edges.txt" % i, "w")
        for (a, b) in undir_edges:
            f.write("%s %s\n" % (a, b))
        f.close()
        f = open("season_%d_undirected_nodes.txt" % i, "w")
        for x in nodes:
            f.write("%s\n" % x)
        f.close()
