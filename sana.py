#!/usr/bin/env python3
import os
import re

# get_sana_alignment()
#
# G1 -- a graph which supports .edges() and .nodes() -- all nodes must have at
#   least one edge
# G2 -- another graph which supports .edges() and .nodes() -- all nodes must
#   have at least one edge
#
# t -- time to run, in minutes -- can be floating point value
#
# G1_name -- allows the program to avoid repeat copying G1 if G1 has already
#   been written to a file in sana's workspace
#
# G2_name -- allows the program to avoid repeat copying G2 if G2 has already
#   been written to a file in sana's workspace
#
# debug -- print info about processing
#
# nodetype -- a type to cast the output file's strings as
#
# multiprocessing -- if true, adds PID to filenames so that 
#
# RETURNS a dict mapping from G1's nodes to G2's nodes.
def get_sana_alignment(G1, G2, G1_name="", G2_name="", t=10.0, debug=True, \
        nodetype=int, multiprocessing=False):
    assert t > 0

    pid = os.getpid()

    if G1_name == "":
        G1_name = "unnamed_input_1"
    if G2_name == "":
        G2_name = "unnamed_input_2"

    if multiprocessing:
        G1_name += "_%d" % pid
        G2_name += "_%d" % pid

    sana_input_dir = "networks/custom_inputs"
    writing_input_dir = "SANA/" + sana_input_dir
    if not os.path.isdir(writing_input_dir):
        os.mkdir(writing_input_dir)

    sana_filenames = []
    for (G, G_name) in [(G1, G1_name), (G2, G2_name)]:
        filename = "%s/%s.el" % (writing_input_dir, G_name)
        sana_filename = "%s/%s.el" % (sana_input_dir, G_name)

        # Protect against malicious "filenames."
        if not bool(re.match('^[a-zA-Z0-9\-\/_\.]+$', filename)):
            raise Exception("Error! Filename must consist of alphanumeric " + \
                "characters and '/', '-', '_', '.' ." + \
                " Filename '%s' does not meet this requirement." % filename)

        sana_filenames.append(sana_filename)
        if not os.path.isfile(filename):
            if debug:
                print("Writing graph '%s' to %s" % (G_name, filename))
            write_sana_edge_list(G, filename)
        elif debug:
            print("File %s already exists, so skipping writing." % filename)

    if debug:
        print("Running SANA for %d minutes..." % t)

    output_name = "sana_output"
    if multiprocessing:
        output_name += "_%d" % pid

    if debug:
        command = "cd SANA; ./sana -fg1 %s -fg2 %s -ec 1 -t %d -o %s" % \
                  (sana_filenames[0], sana_filenames[1], t, output_name)
    else:
        command = "cd SANA; ./sana -fg1 %s -fg2 %s -ec 1 -t %d -o %s >/dev/null" % \
                  (sana_filenames[0], sana_filenames[1], t, output_name)

    if debug:
        print(command)

    os.system(command)

    alignment = {}
    with open("SANA/%s.align" % output_name, "r") as f:
        lines = f.readlines()
        lines = [line.split() for line in lines]
        for line in lines:
            node_a = nodetype(line[0])
            node_b = nodetype(line[1])
            alignment[node_a] = node_b

    return alignment

def write_sana_edge_list(G, filename):
    edges = list(G.edges())
    nodes = list(G.nodes())

    observed_nodes_set = set()
    for (a, b) in edges:
        observed_nodes_set.add(a)
        observed_nodes_set.add(b)

    if len(observed_nodes_set) < len(nodes):
        raise Exception("Error! Cannot align a graph in which a node has " + \
            "zero edges (%s)" % filename)

    with open(filename, "w") as f:
        for i in range(0, len(edges)):
            (a, b) = edges[i]
            f.write("%s %s" % (a, b))
            if i < len(edges) - 1:
                f.write("\n")

# Deletes all files created for sana purposes.
def clear_sana_input_files():
    for root, dirs, files in os.walk("SANA/networks/custom_inputs/"):
        for f in files:
            os.remove(os.path.join(root, f))

def clear_sana_output_files():
    for root, dirs, files in os.walk("SANA/"):
        for f in files:
            if f.startswith("sana_output"):
                os.remove(os.path.join(root, f))

if __name__ == "__main__":
    import networkx as nx
    G1 = nx.gnm_random_graph(100, 1000)
    G2 = nx.gnm_random_graph(100, 1000)

    alignment = get_sana_alignment(G1, G2, G1_name="", G2_name="", t=1.5, \
        debug=False, nodetype=int, multiprocessing=True)
    print(alignment)

    clear_sana_input_files()
    clear_sana_output_files()
