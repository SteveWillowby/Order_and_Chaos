# SCHENO, the SCHEma NOise scoring function

Also includes a pattern-finding genetic algorithm guided by SCHENO as its fitness function.

## About

Given two graphs, SCHENO measures how well the first graph represents the underlying pattern(s) in the second graph.

This code also contains some useful C++ wrappers around the C implementation of Nauty and Traces. In addition to offering a more convenient interface, the wrappers offer the following features: Adding edge colorings to the graph before running Nauty/Traces, running Traces efficiently on directed graphs by creating undirected graphs equivalent to the directed ones behind the scenes, amortized constant-time edits to a graph stored in the Nauty/Traces input format (a non-trivial feature), and a slightly modified version of the Nauty/Traces code that can compile without error for C++ multithreaded use.

We have a separate respository containing only the relevant Nauty and Traces code here (TODO: Link).

The complete list of files and folders needed for the wrappers is:

```
nauty27r4
nauty27r4_modified
nauty_traces_setup.sh
nauty_traces.h
nauty_traces.cpp
graph.h
graph.cpp
sparse_graph.h
sparse_graph.cpp
nt_sparse_graph.h
nt_sparse_graph.cpp
nt_partition.h
nt_partition.cpp
edge.h
coloring.h
```

This code is associated with the research paper (TODO: Title) currently on arXiv at (TODO: link).

## Installation

### Main Binaries

The code for the main binaries is all contained in this repository, so no outside libraries or packages should be necessary.

To compile, run:

`./nauty_traces_setup.sh`

and then run:

`./compile.sh`

That's it! The programs should be ready to run.


### Other Experiments

Much of the experimental code that uses the main binary is coded in Python.

Necessary Python packages are: (TODO: Fill in)

You can install these with the commands: (TODO: Fill in) 


## Running the Main Binary

The main binary is located at `executables/main`.

To see the command-line options, run `executables/main -h`


## Running the Experiments

In the example below, we get the modified version of a college football season team-plays-team graph. Then we extract the modified graph from the result file and re-run the algorithm on the modified graph, this time more slowly (using the heuristic).

```
nice -2 executables/main -graph experiments/real_world_graphs/season_4_undirected_edges.txt -nodes experiments/real_world_graphs/season_4_undirected_nodes.txt -n_itr=200 -o experiments/test_results/season_4_undir >> experiments/test_results/season_4_undir.txt

cd experiments

python3 make_node_list.py test_results/season_4_undir_graph.txt

cd ..

nice -2 executables/main -graph experiments/test_results/season_4_undir_graph.txt -nodes experiments/test_results/season_4_undir_graph_nodes.txt -n_itr=50 -use_heuristic -o experiments/test_results/season_4_core
```


In the example below, we take a 120-node Johnson graph, randomly modify 1% of its connections (measured in terms of the number of its edges), and then see if the algorithm can find the original graph.

```
time nice -2 executables/main -graph experiments/structure_recovery/johnson_10_3_120_edges.txt -nodes experiments/structure_recovery/johnson_10_3_120_nodes.txt -noise- 0.005 -noise+ 0.005 -n_itr 140 -o experiments/test_results/johnson_120
```

(TODO: Fill in)
