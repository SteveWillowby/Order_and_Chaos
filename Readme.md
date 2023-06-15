# Order and Chaos

## About

This code is associated with the research paper (TODO: Title) currently on arXiv at (TODO: link).

It allows for taking a directed or undirected graph and decomposing it into order and chaos, structure and noise.

This code also contains some useful C++ wrappers around the C implementation of Nauty and Traces. In addition to offering a more convenient interface, the wrappers offer the following features: Adding edge colorings to the graph before running Nauty/Traces, running Traces efficiently on directed graphs by creating undirected graphs equivalent to the directed ones behind the scenes, constant-time edits to a graph stored in the Nauty/Traces input format (a non-trivial feature), and a slightly modified version of the Nauty/Traces code that can compile without error for C++ multithreaded use.

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



## Installation

### Main Binary

The code for the main binary is all self-contained in this repository, so no outside libraries or packages should be necessary.

To compile, run:

`./nauty_traces_setup.sh`

and then run:

`./compile.sh`

That's it! The main program should be ready to run.


### Other Experiments

Much of the experimental code that uses the main binary is coded in Python.

Necessary Python packages are: (TODO: Fill in)

You can install these with the commands: (TODO: Fill in) 


## Running the Main Binary

The main binary is located at `executables/main`.

To see the command-line options, run `executables/main -h`


## Running the Experiments

In the below example, we get the modified version of a college football season team-plays-team graph.

```
nice -2 executables/main -graph experiments/real_world_graphs/season_4_undirected_edges.txt -nodes experiments/real_world_graphs/season_4_undirected_nodes.txt -n_itr=200 >> experiments/test_results/season_4_undir.txt

cd experiments/plotting_code

python3 build_result_graph.py ../real_world_graphs/season_4_undirected_nodes.txt ../real_world_graphs/season_4_undirected_edges.txt ../test_results/season_4_undir.txt

cat structure_graph.csv
```


In the below example, we take a 120-node Johnson graph, randomly modify 1% of its connections (measured in terms of the number of its edges), and then see if the algorithm can find the original graph.

```
time nice -2 executables/main -graph experiments/simple_test_graphs/johnson_10_3_120_edges.txt -nodes experiments/simple_test_graphs/johnson_10_3_120_nodes.txt -noise- 0.005 -noise+ 0.005 -n_itr 140
```

(TODO: Fill in)
