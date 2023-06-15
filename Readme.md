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

(TODO: Fill in)
