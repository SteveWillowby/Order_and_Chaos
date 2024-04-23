# `SCHENO`, the SCHEma-NOise scoring function

### ...Also includes a pattern-finding genetic algorithm guided by SCHENO as its fitness function.

## About

The core idea behind `SCHENO` is that real-world graphs are messy manifestations of underlying patterns ("schemas"). `SCHENO` offers a principled way to measure how well you've done at uncovering those patterns.

Given a graph and a "noise set" (a list of node-pairs), `SCHENO` measures how well the graph "minus" the noise set represents the underlying patterns (i.e. schema) of the graph. We put "minus" in scarequotes because if an edge is in the noise set, `SCHENO` considers the graph without that edge, but if a non-edge is in the noise set, `SCHENO` considers the graph with that edge added.

`SCHENO` factors three things into a single numeric score:

- How patterned or structural is the schema?
- How random or noisy is the noise?
- How different is the schema from the original graph?

Much more information is available in the paper [here](http://arxiv.org/abs/2404.13489).

### Isomorphism and Automorphism Calculations

If you are simply interested in this project's handy `C++` wrappers surrounding the `nauty` and `traces` isomorphism code, those wrappers are available in a standalone repository [here](https://github.com/schemanoise/Nauty_and_Traces).

## Setup

The code for the main binaries is all contained in this repository, so no outside libraries or packages should be necessary.

Run `setup.sh` to compile the utilities. Then run `compile.sh` to compile the executables.

The binaries will be located in the `executables` folder.

## Using the Two Programs

You can find a PDF manual describing the use of this repository's code inside the `documentation` folder.

For a full theoretical description of `SCHENO`, see the paper located [here](http://arxiv.org/abs/2404.13489).

### `SCHENO_score`

This program takes a graph (expressed as an edgelist) and a noise set (also expressed as an edgelist). It reports how well the graph without the noise represents the underlying structure of the graph.

The graph can be directed or undirected.

To see the full help menu, run `executables/SCHENO_score -h`


### `SCHENO_ga`

This program takes a graph (expressed as an edgelist) and tries to find the best schema-noise decomposition for that graph using a simple genetic algorithm. The genetic algorithm uses the SCHENO score as its fitness function.

The graph can be directed or undirected.

`SCHENO_ga` will output three files:

- The `graph` file, which represents the pattern (aka "schema") -- given as an edgelist
- The `noise` file, which represents the noise -- given as an edgelist
- the `nodes` file (this one is unimportant -- just the list of the graph's nodes)

To see the full help menu, run `executables/SCHENO_ga -h`

Example:

```
~/Documents/SCHENO$ cat example_graphs/almost_7_cycle.txt
0 1
1 2
2 3
3 4
4 5
5 6
~/Documents/SCHENO$ executables/SCHENO_ga \
                        -graph example_graphs/almost_7_cycle.txt \
                        -o example_graphs/almost_7_cycle \
                        -n_itr 4 -topk 3
Loading graph from files:
    
    example_graphs/almost_7_cycle.txt
  ...graph loaded. It has 7 nodes and 6 edges, 0 of which are self-loops.

The original graph has log2_aut = 1
Running for 4 iterations per trial...
Trial 0 for (noise-, noise+) = (0, 0)

log2_p_plus:          -3.95424
log2_1_minus_p_plus:  -0.0962126
log2_p_minus:         -3.95424
log2_1_minus_p_minus: -0.0962126
p_plus:  0.0645144
p_minus: 0.0645144

  Beginning edge heuristic pre-computation...
  ...Finished edge heuristic pre-computation.
The seed edge set gets a score of 1

Beginning Iteration 1...
	Mutating
	Mating
	Scoring
...Finished iteration with a best score of 2.75669

Beginning Iteration 2...
	Mutating
	Mating
	Scoring
...Finished iteration with a best score of 2.75669

Beginning Iteration 3...
	Mutating
	Mating
	Scoring
...Finished iteration with a best score of 2.75669

Beginning Iteration 4...
	Mutating
	Mating
	Scoring
...Finished iteration with a best score of 2.75669

With a score of 2.75669 we have log2(|Aut(G_H)|) of 3.80735
With edges: 
(0, 6), 

With a score of 1 we have log2(|Aut(G_H)|) of 1
With edges: 


With a score of 0.927809 we have log2(|Aut(G_H)|) of 4.32193
With edges: 
(1, 2), (2, 6),
~/Documents/SCHENO$ cat example_graphs/almost_7_cycle_graph.txt
0 6
0 1
1 2
2 3
3 4
4 5
5 6
```


## Including SCHENO Scores in Your Own Program

You can find a PDF manual describing the interface for this repository's code inside the `documentation` folder.

If you want to call the code directly, you can include all the classes through the file `scheno.h` contained in the `scheno` subfolder.

To compile, include the `scheno/scheno.a` static library. For example:

```
g++ my_program.cpp scheno/scheno.a
```

The file `SCHENO_score.cpp` illustrates the use of most of the classes.

For scoring a bunch of candidate noise sets in parallel, consider using the `ThreadPoolScorer` class.

## License

This repository contains code from several sources, some of which have licenses or copyrights of their own. In particular, see `scheno/Jonker_Volgenant/LICENSE.md`, `scheno/nt_wrappers/nauty27r4/COPYRIGHT`, and `scheno/nt_wrappers/nauty27r4_modified/COPYRIGHT`.

.

The rest of the code is licensed as follows:

Copyright (c) 2024 Justus Hibshman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
