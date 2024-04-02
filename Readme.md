# SCHENO, the SCHEma-NOise scoring function

###...Also includes a pattern-finding genetic algorithm guided by SCHENO as its fitness function.

## About

Given two graphs, SCHENO measures how well the first graph represents the underlying pattern(s) in the second graph.

### Isomorphism and Automorphism Calculations

If you are simply interested in this project's handy `C++` wrappers surrounding the `nauty` and `traces` isomorphism code, those wrappers are available in a standalone repository [here](https://github.com/schemanoise/Nauty_and_Traces).

## Setup

The code for the main binaries is all contained in this repository, so no outside libraries or packages should be necessary.

Run `setup.sh` to compile the utilities. Then run `compile.sh` to compile the executables.

The binaries will be located in the `executables` folder.

## Using the Two Programs

You can find a PDF manual describing the use of this repository's code inside the `documentation` folder.

For a full theoretical description of SCHENO, see the research paper located [here (TODO: Add Link)](https://github.com/schemanoise/SCHENO).

### `SCHENO_score`

This program takes a graph (expressed as an edgelist) and a noise set (also expressed as an edgelist). It reports how well the graph without the noise represents the underlying structure of the graph.

The graph can be directed or undirected.

To see the full help menu, run `executables/SCHENO_score -h`


### `SCHENO_ga`

This program takes a graph (expressed as an edgelist) and tries to find the best schema-noise decomposition for that graph using a simple genetic algorithm. The genetic algorithm uses the SCHENO score as its fitness function.

The graph can be directed or undirected.

To see the full help menu, run `executables/SCHENO_ga -h`


## Including SCHENO Scores in Your Own Program

You can find a PDF manual describing the interface for this repository's code inside the `documentation` folder.

If you want to call the code directly, you can include all the classes through the file `scheno.h` contained in the `scheno` subfolder.

The file `SCHENO_score.cpp` illustrates the use of most of the classes.

For scoring a bunch of candidate noise sets in parallel, consider using the `ThreadPoolScorer` class.

## License

This repository contains code from several sources, some of which have licenses or copyrights of their own. In particular, see `scheno/Jonker_Volgenant/LICENSE.md`, `scheno/nt_wrappers/nauty27r4/COPYRIGHT`, and `scheno/nt_wrappers/nauty27r4_modified/COPYRIGHT`.

The rest of the code was written by Justus Hibshman and is licensed as follows:

*Insert Standard MIT License*
