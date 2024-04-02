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
