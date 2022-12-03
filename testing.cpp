// NOTE: To work, these tests require that the protected members of the
//  NTSparseGraph class be made public. To test, just comment out the
//  "protected:" line in nt_sparse_graph.h

#include<iostream>
#include<random>
#include<string>

#include "debugging.h"
#include "file_utils.h"
#include "nt_sparse_graph.h"
#include "nauty_traces.h"


int main(void) {

#ifdef SYM__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
    std::cout<<"###################################################"<<std::endl
             <<"Checking internal consistency of NT representation."<<std::endl
             <<"###################################################"<<std::endl;

    for (bool directed : {false, true}) {
        std::cout<<std::endl<<"Directed: "<<directed<<std::endl<<std::endl;
        rand_test(0.02, 0.01,
                  0.57, 0.4,
                  directed, 8555, 2000);

        // This means that once we hit a stable state I expect there to be 0
        //  edges roughly ((.57 - .41)/.98) / (1 + ((.57 - .41)/.98)) ~= 13% of
        //  the time.
        rand_test(0.02, 0.0,
                  0.41, 0.57,
                  directed, 8555, 2000);

        rand_test(0.005, 0.0,
                  0.695, 0.3,
                  directed, 8555, 2000);
    }
#endif

    std::cout<<std::endl<<std::endl
             <<"###################################################"<<std::endl
             <<"Checking basic Traces correctness"<<std::endl
             <<"###################################################"<<std::endl
             <<std::endl;

    SYMNautyTracesOptions options;
    options.get_node_orbits = true;
    options.get_edge_orbits = true;
    options.get_canonical_node_order = true;

    Coloring<int> node_coloring = Coloring<int>();

    NTSparseGraph g_undir = NTSparseGraph(false, 7);

    SYMNautyTracesResults result = traces(g_undir, options);
    std::cout<<std::endl;
    std::cout<<"// Undirected Graph on 7 nodes with no edges"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_undir.add_edge(0, 1);
    g_undir.add_edge(2, 1);

    result = traces(g_undir, options);
    std::cout<<std::endl;
    std::cout<<"// Undirected Graph on 7 nodes with edges (0, 1), (1, 2)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    NTSparseGraph g_dir = NTSparseGraph(true, 7);
    g_dir.add_edge(0, 1);
    g_dir.add_edge(2, 1);
    result = traces(g_dir, options);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 7 nodes with edges (0, 1), (2, 1)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_dir.add_edge(0, 0);
    result = traces(g_dir, options);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 7 nodes with edges (0, 1), (2, 1), (0, 0)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_dir.delete_edge(0, 0);
    g_dir.delete_edge(2, 1);
    g_dir.add_edge(1, 2);
    result = traces(g_dir, options);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 7 nodes with edges (0, 1), (1, 2)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_dir.flip_edge(1, 2);
    g_dir.flip_edge(2, 1);
    node_coloring.set(0, 1);
    node_coloring.set(1, 1);
    node_coloring.set(2, 0);
    node_coloring.set(3, 0);
    node_coloring.set(4, 0);
    node_coloring.set(5, 0);
    node_coloring.set(6, 0);
    NTPartition partition = g_dir.nauty_traces_coloring(node_coloring);
    result = traces(g_dir, options, partition);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 7 nodes with edges (0, 1), (2, 1)"<<std::endl;
    std::cout<<"// Nodes 0 and 1 have been given color 1 while all else is 0."<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_dir.add_edge(0, 0);
    result = traces(g_dir, options, partition);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 7 nodes with edges (0, 1), (2, 1), (0, 0)"<<std::endl;
    std::cout<<"// Node 0 has been given a unique color -- (uniqueness overdetermined)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    // Edge colorings
    node_coloring = Coloring<int>();
    Coloring<Edge, EdgeHash> edge_coloring = Coloring<Edge, EdgeHash>();

    g_dir = NTSparseGraph(true, 4);
    g_dir.add_edge(0, 0);
    g_dir.add_edge(0, 1);
    g_dir.add_edge(0, 2);
    g_dir.add_edge(1, 0);
    g_dir.add_edge(1, 1);
    g_dir.add_edge(1, 2);
    g_dir.add_edge(2, 0);
    g_dir.add_edge(2, 1);
    g_dir.add_edge(2, 2);

    g_dir.add_edge(0, 3);
    g_dir.add_edge(1, 3);
    g_dir.add_edge(2, 3);

    edge_coloring.set(EDGE(0, 0, true), 7);
    edge_coloring.set(EDGE(0, 1, true), 0);
    edge_coloring.set(EDGE(0, 2, true), 0);
    edge_coloring.set(EDGE(1, 0, true), 0);
    edge_coloring.set(EDGE(1, 1, true), 7);
    edge_coloring.set(EDGE(1, 2, true), 0);
    edge_coloring.set(EDGE(2, 0, true), 0);
    edge_coloring.set(EDGE(2, 1, true), 0);
    edge_coloring.set(EDGE(2, 2, true), 0);

    edge_coloring.set(EDGE(0, 3, true), 0);
    edge_coloring.set(EDGE(1, 3, true), 0);
    edge_coloring.set(EDGE(2, 3, true), 0);

    partition = g_dir.nauty_traces_coloring(edge_coloring);
    std::cout<<"\n\nStarting Partition: "<<vec_as_string(partition.get_cell_list())<<std::endl;
    std::cout<<vec_as_string(std::vector<int>(partition.get_node_ids(), partition.get_node_ids() + partition.size()))<<std::endl;
    std::cout<<vec_as_string(std::vector<int>(partition.get_partition_ints(), partition.get_partition_ints() + partition.size()))<<std::endl;
    result = traces(g_dir, options, partition);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 4 nodes with ALL possible edges on the first "<<std::endl
             <<"    three (including self-loops) and the first 3 pointing to the 4th"<<std::endl;
    std::cout<<"//  Edges (0, 0) and (1, 1) have been highlighted."<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    edge_coloring.set(EDGE(2, 2, true), 7);
    edge_coloring.set(EDGE(0, 1, true), 7);
    edge_coloring.set(EDGE(1, 2, true), 7);
    edge_coloring.set(EDGE(2, 0, true), 7);
    partition = g_dir.nauty_traces_coloring(edge_coloring);
    result = traces(g_dir, options, partition);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 4 nodes with ALL possible edges on the first "<<std::endl
             <<"    three (including self-loops) and the first 3 pointing to the 4th"<<std::endl;
    std::cout<<"//  Edges (0, 0), (1, 1), (2, 2), (0, 1), (1, 2), and (2, 0) have been highlighted."<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_dir = NTSparseGraph(true, 8);
    g_dir.add_edge(0, 0);
    g_dir.add_edge(1, 1);
    g_dir.add_edge(2, 2);
    g_dir.add_edge(3, 3);
    g_dir.add_edge(4, 4);
    g_dir.add_edge(5, 5);
    g_dir.add_edge(6, 6);
    g_dir.add_edge(7, 7);
    node_coloring = Coloring<int>();
    edge_coloring = Coloring<Edge, EdgeHash>();

    node_coloring.set(0, 0);
    node_coloring.set(1, 0);
    node_coloring.set(2, 1);
    node_coloring.set(3, 1);
    node_coloring.set(4, 1);
    node_coloring.set(5, 1);
    node_coloring.set(6, 1);
    node_coloring.set(7, 2);

    edge_coloring.set(EDGE(0,0,true), 0);
    edge_coloring.set(EDGE(1,1,true), 1);
    edge_coloring.set(EDGE(2,2,true), 1);
    edge_coloring.set(EDGE(3,3,true), 1);
    edge_coloring.set(EDGE(4,4,true), 2);
    edge_coloring.set(EDGE(5,5,true), 2);
    edge_coloring.set(EDGE(6,6,true), 3);
    edge_coloring.set(EDGE(7,7,true), 3);

    partition = g_dir.nauty_traces_coloring(node_coloring, edge_coloring);

    result = traces(g_dir, options, partition);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 8 nodes with self loops only "<<std::endl;
    std::cout<<"//  Node Partition:      0, 1 | 2, 3, 4, 5, 6 | 7"<<std::endl;
    std::cout<<"//  Self-Loop Partition: 0 | 1, 2, 3 | 4, 5 | 6, 7"<<std::endl;
    std::cout<<"// Thus resultant partition should be 0 | 1 | 2, 3 | 4, 5 | 6 | 7"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_dir.add_edge(2, 3);
    g_dir.add_edge(3, 2);
    g_dir.add_edge(4, 5);
    g_dir.add_edge(5, 4);
    edge_coloring.set(EDGE(2,3,true), 0);
    edge_coloring.set(EDGE(3,2,true), 0);
    edge_coloring.set(EDGE(4,5,true), 1);
    edge_coloring.set(EDGE(5,4,true), 1);

    partition = g_dir.nauty_traces_coloring(node_coloring, edge_coloring);
    result = traces(g_dir, options, partition);
    std::cout<<std::endl;
    std::cout<<"// Adding directed edges between the nodes that share initial orbits."<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    std::cout<<std::endl<<std::endl
             <<"###################################################"<<std::endl
             <<"Checking Traces on Real-World Graphs"<<std::endl
             <<"###################################################"<<std::endl
             <<std::endl;

    options.get_node_orbits = true;
    options.get_edge_orbits = true;
    options.get_canonical_node_order = true;

    bool directed = true;

    NTSparseGraph real_graph(directed);
    std::vector<std::string> graph_names = {"not_real.g", "karate.g", "convote.g", "cora.g"};
    std::string graph_name;
    
    for (auto gn_itr = graph_names.begin(); gn_itr != graph_names.end(); gn_itr++) {
        graph_name = *gn_itr;
        real_graph = read_graph(directed, "real_world_graphs/" + graph_name);
        #ifdef SYM__NT_SPARSE_GRAPH_FULL_DEBUG_MODE
        consistency_check(real_graph);
        #endif
        std::cout<<std::endl;
        std::cout<<"// Running on "<<graph_name<<" as a directed graph."<<std::endl;
        std::cout<<"// Num Nodes: "<<real_graph.num_nodes()<<"   Num Edges: "<<real_graph.num_edges()<<std::endl;
        result = traces(real_graph, options);
        std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
        std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
        std::cout<<"Error Status: "<<result.error_status<<std::endl;
    }
    // std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    // std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    return 0;
};
