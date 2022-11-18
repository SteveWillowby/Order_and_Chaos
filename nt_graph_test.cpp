// NOTE: To work, these tests require that the protected members of the
//  NTSparseGraph class be made public. To test, just comment out the
//  "protected:" line in nt_sparse_graph.h

#include<iostream>
#include<random>
#include<string>

#include "debugging.h"
#include "nt_sparse_graph.h"
#include "traces.h"


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

    NTSparseGraph g_undir = NTSparseGraph(false, 7);
    g_undir.add_edge(0, 1);
    g_undir.add_edge(2, 1);
    // 4 * 3 * 2 = 24

    SYMTracesOptions options;
    options.get_node_orbits = true;
    options.get_edge_orbits = true;
    options.get_canonical_node_order = true;

    SYMTracesResults result = traces(g_undir, options);
    std::cout<<"// Undirected Graph on 7 nodes with edges (0, 1), (1, 2)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    NTSparseGraph g_dir = NTSparseGraph(true, 7);
    g_dir.add_edge(0, 1);
    /*
    g_dir.add_edge(2, 1);
    result = traces(g_dir, options);
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 7 nodes with edges (0, 1), (2, 1)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    g_dir.delete_edge(2, 1);
    */
    g_dir.add_edge(1, 2);
    g_dir.add_edge(2, 3);
    std::cout<<"\nWhy does this next line cause a crash???"<<std::endl;
    result = traces(g_dir, options);
    std::cout<<"It doesn't cause a crash. Wat are you talkin about?"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"// Directed Graph on 7 nodes with edge (0, 1), (1, 2)"<<std::endl;
    std::cout<<"|Aut(G)| = "<<result.num_aut_base<<" x 10^"<<result.num_aut_exponent<<std::endl;
    std::cout<<"Num Orbits: "<<result.num_node_orbits<<std::endl;
    std::cout<<"Error Status: "<<result.error_status<<std::endl;
    std::cout<<"Node Orbits: "<<vec_as_string(result.node_orbits)<<std::endl;
    std::cout<<"Canonical Order: "<<vec_as_string(result.canonical_node_order)<<std::endl;

    return 0;
};
