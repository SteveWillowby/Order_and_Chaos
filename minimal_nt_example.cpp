// minimal_nt_example.cpp

#include "scheno/scheno.h"

#include<cmath>
#include<iostream>

int main(void) {

    bool directed = false;
    NTSparseGraph karate(directed);
    karate = read_graph(directed, "nt_test_graphs/karate.txt");

    std::cout<<"# Nodes: "<<karate.num_nodes()<<std::endl;
    std::cout<<"# Edges: "<<karate.num_edges()<<std::endl<<std::endl;

    NautyTracesOptions nto;
    nto.get_node_orbits = true;
    nto.get_edge_orbits = true;
    nto.get_canonical_node_order = true;
    
    NautyTracesResults ntr = traces(karate, nto);
    
    double log10_aut = std::log10(ntr.num_aut_base) + ntr.num_aut_exponent;
    std::cout<<"Log10 of Automorphisms:  "<<log10_aut<<std::endl;
    std::cout<<"Number of Automorphisms: "<<(std::pow(10.0, log10_aut))<<std::endl;
    std::cout<<"Number of Node Orbits:   "<<ntr.node_orbits.colors().size()<<std::endl;
    std::cout<<"Number of Edge Orbits:   "<<ntr.edge_orbits.colors().size()<<std::endl;
    std::cout<<"First Node in Canonical Ordering: "
             <<ntr.canonical_node_order[0]<<std::endl;

    return 0;
}
