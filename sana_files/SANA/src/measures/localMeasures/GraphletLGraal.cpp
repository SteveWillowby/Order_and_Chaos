#include <vector>
#include <iostream>
#include <cmath>
#include "GraphletLGraal.hpp"
using namespace std;

GraphletLGraal::GraphletLGraal(Graph* G1, Graph* G2) : LocalMeasure(G1, G2, "graphletlgraal") {
    string subfolder = autogenMatricesFolder+getName()+"/";
    createFolder(subfolder);
    string fileName = subfolder+G1->getName()+"_"+G2->getName()+"_graphletlgraal.bin";
    loadBinSimMatrix(fileName);
}

GraphletLGraal::~GraphletLGraal() {
}

double GraphletLGraal::gdvSim(uint i, uint j,
    const vector<vector<uint> >& gdvsG1, const vector<vector<uint> >& gdvsG2) {
    double total = 0;
    for (uint k = 0; k < 15; k++) {
        double m = max(gdvsG1[i][k],gdvsG2[j][k]);
        if (m > 0) total += min(gdvsG1[i][k],gdvsG2[j][k])/m;
    }
    return total/15;
}

void GraphletLGraal::initSimMatrix() {
    uint n1 = G1->getNumNodes();
    uint n2 = G2->getNumNodes();
    sims = vector<vector<float> > (n1, vector<float> (n2, 0));
    vector<vector<uint> > gdvs1 = G1->loadGraphletDegreeVectors();
    vector<vector<uint> > gdvs2 = G2->loadGraphletDegreeVectors();

    for (uint i = 0; i < n1; i++) {
        for (uint j = 0; j < n2; j++) {
            sims[i][j] = gdvSim(i,j,gdvs1,gdvs2);
        }
    }
}
