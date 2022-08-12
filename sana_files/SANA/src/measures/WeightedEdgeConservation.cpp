#include "WeightedEdgeConservation.hpp"

WeightedEdgeConservation::WeightedEdgeConservation(Graph* G1, Graph* G2, LocalMeasure* m) : Measure(G1, G2, "wec") {
    nodeSim = m;
}

WeightedEdgeConservation::~WeightedEdgeConservation() {
    //delete nodeSim; not sure if necessary
}

LocalMeasure* WeightedEdgeConservation::getNodeSimMeasure() {
    return nodeSim;
}

double WeightedEdgeConservation::eval(const Alignment& A) {
    vector<vector<float> >* simMatrix = nodeSim->getSimMatrix();
    vector<vector<uint> > edgeListG1;
    G1->getEdgeList(edgeListG1);
    Matrix<MATRIX_UNIT> matrixG2;
    G2->getMatrix(matrixG2);
    double score = 0;
    for (const auto& edge: edgeListG1) {
        uint node1 = edge[0], node2 = edge[1];
        if (matrixG2[A[node1]][A[node2]]) {
            score += (*simMatrix)[node1][A[node1]];
            score += (*simMatrix)[node2][A[node2]];
        }
    }
    //normalization factor to ensure 0 <= score <= 1
    //note that it is not part of the original WEC defined in WAVE
    score /= (2*G1->getNumEdges()); // Commenting it out breaks our duplication of LGRAAL, but may need to go when comparing to WAVE
    return score;
}
