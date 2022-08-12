#ifndef GRAPHLETNORM_HPP
#define GRAPHLETNORM_HPP
#include "LocalMeasure.hpp"

class GraphletNorm: public LocalMeasure {
public:
    GraphletNorm(Graph* G1, Graph* G2);
    virtual ~GraphletNorm();

private:
    void initSimMatrix();
    double magnitude(vector<uint> &vector);
    const uint NUM_ORBITS = 73;
    vector<double> NODV(vector<uint> &v);
    double ODVratio(vector<double> &u, vector<double> &v, uint i);
    double RMS_ODVdiff1(vector<uint> &u, vector<uint> &v);
    double ODVsim(vector<uint> &u, vector<uint> &v);
    vector<uint> reduce(vector<uint> &v);
    
};

#endif
