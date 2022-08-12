#include "Experiment.hpp"
#include "AlphaEstimation.hpp"
#include "../measures/Measure.hpp"
#include "../measures/localMeasures/LocalMeasure.hpp"
#include "../measures/localMeasures/Sequence.hpp"
#include "../measures/localMeasures/Importance.hpp"
#include "../measures/EdgeCorrectness.hpp"
#include "../measures/SymmetricSubstructureScore.hpp"
#include "../measures/WeightedEdgeConservation.hpp"
#include "../measures/localMeasures/GraphletLGraal.hpp"
#include "../utils/utils.hpp"
#include "../Alignment.hpp"
#include "../utils/Timer.hpp"
#include <cassert>

double AlphaEstimation::computeAlpha(Graph& G1, Graph& G2, string methodName, Measure* topMeasure) {
    string g1Name = G1.getName();
    string g2Name = G2.getName();
    methodName = toLowerCase(methodName);
    string aligFileAlpha0 = "experiments/"+methodName+"_alpha0/"+g1Name+"_"+g2Name+"_"+methodName+"_alpha0.txt";
    string aligFileAlpha1 = "experiments/"+methodName+"_alpha1/"+g1Name+"_"+g2Name+"_"+methodName+"_alpha1.txt";
    if (methodName == "hubalign") {
        string aligFileAlpha0 = "experiments/"+methodName+"_alpha0/"+g1Name+"_"+g2Name+"_"+methodName+"_alpha0.txt";
        string aligFileAlpha1 = "experiments/"+methodName+"_alpha0,9999/"+g1Name+"_"+g2Name+"_"+methodName+"_alpha0.9999.txt";
    }
    if (not fileExists(aligFileAlpha0)) return -1;
    if (not fileExists(aligFileAlpha1)) return -1;

    double topScore = topMeasure->eval(Alignment::loadMapping(aligFileAlpha0));
    Sequence seq(&G1, &G2);
    double seqScore = seq.eval(Alignment::loadMapping(aligFileAlpha1));

    return topScore/(topScore+seqScore);
}

void AlphaEstimation::run(ArgumentParser& args) {
    string alphaEstimation = args.strings["-alphaestimation"];
    string experFile = "experiments/"+alphaEstimation;
    assert(fileExists(experFile));
    init(experFile);
    printData(experFile+".out");
}

AlphaEstimation::AlphaEstimation(){}

string AlphaEstimation::getName(void) {
    return "AlphaEstimation";
}

void AlphaEstimation::init(string alphaFile) {
    vector<vector<string> > content = fileToStringsByLines(alphaFile);
    methods = content[0];
    for (uint i = 1; i < content.size(); i++) {
        networkPairs.push_back(content[i]);
    }

    cout << "Loading graphs...";
    Timer T;
    T.start();
    for (auto pair : networkPairs) {
        string g1Name = pair[0];
        string g2Name = pair[1];
        if (graphs.count(g1Name) == 0) {
            Graph::loadGraph(g1Name, graphs[g1Name]);
        }
        if (graphs.count(g2Name) == 0) {
            Graph::loadGraph(g2Name, graphs[g2Name]);
        }
    }
    cout << "graph loading done ("+T.elapsedString()+")" << endl;

    cout << "Computing alphas...";
    T.start();
    computeAlphas();
    cout << "AlphaEstimation::init done ("+T.elapsedString()+")" << endl;
}

double AlphaEstimation::computeAlphaSANA(Graph& G1, Graph& G2, Measure* topMeasure) {
    string g1Name = G1.getName();
    string g2Name = G2.getName();
    string optName = toLowerCase(topMeasure->getName());
    if (optName.substr(0, 6) == "lgraal") optName = optName.substr(0, 6);
    string aligFileAlpha0 = "experiments/sana_"+optName+"_alpha0/"+g1Name+"_"+g2Name+"_sana_alpha0.txt";
    double topScore = 0;
    uint sampleCount = 10;
    for (uint i = 0; i < 10; i++) {
        if (i == 1) aligFileAlpha0 += "_2";
        else if (i > 1 and i < 9) aligFileAlpha0[aligFileAlpha0.size()-1] = '0'+i+1;
        else if (i == 9) {
            aligFileAlpha0[aligFileAlpha0.size()-1] = '1';
            aligFileAlpha0 += "0";
        }
        if (not fileExists(aligFileAlpha0)) {
            sampleCount--;
        }
        else {
            topScore += topMeasure->eval(Alignment::loadMapping(aligFileAlpha0));
        }

    }
    topScore /= sampleCount;

    string aligFileAlpha1 = "experiments/sana_seq/"+g1Name+"_"+g2Name+"_sana_alpha1.txt";
    double seqScore = 0;
    Sequence seq(&G1, &G2);
    sampleCount = 10;
    for (uint i = 0; i < 10; i++) {
        if (i == 1) aligFileAlpha1 += "_2";
        else if (i > 1 and i < 9) aligFileAlpha1[aligFileAlpha1.size()-1] = '0'+i+1;
        else if (i == 9) {
            aligFileAlpha1[aligFileAlpha1.size()-1] = '1';
            aligFileAlpha1 += "0";
        }
        if (not fileExists(aligFileAlpha1)) {
            sampleCount--;
        }
        else {
            seqScore += seq.eval(Alignment::loadMapping(aligFileAlpha1));
        }
    }
    seqScore /= sampleCount;
    return topScore/(topScore+seqScore);
}

void AlphaEstimation::printData(string outputFile) {
    uint nMethods = methods.size();
    uint nPairs = networkPairs.size();

    ofstream fout(outputFile.c_str());
    for (uint i = 0; i < nMethods; i++) {
        for (uint j = 0; j < nPairs; j++) {
            fout << methods[i] << "\t" << networkPairs[j][0] << "\t";
            fout << networkPairs[j][1] << "\t" << alphas[i][j] << endl;
        }
    }
}

//needs to be refactorized
void AlphaEstimation::computeAlphas() {
    uint nMethods = methods.size();
    uint nPairs = networkPairs.size();
    alphas = vector<vector<double> > (nMethods, vector<double> (nPairs, -1));
    for (uint i = 0; i < nMethods; i++) {
        string methodName = methods[i];
        cout << methodName << endl;
        Measure* topMeasure;
        for (uint j = 0; j < nPairs; j++) {
            string g1Name = networkPairs[j][0];
            string g2Name = networkPairs[j][1];

            if (methodName == "LGRAAL") {
                LocalMeasure* m = new GraphletLGraal(&graphs[g1Name], &graphs[g2Name]);
                topMeasure = new WeightedEdgeConservation(&graphs[g1Name], &graphs[g2Name], m);
            }
            else if (methodName == "HubAlign") topMeasure = new Importance(&graphs[g1Name], &graphs[g2Name]);
            else if (methodName == "SANA_EC") topMeasure = new EdgeCorrectness(&graphs[g1Name], &graphs[g2Name]);
            else /*if (methodName == "sanas3")*/ topMeasure = new SymmetricSubstructureScore(&graphs[g1Name], &graphs[g2Name]);

            if (methodName.substr(0,4) == "SANA") {
                alphas[i][j] = computeAlphaSANA(graphs[g1Name], graphs[g2Name], topMeasure);
            }
            else {
                alphas[i][j] = computeAlpha(graphs[g1Name], graphs[g2Name], methodName, topMeasure);
            }
            delete topMeasure;
        }
    }
}

AlphaEstimation::AlphaEstimation(string alphaFile) {
    init(alphaFile);
}

double AlphaEstimation::getAlpha(string alphaFile, string methodName, string G1Name, string G2Name) {
    vector<vector<string> > content = fileToStringsByLines(alphaFile);
    for (uint i = 0; i < content.size(); i++) {
        if (content[i][0] == methodName and content[i][1] == G1Name and content[i][2] == G2Name)
            return stod(content[i][3]);
    }
    throw runtime_error("Could not find the alpha value");
}
