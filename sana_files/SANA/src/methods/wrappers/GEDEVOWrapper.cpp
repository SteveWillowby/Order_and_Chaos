#include "GEDEVOWrapper.hpp"

using namespace std;

const string CONVERTER = "python GWtoNTW.py";
const string PROGRAM   = "gedevo";
const string OUTPUT_CONVERTER = "./GEDEVOoutput.sh";
const string GLOBAL_PARAMETERS = " --undirected --no-prematch --pop 400 ";

//ARGUMENTS: --maxsecs <seconds> --blastpairlist [3 columns] --pop [400] --threads <N> [recommended runtime is "--maxsame 3000"]

GEDEVOWrapper::GEDEVOWrapper(Graph* G1, Graph* G2, string args): WrappedMethod(G1, G2, "GEDEVO", args) {
    wrappedDir = "wrappedAlgorithms/GEDEVO";
}

void GEDEVOWrapper::loadDefaultParameters() {
    parameters = "--maxsame 3000"; // maxsame 3000 is what they recommend, runtime many hours
}

string GEDEVOWrapper::convertAndSaveGraph(Graph* graph, string name) {
    string gwFile  = name + ".gw";
    string ntwFile = name + ".ntw";

    graph->saveInGWFormatWithNames(gwFile);
    exec("mv " + gwFile + " " + wrappedDir + "/" + gwFile);
    exec("cd " + wrappedDir + "; " + CONVERTER + " " + gwFile + " " + ntwFile);
    exec("mv " + wrappedDir + "/" + ntwFile + " " + ntwFile);

    return ntwFile;
}

string GEDEVOWrapper::generateAlignment() {
    string g1Sigs = g1TmpName + ".sigs";
    string g2Sigs = g2TmpName + ".sigs";
    G1->saveGraphletsAsSigs(wrappedDir + "/" + g1Sigs);
    G2->saveGraphletsAsSigs(wrappedDir + "/" + g2Sigs);

    string cmd = GLOBAL_PARAMETERS + " --save " + alignmentTmpName + " --no-save --groups " +
            g1TmpName + " " + g2TmpName +
            " --ntw " + g1File +
            " --ntw " + g2File +
            " --grsig " + g1Sigs + " " + g1TmpName +
            " --grsig " + g2Sigs + " " + g2TmpName +
            " --no-workfiles " +
            " " + parameters;

    cout << "\n\n\nrunning with: \"" + cmd + "\"" << endl << flush;

    // exec("cd " + wrappedDir + "; chmod +x " + PROGRAM);
    execPrintOutput("cd " + wrappedDir + "; " + "./" + PROGRAM + " " + cmd);

    exec("cd " + wrappedDir + ";" + OUTPUT_CONVERTER + " " + alignmentTmpName + ".matching " + alignmentTmpName);

    //exec("cd " + wrappedDir + "; rm " + g1Sigs + " " + g2Sigs);
    //exec("cd " + wrappedDir + "; rm " + alignmentTmpName + ".matching");

    return wrappedDir + "/" + alignmentTmpName;
}

Alignment GEDEVOWrapper::loadAlignment(Graph* G1, Graph* G2, string fileName) {
    return Alignment::loadPartialEdgeList(G1, G2, fileName, true);
}

void GEDEVOWrapper::deleteAuxFiles() {
    //exec("cd " + wrappedDir + "; rm -f " + g1File + " " + g2File + " *.gw" + " *.algn");
}
