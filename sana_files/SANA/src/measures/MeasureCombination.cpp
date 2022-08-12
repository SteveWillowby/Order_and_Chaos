#include "MeasureCombination.hpp"
#include "localMeasures/LocalMeasure.hpp"
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <functional>

MeasureCombination::MeasureCombination(){
}

MeasureCombination::~MeasureCombination() {
  for(auto m: measures){
    delete m;
  }
}

double MeasureCombination::eval(const Alignment& A) const {
    uint n = measures.size();
    double total = 0;
    for (uint i = 0; i < n; i++) {
        if (weights[i] > 0) {
            total += measures[i]->eval(A) * weights[i];
        }
    }
    return total;
}

double MeasureCombination::eval(const string& measureName, const Alignment& A) const {
    uint n = measures.size();
    for (uint i = 0; i < n; i++) {
        if (measures[i]->getName() == measureName) {
            return measures[i]->eval(A);
        }
    }
    throw runtime_error("Measure not found.");
    return 0;
}

void MeasureCombination::addMeasure(Measure* m, double weight) {
    measures.push_back(m);
    weights.push_back(weight);
}

void MeasureCombination::rebalanceWeight(string& input){
    istringstream iss(input);
    vector<string> split{istream_iterator<string>{iss},istream_iterator<string>{}};

    for(uint i = 0; i < split.size(); i++){
        string name = split[i];
        int rebalancePos = -1;
        for(uint j = 0; j < measures.size(); j++){
            if(measures[j]->getName() == name){
                rebalancePos = j;
            }
        }

        if(rebalancePos == -1)
            throw runtime_error(name + " is not a measure that is being optimized.");

        if(!measures[rebalancePos]->isLocal())
            throw runtime_error(name + " is not a local measure.  You can only balance local measures.");

        double newWeight = weights[rebalancePos]*(measures[rebalancePos]->balanceWeight());

        cout << "Rebalancing " << name << " from weight = " << weights[rebalancePos] << " to weight = " << newWeight << ".  This will be renormalized with the other weights before SANA begins." << endl;

        weights[rebalancePos]=newWeight;
    }
    normalize();
}

void MeasureCombination::rebalanceWeight(){
    cout << "Rebalancing all local measures." << endl;
    for(uint i = 0; i < measures.size(); i++){
        if(measures[i]->isLocal()){
            double newWeight = weights[i]*(measures[i]->balanceWeight());

            cout << "Rebalancing " << measures[i]->getName() << " from weight = " << weights[i] << " to weight = " << newWeight << ".  This will be renormalized with the other weights before SANA begins." << endl;

            weights[i] = newWeight;
        }
    }
    normalize();
}

void MeasureCombination::addMeasure(Measure* m) {
    measures.push_back(m);
    weights.push_back(0);
}

void MeasureCombination::printWeights(ostream& ofs) const {
    uint n = measures.size();
    for (uint i = 0; i < n; i++) {
        if (weights[i] > 0) {
            ofs << "weight " << measures[i]->getName();
            ofs << ": " << weights[i] << endl;
        }
    }
}

void MeasureCombination::printMeasures(const Alignment& A, ostream& ofs) const {
    uint n = measures.size();
    for (uint i = 0; i < n; i++) {
        ofs << measures[i]->getName();
        ofs << ": " << measures[i]->eval(A) << endl;
    }
    uint count = 0;
    for (uint i = 0; i < n; i++) {
        if (weights[i] > 0) count++;
    }
    if (count > 1) {
        ofs << "Combined: " << eval(A) << " ( ";
        for (uint i = 0; i < n; i++) {
            if (weights[i] > 0) ofs << measures[i]->getName() << " " << weights[i] << " ";
        }
        ofs << ")" << endl;
    }
}

struct Crit {
    string name;
    double score;
    double weight;
};

bool critComp(Crit a, Crit b) {
    return a.weight > b.weight;
}

string MeasureCombination::toString() const {
    uint n = measures.size();
    vector<Crit> crits(0);
    Crit c;
    for (uint i = 0; i < n; i++) {
        if (weights[i] > 0) {
            c.name=measures[i]->getName();
            c.weight=weights[i];
            crits.push_back(c);
        }
    }
    sort(crits.begin(), crits.end(), critComp);
    string alig = "";
    for (uint i = 0; i < crits.size(); i++) {
        alig += crits[i].name + "_";
        if (crits[i].weight < 1) alig += extractDecimals(crits[i].weight, 2) + "_";
    }
    return alig.substr(0, alig.size()-1);;
}

void MeasureCombination::normalize() {
    normalizeWeights(weights);
}

double MeasureCombination::getWeight(const string& measureName) const {
    uint n = measures.size();
    for (uint i = 0; i < n; i++) {
        if (measures[i]->getName() == measureName) {
            return weights[i];
        }
    }
    throw runtime_error("Measure not found.");
}

Measure* MeasureCombination::getMeasure(const string& measureName) const {
    uint n = measures.size();
    for (uint i = 0; i < n; i++) {
        if (measures[i]->getName() == measureName) {
            return measures[i];
        }
    }
    throw runtime_error("Measure not found.");
}

Measure* MeasureCombination::getMeasure(int i) const {
    return measures[i];
}

bool MeasureCombination::containsMeasure(const string& measureName) {
    uint n = measures.size();
    for (uint i = 0; i < n; i++) {
        if (measures[i]->getName() == measureName) {
            return true;
        }
    }
    return false;
}

uint MeasureCombination::numMeasures() const {
    return measures.size();
}

void MeasureCombination::initn1n2(uint& n1, uint& n2) const {
    if(n1 != 0 && n2 != 0)
      return;
    for (uint i = 0; i < numMeasures(); i++) {
        Measure* m = measures[i];
        if (m->isLocal()) {
            vector<vector<float> >* mSims = ((LocalMeasure*) m)->getSimMatrix();
            n1 = mSims->size();
            n2 = (*mSims)[0].size();
            return;
        }
    }
    throw runtime_error("There are no local measures");
}

typedef vector<vector<float> > SimMatrix;
typedef function<void(SimMatrix &, uint const &, uint const &)> SimMatrixRecipe;

//Returns a reference to the similarity matrix of the weighted sum of local measures.
//Only initializes the matrix on the first call.
vector<vector<float> >& MeasureCombination::getAggregatedLocalSims() {
    //A flag to check if the map has been initialized.
    static bool is_init = false;
    //The "recipe" that describes how to create the sim matrix,
    //namely to combine all locals into a new localdo.
    static function<void(vector<vector<float> > &, uint const &, uint const &)> const initFunc =
      [this] (vector<vector<float> > & sim, uint const & n1, uint const & n2) {
        Measure* m;
        double w;
        for (uint i = 0; i < numMeasures(); i++) {
            m = measures[i];
            w = weights[i];
            if (m->isLocal() and w > 0) {
                vector<vector<float> >* mSims = ((LocalMeasure*) m)->getSimMatrix();
                for (uint i = 0; i < n1; i++) {
                    for (uint j = 0; j < n2; j++) {
                        sim[i][j] += w * (*mSims)[i][j];
                    }
                }
            }
        }
      };
    if(!is_init) {
      is_init = true;
      localAggregatedSim = initSim(initFunc);
    }
    return localAggregatedSim;
}

//Returns a map defined as LocalMeasure (string) -> SimMatrix
//Returns a reference to the map, initializing it on the first call
//and returning the existing map on subsequent calls.
map<string, SimMatrix>& MeasureCombination::getLocalSimMap() {
    Measure* m;
    double w;
    //A flag to check if the map has been initialized.
    static bool is_init = false;
    //The "recipe" that describes how to create the sim matrix.
    static SimMatrixRecipe const initFunc =
        [this, &m, &w] (SimMatrix & sim, uint const & n1, uint const & n2) {
            SimMatrix* mSims = ((LocalMeasure*) m)->getSimMatrix();
            for (uint i = 0; i < n1; i++) {
                for (uint j = 0; j < n2; j++) {
                    sim[i][j] = w * (*mSims)[i][j];
                }
            }
        };
    if(!is_init) {
      is_init = true;
      for (uint i = 0; i < numMeasures(); ++i) {
          m = measures[i];
          w = weights[i];
          if (m->isLocal() and w > 0) {
              localScoreSimMap[m->getName()] = initSim(initFunc);
          }
      }
    }
    return localScoreSimMap;
}

//Abstracts the construction of the similarity matrix. Instead of the get..()
//functions producing possibly different implementations of similarity matrices,
//a common type of similarity matrix is produced in initSim and populated
//by a Recipe function.
SimMatrix MeasureCombination::initSim(SimMatrixRecipe recipe) const {
  static uint n1 = 0, n2 = 0;
  initn1n2(n1, n2);
  vector<vector<float> > sim(n1, vector<float> (n2, 0));
  recipe(sim, n1, n2);
  return sim;
}

int MeasureCombination::getNumberOfLocalMeasures() const {
    int localMeasureCount = 0;
    for (uint i = 0; i < numMeasures(); i++) {
        Measure* m = measures[i];
        double w = weights[i];
        if (w > 0 and m->isLocal()) {
            localMeasureCount++;
        }
    }
    return localMeasureCount;
}

double MeasureCombination::getSumLocalWeight() const {
    double res = 0;
    for (uint i = 0; i < numMeasures(); i++) {
        if (measures[i]->isLocal()) {
            res += weights[i];
        }
    }
    return res;
}

void MeasureCombination::clearWeights() {
    for (uint i = 0; i < weights.size(); i++) {
        weights[i] = 0;
    }
}

void MeasureCombination::setWeight(const string& measureName, double weight) {
    uint n = measures.size();
    for (uint i = 0; i < n; i++) {
        if (measures[i]->getName() == measureName) {
            weights[i] = weight;
            return;
        }
    }
    throw runtime_error("Measure not found: "+measureName);
}

/*Writes out the local scores file in this format (example only of course):
Pairwise Alignment  LocalMeasure1       LocalMeasure2       Weighted Sum
821    723            0.334               0.214               0.548
*/
typedef unordered_map<uint,string> NodeIndexMap;
void MeasureCombination::writeLocalScores(ostream & outFile, Graph const & G1, Graph const & G2, Alignment const & A) const {
  NodeIndexMap mapG1 = G1.getIndexToNodeNameMap();
  NodeIndexMap mapG2 = G2.getIndexToNodeNameMap();
  int const COL_WIDTH = 20,
            PRECISION = 3;
  outFile << setw(COL_WIDTH) << left << "Pairwise Alignment";
  for(auto const & mapping : localScoreSimMap)
      outFile << setw(COL_WIDTH) << left << mapping.first;
  outFile << setw(COL_WIDTH) << left << "Weighted Sum" << endl;
  ostringstream edgeStream;
#define GENERATE_FULL_SIM_FILE 0
#if GENERATE_FULL_SIM_FILE
#define OTHER j
  for(uint i = 0; i < mapG1.size(); i++)
  for(uint j = 0; j < mapG2.size(); j++) {
#else // output only aligned pairs
#define OTHER A[i]
  for(uint i = 0; i < A.size(); ++i) {
#endif
      edgeStream << mapG1[i] << "\t" << mapG2[OTHER];
      outFile << setw(COL_WIDTH) << edgeStream.str();
      edgeStream.str("");
      edgeStream.clear();
      for(auto const & mapping : localScoreSimMap)
          outFile << setw(COL_WIDTH) << left << setprecision(PRECISION) << mapping.second[i][OTHER];
      outFile << setw(COL_WIDTH) << left << setprecision(PRECISION) << localAggregatedSim[i][OTHER] << endl;
  }
}
