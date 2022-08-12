#ifndef REPORT_HPP_
#define REPORT_HPP_

#include <string>

#include "Graph.hpp"
#include "Alignment.hpp"
#include "measures/MeasureCombination.hpp"
#include "methods/Method.hpp"

void saveReport(const Graph& G1, Graph& G2, const Alignment& A,
  const MeasureCombination& M, Method* method, string reportFile);

  void saveReport(const Graph& G1, Graph& G2, const Alignment& A,
    const MeasureCombination& M, Method* method, string reportFile, bool multiPairwiseIteration);

void saveLocalMeasures(Graph const & G1, Graph const & G2, Alignment const & A,
  MeasureCombination const & M, Method * const method, string & localMeasureFile);

void makeReport(const Graph& G1, Graph& G2, const Alignment& A,
  const MeasureCombination& M, Method* method, ofstream& stream);

void makeReport(const Graph& G1, Graph& G2, const Alignment& A,
    const MeasureCombination& M, Method* method, ofstream& stream, bool multiPairwiseIteration);

string ensureFileNameExistsAndOpenOutFile(string const & fileType, string outFileName, ofstream & outfile, string && G1Name, string && G2Name, Method * const & method, Alignment const & A);

#endif /* REPORT_HPP_ */
