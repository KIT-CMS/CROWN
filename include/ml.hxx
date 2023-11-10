#ifndef GUARD_ML_H
#define GUARD_ML_H

#include "basefunctions.hxx"
#include "defaults.hxx"
#include "utility/Logger.hxx"
#include "vectoroperations.hxx"
#include <Math/Vector4D.h>


std::vector<std::vector<double>> readCSV(const std::string& filename);

namespace ml {

ROOT::RDF::RNode GaussianTransform(ROOT::RDF::RNode df,
				   const std::string &inputname,
				   const std::string &outputname,
				   const std::string &paramfile,
				   const unsigned position,
				   const std::string &vartype);


namespace sofie {

  void CompileModelForRDF(const std::string & headerModelFile, unsigned int ninputs, unsigned int nslots=0);

  ROOT::RDF::RNode KerasEvaluate(ROOT::RDF::RNode df, const std::vector<std::string> &input_vec,
				 const std::string &outputname,
				 const std::string &modelfile);

} // end namespace sofie
} // end namespace ml
#endif /* GUARD_ML_H */
