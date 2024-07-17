#ifndef GUARDHADRECOIL_H
#define GUARDHADRECOIL_H

#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>

/// Namespace used for calculating quantities related to the hadronic recoil

namespace hadrecoil {

ROOT::RDF::RNode scalar_ht(ROOT::RDF::RNode df, const std::string &outputname, const std::string &jet_pt, const std::string &jet_collection);

ROOT::RDF::RNode vectorial_mht(ROOT::RDF::RNode df, const std::string &outputname, const std::string &jet_pt, const std::string &jet_eta, const std::string &jet_collection);

} // namespace hadrecoil

#endif /* GUARDHADRECOIL_H */