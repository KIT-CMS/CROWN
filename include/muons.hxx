#ifndef GUARD_MUONS_H
#define GUARD_MUONS_H

namespace physicsobject {
namespace muon {

ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df, const std::string &outputname,
               const std::string &charge, const std::string &pt, 
               const std::string &eta, const std::string &phi, 
               const std::string &gen_pt, const std::string &n_tracker_layers,
               const std::string &rndm_values, const std::string &index_vector,
               const int &position, const std::string &filename, 
               const int error_set, const int error_member);
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df, const std::string &outputname,
                 const std::string &charge, const std::string &pt,
                 const std::string &eta, const std::string &phi,
                 const std::string &index_vector, const int &position, 
                 const std::string &filename, int error_set, int error_member);
} // namespace muon
} // namespace physicsobject
#endif /* GUARD_MUONS_H */