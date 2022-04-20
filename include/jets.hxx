#ifndef GUARDJETS_H
#define GUARDJETS_H

namespace jet {
ROOT::RDF::RNode VetoOverlappingJets(auto &df, const std::string &output_col,
                                     const std::string &jet_eta,
                                     const std::string &jet_phi,
                                     const std::string &p4_1,
                                     const std::string &p4_2,
                                     const float &deltaRmin);
ROOT::RDF::RNode OrderJetsByPt(auto &df, const std::string &output_col,
                               const std::string &jet_pt,
                               const std::string &jetmask_name);
} // end namespace jet

namespace physicsobject {
namespace jet {

ROOT::RDF::RNode CutID(auto &df, const std::string &maskname,
                       const std::string &nameID, const int &idxID);
ROOT::RDF::RNode CutPUID(auto &df, const std::string &maskname,
                         const std::string &nameID, const std::string &jet_pt,
                         const int &idxID, const float &jet_pt_cut);
ROOT::RDF::RNode
JetPtCorrection(auto &df, const std::string &corrected_jet_pt,
                const std::string &jet_pt, const std::string &jet_eta,
                const std::string &jet_phi, const std::string &gen_jet_pt,
                const std::string &gen_jet_eta, const std::string &gen_jet_phi,
                const std::string &rho,
                const std::vector<std::string> &energy_shift_sources,
                const int &energy_shift_state, const int &energy_reso_shift);
ROOT::RDF::RNode CutRawID(auto &df, const std::string &quantity,
                          const std::string &maskname,
                          const float &idThreshold);
} // end namespace jet
} // end namespace physicsobject

namespace quantities {
namespace jet {
ROOT::RDF::RNode NumberOfJets(auto &df, const std::string &outputname,
                              const std::string &jetcollection);
ROOT::RDF::RNode btagValue(auto &df, const std::string &outputname,
                           const std::string &btagcolumn,
                           const std::string &jetcollection,
                           const int &position);
} // end namespace jet
} // end namespace quantities
#endif /* GUARDJETS_H */