#ifndef GUARDJETS_H
#define GUARDJETS_H

namespace jet {
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &output_col,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &p4_1, const std::string &p4_2,
                    const float &deltaRmin);
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &output_col,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &p4_1, const float &deltaRmin);
ROOT::RDF::RNode VetoOverlappingJetsIsoLepOnly(ROOT::RDF::RNode df,
                                               const std::string &output_col,
                                               const std::string &jet_eta,
                                               const std::string &jet_phi,
                                               const std::string &p4_1,
                                               const std::string &lep_is_iso,
                                               const float &deltaRmin);
ROOT::RDF::RNode OrderJetsByPt(ROOT::RDF::RNode df,
                               const std::string &output_col,
                               const std::string &jet_pt,
                               const std::string &jetmask_name);
} // end namespace jet

namespace physicsobject {
namespace jet {

ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID, const int &idxID);
ROOT::RDF::RNode CutPUID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const std::string &jet_pt,
                         const int &idxID, const float &jet_pt_cut);
ROOT::RDF::RNode
JetPtCorrection(ROOT::RDF::RNode df, const std::string &corrected_jet_pt,
                const std::string &jet_pt, const std::string &jet_eta,
                const std::string &jet_phi, const std::string &jet_area,
                const std::string &jet_rawFactor, const std::string &jet_ID,
                const std::string &gen_jet_pt, const std::string &gen_jet_eta,
                const std::string &gen_jet_phi, const std::string &rho,
                bool reapplyJES,
                const std::vector<std::string> &jes_shift_sources,
                const int &jes_shift, const std::string &jer_shift,
                const std::string &jec_file, const std::string &jer_tag,
                const std::string &jes_tag, const std::string &jec_algo);
ROOT::RDF::RNode
JetPtCorrection_data(ROOT::RDF::RNode df, const std::string &corrected_jet_pt,
                     const std::string &jet_pt, const std::string &jet_eta,
                     const std::string &jet_area,
                     const std::string &jet_rawFactor, const std::string &rho,
                     const std::string &jec_file, const std::string &jes_tag,
                     const std::string &jec_algo);
ROOT::RDF::RNode
BJetPtCorrection(ROOT::RDF::RNode df, const std::string &corrected_bjet_pt,
                           const std::string &jet_pt, const std::string &good_bjet_mask,
                           const std::string &corr_factor);
ROOT::RDF::RNode CutRawID(ROOT::RDF::RNode df, const std::string &quantity,
                          const std::string &maskname,
                          const float &idThreshold);
ROOT::RDF::RNode AntiCutRawID(ROOT::RDF::RNode df, const std::string &quantity,
                              const std::string &maskname,
                              const float &idThreshold);
} // end namespace jet
} // end namespace physicsobject

namespace quantities {
namespace jet {
ROOT::RDF::RNode NumberOfJets(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &jetcollection);
ROOT::RDF::RNode btagValue(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &btagcolumn,
                           const std::string &jetcollection,
                           const int &position);
ROOT::RDF::RNode flavor(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &flavorcolumn,
                        const std::string &jetcollection, const int &position);
ROOT::RDF::RNode gen_flavor(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &flavorcolumn,
                        const std::string &jetcollection, const int &position);
} // end namespace jet
} // end namespace quantities
#endif /* GUARDJETS_H */
