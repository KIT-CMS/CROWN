#ifndef GUARD_JETS_H
#define GUARD_JETS_H

namespace physicsobject {
namespace jet {

ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &jet_pt,
               const std::string &jet_eta, const std::string &jet_phi,
               const std::string &jet_area, const std::string &jet_raw_factor,
               const std::string &jet_id, const std::string &gen_jet_pt,
               const std::string &gen_jet_eta, const std::string &gen_jet_phi,
               const std::string &rho, const std::string &jec_file,
               const std::string &jec_algo, const std::string &jes_tag,
               const std::vector<std::string> &jes_shift_sources,
               const std::string &jer_tag, bool reapply_jes,
               const int &jes_shift, const std::string &jer_shift,
               const int jer_seed = 42, const int& lhc_run = 2);
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correction_manager,
                 const std::string &outputname, const std::string &jet_pt,
                 const std::string &jet_eta, const std::string &jet_area,
                 const std::string &jet_raw_factor, const std::string &rho,
                 const std::string &jec_file, const std::string &jec_algo,
                 const std::string &jes_tag);
ROOT::RDF::RNode
PtCorrectionBJets(ROOT::RDF::RNode df,
                  const std::string &outputname,
                  const std::string &jet_pt,
                  const std::string &scale_factor,
                  const std::string &bjet_mask);
ROOT::RDF::RNode CutPileupID(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &jet_pu_id,
                             const std::string &jet_pt, const int &pu_id_cut,
                             const float &pt_cut);
ROOT::RDF::RNode
ApplyVetoMap(ROOT::RDF::RNode df, const std::string &outputname,
             const std::string &jet_eta, const std::string &jet_phi,
             const std::string &vetomap_file, const std::string &vetomap_name,
             const std::string &vetomap_type);
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &target_p4_1,
                    const std::string &target_p4_2, const float &min_delta_r);
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &target_p4, const float &min_delta_r);
ROOT::RDF::RNode VetoOverlappingJetsWithIsoLepton(ROOT::RDF::RNode df,
                                                  const std::string &outputname,
                                                  const std::string &jet_eta,
                                                  const std::string &jet_phi,
                                                  const std::string &lepton_p4,
                                                  const std::string &lepton_iso,
                                                  const float &min_delta_r);

namespace quantity {
ROOT::RDF::RNode 
ID(ROOT::RDF::RNode df,
              correctionManager::CorrectionManager &correction_manager,
              const std::string &outputname,
              const std::string &jet_eta,
              const std::string &jet_chHEF,
              const std::string &jet_neHEF,
              const std::string &jet_chEmEF,
              const std::string &jet_neEmEF,
              const std::string &jet_muEF,
              const std::string &jet_chMult,
              const std::string &jet_neMult,
              const std::string &jet_id_file,
              const std::string &jet_name);
} // end namespace quantity

namespace scalefactor {

ROOT::RDF::RNode
Btagging(ROOT::RDF::RNode df,
       correctionManager::CorrectionManager &correction_manager,
       const std::string &outputname, const std::string &pt, 
       const std::string &eta, const std::string &btag_value, 
       const std::string &flavor, const std::string &jet_mask, 
       const std::string &bjet_mask, const std::string &jet_veto_mask, 
       const std::string &sf_file, const std::string &sf_name, 
       const std::string &variation);  
} // end namespace scalefactor
} // end namespace jet
} // end namespace physicsobject
#endif /* GUARD_JETS_H */
