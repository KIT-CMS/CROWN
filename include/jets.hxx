#ifndef GUARD_JETS_H
#define GUARD_JETS_H

#include "../include/utility/CorrectionManager.hxx"
#include "TRandom3.h"
#include "correction.h"


namespace physicsobject {
namespace jet {
namespace jec {
typedef struct jec_result_t {
    float jet_pt_l1;
    float jet_pt_l2rel;
    float jet_pt_l2l3res;
    float jet_pt_syst;
    float jet_pt_corr;
} JECResult;
const correction::Correction* load_nominal_jes_correction(
    correctionManager::CorrectionManager &correction_manager,
    const std::string &jec_file,
    const std::string &jes_tag,
    const std::string &type_tag,
    const std::string &jes_level,
    const std::string &jec_algo
);
const correction::Correction* load_shifted_jes_correction(
    correctionManager::CorrectionManager &correction_manager,
    const std::string &jec_file,
    const std::string &jer_tag,
    const std::string &type_tag,
    const std::string &jes_shift,
    const std::string &jec_algo
);
const correction::Correction* load_jer_correction(
    correctionManager::CorrectionManager &correction_manager,
    const std::string &jec_file,
    const std::string &jer_tag,
    const std::string &type_tag,
    const std::string &jer_parameter,
    const std::string &jec_algo
);
float apply_jes_l1 (
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_area,
    const float &rho,
    const correction::Correction *jes_l1_evaluator
);
float apply_jes_l2rel (
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const correction::Correction *jes_l2rel_evaluator
);
float apply_jes_l2l3res (
    const float &jet_pt,
    const float &jet_eta,
    const correction::Correction *jes_l2l3res_evaluator
);
float apply_jes_shifts (
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const UChar_t &jet_id,
    const std::vector<std::string> &jes_shift_sources,
    const int &jes_shift_factor,
    const std::vector<correction::Correction*> &jes_shift_evaluators
);
float apply_jer(
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const float &rho,
    const ROOT::RVec<float> &genjet_pt,
    const ROOT::RVec<float> &genjet_eta,
    const ROOT::RVec<float> &genjet_phi,
    const correction::Correction *jer_resolution_evaluator,
    const correction::Correction *jer_scalefactor_evaluator,
    const std::string &jer_shift,
    const float &jet_radius,
    const std::string &era,
    TRandom3 randgen
);
JECResult apply_full_jec_mc(
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const UChar_t &jet_id,
    const float &jet_area,
    const float &rho,
    const ROOT::RVec<float> &genjet_pt,
    const ROOT::RVec<float> &genjet_eta,
    const ROOT::RVec<float> &genjet_phi,
    const std::vector<std::string> &jes_shift_sources,
    const int &jes_shift_factor,
    const std::string &jer_shift,
    const float &jet_radius,
    const std::string &era,
    TRandom3 randgen,
    const correction::Correction* jes_l1_evaluator,
    const correction::Correction* jes_l2rel_evaluator,
    const std::vector<correction::Correction*> &jes_shift_evaluators,
    const correction::Correction *jer_resolution_evaluator,
    const correction::Correction *jer_scalefactor_evaluator
);
JECResult apply_jes_shifts_and_jer_mc(
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const UChar_t &jet_id,
    const float &rho,
    const ROOT::RVec<float> &genjet_pt,
    const ROOT::RVec<float> &genjet_eta,
    const ROOT::RVec<float> &genjet_phi,
    const std::vector<std::string> &jes_shift_sources,
    const int &jes_shift_factor,
    const std::string &jer_shift,
    const float &jet_radius,
    const std::string &era,
    TRandom3 randgen,
    const std::vector<correction::Correction*> &jes_shift_evaluators,
    const correction::Correction *jer_resolution_evaluator,
    const correction::Correction *jer_scalefactor_evaluator
);
JECResult apply_full_jec_data(
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const float &jet_area,
    const float &rho,
    const correction::Correction* jes_l1_evaluator,
    const correction::Correction* jes_l2rel_evaluator,
    const correction::Correction* jes_l2l3res_evaluator
);
ROOT::RDF::RNode Raw(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &jet_quantity,
    const std::string &jet_raw_factor
);
ROOT::RDF::RNode RawMuonSubtr(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &jet_quantity,
    const std::string &jet_raw_factor,
    const std::string &jet_muon_subtr_factor
);
ROOT::RDF::RNode PtCorrectionMC(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &output_jec_result,
    const std::string &output_l1,
    const std::string &output_l2rel,
    const std::string &output_l2l3res,
    const std::string &output_full,
    const std::string &jet_pt_raw,
    const std::string &jet_eta,
    const std::string &jet_phi,
    const std::string &jet_area,
    const std::string &jet_id,
    const std::string &genjet_pt,
    const std::string &genjet_eta,
    const std::string &genjet_phi,
    const std::string &rho,
    const std::string &jer_seed,
    const std::string &jec_file,
    const std::string &jec_algo,
    const std::string &jes_tag,
    const std::string &jer_tag,
    const std::vector<std::string> &jes_shift_sources,
    const int &jes_shift_factor,
    const std::string &jer_shift,
    const bool &reapply_jes,
    const std::string &era
);
ROOT::RDF::RNode PtCorrectionData(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &output_jec_result,
    const std::string &output_l1,
    const std::string &output_l2rel,
    const std::string &output_l2l3res,
    const std::string &output_full,
    const std::string &jet_pt_raw,
    const std::string &jet_eta,
    const std::string &jet_phi,
    const std::string &jet_area,
    const std::string &rho,
    const std::string &jec_file,
    const std::string &jec_algo,
    const std::string &jes_tag,
    const bool &reapply_jes,
    const std::string &era
);
ROOT::RDF::RNode MassCorrectionFromPt(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &jet_mass_raw,
    const std::string &jet_pt_raw,
    const std::string &jet_pt_corrected
);
}

ROOT::RDF::RNode RawPt(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &pts,
                        const std::string &jet_raw_factor);

ROOT::RDF::RNode
PtCorrectionL1(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &jet_pt,
        const std::string &jet_eta, 
        const std::string &jet_phi,
        const std::string &jet_area, 
        const std::string &jet_raw_factor,
        const std::string &jet_raw_muonfactor,
        const std::string &corrjet_pt,
        const std::string &corrjet_eta, 
        const std::string &corrjet_phi,
        const std::string &corrjet_area, 
        const std::string &corrjet_raw_muonfactor,
        const std::string &rho, 
        const std::string &jec_file, 
        const std::string &jec_algo,
        const std::string &jes_tag_mc, 
        const std::string &jes_tag_data, 
        const std::string &era,
        const bool &is_data,
        const bool &is_embedding);

ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &jet_pts,
        const std::string &jet_eta, 
        const std::string &jet_phi,
        const std::string &jet_area, 
        const std::string &jet_id,
        const std::string &corrjet_eta, 
        const std::string &corrjet_phi,
        const std::string &corrjet_area, 
        const std::string &gen_jet_pt,
        const std::string &gen_jet_eta, 
        const std::string &gen_jet_phi,
        const std::string &rho, 
        const std::string &jer_seed,
        const std::string &run, 
        const std::string &jec_file, 
        const std::string &jec_algo,
        const std::string &jes_tag_mc, 
        const std::string &jes_tag_data, 
        const std::vector<std::string> &jes_shift_sources,
        const std::string &jer_tag,
        const int &jes_shift, const std::string &jer_shift,
        const std::string &era, const bool &is_data,
        const bool &is_embedding);

ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &jet_pt,
               const std::string &jet_eta, const std::string &jet_phi,
               const std::string &jet_area, const std::string &jet_raw_factor,
               const std::string &jet_id, const std::string &gen_jet_pt,
               const std::string &gen_jet_eta, const std::string &gen_jet_phi,
               const std::string &rho, const std::string &jer_seed,
               const std::string &jec_file, const std::string &jec_algo,
               const std::string &jes_tag, const std::vector<std::string> &jes_shift_sources,
               const std::string &jer_tag, bool reapply_jes,
               const int &jes_shift, const std::string &jer_shift,
               const std::string &era, const bool &no_jer_for_unmatched_forward_jets = false);
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correction_manager,
                 const std::string &outputname, const std::string &jet_pt,
                 const std::string &jet_eta, const std::string &jet_phi,
                 const std::string &jet_area, const std::string &jet_raw_factor,
                 const std::string &rho, const std::string &run,
                 const std::string &jec_file, const std::string &jec_algo,
                 const std::string &jes_tag, const std::string &era);
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
BtaggingShape(ROOT::RDF::RNode df,
       correctionManager::CorrectionManager &correction_manager,
       const std::string &outputname, const std::string &pt, 
       const std::string &eta, const std::string &btag_value, 
       const std::string &flavor, const std::string &jet_mask, 
       const std::string &bjet_mask, const std::string &jet_veto_mask, 
       const std::string &sf_file, const std::string &sf_name, 
       const std::string &variation);  
ROOT::RDF::RNode
BtaggingWP(ROOT::RDF::RNode df,
       correctionManager::CorrectionManager &correction_manager,
       const std::string &outputname, const std::string &pt, 
       const std::string &eta, const std::string &flavor,
       const std::string &jet_mask, const std::string &bjet_mask,
       const std::string &jet_veto_mask, const std::string &sf_file,
       const std::string &sf_name, const std::string &variation,
       const std::string &btag_wp);  
} // end namespace scalefactor
} // end namespace jet
} // end namespace physicsobject
#endif /* GUARD_JETS_H */
