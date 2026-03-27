#include "../include/defaults.hxx"
#include "../include/jets.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <algorithm>

namespace physicsobject {
namespace jet {
namespace jec {

const correction::Correction* load_nominal_jes_correction(
    correctionManager::CorrectionManager &correction_manager,
    const std::string &jec_file,
    const std::string &jes_tag,
    const std::string &type_tag,
    const std::string &jes_level,
    const std::string &jec_algo
) {
    return correction_manager.loadCorrection(
        jec_file,
        jes_tag + "_" + type_tag + "_" + jes_level + "_" + jec_algo
    );
}

const correction::Correction* load_shifted_jes_correction(
    correctionManager::CorrectionManager &correction_manager,
    const std::string &jec_file,
    const std::string &jer_tag,
    const std::string &type_tag,
    const std::string &jes_shift,
    const std::string &jec_algo
) {
    return correction_manager.loadCorrection(
        jec_file,
        jer_tag + "_" + type_tag + "_" + jes_shift + "_" + jec_algo
    );
}

const correction::Correction* load_jer_correction(
    correctionManager::CorrectionManager &correction_manager,
    const std::string &jec_file,
    const std::string &jer_tag,
    const std::string &type_tag,
    const std::string &jer_parameter,
    const std::string &jec_algo
) {
    return correction_manager.loadCorrection(
        jec_file,
        jer_tag + "_" + type_tag + "_" + jer_parameter + "_" + jec_algo
    );
}

float apply_jes_l1 (
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_area,
    const float &rho,
    const correction::Correction *jes_l1_evaluator
) {
    // Calculate L1FastJet-corrected pt
    return jet_pt * jes_l1_evaluator->evaluate({
        jet_area, jet_eta, jet_pt, rho
    });
}

float apply_jes_l2rel (
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const correction::Correction *jes_l2rel_evaluator
) {
    // Calculate the L2rel-corrected pt
    return jet_pt * jes_l2rel_evaluator->evaluate({
        jet_eta, jet_phi, jet_pt 
    });
}


float apply_jes_l2l3res (
    const float &jet_pt,
    const float &jet_eta,
    const float &run,
    const correction::Correction *jes_l2l3res_evaluator
) {
    // Calculate the L2rel-corrected pt
    return jet_pt * jes_l2l3res_evaluator->evaluate({
        run, jet_eta, jet_pt 
    });
}

float apply_jes_shifts (
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const UChar_t &jet_id,
    const std::vector<std::string> &jes_shift_sources,
    const int &jes_shift_factor,
    const std::vector<correction::Correction*> &jes_shift_evaluators
) {
    float jet_pt_corr;
    if (jes_shift_sources.at(0) == "HEMIssue") {
        // To assign an uncertainty to the HEM issue, the jet pt needs to be
        // manually scaled in a specific phase space region.
        float sf = 1.0;
        if (
            jes_shift_factor == -1
            && jet_pt > 15.0
            && jet_phi > -1.57
            && jet_phi < -0.87
            && jet_id == 2
        ) {
            if (jet_eta > -2.5 && jet_eta < -1.3) {
                sf = 0.8;
            } else if (jet_eta > -3.0 && jet_eta <= -2.5) {
                sf = 0.65;
            }
        }
        jet_pt_corr = sf * jet_pt;
    } else {
        // Calculate the relative difference to the nominal corrected pt.
        // If multiple sources are combined, the squared sum of the
        // differences is computed.
        float delta_squared = 0.;
        for (const auto &evaluator : jes_shift_evaluators) {
            delta_squared += std::pow(
                evaluator->evaluate({jet_eta, jet_pt}),
                2
            );
        }

        jet_pt_corr = jet_pt * (
            1 + jes_shift_factor * std::sqrt(delta_squared)
        );
    }
    return jet_pt_corr;
}

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
) {
    // Get the JER MC resolution and data-MC scale factor for the smearing
    auto resol = jer_resolution_evaluator->evaluate({
        jet_eta, jet_pt, rho
    });
    auto sf = 1.0;
    if (std::stoi(era.substr(0, 4)) <= 2018) { // with run 2 inputs 
        sf = jer_scalefactor_evaluator->evaluate({
            jet_eta, jer_shift
        });
    } else {
        sf = jer_scalefactor_evaluator->evaluate({ // with run 3 inputs
            jet_eta, jet_pt, jer_shift
        });
    }

    // Match the jet to a generator-level jet
    float min_delta_r = 0;
    float gen_pt = -10.0;
    for (int i = 0; i < genjet_pt.size(); ++i) {
        float delta_r = ROOT::VecOps::DeltaR(
            jet_eta, jet_phi, genjet_eta[i], genjet_phi[i]
        );
        if (delta_r > min_delta_r) {
            continue;
        } 
        if (
            delta_r < (jet_radius / 2.)
            && std::abs(jet_pt - genjet_pt[i]) < (3.0 * resol * jet_pt)
        ) {
            min_delta_r = delta_r;
            gen_pt = genjet_pt[i];
        }
    }

    // Smear the JER by either using the rescaling or the random smearing
    // method
    float delta_jer = 0.0;
    if (gen_pt >= 0.0) {
        delta_jer = (sf - 1.0) * (jet_pt - gen_pt) / jet_pt;
    } else {
        // Jet horn mitigation for run 3
        // If no generator-level jet is found for a reconstructed jet in
        // 2.5 < eta < 3.0, no smearing shall be applied.
        if (
            std::stoi(era.substr(0, 4)) > 2022
            && abs(jet_eta) > 2.5
            && abs(jet_eta) < 3.0
        ) {
            delta_jer = 0.0;
        } else {
            delta_jer = randgen.Gaus(0, resol) * (
                std::sqrt(std::max(std::pow(sf, 2) - 1.0, 0.0))
            );
        }
    }

    return jet_pt * std::max(0.f, 1.f + delta_jer);
}

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
) {
    // Apply the consecutive steps of the jet energy calibration
    auto jet_pt_l1 = apply_jes_l1(
        jet_pt,
        jet_eta,
        jet_area,
        rho,
        jes_l1_evaluator
    );
    auto jet_pt_l2rel = apply_jes_l2rel (
        jet_pt_l1,
        jet_eta,
        jet_phi,
        jes_l2rel_evaluator
    );
    auto jet_pt_syst = apply_jes_shifts(
        jet_pt_l2rel,
        jet_eta,
        jet_phi,
        jet_id,
        jes_shift_sources,
        jes_shift_factor,
        jes_shift_evaluators
    );
    auto jet_pt_jer = apply_jer(
        jet_pt_syst,
        jet_eta,
        jet_phi,
        rho,
        genjet_pt,
        genjet_eta,
        genjet_phi,
        jer_resolution_evaluator,
        jer_scalefactor_evaluator,
        jer_shift,
        jet_radius,
        era,
        randgen
    );
     
    // Create the JECResult which also contains intermediate results of the
    // calibration
    auto result = physicsobject::jet::jec::JECResult{
        jet_pt_l1,
        jet_pt_l2rel,
        jet_pt_l2rel, // not applied to MC, so this is just the correction after L2rel
        jet_pt_syst,
        jet_pt_jer
    };

    return result;
}

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
) {
    // Apply the jet energy scale shifts and the resolution smearing
    auto jet_pt_syst = apply_jes_shifts(
        jet_pt,
        jet_eta,
        jet_phi,
        jet_id,
        jes_shift_sources,
        jes_shift_factor,
        jes_shift_evaluators
    );
    auto jet_pt_jer = apply_jer(
        jet_pt_syst,
        jet_eta,
        jet_phi,
        rho,
        genjet_pt,
        genjet_eta,
        genjet_phi,
        jer_resolution_evaluator,
        jer_scalefactor_evaluator,
        jer_shift,
        jet_radius,
        era,
        randgen
    );

    // Create the JECResult which also contains intermediate results of the
    // calibration
    auto result = JECResult{
        jet_pt,  // the input is the already JES-corrected pt
        jet_pt,  // the input is the already JES-corrected pt
        jet_pt,  // the input is the already JES-corrected pt
        jet_pt_syst,
        jet_pt_jer
    };

    return result;
}

JECResult apply_full_jec_data(
    const float &jet_pt,
    const float &jet_eta,
    const float &jet_phi,
    const float &jet_area,
    const float &rho,
    const unsigned int &run,
    const correction::Correction* jes_l1_evaluator,
    const correction::Correction* jes_l2rel_evaluator,
    const correction::Correction* jes_l2l3res_evaluator
) {
    // Apply the consecutive steps of the jet energy calibration
    auto jet_pt_l1 = apply_jes_l1(
        jet_pt,
        jet_eta,
        jet_area,
        rho,
        jes_l1_evaluator
    );
    auto jet_pt_l2rel = apply_jes_l2rel (
        jet_pt_l1,
        jet_eta,
        jet_phi,
        jes_l2rel_evaluator
    );
    auto jet_pt_l2l3res = apply_jes_l2l3res(
        jet_pt_l2rel,
        jet_eta,
        static_cast<float>(run),
        jes_l2l3res_evaluator
    );
     
    // Create the JECResult which also contains intermediate results of the
    // calibration
    auto result = JECResult{
        jet_pt_l1,
        jet_pt_l2rel,
        jet_pt_l2l3res, 
        jet_pt_l2l3res, // not applied to data, so this is just the correction after L2L3res
        jet_pt_l2l3res  // no JER smearing applied to data, so this is just the correction after L2L3res
    };

    return result;
}

ROOT::RDF::RNode Raw(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &jet_quantity,
    const std::string &jet_raw_factor
) {
    // Event-level lambda that computes the jet pt before JEC for each jet in
    // the event
    auto func = [] (
        const ROOT::RVec<float> &jet_quantity,
        const ROOT::RVec<float> &jet_raw_factor
    ) {
        return jet_quantity * (1 - jet_raw_factor);
    };

    return df.Define(
        outputname, 
        func,
        {
            jet_quantity,
            jet_raw_factor
        }
    );
}

ROOT::RDF::RNode RawMuonSubtr(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &jet_quantity,
    const std::string &jet_raw_factor,
    const std::string &jet_muon_subtr_factor
) {
    // Event-level lambda that computes the jet pt before JEC and with 
    // subtracting the muon component for each jet in the event
    auto callable = [] (
        const ROOT::RVec<float> &jet_quantity,
        const ROOT::RVec<float> &jet_raw_factor,
        const ROOT::RVec<float> &jet_muon_subtr_factor
    ) {
        return jet_quantity * (1 - jet_raw_factor) * (1 - jet_muon_subtr_factor);
    };

    return df.Define(
        outputname, 
        callable,
        {
            jet_quantity,
            jet_raw_factor,
            jet_muon_subtr_factor
        }
    );
}

ROOT::RDF::RNode RawMuonSubtr(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &jet_quantity,
    const std::string &jet_muon_subtr_factor
) {
    // Event-level lambda that computes the jet pt before JEC and with 
    // subtracting the muon component for each jet in the event
    auto callable = [] (
        const ROOT::RVec<float> &jet_quantity,
        const ROOT::RVec<float> &jet_muon_subtr_factor
    ) {
        return jet_quantity * ( 1 - jet_muon_subtr_factor);
    };

    return df.Define(
        outputname, 
        callable,
        {
            jet_quantity,
            jet_muon_subtr_factor
        }
    );
}

/**
 * @brief This function applies the full jet energy calibration (JEC) procedure
 * to MC according to the recommendations of the JME POG. The corrections are
 * implemented as a multi-step procedure, where the corrected jet \f$p_T\f$ of
 * the previous step serves as input for the next step.
 *
 * The jet energy scale (JES) in simulation is corrected in two steps
 * ([JES recipe](https://cms-jerc.web.cern.ch/JES/)):
 * - `L1FastJet`: Correction for the offset of the energy measurement due to
 *   pileup.
 * - `L2Relative`: Correction to bring the response of reconstructed-level jets
 *   relative to particle-level jets to 1.
 * These steps are only processed by this function if the parameter `reapply_jes` is set to `true`. Otherwise, make sure that the `jet_pt_raw` function contains the already JES-corrected jet \f$p_T\f$ values, e.g., from the input NANOAOD file.
 * 
 * Systematic shifts encoding different sources of uncertainties are
 * applied on top of the outcome of these corrections. Recommendations for the
 * uncertainty scheme to use can be found in the
 * [JES uncertainty recipe](https://cms-jerc.web.cern.ch/JECUncertaintySources/).
 *
 * On top of the outcome of the previous steps, a jet energy resolution smearing
 * is performed, details can be found in the 
 * [JER recipe](https://cms-jerc.web.cern.ch/JER/). This function uses the
 * hybrid method to smear the resolution of reconstructed jets in simulation.
 * The JER corrections are also associated with systematic uncertainties.
 *
 * The corrections are provided as JSON files provided by the JME POG. The
 * most recent recommendations on the use of these files can be found in
 * the
 * [JERC documentation](https://cms-jerc.web.cern.ch/Recommendations/#jet-energy-scale).
 * Specifications of the data format for corresponding eras can be found in the
 * [CMS analysis corrections documentation](https://cms-analysis-corrections.docs.cern.ch/corrections/JME/).
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the jet energy correction file
 * @param output_jec_result name of the output column for storing a `JECResult` object containing intermediate results of the calibration of the jet \f$p_T\f$
 * @param output_l1 name of the output column for corrected jet \f$p_T\f$ after the `L1FastJet` correction level
 * @param output_l2rel name of the output column for corrected jet \f$p_T\f$ after the `L2Relative` correction level
 * @param output_l2l3res name of the output column for corrected jet \f$p_T\f$ after the `L2L3Residual` correction level
 * @param output_full name of the output column for corrected jet \f$p_T\f$ after the full procedure
 * @param jet_pt_raw collection column of raw jet \f$p_T\f$ before the calibration; if `reapply_jes` is set to `true`, this column is interpreted as the already JES-corrected jet \f$p_T\f$.
 * @param jet_eta collection column of jet \f$\eta\f$
 * @param jet_phi collection column of jet \f$\phi\f$
 * @param jet_area collection column of area that the clustered object covers in \f$\eta\f$-\f$\phi\f$ plane
 * @param jet_id collection column of jet ID
 * @param genjet_pt collection column of particle-level jet \f$p_T\f$
 * @param genjet_eta collection column of particle-level jet \f$\eta\f$
 * @param genjet_phi collection column of particle-level jet \f$\phi\f$
 * @param rho column containing the average event energy density
 * @param jer_seed column with eventwise seed value for the random number generator that is used for the jet energy resolution smearing
 * @param jec_file path to the JEC correction file
 * @param jec_algo name of the jet reconstruction algorithm (e.g., "AK4PFchs" or "AK8PFPuppi")
 * @param jes_tag tag of the JES correction campaign (e.g., "Summer19UL18_V5")
 * @param jer_tag tag of the JER correction campaign (e.g., "Summer19UL18_JRV2")
 * @param jes_shift_sources list of JES shift sources for the systematic shift to be applied
 * @param jes_shift_factor factor of the JES shift variation (0 = nominal, +/-1 = up/down)
 * @param jer_shift name of the JER shift variation ("nom", "up", or "down")
 * @param reapply_jes flag to reapply the jet energy calibration, otherwise `jet_pt_raw` is taken as the already JES-corrected \f$p_T\f$
 * @param era string defining the currently processed era, needed due to different kind of recommendations from JME POG for different eras
 * 
 * @return A dataframe with a new column of corrected jet \f$p_T\f$'s
 */
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
) {
    // In nanoAODv12 the type of jet/fatjet ID was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, jet_id_v12] =
        utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, jet_id+"_v12", "ROOT::VecOps::RVec<UChar_t>", jet_id
        )
    ;

    // Identify jet radius from algorithm
    float jet_radius = 0.4;
    if (jec_algo.find("AK8") != std::string::npos) {
        jet_radius = 0.8;
    }

    // Set the type tag to "MC"
    const std::string type_tag = "MC";
    
    // Load the nominal jet energy scale evaluators
    auto jes_l1_evaluator = load_nominal_jes_correction(
        correction_manager,
        jec_file,
        jes_tag,
        type_tag,
        "L1FastJet",
        jec_algo
    );
    auto jes_l2rel_evaluator = load_nominal_jes_correction(
        correction_manager,
        jec_file,
        jes_tag,
        type_tag,
        "L2Relative",
        jec_algo
    );

    // Load the jet energy scale variation evaluators
    std::vector<correction::Correction *> jes_shift_evaluators;
    for (const auto &source : jes_shift_sources) {
        if (source != "" && source != "HEMIssue") {
            auto evaluator = const_cast<correction::Correction*>(
                load_shifted_jes_correction(
                    correction_manager,
                    jec_file,
                    jes_tag,
                    type_tag,
                    source,
                    jec_algo
                )
            );
            jes_shift_evaluators.push_back(evaluator);
        }
    }

    // Load the jet energy resolution evaluators
    auto jer_resolution_evaluator = load_jer_correction(
        correction_manager,
        jec_file,
        jer_tag,
        type_tag,
        "PtResolution",
        jec_algo
    );
    auto jer_scalefactor_evaluator = load_jer_correction(
        correction_manager,
        jec_file,
        jer_tag,
        type_tag,
        "ScaleFactor",
        jec_algo
    );

    // Function to retrieve the JEC result with intermediate steps
    auto func_jec_result = [
        jes_shift_sources,
        jes_shift_factor,
        jer_shift,
        jet_radius, 
        reapply_jes,
        era,
        jes_l1_evaluator,
        jes_l2rel_evaluator,
        jes_shift_evaluators,
        jer_resolution_evaluator,
        jer_scalefactor_evaluator
    ] (
        const ROOT::RVec<float> &jet_pt_raw,
        const ROOT::RVec<float> &jet_eta,
        const ROOT::RVec<float> &jet_phi,
        const ROOT::RVec<UChar_t> &jet_id,
        const ROOT::RVec<float> &jet_area,
        const ROOT::RVec<float> &genjet_pt,
        const ROOT::RVec<float> &genjet_eta,
        const ROOT::RVec<float> &genjet_phi,
        const float &rho, 
        const unsigned int &seed
    ) {
        // Random value generator for jet energy resolution smearing
        TRandom3 randgen = TRandom3(seed);

        ROOT::RVec<JECResult> jet_jec_result;
        if (reapply_jes) {
            // Apply the jet energy scale and resolution corrections to MC
            // events. This is done by using the corresponding helper function 
            // for single jets and wrap it with ROOT::VecOps::Map to retrieve 
            // the calibrated momenta for the full collection.
            jet_jec_result = ROOT::VecOps::Map(
                jet_pt_raw,
                jet_eta,
                jet_phi,
                jet_id,
                jet_area,
                [
                    rho,
                    genjet_pt,
                    genjet_eta,
                    genjet_phi,
                    jes_shift_sources,
                    jes_shift_factor,
                    jer_shift,
                    jet_radius,
                    era,
                    randgen,
                    jes_l1_evaluator,
                    jes_l2rel_evaluator,
                    jes_shift_evaluators,
                    jer_resolution_evaluator,
                    jer_scalefactor_evaluator
                ] (
                    const float &jet_pt,
                    const float &jet_eta,
                    const float &jet_phi,
                    const float &jet_id,
                    const float &jet_area
                ) {
                    return apply_full_jec_mc(
                        jet_pt,
                        jet_eta,
                        jet_phi,
                        jet_id,
                        jet_area,
                        rho,
                        genjet_pt,
                        genjet_eta,
                        genjet_phi,
                        jes_shift_sources,
                        jes_shift_factor,
                        jer_shift,
                        jet_radius,
                        era,
                        randgen,
                        jes_l1_evaluator,
                        jes_l2rel_evaluator,
                        jes_shift_evaluators,
                        jer_resolution_evaluator,
                        jer_scalefactor_evaluator
                    );
                }
            );
        } else {
            // The jet_pt_raw is already the jet energy scale-corrected pt.
            // Only shifts and and resolution corrections need to be applied
            // on top of it. This is done by using the corresponding helper
            // function for single jets and wrap it with ROOT::VecOps::Map to 
            // retrieve the calibrated momenta for the full collection.
            jet_jec_result = ROOT::VecOps::Map(
                jet_pt_raw,
                jet_eta,
                jet_phi,
                jet_id,
                jet_area,
                [
                    rho,
                    genjet_pt,
                    genjet_eta,
                    genjet_phi,
                    jes_shift_sources,
                    jes_shift_factor,
                    jer_shift,
                    jet_radius,
                    era,
                    randgen,
                    jes_shift_evaluators,
                    jer_resolution_evaluator,
                    jer_scalefactor_evaluator
                ] (
                    const float &jet_pt,
                    const float &jet_eta,
                    const float &jet_phi,
                    const float &jet_id,
                    const float &jet_area
                ) {
                    return apply_jes_shifts_and_jer_mc(
                        jet_pt,
                        jet_eta,
                        jet_phi,
                        jet_id,
                        rho,
                        genjet_pt,
                        genjet_eta,
                        genjet_phi,
                        jes_shift_sources,
                        jes_shift_factor,
                        jer_shift,
                        jet_radius,
                        era,
                        randgen,
                        jes_shift_evaluators,
                        jer_resolution_evaluator,
                        jer_scalefactor_evaluator
                    );
                }
            );
        }
        return jet_jec_result;
    };

    // Function to store the L1FastJet step outcome in a column
    auto func_pt_l1 = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_l1 = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_l1;
            }
        );
        return jet_pt_l1;
    };

    // Function to store the L2Rel step outcome in a column
    auto func_pt_l2rel = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_l2rel = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_l2rel;
            }
        );
        return jet_pt_l2rel;
    };

    // Function to store the L2L3Residual step outcome in a column
    auto func_pt_l2l3res = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_l2l3res = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_l2l3res;
            }
        );
        return jet_pt_l2l3res;
    };

    // Function to store the full procedure outcome in a column
    auto func_pt_full = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_final = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_corr;
            }
        );
        return jet_pt_final;
    };

    // Store the JECResult
    auto df2 = df1.Define(
        output_jec_result,
        func_jec_result,
        {
            jet_pt_raw,
            jet_eta,
            jet_phi,
            jet_id_v12,
            jet_area,
            genjet_pt,
            genjet_eta,
            genjet_phi,
            rho, 
            jer_seed
        }
    );

    // Store the L1FastJet-corrected pt
    auto df3 = df2.Define(output_l1, func_pt_l1, { output_jec_result });
 
    // Store the L2Relative-corrected pt
    auto df4 = df3.Define(output_l2rel, func_pt_l2rel, { output_jec_result });

    // Store the L2L3Residual-corrected pt
    auto df5 = df4.Define(output_l2l3res, func_pt_l2l3res, { output_jec_result });

    // Store the corrected pt after the full JEC procedure
    auto df6 = df5.Define(output_full, func_pt_full, { output_jec_result });

    return df6;
}

/**
 * @brief This function applies the full jet energy calibration (JEC) procedure
 * to data according to the recommendations of the JME POG. The corrections are
 * implemented as a multi-step procedure, where the corrected jet \f$p_T\f$ of
 * the previous step serves as input for the next step.
 *
 * The jet energy scale (JES) in simulation is corrected in two steps
 * ([JES recipe](https://cms-jerc.web.cern.ch/JES/)):
 * - `L1FastJet`: Correction for the offset of the energy measurement due to
 *   pileup.
 * - `L2Relative`: Correction to bring the response of reconstructed-level jets
 *   relative to particle-level jets to 1.
 * - `L2L3Residual`: Relative (\f$\eta\f$-dependent) and absolute
 *   (\f$p_T\f$-dependent) corrections to the jet \f$p_T\f to correct for
 *   residual differences between data and simulation, introduced in the
 *   `L2Relative` step.
 * These steps are only processed by this function if the parameter
 * `reapply_jes` is set to `true`. Otherwise, make sure that the `jet_pt_raw`
 * function contains the already JES-corrected jet \f$p_T\f$ values, e.g., from
 * the input NANOAOD file.
 *
 * The corrections are provided as JSON files provided by the JME POG. The
 * most recent recommendations on the use of these files can be found in
 * the
 * [JERC documentation](https://cms-jerc.web.cern.ch/Recommendations/#jet-energy-scale).
 * Specifications of the data format for corresponding eras can be found in the
 * [CMS analysis corrections documentation](https://cms-analysis-corrections.docs.cern.ch/corrections/JME/).
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the jet energy correction file
 * @param output_jec_result name of the output column for storing a `JECResult` object containing intermediate results of the calibration of the jet \f$p_T\f$
 * @param output_l1 name of the output column for corrected jet \f$p_T\f$ after the `L1FastJet` correction level
 * @param output_l2rel name of the output column for corrected jet \f$p_T\f$ after the `L2Relative` correction level
 * @param output_l2l3res name of the output column for corrected jet \f$p_T\f$ after the `L2L3Residual` correction level
 * @param output_full name of the output column for corrected jet \f$p_T\f$ after the full procedure
 * @param jet_pt_raw collection column of raw jet \f$p_T\f$ before the calibration; if `reapply_jes` is set to `true`, this column is interpreted as the already JES-corrected jet \f$p_T\f$.
 * @param jet_eta collection column of jet \f$\eta\f$
 * @param jet_phi collection column of jet \f$\phi\f$
 * @param jet_area collection column of area that the clustered object covers in \f$\eta\f$-\f$\phi\f$ plane
 * @param rho column containing the average event energy density
 * @param run run index of the event
 * @param jec_file path to the JEC correction file
 * @param jec_algo name of the jet reconstruction algorithm (e.g., "AK4PFchs" or "AK8PFPuppi")
 * @param jes_tag tag of the JES correction campaign (e.g., "Summer19UL18_V5")
 * @param jer_shift name of the JER shift variation ("nom", "up", or "down")
 * @param reapply_jes flag to reapply the jet energy calibration, otherwise `jet_pt_raw` is taken as the already JES-corrected \f$p_T\f$
 * @param era string defining the currently processed era, needed due to different kind of recommendations from JME POG for different eras
 *
 * @return A dataframe with a new column of corrected jet \f$p_T\f$'s
 */
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
    const std::string &run,
    const std::string &jec_file,
    const std::string &jec_algo,
    const std::string &jes_tag,
    const bool &reapply_jes,
    const std::string &era
) {
    // Identify jet radius from algorithm
    float jet_radius = 0.4;
    if (jec_algo.find("AK8") != std::string::npos) {
        jet_radius = 0.8;
    }

    // Set the type tag to "DATA"
    const std::string type_tag = "DATA";
    
    // Load the nominal jet energy scale evaluators
    auto jes_l1_evaluator = load_nominal_jes_correction(
        correction_manager,
        jec_file,
        jes_tag,
        type_tag,
        "L1FastJet",
        jec_algo
    );
    auto jes_l2rel_evaluator = load_nominal_jes_correction(
        correction_manager,
        jec_file,
        jes_tag,
        type_tag,
        "L2Relative",
        jec_algo
    );
    auto jes_l2l3res_evaluator = load_nominal_jes_correction(
        correction_manager,
        jec_file,
        jes_tag,
        type_tag,
        "L2L3Residual",
        jec_algo
    );

    // Function to retrieve the JEC result with intermediate steps
    auto func_jec_result = [
        reapply_jes,
        jes_l1_evaluator,
        jes_l2rel_evaluator,
        jes_l2l3res_evaluator
    ] (
        const ROOT::RVec<float> &jet_pt_raw,
        const ROOT::RVec<float> &jet_eta,
        const ROOT::RVec<float> &jet_phi,
        const ROOT::RVec<float> &jet_area,
        const float &rho,
        const unsigned int &run
    ) {
        ROOT::RVec<JECResult> jet_jec_result;
        if (reapply_jes) {
            // Apply the jet energy scale corrections to data. This is done by
            // using the corresponding helper function for single jets and wrap 
            // it with ROOT::VecOps::Map to retrieve the calibrated momenta for
            // the full collection.
            jet_jec_result = ROOT::VecOps::Map(
                jet_pt_raw,
                jet_eta,
                jet_phi,
                jet_area,
                [
                    rho,
                    run,
                    jes_l1_evaluator,
                    jes_l2rel_evaluator,
                    jes_l2l3res_evaluator
                ] (
                    const float &jet_pt,
                    const float &jet_eta,
                    const float &jet_phi,
                    const float &jet_area
                ) {
                    return apply_full_jec_data(
                        jet_pt,
                        jet_eta,
                        jet_phi,
                        jet_area,
                        rho,
                        run,
                        jes_l1_evaluator,
                        jes_l2rel_evaluator,
                        jes_l2l3res_evaluator
                    );
                }
            );
        } else {
            // Nothing needs to be done here as the input jet_pt_raw is
            // interpreted as the already corrected jet pt. ROOT::VecOps::Map
            // is used to create dummy JECResult objects in order to allow
            // consistent processing with the other branch of this if else
            // query.
            jet_jec_result = ROOT::VecOps::Map(
                jet_pt_raw,
                [] (const float &jet_pt) {
                    return JECResult{
                        jet_pt,  // the input is the already JES-corrected pt
                        jet_pt,  // the input is the already JES-corrected pt
                        jet_pt,  // the input is the already JES-corrected pt
                        jet_pt,  // the input is the already JES-corrected pt
                        jet_pt,  // the input is the already JES-corrected pt
                    };
                }
            );
        }

        return jet_jec_result;
    };

    // Function to store the L1FastJet step outcome in a column
    auto func_pt_l1 = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_l1 = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_l1;
            }
        );
        return jet_pt_l1;
    };

    // Function to store the L2Rel step outcome in a column
    auto func_pt_l2rel = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_l2rel = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_l2rel;
            }
        );
        return jet_pt_l2rel;
    };

    // Function to store the L2L3Residual step outcome in a column
    auto func_pt_l2l3res = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_l2l3res = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_l2l3res;
            }
        );
        return jet_pt_l2l3res;
    };

    // Function to store the full procedure outcome in a column
    auto func_pt_full = [] (
        const ROOT::RVec<JECResult> &jec_result
    ) {
        // Retrieve the result from the JECResult struct eventwise and wrap
        // with ROOT::VecOps::Map to get the result collection.
        auto jet_pt_final = ROOT::VecOps::Map(
            jec_result,
            [] (const JECResult &jec_result) {
                return jec_result.jet_pt_corr;
            }
        );
        return jet_pt_final;
    };

    // Store the JECResult
    auto df1 = df.Define(
        output_jec_result,
        func_jec_result,
        {
            jet_pt_raw,
            jet_eta,
            jet_phi,
            jet_area,
            rho,
            run
        }
    );

    // Store the L1FastJet-corrected pt
    auto df2 = df1.Define(output_l1, func_pt_l1, { output_jec_result });
 
    // Store the L2Relative-corrected pt
    auto df3 = df2.Define(output_l2rel, func_pt_l2rel, { output_jec_result });

    // Store the L2L3Residual-corrected pt
    auto df4 = df3.Define(output_l2l3res, func_pt_l2l3res, { output_jec_result });

    // Store the corrected pt after the full JEC procedure
    auto df5 = df4.Define(output_full, func_pt_full, { output_jec_result });

    return df5;
}

/**
 * @brief This function modifies the mass of objects in an event using the
 * formula \f[ M_{\text{corrected},i} = M_{\text{raw},i} \times
 * \frac{p_{T,\text{corrected},i}}{p_{T,\text{raw},i}} \f] for each object of an
 * object collection in the event. The correction is applied element-wise to the
 * mass vector and is needed as part of for example energy scale corrections
 * that were before to the transverse momenta.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the corrected masses
 * @param raw_mass name of the column containing raw object masses
 * @param raw_pt name of the column containing raw object transverse momenta
 * @param corrected_pt name of the column containing corrected transverse
 * momenta
 *
 * @return a dataframe with a new column
 */
ROOT::RDF::RNode MassCorrectionFromPt(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &jet_mass_raw,
    const std::string &jet_pt_raw,
    const std::string &jet_pt_corrected
) {
    // Function to align the jet mass to the corrected jet pt
    auto func = [](
        const ROOT::RVec<float> &jet_mass_raw,
        const ROOT::RVec<float> &jet_pt_raw,
        const ROOT::RVec<float> &jet_pt_corrected
    ) {
        return jet_mass_raw * jet_pt_corrected / jet_pt_raw;
    };

    auto df1 = df.Define(
        outputname,
        func,
        {
            jet_mass_raw,
            jet_pt_raw,
            jet_pt_corrected
        }
    );
    return df1;
}

} // end namespace jec

/**
 * @brief This function applies jet energy scale corrections (JES) and jet
 * energy resolution (JER) smearing to simulated jets using correction factors.
 * In nanoAOD the JES corrections are already applied to jets, however, if new
 * corrections are measured it is possible to reapply them using the newest
 * correction files by setting `reapply_jes` to `true`.
 *
 * The second part of the jet energy corrections (JEC) is the resolution
 * correction. For that this function follows CMS recommendations on how to
 * apply the jet energy smearing using the hybrid method.
 * (Ref. https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution)
 *
 * This function is able to process the JEC uncertainties, which includes the
 * induvidual jes uncertainties as well as the merged uncertainty scheme.
 * Additionally, the HEM issue (2018) can be included as an uncertainty
 * based on https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * jet energy correction file
 * @param outputname name of the output column for corrected jet \f$p_T\f$'s
 * @param jet_pt name of the column containing the jet transverse momenta
 * @param jet_eta name of the column containing the jet pseudorapidities
 * @param jet_phi name of the column containing the jet azimuthal angles
 * @param jet_area name of the column containing the jet catchment area
 * @param jet_raw_factor name of the column containing the jet factors to
 * calculate back to the raw jet \f$p_T\f$'s
 * @param jet_id name of the column containing the jet IDs
 * @param gen_jet_pt name of the column containing the generator-level jet
 * \f$p_T\f$'s
 * @param gen_jet_eta name of the column containing generator-level jet etas
 * @param gen_jet_phi name of the column containing generator-level jet phis
 * @param rho name of the column containing the event energy density
 * @param jer_seed seed value for the random number generator that is used for
 * the jet energy resolution smearing
 * @param jec_file path to the JEC correction file
 * @param jec_algo name of the jet reconstruction algorithm (e.g., "AK4PFchs" or
 * "AK8PFPuppi")
 * @param jes_tag tag of the JES correction campaign (e.g.,
 * "Summer19UL18_V5_MC")
 * @param jes_shift_sources list of JES shift sources for systematic
 * uncertainties
 * @param jer_tag tag of the JER correction campaign (e.g.,
 * "Summer19UL18_JRV2_MC")
 * @param reapply_jes boolean flag to indicate whether to the JES correction
 * should be reapplied
 * @param jes_shift JES shift variation (0 = nominal, +/-1 = up/down)
 * @param jer_shift JER shift variation ("nom", "up", or "down")
 * @param era string defining the currently processed era, needed due to different
 * kind of recommendations from JME POG for different eras
 * @param no_jer_for_unmatched_forward_jets if true, no jet energy resolution
 * smearing is applied to unmatched jets in the forward region
 * (\f$|\eta| > 2.5\f$).
 * 
 * @return a dataframe with a new column of corrected jet \f$p_T\f$'s
 *
 * @note If jets with \f$p_T\f$ > 15 GeV are corrected, this change should be
 * propagated to the missing transverse momentum.
 * (see `physicsobject::PropagateToMET`)
 *
 * @note The option `no_jer_for_unmatched_forward_jets` is introduced to
 * mitigate jet horns appearing in the eta distribution of jets at
 * \f$|\eta| \sim 2-3\f$ in Run 3 analyses. If the option is set to `true`,
 * no jet energy resolution smearing is applied to jets without a matching
 * generator-level jet for \f$|\eta| > 2.5\f$. Only use this option set to
 * `true` if the recipe to mitigate jet horns is implemented in the analysis.
 */
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
               const std::string &era, const bool &no_jer_for_unmatched_forward_jets) {
    // In nanoAODv12 the type of jet/fatjet ID was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, jet_id_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, jet_id+"_v12", "ROOT::VecOps::RVec<UChar_t>", jet_id);

    // identifying jet radius from algorithm
    float jet_radius = 0.4;
    if (jec_algo.find("AK8") != std::string::npos) {
        jet_radius = 0.8;
    }
    // loading JES variations
    std::vector<correction::Correction *> jet_energy_scale_shifts;
    for (const auto &source : jes_shift_sources) {
        if (source != "" && source != "HEMIssue") {
            auto jes_source_evaluator = const_cast<correction::Correction *>(
                correction_manager.loadCorrection(
                    jec_file, jes_tag + "_" + source + "_" + jec_algo));
            jet_energy_scale_shifts.push_back(jes_source_evaluator);
        }
    };
    // loading jet energy correction scale factor function
    auto jes_evaluator = correction_manager.loadCompoundCorrection(
        jec_file, jes_tag + "_L1L2L3Res_" + jec_algo);
    
    // Create a unified lambda that handles both era cases
    auto jet_energy_scale_sf = [jes_evaluator](const float area, const float eta,
                                               const float pt, const float rho,
                                               const float phi, const std::string era) {
        if (std::stoi(era.substr(0, 4)) <= 2022 || era == "2023preBPix") {
            // run 2 and 2022 to 2023preBPix cases
            return jes_evaluator->evaluate({area, eta, pt, rho});
        } else {
            // run 3 from 2023postBpix onwards
            return jes_evaluator->evaluate({area, eta, pt, rho, phi});
        }
    };
    // loading relative pT resolution function
    auto jer_resolution_evaluator = correction_manager.loadCorrection(
        jec_file, jer_tag + "_PtResolution_" + jec_algo);
    auto jet_energy_resolution = [jer_resolution_evaluator](const float eta,
                                                            const float pt,
                                                            const float rho) {
        return jer_resolution_evaluator->evaluate({eta, pt, rho});
    };
    // loading JER scale factor function
    auto jer_sf_evaluator = correction_manager.loadCorrection(
        jec_file, jer_tag + "_ScaleFactor_" + jec_algo);

    // lambda run with dataframe
    auto correction_lambda = [reapply_jes, jet_energy_scale_shifts,
                              jet_energy_scale_sf, jet_energy_resolution,
                              jer_sf_evaluator, jes_shift_sources,
                              jes_shift, jer_shift, jet_radius, 
                              era, no_jer_for_unmatched_forward_jets](
                                       const ROOT::RVec<float> &pts,
                                       const ROOT::RVec<float> &etas,
                                       const ROOT::RVec<float> &phis,
                                       const ROOT::RVec<float> &area,
                                       const ROOT::RVec<float> &raw_factors,
                                       const ROOT::RVec<UChar_t> &ids_v12,
                                       const ROOT::RVec<float> &gen_pts,
                                       const ROOT::RVec<float> &gen_etas,
                                       const ROOT::RVec<float> &gen_phis,
                                       const float &rho, const unsigned int &seed) {
        // random value generator for jet smearing
        TRandom3 randm = TRandom3(seed);
        
        auto ids = static_cast<ROOT::RVec<int>>(ids_v12);
        ROOT::RVec<float> corrected_pts;
        for (int i = 0; i < pts.size(); i++) {
            float corr_pt = pts.at(i);
            if (reapply_jes) {
                // reapplying the JES correction
                float raw_pt = pts.at(i) * (1 - raw_factors.at(i));
                float corr =
                    jet_energy_scale_sf(area.at(i), etas.at(i), raw_pt, rho, phis.at(i), era);
                corr_pt = raw_pt * corr;
                Logger::get("physicsobject::jet::PtCorrectionMC")
                    ->debug("reapplying JE scale: orig. jet pt {} to raw "
                            "jet pt {} to recorr. jet pt {}",
                            pts.at(i), raw_pt, corr_pt);
            }
            corrected_pts.push_back(corr_pt);

            // apply jet energy smearing
            float reso =
                jet_energy_resolution(etas.at(i), corrected_pts.at(i), rho);
            float reso_sf = 1.0;
            if (std::stoi(era.substr(0, 4)) <= 2018) {
                // run 2 case
                reso_sf = jer_sf_evaluator->evaluate({etas.at(i), jer_shift});
            } else {
                // run 3 case
                reso_sf = jer_sf_evaluator->evaluate({etas.at(i), corrected_pts.at(i), jer_shift});
            }
            Logger::get("physicsobject::jet::PtCorrectionMC")
                ->debug("Calculate JER {}: SF: {} resolution: {} ", jer_shift,
                        reso_sf, reso);
            // gen jet matching algorithm for JER
            ROOT::Math::RhoEtaPhiVectorF jet(corrected_pts.at(i), etas.at(i),
                                             phis.at(i));
            float gen_jet_pt = default_float;
            Logger::get("physicsobject::jet::PtCorrectionMC")
                ->debug("Going to smear jet: Eta: {} Phi: {} ", jet.Eta(),
                        jet.Phi());
            double min_delta_r = std::numeric_limits<double>::infinity();
            for (int j = 0; j < gen_pts.size(); j++) {
                ROOT::Math::RhoEtaPhiVectorF gen_jet(
                    gen_pts.at(j), gen_etas.at(j), gen_phis.at(j));
                Logger::get("physicsobject::jet::PtCorrectionMC")
                    ->debug("Checking gen Jet: Eta: {} Phi: {}", gen_jet.Eta(),
                            gen_jet.Phi());
                auto delta_r = ROOT::Math::VectorUtil::DeltaR(jet, gen_jet);
                if (delta_r > min_delta_r)
                    continue;
                if (delta_r < (jet_radius / 2.) &&
                    std::abs(corrected_pts.at(i) - gen_pts.at(j)) <
                        (3.0 * reso * corrected_pts.at(i))) {
                    min_delta_r = delta_r;
                    gen_jet_pt = gen_pts.at(j);
                }
            }
            // if jet matches a gen. jet the scaling method is applied,
            // otherwise the stochastic method
            if (gen_jet_pt > 0.0) {
                Logger::get("physicsobject::jet::PtCorrectionMC")
                    ->debug("Found gen jet for hybrid smearing method");
                double shift = (reso_sf - 1.0) *
                               (corrected_pts.at(i) - gen_jet_pt) /
                               corrected_pts.at(i);
                corrected_pts.at(i) *= std::max(0.0, 1.0 + shift);
            } else {
                Logger::get("physicsobject::jet::PtCorrectionMC")
                    ->debug("No gen jet found. Applying stochastic smearing.");
                double shift = 0.0;
                if (no_jer_for_unmatched_forward_jets && abs(etas.at(i)) > 2.5) {
                    Logger::get("physicsobject::jet::PtCorrectionMC")
                        ->debug(
                            "Jet has eta > 2.5, and no JER applied to "
                            "unmatched forward jets turned on."
                        );
                    shift = 0.0;
                } else {
                    shift = randm.Gaus(0, reso) *
                            std::sqrt(std::max(reso_sf * reso_sf - 1., 0.0));
                }
                corrected_pts.at(i) *= std::max(0.0, 1.0 + shift);
            }
            Logger::get("physicsobject::jet::PtCorrectionMC")
                ->debug("JER shift of jet pt from {} to {} ", corr_pt,
                        corrected_pts.at(i));

            float pt_scale_sf = 1.0;
            if (jes_shift != 0.0) {
                if (jes_shift_sources.at(0) != "HEMIssue") {
                    // Differentiate between single source and combined source
                    // for reduced scheme
                    if (jet_energy_scale_shifts.size() == 1) {
                        pt_scale_sf =
                            1. +
                            jes_shift * jet_energy_scale_shifts.at(0)->evaluate(
                                            {etas.at(i), corrected_pts.at(i)});
                        Logger::get("physicsobject::jet::PtCorrectionMC")
                            ->debug("JES shift of jet pt by {} for single source "
                                    "with SF {}",
                                    jes_shift, pt_scale_sf);
                    } else {
                        float quad_sum = 0.;
                        for (const auto &evaluator : jet_energy_scale_shifts) {
                            quad_sum +=
                                std::pow(evaluator->evaluate(
                                             {etas.at(i), corrected_pts.at(i)}),
                                         2.0);
                        }
                        pt_scale_sf = 1. + jes_shift * std::sqrt(quad_sum);
                        Logger::get("physicsobject::jet::PtCorrectionMC")
                            ->debug("JES shift of jet pt by {} for multiple "
                                    "sources with SF {}",
                                    jes_shift, pt_scale_sf);
                    }
                } else if (jes_shift_sources.at(0) == "HEMIssue") {
                    if (jes_shift == (-1.) && corrected_pts.at(i) > 15. &&
                        phis.at(i) > (-1.57) && phis.at(i) < (-0.87) &&
                        ids.at(i) == 2) {
                        if (etas.at(i) > (-2.5) && etas.at(i) < (-1.3))
                            pt_scale_sf = 0.8;
                        else if (etas.at(i) > (-3.) && etas.at(i) <= (-2.5))
                            pt_scale_sf = 0.65;
                    }
                }
            }
            corrected_pts.at(i) *= pt_scale_sf;
            Logger::get("physicsobject::jet::PtCorrectionMC")
                ->debug("JES shift of jet pt from {} to {} ",
                        corrected_pts.at(i) / pt_scale_sf, corrected_pts.at(i));
        }
        return corrected_pts;
    };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {jet_pt, jet_eta, jet_phi, jet_area, 
                          jet_raw_factor, jet_id_column, gen_jet_pt, 
                          gen_jet_eta, gen_jet_phi, rho, jer_seed});
    return df2;
}

/**
 * @brief This function applies jet energy corrections (JEC) to jets in real
 * data using correction factors. In nanoAOD the JEC corrections are already
 * applied to jets, however, if new corrections are measured it is possible to
 * reapply them using the newest correction files. The JEC is reapplied to data
 * if a `jes_tag` is specified.
 *
 * Unlike in Monte Carlo (MC), no smearing is applied, as resolution corrections
 * are not necessary for data.
 *
 * @warning It is not recommended to use this function because CROWN does not
 * yet have a functionality to differenciate between individual runs in eras.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * jet energy correction file
 * @param outputname name of the output column for corrected jet \f$p_T\f$'s
 * @param jet_pt name of the column containing the jet transverse momenta
 * @param jet_eta name of the column containing the jet pseudorapidities
 * @param jet_phi name of the column containing the jet azimuthal angles
 * @param jet_area name of the column containing the jet catchment area
 * @param jet_raw_factor name of the column containing the jet factors to
 * calculate back to the raw jet \f$p_T\f$'s
 * @param rho name of the column containing the event energy density
 * @param run name of the column containing the run number
 * @param jec_file path to the JEC correction file
 * @param jec_algo name of the jet reconstruction algorithm (e.g., "AK4PFchs" or
 * "AK8PFPuppi")
 * @param jes_tag tag of the JES correction campaign which is run dependent
 * (e.g., "Summer19UL18_RunA_V5_DATA")
 * @param era string defining the currently processed era, needed due to different
 * kind of recommendations from JME POG for different eras
 *
 * @return a dataframe with a new column of corrected jet \f$p_T\f$'s
 *
 * @note If jets with \f$p_T\f$ > 15 GeV are corrected, this change should be
 * propagated to the missing transverse momentum.
 * (see `physicsobject::PropagateToMET`)
 */
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correction_manager,
                 const std::string &outputname, const std::string &jet_pt,
                 const std::string &jet_eta, const std::string &jet_phi,
                 const std::string &jet_area, const std::string &jet_raw_factor,
                 const std::string &rho, const std::string &run,
                 const std::string &jec_file, const std::string &jec_algo,
                 const std::string &jes_tag, const std::string &era) {
    if (jes_tag != "") {
        // loading jet energy correction scale factor function
        auto jes_evaluator = correction_manager.loadCompoundCorrection(
            jec_file, jes_tag + "_L1L2L3Res_" + jec_algo);
        Logger::get("physicsobject::jet::PtCorrectionData")
            ->debug("file: {}, function {}", jec_file,
                    (jes_tag + "_L1L2L3Res_" + jec_algo));
        auto jet_energy_scale_sf =
            [jes_evaluator](const float area, const float eta, const float pt,
                            const float rho, const float phi, const unsigned int run,
                            const std::string era) {
                if (std::stoi(era.substr(0, 4)) <= 2022) {
                    // run 2 and 2022 cases
                    return jes_evaluator->evaluate({area, eta, pt, rho});
                } else if (era == "2023preBPix") {
                    // run 3 2023preBPix
                    return jes_evaluator->evaluate({area, eta, pt, rho, (float)run});
                } else {
                    // run 3 from 2023postBpix onwards
                    return jes_evaluator->evaluate({area, eta, pt, rho, phi, (float)run});
                }
            };

        auto correction_lambda =
            [jes_tag, jet_energy_scale_sf, era](
                const ROOT::RVec<float> &pts, const ROOT::RVec<float> &etas,
                const ROOT::RVec<float> &phis, const ROOT::RVec<float> &area,
                const ROOT::RVec<float> &raw_factors, const float &rho,
                const unsigned int &run) {
                ROOT::RVec<float> corrected_pts;
                for (int i = 0; i < pts.size(); i++) {
                    float corr_pt = pts.at(i);

                    // reapplying the JES correction
                    float raw_pt = pts.at(i) * (1 - raw_factors.at(i));
                    float corr = jet_energy_scale_sf(area.at(i), etas.at(i), raw_pt, rho, phis.at(i), run, era);

                    corr_pt = raw_pt * corr;
                    Logger::get("physicsobject::jet::PtCorrectionData")
                        ->debug("reapplying JE scale for data: orig. jet "
                                "pt {} to raw "
                                "jet pt {} to recorr. jet pt {}",
                                pts.at(i), raw_pt, corr_pt);
    
                    corrected_pts.push_back(corr_pt);
                }
                return corrected_pts;
            };
            auto df1 = df.Define(outputname, correction_lambda,
                                 {jet_pt, jet_eta, jet_phi, jet_area, jet_raw_factor, rho, run});
            return df1;
    } else {
        auto df1 = df.Define(outputname,
                             [](const ROOT::RVec<float> &pts) { return pts; },
                             {jet_pt});
        return df1;
    }
}

/** 
 * @brief This function applies a b-jet energy regression. The \f$p_T\f$ correction 
 * is applied to all b-jets identified with a b-tagging algorithm, e.g. DeepJet.
 * The goal of the regression is to better estimate the energy of b-jets because 
 * compared to other jet flavors, b-jets have a significantly higher rate of leptons
 * and therefore also neutrinos in the decay, which leads to a lower reconstructed energy.
 * The correction is determined with a neural network that was trained to simultaneously
 * estimate the b-jet energy and resolution. The application can be done on top of the
 * general jet energy scale corrections. Ref. http://cds.cern.ch/record/2690804
 *
 * @note This function should only be used for Run2 since the regression was not further
 * developed for Run3 and is also not present in the nanoAODs anymore.
 *
 * @param df input dataframe
 * @param outputname name of the output column for corrected b-jet \f$p_T\f$'s
 * @param jet_pt name of the column containing the jet \f$p_T\f$'s
 * @param scale_factor name of the column containing the scale factors for the 
 * b-jet \f$p_T\f$
 * @param bjet_mask name of the column containing the jet mask with identified 
 * b-jets
 *
 * @return a dataframe with a new column
 */
ROOT::RDF::RNode PtCorrectionBJets(ROOT::RDF::RNode df, 
                            const std::string &outputname,
                            const std::string &jet_pt, 
                            const std::string &scale_factor,
                            const std::string &bjet_mask) {
    auto bjet_pt_correction = [](const ROOT::RVec<float> &pts, 
                                 const ROOT::RVec<float> &scale_factors,
                                 const ROOT::RVec<int> &bjet_mask) {
            ROOT::RVec<float> corrected_pts;
            for (int i = 0; i < pts.size(); i++) {
                float corr_pt = pts.at(i);
                if (bjet_mask.at(i)) {
                    // applying b jet energy correction
                    corr_pt = pts.at(i) * scale_factors.at(i);
                    Logger::get("physicsobject::jet::PtCorrectionBJets")
                        ->debug("applying b jet energy correction: orig. jet "
                                "pt {} to corrected "
                                "jet pt {} with correction factor {}",
                                pts.at(i), corr_pt, scale_factors.at(i));
                }
                corrected_pts.push_back(corr_pt);
            }
            return corrected_pts;
            };

    auto df1 = df.Define(outputname, bjet_pt_correction,
                         {jet_pt, scale_factor, bjet_mask});
    return df1;
}

/**
 * @brief This function applies a jet pileup ID cut to jets. The pileup ID
 * is recommended to be applied in addition to the usual jet ID and only for
 * jets with \f$p_T\f$ < 50 GeV.
 * (Run2 ref. https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL)

 * The jet pileup ID has four possible values:
 * - 0: Fail
 * - 4: Pass loose
 * - 6: Pass loose and medium
 * - 7: Pass loose, medium, and tight
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the selection mask
 * @param jet_pu_id name of the column containing jet pileup ID values
 * @param jet_pt name of the column containing jet transverse momenta
 * @param pu_id_cut minimum pileup ID value required for a jet to pass
 * @param pt_cut minimum \f$p_T\f$ threshold for a jet to bypass the pileup ID
 cut
 *
 * @return a dataframe containing the new mask as a column
 */
ROOT::RDF::RNode CutPileupID(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &jet_pu_id,
                             const std::string &jet_pt, const int &pu_id_cut,
                             const float &pt_cut) {
    // In nanoAODv12 the type of jet PU ID was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, jet_pu_id_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, jet_pu_id+"_v12", "ROOT::VecOps::RVec<UChar_t>", jet_pu_id);

    auto pass_pu_id = [pu_id_cut, pt_cut](const ROOT::RVec<UChar_t> &pu_ids_v12,
                                          const ROOT::RVec<float> &jet_pts) {
        auto pu_ids = static_cast<ROOT::RVec<int>>(pu_ids_v12);
        ROOT::RVec<int> id_mask = pu_ids >= pu_id_cut;
        ROOT::RVec<int> pt_mask = jet_pts >= pt_cut;
        ROOT::RVec<int> mask = (id_mask + pt_mask) > 0;
        return mask;
    };
    auto df2 = df1.Define(outputname, pass_pu_id, {jet_pu_id_column, jet_pt});
    return df2;
}

/**
 * @brief This function applies a veto map to jets based on their
 * pseudorapidity and azimuthal angle. This function checks if jets
 * fall within vetoed regions defined in an external correction file
 * from JetMET POG and marks them accordingly in the output mask.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * jet veto map file
 * @param outputname name of the output column storing the veto mask
 * @param jet_eta name of the column containing jet pseudorapidities
 * @param jet_phi name of the column containing jet azimuthal angles
 * @param vetomap_file path to the veto map file
 * @param vetomap_name name of the veto map correction within the file
 * @param vetomap_type type of veto map to apply, recommented is "jetvetomap"
 *
 * @return a new dataframe containing the veto mask
 */
ROOT::RDF::RNode
ApplyVetoMap(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correction_manager,
             const std::string &outputname, const std::string &jet_eta,
             const std::string &jet_phi, const std::string &vetomap_file,
             const std::string &vetomap_name, const std::string &vetomap_type) {
    auto vetomap_evaluator =
        correction_manager.loadCorrection(vetomap_file, vetomap_name);

    auto lambda = [vetomap_evaluator,
                   vetomap_type](const ROOT::RVec<float> &etas,
                                 const ROOT::RVec<float> &phis) {
        ROOT::RVec<int> mask(etas.size(), 1.);
        bool veto = false;

        for (int i = 0; i < etas.size(); i++) {
            float phi_tmp = std::clamp(phis.at(i), (float)-3.14, (float)3.14);
            float eta_tmp = std::clamp(etas.at(i), (float)-5.1, (float)5.1);
            veto = bool(
                vetomap_evaluator->evaluate({vetomap_type, eta_tmp, phi_tmp}));

            Logger::get("physicsobject::jet::ApplyVetoMap")
                ->debug("checking object with eta {} and phi {}: should object "
                        "get vetoed? -> {}",
                        etas.at(i), phis.at(i), veto);
            mask[i] = 1 - veto;
        }
        Logger::get("physicsobject::jet::ApplyVetoMap")
            ->debug("final object veto mask {}", mask);

        return mask;
    };
    auto df1 = df.Define(outputname, lambda, {jet_eta, jet_phi});
    return df1;
}

/**
 * @brief This function checks the separation (deltaR) between each jet and
 * the given four-momentum vectors which can be for example a selected lepton
 * pair. If the deltaR is below `min_delta_r`, the jet is vetoed.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the veto mask
 * @param jet_eta name of the column containing the jet pseudorapidities
 * @param jet_phi name of the column containing the jet azimuthal angles
 * @param target_p4_1 first four-momentum vector of a target object pair
 * (e.g., a selected lepton pair)
 * @param target_p4_2 second four-momentum vector of the target object pair
 * (e.g., a selected lepton pair)
 * @param min_delta_r minimal deltaR distance allowed between jets and targets
 * to count as not overlapping
 *
 * @return a new dataframe containing the veto mask
 */
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &target_p4_1,
                    const std::string &target_p4_2, const float &min_delta_r) {
    auto df1 = df.Define(
        outputname,
        [min_delta_r](const ROOT::RVec<float> &jet_etas,
                      const ROOT::RVec<float> &jet_phis,
                      const ROOT::Math::PtEtaPhiMVector &p4_1,
                      const ROOT::Math::PtEtaPhiMVector &p4_2) {
            Logger::get("physicsobject::jet::VetoOverlappingJets")
                ->debug("Checking jets");
            ROOT::RVec<int> mask(jet_etas.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_etas.at(idx),
                                                 jet_phis.at(idx));
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("Jet: Eta: {} Phi: {}", jet.Eta(), jet.Phi());
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("Target 1 {}:  Eta: {} Phi: {}, Pt{}", p4_1,
                            p4_1.Eta(), p4_1.Phi(), p4_1.Pt());
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("Target 2 {}:  Eta: {} Phi: {}, Pt{}", p4_2,
                            p4_2.Eta(), p4_2.Phi(), p4_2.Pt());
                auto delta_r_1 = ROOT::Math::VectorUtil::DeltaR(jet, p4_1);
                auto delta_r_2 = ROOT::Math::VectorUtil::DeltaR(jet, p4_2);
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("DeltaR 1 {}", delta_r_1);
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("DeltaR 2 {}", delta_r_2);
                mask[idx] =
                    (delta_r_1 > min_delta_r && delta_r_2 > min_delta_r);
            }
            Logger::get("physicsobject::jet::VetoOverlappingJets")
                ->debug("vetomask due to overlap: {}", mask);
            return mask;
        },
        {jet_eta, jet_phi, target_p4_1, target_p4_2});
    return df1;
}

/**
 * @brief This function checks the separation (deltaR) between each jet and
 * the given four-momentum vector which can be for example a selected lepton.
 * If the deltaR is smaller than `min_delta_r`, the jet is vetoed.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the veto mask
 * @param jet_eta name of the column containing the jet pseudorapidities
 * @param jet_phi name of the column containing the jet azimuthal angles
 * @param target_p4 four-momentum vector of the target object (e.g., a
 * selected lepton)
 * @param min_delta_r minimal deltaR distance allowed between jets and target
 * to count as not overlapping
 *
 * @return a new dataframe containing the veto mask
 */
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &target_p4, const float &min_delta_r) {
    auto df1 = df.Define(
        outputname,
        [min_delta_r](const ROOT::RVec<float> &jet_etas,
                      const ROOT::RVec<float> &jet_phis,
                      const ROOT::Math::PtEtaPhiMVector &p4) {
            Logger::get("physicsobject::jet::VetoOverlappingJets")
                ->debug("Checking jets");
            ROOT::RVec<int> mask(jet_etas.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_etas.at(idx),
                                                 jet_phis.at(idx));
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("Jet: Eta: {} Phi: {}", jet.Eta(), jet.Phi());
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("Target {}:  Eta: {} Phi: {}, Pt{}", p4, p4.Eta(),
                            p4.Phi(), p4.Pt());
                auto delta_r = ROOT::Math::VectorUtil::DeltaR(jet, p4);
                Logger::get("physicsobject::jet::VetoOverlappingJets")
                    ->debug("DeltaR {}", delta_r);
                mask[idx] = (delta_r > min_delta_r);
            }
            Logger::get("physicsobject::jet::VetoOverlappingJets")
                ->debug("vetomask due to overlap: {}", mask);
            return mask;
        },
        {jet_eta, jet_phi, target_p4});
    return df1;
}

/**
 * @brief This function vetos jets that overlap with an isolated lepton
 * with a given deltaR. A jet is only vetoed if the lepton is flagged as
 * isolated (`lepton_iso` of `+1`) and deltaR is smaller than `min_delta_r`.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the veto mask
 * @param jet_eta name of the column containing the jet pseudorapidities
 * @param jet_phi name of the column containing the jet azimuthal angles
 * @param lepton_p4 four-momentum vector of the target lepton
 * @param lepton_iso name of the column name containing the lepton isolation
 * flag with values of `+1`, `-1` or `0`
 * @param min_delta_r minimal deltaR distance allowed between jets and target
 * to count as not overlapping
 *
 * @return a new dataframe containing the veto mask
 */
ROOT::RDF::RNode VetoOverlappingJetsWithIsoLepton(ROOT::RDF::RNode df,
                                                  const std::string &outputname,
                                                  const std::string &jet_eta,
                                                  const std::string &jet_phi,
                                                  const std::string &lepton_p4,
                                                  const std::string &lepton_iso,
                                                  const float &min_delta_r) {
    auto df1 = df.Define(
        outputname,
        [min_delta_r](const ROOT::RVec<float> &jet_etas,
                      const ROOT::RVec<float> &jet_phis,
                      const ROOT::Math::PtEtaPhiMVector &lep_p4,
                      const int &lep_iso) {
            Logger::get("physicsobject::jet::VetoOverlappingJetsWithIsoLepton")
                ->debug("Checking jets");
            ROOT::RVec<int> mask(jet_etas.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_etas.at(idx),
                                                 jet_phis.at(idx));
                Logger::get(
                    "physicsobject::jet::VetoOverlappingJetsWithIsoLepton")
                    ->debug("Jet:  Eta: {} Phi: {} ", jet.Eta(), jet.Phi());
                Logger::get(
                    "physicsobject::jet::VetoOverlappingJetsWithIsoLepton")
                    ->debug("Letpon {}:  Eta: {} Phi: {}, Pt{}", lep_p4,
                            lep_p4.Eta(), lep_p4.Phi(), lep_p4.Pt());
                auto delta_r = ROOT::Math::VectorUtil::DeltaR(jet, lep_p4);
                Logger::get(
                    "physicsobject::jet::VetoOverlappingJetsWithIsoLepton")
                    ->debug("DeltaR {}", delta_r);
                if (lep_iso == +1)
                    mask[idx] = (delta_r > min_delta_r);
            }
            Logger::get("physicsobject::jet::VetoOverlappingJetsWithIsoLepton")
                ->debug("vetomask due to overlap: {}", mask);
            return mask;
        },
        {jet_eta, jet_phi, lepton_p4, lepton_iso});
    return df1;
}

namespace quantity {

/**
 * @brief Applies jet identification criteria based on JSON-defined jet ID corrections.
 *
 * This function loads jet ID definitions from correctionlib JSON files for the specified
 * jet collection and evaluates two sets of criteria:
 *  - Tight ID
 *  - Tight Lepton Veto ID
 *
 * It uses these evaluations to assign a jet ID code to each jet in the input dataframe:
 *  - 6 : passes Tight and Tight Lepton Veto IDs
 *  - 2 : passes Tight ID but fails Tight Lepton Veto ID
 *  - 0 : fails Tight ID
 * (Ref. https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV#Recommendations_for_the_13_6_AN1)
 *
 * The jet ID is returned as a vector of int, compatible with NanoAOD v9 conventions.
 *
 * @param df Input ROOT RDataFrame containing jet variables
 * @param correction_manager correction manager responsible for loading the
 * correction scale uncertainty patch file
 * @param outputname Name of the new column to hold the computed jet ID flags
 * @param jet_eta Name of the branch for jet pseudorapidity
 * @param jet_chHEF Name of the branch for charged hadron energy fraction
 * @param jet_neHEF Name of the branch for neutral hadron energy fraction
 * @param jet_chEmEF Name of the branch for charged electromagnetic energy fraction
 * @param jet_neEmEF Name of the branch for neutral electromagnetic energy fraction
 * @param jet_muEF Name of the branch for muon energy fraction
 * @param jet_chMult Name of the branch for charged multiplicity
 * @param jet_neMult Name of the branch for neutral multiplicity
 * @param jet_id_file Path to the jet ID JSON file containing correction definitions
 * @param jet_name Prefix of the jet collection used to select the appropriate corrections
 *
 * @return a RDataFrame with the new jet ID column appended
 */

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
        const std::string &jet_name) {

    // Load the jet ID correction file with Tight criteria
    auto tightID = 
        correction_manager.loadCorrection(jet_id_file, jet_name + "_Tight");  
    
    // Load the jet ID correction file with TightLeptonVeto criteria
    auto tightLepVetoID = 
        correction_manager.loadCorrection(jet_id_file, jet_name + "_TightLeptonVeto"); 

    auto compute_jet_id = [tightID, tightLepVetoID](const ROOT::RVec<float> &eta,
                                                    const ROOT::RVec<float> &chHEF,
                                                    const ROOT::RVec<float> &neHEF,
                                                    const ROOT::RVec<float> &chEmEF,
                                                    const ROOT::RVec<float> &neEmEF,
                                                    const ROOT::RVec<float> &muEF,
                                                    const ROOT::RVec<UChar_t> &chMult,
                                                    const ROOT::RVec<UChar_t> &neMult) {

        size_t nJets = eta.size();
        ROOT::RVec<int> jetId(nJets); 
        for (size_t i = 0; i < nJets; ++i) {
            UChar_t mult = chMult.at(i) + neMult.at(i);
            bool passTight = false, passTightLepVeto = false;

            passTight = (tightID->evaluate(
                {eta.at(i), chHEF.at(i), neHEF.at(i), chEmEF.at(i),
                neEmEF.at(i), muEF.at(i), chMult.at(i), neMult.at(i), mult}
            ) > 0.5);

            passTightLepVeto = (tightLepVetoID->evaluate(
                {eta.at(i), chHEF.at(i), neHEF.at(i), chEmEF.at(i),
                neEmEF.at(i), muEF.at(i), chMult.at(i), neMult.at(i), mult}
            ) > 0.5);

            if (passTight && passTightLepVeto) jetId[i] = 6;
            else if (passTight && !passTightLepVeto) jetId[i] = 2;
            else jetId[i] = 0;
        }
        return jetId;
    };

    auto df1 = df.Define(outputname, compute_jet_id,
                     {jet_eta, jet_chHEF, jet_neHEF,
                      jet_chEmEF, jet_neEmEF, jet_muEF,
                      jet_chMult, jet_neMult});
    return df1;
}
} // end namespace quantity

namespace scalefactor {

/**
 * @brief This function calculates the b-tagging scale factor. The scale 
 * factor corrects inconsistencies in the b-tagging efficiency between data and
 * simulation. The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation.
 *
 * This producer is optimized for evaluating b-tagging scale factors for a 
 * full shape correction of the b-tagging discriminant (DeepJet). Working point 
 * based scale factors have different dependencies.
 *
 * More information from BTV POG can be found here https://btv-wiki.docs.cern.ch/ScaleFactors/
 *
 * and about the correctionlib files:
 * - [UL2018 b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/BTV_2018_UL_btagging.html)
 * - [UL2017 b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/BTV_2017_UL_btagging.html)
 * - [UL2016postVFP b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/BTV_2016postVFP_UL_btagging.html)
 * - [UL2016preVFP b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/BTV_2016preVFP_UL_btagging.html)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column containing the b-tagging scale factor
 * @param pt name of the column containing the transverse momenta of jets
 * @param eta name of the column containing the pseudorapidity of jets
 * @param btag_value name of the column containing the btag values of jets
 * based on a b-jet tagger (e.g. DeepJet)
 * @param flavor name of the column containing the flavors of jets, usually used
 * flavors are: 5=b-jet, 4=c-jet, 0=light jet (g, u, d, s)
 * @param jet_mask name of the column containing the mask for good/selected jets
 * @param bjet_mask name of the column containing themask for good/selected b-jets
 * @param jet_veto_mask name of the column containing the veto mask for 
 * overlapping jets (e.g. with selected lepton pairs)
 * @param sf_file path to the file with the b-tagging scale factors
 * @param sf_name name of the b-tagging scale factor correction e.g. "deepJet_shape"
 * @param variation name the scale factor variation, available values:
 * central, down_*, up_* (* name of specific variation)
 *
 * @return a new dataframe containing the new column
 *
 * @note TODO: Add Run3 support and additional producers for working point
 * based scale factors.
 */
ROOT::RDF::RNode
BtaggingShape(ROOT::RDF::RNode df,
       correctionManager::CorrectionManager &correction_manager,
       const std::string &outputname, const std::string &pt, 
       const std::string &eta, const std::string &btag_value, 
       const std::string &flavor, const std::string &jet_mask, 
       const std::string &bjet_mask, const std::string &jet_veto_mask, 
       const std::string &sf_file, const std::string &sf_name,
       const std::string &variation) {
    Logger::get("physicsobject::jet::scalefactor::BtaggingShape")->debug(
        "Setting up functions for b-tag sf with correctionlib");
    Logger::get("physicsobject::jet::scalefactor::BtaggingShape")->debug("Correction algorithm - Name {}",
                                 sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    
    // In nanoAODv12 the type of jet flavor was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, flavor_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, flavor+"_v12", "ROOT::VecOps::RVec<UChar_t>", flavor);

    auto btagSF_lambda = [evaluator,
                          variation](const ROOT::RVec<float> &pts,
                                     const ROOT::RVec<float> &etas,
                                     const ROOT::RVec<float> &btag_values,
                                     const ROOT::RVec<UChar_t> &flavors_v12,
                                     const ROOT::RVec<int> &jet_mask,
                                     const ROOT::RVec<int> &bjet_mask,
                                     const ROOT::RVec<int> &jet_veto_mask) {
        auto flavors = static_cast<ROOT::RVec<int>>(flavors_v12);
        Logger::get("physicsobject::jet::scalefactor::BtaggingShape")->debug("Variation - Name {}", variation);
        float sf = 1.;
        for (int i = 0; i < pts.size(); i++) {
            Logger::get("physicsobject::jet::scalefactor::BtaggingShape")->debug(
                "jet masks - jet {}, b-jet {}, jet veto {}", jet_mask.at(i),
                bjet_mask.at(i), jet_veto_mask.at(i));
            // considering only good jets/b-jets, this is needed since jets and
            // bjets might have different quality cuts depending on the analysis
            if ((jet_mask.at(i) || bjet_mask.at(i)) && jet_veto_mask.at(i)) {
                Logger::get("physicsobject::jet::scalefactor::BtaggingShape")->debug(
                    "SF - pt {}, eta {}, btag value {}, flavor {}",
                    pts.at(i), etas.at(i), btag_values.at(i),
                    flavors.at(i));
                float jet_sf = 1.;
                // considering only phase space where the scale factors are
                // defined
                if (
                    pts.at(i) >= 20.0
                    && pts.at(i) < 10000.0
                    && std::abs(etas.at(i)) < 2.5
                    && !(std::isnan(btag_values.at(i)) || btag_values.at(i) == -1.0)
                ) {
                    // for c-jet related uncertainties only scale factors of
                    // c-jets are varied, the rest is nominal/central
                    if (variation.find("cferr") != std::string::npos) {
                        // flavor=4 means c-flavor
                        if (flavors.at(i) == 4) {
                            jet_sf = evaluator->evaluate(
                                {variation, flavors.at(i),
                                 std::abs(etas.at(i)), pts.at(i),
                                 btag_values.at(i)});
                        } else {
                            jet_sf = evaluator->evaluate(
                                {"central", flavors.at(i),
                                 std::abs(etas.at(i)), pts.at(i),
                                 btag_values.at(i)});
                        }
                    }
                    // for nominal/central and all other uncertainties c-jets
                    // have a scale factor of 1 (only for central defined in
                    // json file from BTV)
                    else {
                        if (flavors.at(i) != 4) {
                            jet_sf = evaluator->evaluate(
                                {variation, flavors.at(i),
                                 std::abs(etas.at(i)), pts.at(i),
                                 btag_values.at(i)});
                        } else {
                            jet_sf = evaluator->evaluate(
                                {"central", flavors.at(i),
                                 std::abs(etas.at(i)), pts.at(i),
                                 btag_values.at(i)});
                        }
                    }
                }
                Logger::get("physicsobject::jet::scalefactor::BtaggingShape")->debug("Jet Scale Factor {}", jet_sf);
                sf *= jet_sf;
            }
        };
        Logger::get("physicsobject::jet::scalefactor::BtaggingShape")->debug("Event Scale Factor {}", sf);
        return sf;
    };
    auto df2 = df1.Define(
        outputname, btagSF_lambda,
        {pt, eta, btag_value, flavor_column, jet_mask, bjet_mask, jet_veto_mask});
    return df2;
}

/**
 * @brief This function calculates the b-tagging scale factor. The scale 
 * factor corrects inconsistencies in the b-tagging efficiency between data and
 * simulation. The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation.
 *
 * This producer can be used to evaluate working point based scale factors. It is
 * defined based on scale factors provided by BTV POG for Run3 2024.
 *
 * More information from BTV POG can be found here https://btv-wiki.docs.cern.ch/ScaleFactors/
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column containing the b-tagging scale factor
 * @param pt name of the column containing the transverse momenta of jets
 * @param eta name of the column containing the pseudorapidity of jets
 * @param flavor name of the column containing the flavors of jets, usually used
 * flavors are: 5=b-jet, 4=c-jet, 0=light jet (g, u, d, s)
 * @param jet_mask name of the column containing the mask for good/selected jets
 * @param bjet_mask name of the column containing the mask for good/selected b-jets
 * @param jet_veto_mask name of the column containing the veto mask for 
 * overlapping jets (e.g. with selected lepton pairs)
 * @param sf_file path to the file with the b-tagging scale factors
 * @param sf_name name of the b-tagging scale factor correction e.g. "deepJet_shape"
 * @param variation name the scale factor variation, available values:
 * central, down_*, up_* (* name of specific variation)
 * @param btag_wp string that specifies the b-tagging working point used in an
 * analysis e.g. "L", "M", "T", ...
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
BtaggingWP(ROOT::RDF::RNode df,
       correctionManager::CorrectionManager &correction_manager,
       const std::string &outputname, const std::string &pt, 
       const std::string &eta, const std::string &flavor,
       const std::string &jet_mask, 
       const std::string &bjet_mask, const std::string &jet_veto_mask, 
       const std::string &sf_file, const std::string &sf_name,
       const std::string &variation, const std::string &btag_wp) {
    Logger::get("physicsobject::jet::scalefactor::BtaggingWP")->debug(
        "Setting up functions for wp based b-tag sf with correctionlib");
    Logger::get("physicsobject::jet::scalefactor::BtaggingWP")->debug("Correction algorithm - Name {}",
                                 sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);

    // In nanoAODv12 the type of jet flavor was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, flavor_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, flavor+"_v12", "ROOT::VecOps::RVec<UChar_t>", flavor);

    auto btagSF_lambda = [evaluator,
                          variation, btag_wp](
                                    const ROOT::RVec<float> &etas,
                                    const ROOT::RVec<float> &pts,
                                    const ROOT::RVec<UChar_t> &flavors_v12,
                                    const ROOT::RVec<int> &jet_mask,
                                    const ROOT::RVec<int> &bjet_mask,
                                    const ROOT::RVec<int> &jet_veto_mask) {
        auto flavors = static_cast<ROOT::RVec<int>>(flavors_v12);
        Logger::get("physicsobject::jet::scalefactor::BtaggingWP")->debug("Variation - Name {}", variation);
        float sf = 1.;
        for (int i = 0; i < pts.size(); i++) {
            Logger::get("physicsobject::jet::scalefactor::BtaggingWP")->debug(
                "jet masks - jet {}, b-jet {}, jet veto {}", jet_mask.at(i),
                bjet_mask.at(i), jet_veto_mask.at(i));
            // considering only good jets/b-jets, this is needed since jets and
            // bjets might have different quality cuts depending on the analysis
            if ((jet_mask.at(i) || bjet_mask.at(i)) && jet_veto_mask.at(i)) {
                Logger::get("physicsobject::jet::scalefactor::BtaggingWP")->debug(
                    "SF - pt {}, eta {}, btag wp {}, flavor {}",
                    pts.at(i), etas.at(i), btag_wp,
                    flavors.at(i));
                float jet_sf = 1.;
                // considering only phase space where the scale factors are
                // defined
                if (pts.at(i) >= 20.0 && pts.at(i) < 10000.0 &&
                    std::abs(etas.at(i)) < 2.5) {
                        jet_sf = evaluator->evaluate(
                                {variation, btag_wp, flavors.at(i),
                                 std::abs(etas.at(i)), pts.at(i)});
                }
                Logger::get("physicsobject::jet::scalefactor::BtaggingWP")->debug("Jet Scale Factor {}", jet_sf);
                sf *= jet_sf;
            }
        };
        Logger::get("physicsobject::jet::scalefactor::BtaggingWP")->debug("Event Scale Factor {}", sf);
        return sf;
    };
    auto df2 = df1.Define(
        outputname, btagSF_lambda,
        {pt, eta, flavor_column, jet_mask, bjet_mask, jet_veto_mask});
    return df2;
}

} // end namespace scalefactor
} // end namespace jet
} // end namespace physicsobject
