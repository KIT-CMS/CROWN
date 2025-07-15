#ifndef GUARD_JETS_H
#define GUARD_JETS_H

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
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
 * @param jer_seed seed value for the random number generator that is used for 
 * the jet energy resolution smearing, if not set the answer to everything is 
 * used as default `42`
 *
 * @return a dataframe with a new column of corrected jet \f$p_T\f$'s
 *
 * @note If jets with \f$p_T\f$ > 15 GeV are corrected, this change should be
 * propagated to the missing transverse momentum.
 */
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
               const int jer_seed) {
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
    auto jet_energy_scale_sf = [jes_evaluator](const float area,
                                               const float eta, const float pt,
                                               const float rho) {
        return jes_evaluator->evaluate({area, eta, pt, rho});
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
    auto jet_energy_resolution_sf =
        [jer_sf_evaluator](const float eta, const std::string jer_shift) {
            return jer_sf_evaluator->evaluate({eta, jer_shift});
        };
    // lambda run with dataframe
    auto correction_lambda = [reapply_jes, jet_energy_scale_shifts,
                              jet_energy_scale_sf, jet_energy_resolution,
                              jet_energy_resolution_sf, jes_shift_sources,
                              jes_shift, jer_shift, jer_seed,
                              jet_radius](const ROOT::RVec<float> &pts,
                                          const ROOT::RVec<float> &etas,
                                          const ROOT::RVec<float> &phis,
                                          const ROOT::RVec<float> &area,
                                          const ROOT::RVec<float> &raw_factors,
                                          const ROOT::RVec<int> &ids,
                                          const ROOT::RVec<float> &gen_pts,
                                          const ROOT::RVec<float> &gen_etas,
                                          const ROOT::RVec<float> &gen_phis,
                                          const float &rho) {
        // random value generator for jet smearing
        TRandom3 randm = TRandom3(jer_seed);

        ROOT::RVec<float> corrected_pts;
        for (int i = 0; i < pts.size(); i++) {
            float corr_pt = pts.at(i);
            if (reapply_jes) {
                // reapplying the JES correction
                float raw_pt = pts.at(i) * (1 - raw_factors.at(i));
                float corr =
                    jet_energy_scale_sf(area.at(i), etas.at(i), raw_pt, rho);
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
            float reso_sf = jet_energy_resolution_sf(etas.at(i), jer_shift);
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
                double shift = randm.Gaus(0, reso) *
                               std::sqrt(std::max(reso_sf * reso_sf - 1., 0.0));
                corrected_pts.at(i) *= std::max(0.0, 1.0 + shift);
            }
            Logger::get("physicsobject::jet::PtCorrectionMC")
                ->debug("Shifting jet pt from {} to {} ", corr_pt,
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
                            ->debug("Shifting jet pt by {} for single source "
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
                            ->debug("Shifting jet pt by {} for multiple "
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
                ->debug("Shifting jet pt from {} to {} ",
                        corrected_pts.at(i) / pt_scale_sf, corrected_pts.at(i));
        }
        return corrected_pts;
    };
    auto df1 = df.Define(outputname, correction_lambda,
                         {jet_pt, jet_eta, jet_phi, jet_area, jet_raw_factor,
                          jet_id, gen_jet_pt, gen_jet_eta, gen_jet_phi, rho});
    return df1;
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
 * @param jet_area name of the column containing the jet catchment area
 * @param jet_raw_factor name of the column containing the jet factors to
 * calculate back to the raw jet \f$p_T\f$'s
 * @param rho name of the column containing the event energy density
 * @param jec_file path to the JEC correction file
 * @param jec_algo name of the jet reconstruction algorithm (e.g., "AK4PFchs" or
 * "AK8PFPuppi")
 * @param jes_tag tag of the JES correction campaign which is run dependent
 * (e.g., "Summer19UL18_RunA_V5_DATA")
 *
 * @return a dataframe with a new column of corrected jet \f$p_T\f$'s
 *
 * @note If jets with \f$p_T\f$ > 15 GeV are corrected, this change should be
 * propagated to the missing transverse momentum.
 */
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correction_manager,
                 const std::string &outputname, const std::string &jet_pt,
                 const std::string &jet_eta, const std::string &jet_area,
                 const std::string &jet_raw_factor, const std::string &rho,
                 const std::string &jec_file, const std::string &jec_algo,
                 const std::string &jes_tag) {
    if (jes_tag != "") {
        // loading jet energy correction scale factor function
        auto jes_evaluator = correction_manager.loadCompoundCorrection(
            jec_file, jes_tag + "_L1L2L3Res_" + jec_algo);
        Logger::get("physicsobject::jet::PtCorrectionData")
            ->debug("file: {}, function {}", jec_file,
                    (jes_tag + "_L1L2L3Res_" + jec_algo));
        auto jet_energy_scale_sf =
            [jes_evaluator](const float area, const float eta, const float pt,
                            const float rho) {
                return jes_evaluator->evaluate({area, eta, pt, rho});
            };

        auto correction_lambda =
            [jes_tag, jet_energy_scale_sf](
                const ROOT::RVec<float> &pts, const ROOT::RVec<float> &etas,
                const ROOT::RVec<float> &area,
                const ROOT::RVec<float> &raw_factors, const float &rho) {
                ROOT::RVec<float> corrected_pts;
                for (int i = 0; i < pts.size(); i++) {
                    float corr_pt = pts.at(i);
                    if (jes_tag != "") {
                        // reapplying the JES correction
                        float raw_pt = pts.at(i) * (1 - raw_factors.at(i));
                        float corr = jet_energy_scale_sf(area.at(i), etas.at(i),
                                                         raw_pt, rho);
                        corr_pt = raw_pt * corr;
                        Logger::get("physicsobject::jet::PtCorrectionData")
                            ->debug("reapplying JE scale for data: orig. jet "
                                    "pt {} to raw "
                                    "jet pt {} to recorr. jet pt {}",
                                    pts.at(i), raw_pt, corr_pt);
                    }
                    corrected_pts.push_back(corr_pt);
                }
                return corrected_pts;
            };
        auto df1 = df.Define(outputname, correction_lambda,
                             {jet_pt, jet_eta, jet_area, jet_raw_factor, rho});
        return df1;
    } else {
        auto df1 = df.Define(outputname,
                             [](const ROOT::RVec<float> &pts) { return pts; },
                             {jet_pt});
        return df1;
    }
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
    auto jet_pu_id_int = jet_pu_id;
    auto df_int = df;
    // Recasting to int if it is UChar_t
    if (df.GetColumnType(jet_pu_id) == "ROOT::VecOps::RVec<UChar_t>") {
        jet_pu_id_int = jet_pu_id + "_int";
        df_int = df_int.Define(jet_pu_id_int, [](const ROOT::RVec<UChar_t>& v) {
            return ROOT::RVec<int>(v.begin(), v.end());
        }, {jet_pu_id});
    }

    auto pass_pu_id = [pu_id_cut, pt_cut](const ROOT::RVec<int> &pu_ids,
                                          const ROOT::RVec<float> &jet_pts) {
        ROOT::RVec<int> id_mask = pu_ids >= pu_id_cut;
        ROOT::RVec<int> pt_mask = jet_pts >= pt_cut;
        ROOT::RVec<int> mask = (id_mask + pt_mask) > 0;
        return mask;
    };
    auto df1 = df.Define(outputname, pass_pu_id, {jet_pu_id_int, jet_pt});
    return df1;
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
 * @param btag_value name of the column containing the btag values of jetd 
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
Btagging(ROOT::RDF::RNode df,
       correctionManager::CorrectionManager &correction_manager,
       const std::string &outputname, const std::string &pt, 
       const std::string &eta, const std::string &btag_value, 
       const std::string &flavor, const std::string &jet_mask, 
       const std::string &bjet_mask, const std::string &jet_veto_mask, 
       const std::string &sf_file, const std::string &sf_name,
       const std::string &variation) {
    Logger::get("physicsobject::jet::scalefactor::Btagging")->debug(
        "Setting up functions for b-tag sf with correctionlib");
    Logger::get("physicsobject::jet::scalefactor::Btagging")->debug("Correction algorithm - Name {}",
                                 sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);

    // In nanoAODv12 the type of jet flavor was changed to UChar_t
    auto flavor_int = flavor;
    auto df_int = df;
    // Recasting to int if it is UChar_t
    if (df.GetColumnType(flavor) == "ROOT::VecOps::RVec<UChar_t>") {
        flavor_int = flavor + "_int";
        df_int = df_int.Define(flavor_int, [](const ROOT::RVec<UChar_t>& v) {
            return ROOT::RVec<int>(v.begin(), v.end());
        }, {flavor});
    }

    auto btagSF_lambda = [evaluator,
                          variation](const ROOT::RVec<float> &pts,
                                     const ROOT::RVec<float> &etas,
                                     const ROOT::RVec<float> &btag_values,
                                     const ROOT::RVec<int> &flavors,
                                     const ROOT::RVec<int> &jet_mask,
                                     const ROOT::RVec<int> &bjet_mask,
                                     const ROOT::RVec<int> &jet_veto_mask) {
        Logger::get("physicsobject::jet::scalefactor::Btagging")->debug("Vatiation - Name {}", variation);
        float sf = 1.;
        for (int i = 0; i < pts.size(); i++) {
            Logger::get("physicsobject::jet::scalefactor::Btagging")->debug(
                "jet masks - jet {}, b-jet {}, jet veto {}", jet_mask.at(i),
                bjet_mask.at(i), jet_veto_mask.at(i));
            // considering only good jets/b-jets, this is needed since jets and
            // bjets might have different quality cuts depending on the analysis
            if ((jet_mask.at(i) || bjet_mask.at(i)) && jet_veto_mask.at(i)) {
                Logger::get("physicsobject::jet::scalefactor::Btagging")->debug(
                    "SF - pt {}, eta {}, btag value {}, flavor {}",
                    pts.at(i), etas.at(i), btag_values.at(i),
                    flavors.at(i));
                float jet_sf = 1.;
                // considering only phase space where the scale factors are
                // defined
                if (pts.at(i) >= 20.0 && pts.at(i) < 10000.0 &&
                    std::abs(etas.at(i)) < 2.5) {
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
                Logger::get("physicsobject::jet::scalefactor::Btagging")->debug("Jet Scale Factor {}", jet_sf);
                sf *= jet_sf;
            }
        };
        Logger::get("physicsobject::jet::scalefactor::Btagging")->debug("Event Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df_int.Define(
        outputname, btagSF_lambda,
        {pt, eta, btag_value, flavor_int, jet_mask, bjet_mask, jet_veto_mask});
    return df1;
}
} // end namespace scalefactor
} // end namespace jet
} // end namespace physicsobject
#endif /* GUARD_JETS_H */
