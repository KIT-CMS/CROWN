#ifndef GUARD_TAUS_H
#define GUARD_TAUS_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "../include/defaults.hxx"
#include "ROOT/RDataFrame.hxx"
#include "correction.h"
#include <stdexcept>
#include <unordered_map>

namespace physicsobject {
namespace tau {

/**
 * @brief Helper to get the variation of the tau pt scale correction to use in
 * the correctionlib evaluator, depending on the absolute pseudorapidity, decay
 * mode, and gen match of the tau. If no criterion for any of the variations is
 * matched, the nominal shift "nom" is returned.
 * 
 * @param abs_eta absolute pseudorapidity of the tau
 * @param decay_mode decay mode of the tau
 * @param gen_match gen match of the tau
 * @param variation_efake_dm0_barrel variation for electron faking a tau with
 * decay mode 0 in the barrel region
 * @param variation_efake_dm1_barrel variation for electron faking a tau with
 * decay mode 1 in the barrel region
 * @param variation_efake_dm0_endcap variation for electron faking a tau with
 * decay mode 0 in the endcap region
 * @param variation_efake_dm1_endcap variation for electron faking a tau with
 * decay mode 1 in the endcap region
 * @param variation_mufake variation for muon faking a tau
 * @param variation_gentau_dm0 variation for genuine tau with decay mode 0
 * @param variation_gentau_dm1 variation for genuine tau with decay mode 1
 * @param variation_gentau_dm10 variation for genuine tau with decay mode 10
 * @param variation_gentau_dm11 variation for genuine tau with decay mode 11
 * 
 * @return the variation to use in the correctionlib evaluator
 */
std::string get_tes_variation(
    const float &abs_eta,
    const int &decay_mode,
    const int &gen_match,
    const std::string &variation_efake_dm0_barrel,
    const std::string &variation_efake_dm1_barrel,
    const std::string &variation_efake_dm0_endcap,
    const std::string &variation_efake_dm1_endcap,
    const std::string &variation_mufake,
    const std::string &variation_gentau_dm0,
    const std::string &variation_gentau_dm1,
    const std::string &variation_gentau_dm10,
    const std::string &variation_gentau_dm11
) {
    // values for pseudorapidity cuts to define barrel and endcap regions
    const float barrel_end_cut = 1.5;
    const float endcap_end_cut = 2.5;

    // set the variation depending on the gen match, decay mode, and eta
    // for uncovered cases, "nom" is returned
    std::string variation = "nom";
    if ((gen_match == 1) || (gen_match == 3)) {

        // energy scale correction for an electron faking a tau
        if (decay_mode == 0) {
            if (abs_eta <= barrel_end_cut) {
                variation = variation_efake_dm0_barrel;
            } else if (abs_eta > barrel_end_cut and abs_eta <= endcap_end_cut) {
                variation = variation_efake_dm0_endcap;
            }
        } else if (decay_mode == 1) {
            if (abs_eta <= barrel_end_cut) {
                variation = variation_efake_dm1_barrel;
            } else if (abs_eta > barrel_end_cut and abs_eta <= endcap_end_cut) {
                variation = variation_efake_dm1_endcap;
            }
        }

    } else if ((gen_match == 2) || (gen_match == 4)) {

        // energy scale correction for a muon faking a tau
        variation = variation_mufake;

    } else if (gen_match == 5) {

        // energy scale correction for a genuine tau
        if (decay_mode == 0) {
            variation = variation_gentau_dm0;
        } else if (decay_mode == 1) {
            variation = variation_gentau_dm1;
        } else if (decay_mode == 10) {
            variation = variation_gentau_dm10;
        } else if (decay_mode == 11) {
            variation = variation_gentau_dm11;
        }

    }

    // debug information
    Logger::get("physicsobject::tau::get_tes_variation")
        ->debug(
            "variation for tau abs eta {}, decaymode {}, gen match {}: "
            "variation {}",
            abs_eta,
            decay_mode,
            gen_match,
            variation
        );

    return variation;
}

/**
 * @brief This function applies a transverse momentum (\f$p_T\f$) correction to
 * hadronic taus in MC simulations.
 * 
 * The correction depends on the physical origin of the tau (electron fake, muon
 * fake, or genuine tau), the decay mode, the \f$p_T\f$, and the pseudorapidity.
 * 
 * For Run 3 analyses, the corrections are calculated for different working
 * points of the `DeepTau` algorithm, regarding the identification against
 * jets and against electrons. This is not the case for Run 2 analyses. This
 * function can be used for both Run 2 and Run 3 analyses. For Run 2 analyses,
 * the values of `id_vs_jet_wp` and `id_vs_ele_wp` can be set to `""`
 * to obtain the corrections.
 * 
 * The uncertainty scheme is split into nine different uncertainty sources:
 * - For electrons faking taus, four uncertainty sources are considered, split
 *   by the decay modes 0 and 1, and by whether the tau is found in the barrel
 *   or the endcap region.
 * - For muons faking taus, one inclusive uncertainty source is considered.
 * - For genuine taus, four uncertainty sources are considered, split by the
 *   decay modes 0, 1, 10, and 11.
 * For each source, the variation can be set individually. The variations can
 * take the values ``nom``, ``up``, or ``down``.
 * 
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * The implementation of this function is based on the TAU POG
 * [recommendations for Run 2](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2)
 * and [recommendations for Run 3](https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3).
 * 
 * The specification of the correctionlib files used to evaluate the corrections can be found here:
 * - [2016preVFP](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html)
 * - [2016postVFP](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html)
 * - [2017](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html)
 * - [2018](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html)
 * - [2022preEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2022_Summer22_tau_DeepTau2018v2p5_2022_preEE.html)
 * - [2022postEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2022_Summer22EE_tau_DeepTau2018v2p5_2022_postEE.html)
 * - [2023preBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2023_Summer23_tau_DeepTau2018v2p5_2023_preBPix.html)
 * - [2023postBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2023_Summer23BPix_tau_DeepTau2018v2p5_2023_postBPix.html)
 * 
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param id_vs_jet_wp working point for the identification against jets; set to `""` if the corrections do not depend on this parameter
 * @param id_vs_ele_wp working point for the identification against electrons; set to `""` if the corrections do not depend on this parameter
 * @param variation_efake_dm0_barrel variation for electron faking a tau for decay mode 0 in the barrel region,
 * options are "nom", "up", "down"
 * @param variation_efake_dm1_barrel variation for electron faking a tau for decay mode 1 in the barrel region,
 * options are "nom", "up", "down"
 * @param variation_efake_dm0_endcap variation for electron faking a tau for decay mode 0 in the endcap regionefake_,
 * options are "nom", "up", "down"
 * @param variation_efake_dm1_endcap variation for electron faking a tau for decay mode 1 in the endcap region,
 * options are "nom", "up", "down"
 * @param variation_mufake variation for muon faking a tau, options are "nom", "up", "down"
 * @param variation_gentau_dm0 variation for genuine tau for decay mode 0, options are "nom", "up",
 * "down"
 * @param variation_gentau_dm1 variation for genuine tau for decay mode 1, options are "nom", "up",
 * "down"
 * @param variation_gentau_dm10 variation for genuine tau for decay mode 10, options are "nom", "up",
 * "down"
 * @param variation_gentau_dm11 variation for genuine tau for decay mode 11, options are "nom", "up",
 * "down"
 * 
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This function is intended to be used for Run 3 analyses. In Run 3,
 * the tau energy scale corrections also depend on the DeepTau working points
 * for ID vs. electrons and vs. jets. An overloaded version of this function
 * exists for this purpose.
 */
ROOT::RDF::RNode
PtCorrectionMC(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname,
    const std::string &pt,
    const std::string &eta,
    const std::string &decay_mode,
    const std::string &gen_match,
    const std::string &es_file,
    const std::string &correction_name,
    const std::string &id_algorithm,
    const std::string &variation_efake_dm0_barrel,
    const std::string &variation_efake_dm1_barrel,
    const std::string &variation_efake_dm0_endcap,
    const std::string &variation_efake_dm1_endcap,
    const std::string &variation_mufake,
    const std::string &variation_gentau_dm0,
    const std::string &variation_gentau_dm1,
    const std::string &variation_gentau_dm10,
    const std::string &variation_gentau_dm11,
    const std::string &id_vs_jet_wp = "",
    const std::string &id_vs_ele_wp = ""
) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::VecOps::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);

    auto correction_lambda =
        [
            evaluator,
            id_algorithm,
            id_vs_jet_wp,
            id_vs_ele_wp,
            variation_efake_dm0_barrel,
            variation_efake_dm1_barrel,
            variation_efake_dm0_endcap,
            variation_efake_dm1_endcap,
            variation_mufake,
            variation_gentau_dm0,
            variation_gentau_dm1,
            variation_gentau_dm10,
            variation_gentau_dm11
        ] (
            const ROOT::RVec<float> &pts,
            const ROOT::RVec<float> &etas,
            const ROOT::RVec<UChar_t> &decay_modes_v12,
            const ROOT::RVec<UChar_t> &gen_matches_char
        ) {
            // convert decay modes and gen matches to integers
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);

            // container for corrected pts
            ROOT::RVec<float> corrected_pts(pts.size());

            for (int i = 0; i < pts.size(); i++) {
                // get tau variables that we need for scale factor evaluation
                auto pt = pts.at(i);
                auto abs_eta = etas.at(i);
                auto decay_mode = decay_modes.at(i);
                auto gen_match = gen_matches.at(i);

                // set the variation depending on the gen match, decay mode, and barrel/endcap region
                std::string variation = get_tes_variation(
                    abs_eta,
                    decay_mode,
                    gen_match,
                    variation_efake_dm0_barrel,
                    variation_efake_dm1_barrel,
                    variation_efake_dm0_endcap,
                    variation_efake_dm1_endcap,
                    variation_mufake,
                    variation_gentau_dm0,
                    variation_gentau_dm1,
                    variation_gentau_dm10,
                    variation_gentau_dm11
                );

                // evaluate the correction factor
                // ensure that the tau fulfills the selection criteria for application of the correction,
                // set the correction factor to 1 otherwise
                float correction_factor = 1.0;
                const std::unordered_set<int> valid_modes = {0, 1, 10, 11};
                if (valid_modes.count(decay_mode)) {
                    if ((id_vs_jet_wp == "") && (id_vs_ele_wp == "")) {
                        correction_factor = evaluator->evaluate(
                            {
                                pt,
                                abs_eta,
                                decay_mode,
                                gen_match,
                                id_algorithm,
                                variation
                            }
                        );
                    } else {
                        correction_factor = evaluator->evaluate(
                            {
                                pt,
                                abs_eta,
                                decay_mode,
                                gen_match,
                                id_algorithm,
                                id_vs_jet_wp,
                                id_vs_ele_wp,
                                variation
                            }
                        );
                    }
                } else {
                    correction_factor = 1.0;
                }

                // calculate the corrected pt
                corrected_pts[i] = pt * correction_factor;

                // debug information
                Logger::get("physicsobject::tau::PtCorrectionMC")
                    ->debug(
                        "apply tau pt correction to tau pt {}, decaymode {}, "
                        "gen match {}, variation {} --> corrected pt {}, "
                        "correction factor {}",
                        pt,
                        decay_mode,
                        gen_match,
                        variation,
                        corrected_pts.at(i),
                        correction_factor
                    );
                }

           return corrected_pts;
        };

    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}

/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of hadronic taus that originate from electrons that are
 * misidentified. The energy scale correction for these objects is measured for
 * two tau decay modes (dm0 and dm1) and depends on the barrel and endcap region
 * of the detector. This configuration corresponds to the officially provided
 * scale factors in Run 2.
 * 
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param variation_dm0_barrel variation for decay mode 0 in the barrel region,
 * options are "nom", "up", "down"
 * @param variation_dm1_barrel variation for decay mode 1 in the barrel region,
 * options are "nom", "up", "down"
 * @param variation_dm0_endcap variation for decay mode 0 in the endcap region,
 * options are "nom", "up", "down"
 * @param variation_dm1_endcap variation for decay mode 1 in the endcap region,
 * options are "nom", "up", "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to misidentified hadronic taus
 * originating from prompt electrons (`gen_match=1`) and electrons that decayed
 * from a tau lepton
 * (`gem_match=3`).
 * 
 * @note This function is intended to be used for Run 2 analyses. In Run 3,
 * the tau energy scale corrections also depend on the DeepTau working points
 * for ID vs. electrons and vs. jets. Use
 * `physicsobject::tau::TauCorrectionMC` instead.
 */
ROOT::RDF::RNode
PtCorrectionMC_eleFake(ROOT::RDF::RNode df,
                       correctionManager::CorrectionManager &correction_manager,
                       const std::string &outputname, const std::string &pt,
                       const std::string &eta, const std::string &decay_mode,
                       const std::string &gen_match, const std::string &es_file,
                       const std::string &correction_name,
                       const std::string &id_algorithm,
                       const std::string &variation_dm0_barrel,
                       const std::string &variation_dm1_barrel,
                       const std::string &variation_dm0_endcap,
                       const std::string &variation_dm1_endcap) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::VecOps::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda =
        [evaluator, id_algorithm, variation_dm0_barrel, variation_dm1_barrel,
         variation_dm0_endcap, variation_dm1_endcap](
            const ROOT::RVec<float> &pts, const ROOT::RVec<float> &etas,
            const ROOT::RVec<UChar_t> &decay_modes_v12,
            const ROOT::RVec<UChar_t> &gen_matches_char) {
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
            ROOT::RVec<float> corrected_pts(pts.size());
            const float barrel_end_cut = 1.5;
            const float endcap_end_cut = 2.5;
            const int dm0 = 0;
            const int dm1 = 1;

            for (int i = 0; i < pts.size(); i++) {
                if (gen_matches.at(i) == 1 || gen_matches.at(i) == 3) {
                    if (decay_modes.at(i) == dm0 &&
                        std::abs(etas.at(i)) <= barrel_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm0_barrel});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == dm0 &&
                               std::abs(etas.at(i)) > barrel_end_cut &&
                               std::abs(etas.at(i)) <= endcap_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm0_endcap});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == dm1 &&
                               std::abs(etas.at(i)) <= barrel_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm1_barrel});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == dm1 &&
                               std::abs(etas.at(i)) > barrel_end_cut &&
                               std::abs(etas.at(i)) <= endcap_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm1_endcap});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    }
                } else {
                    corrected_pts[i] = pts.at(i);
                }
                Logger::get("physicsobject::tau::PtCorrectionMC_eleFake")
                    ->debug("tau pt before {}, tau pt after {}", pts.at(i),
                            corrected_pts.at(i));
            }
            return corrected_pts;
        };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}

/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of hadronic taus that originate from muons that are
 * misidentified. The energy scale correction is always set to `1` and for the
 * uncertainty varied by `1%`.
 *
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param variation variation of the correction, options are "nom", "up", "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to misidentified hadronic taus
 * originating from prompt muons (`gen_match=2`) and muons that decayed from a
 * tau lepton (`gen_match=4`).
 *
 * @note This function is intended to be used for Run 2 analyses. In Run 3,
 * the tau energy scale corrections also depend on the DeepTau working points
 * for ID vs. electrons and vs. jets. Use 
 * `physicsobject::tau::TauCorrectionMC` instead.
 */
ROOT::RDF::RNode
PtCorrectionMC_muFake(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correction_manager,
                      const std::string &outputname, const std::string &pt,
                      const std::string &eta, const std::string &decay_mode,
                      const std::string &gen_match, const std::string &es_file,
                      const std::string &correction_name,
                      const std::string &id_algorithm,
                      const std::string &variation) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::VecOps::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda = [evaluator, id_algorithm, variation](
                                 const ROOT::RVec<float> &pts,
                                 const ROOT::RVec<float> &etas,
                                 const ROOT::RVec<UChar_t> &decay_modes_v12,
                                 const ROOT::RVec<UChar_t> &gen_matches_char) {
        auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
        auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
        ROOT::RVec<float> corrected_pts(pts.size());
        for (int i = 0; i < pts.size(); i++) {
            if (gen_matches.at(i) == 2 || gen_matches.at(i) == 4) {
                auto correction_factor = evaluator->evaluate(
                    {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                     gen_matches.at(i), id_algorithm,
                     variation});
                corrected_pts[i] = pts.at(i) * correction_factor;
            } else {
                corrected_pts[i] = pts.at(i);
            }
            Logger::get("physicsobject::tau::PtCorrectionMC_muFake")
                ->debug("tau pt before {}, tau pt after {}", pts.at(i),
                        corrected_pts.at(i));
        }
        return corrected_pts;
    };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}

/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of genuine hadronic taus. The energy scale correction for these
 * objects is measured for four tau decay modes (dm0, dm1, dm10 and dm11) and
 * depends on the transverse momentum of the hadronic tau.
 *
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param variation_dm0 variation for decay mode 0, options are "nom", "up",
 * "down"
 * @param variation_dm1 variation for decay mode 1, options are "nom", "up",
 * "down"
 * @param variation_dm10 variation for decay mode 10, options are "nom", "up",
 * "down"
 * @param variation_dm11 variation for decay mode 11, options are "nom", "up",
 * "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to genuine hadronic taus
 * (`gen_match=5`).
 */
ROOT::RDF::RNode PtCorrectionMC_genuineTau(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &pt,
    const std::string &eta, const std::string &decay_mode,
    const std::string &gen_match, const std::string &es_file,
    const std::string &correction_name, const std::string &id_algorithm,
    const std::string &variation_dm0, const std::string &variation_dm1,
    const std::string &variation_dm10, const std::string &variation_dm11) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::VecOps::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda =
        [evaluator, id_algorithm, variation_dm0, variation_dm1, variation_dm10,
         variation_dm11](const ROOT::RVec<float> &pts,
                         const ROOT::RVec<float> &etas,
                         const ROOT::RVec<UChar_t> &decay_modes_v12,
                         const ROOT::RVec<UChar_t> &gen_matches_char) {
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (gen_matches.at(i) == 5) {
                    if (decay_modes.at(i) == 0) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm0});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 1) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm1});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 10) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm10});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 11) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm11});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    }
                } else {
                    corrected_pts[i] = pts.at(i);
                }
                Logger::get("physicsobject::tau::PtCorrection_genuineTau")
                    ->debug("tau pt before {}, tau pt after {}, decaymode {}, gen match {}",
                            pts.at(i), corrected_pts.at(i), decay_modes.at(i), gen_matches.at(i));
            }
            return corrected_pts;
        };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}


ROOT::RDF::RNode PtCorrectionMC_genuineTau(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &pt,
    const std::string &eta, const std::string &decay_mode,
    const std::string &gen_match, const std::string &es_file,
    const std::string &correction_name, const std::string &id_algorithm,
    // New: 8 arguments for DM/Pt splitting
    const std::string &variation_dm0_pt20to40,
    const std::string &variation_dm0_pt40toInf,
    const std::string &variation_dm1_pt20to40,
    const std::string &variation_dm1_pt40toInf,
    const std::string &variation_dm10_pt20to40,
    const std::string &variation_dm10_pt40toInf,
    const std::string &variation_dm11_pt20to40,
    const std::string &variation_dm11_pt40toInf) {

    // 1. Define the variation map structure
    // This map acts as both the variation lookup AND the "Allowed Decay Modes" list
    const std::unordered_map<int, std::map<float, std::string>> variations = {
        {0,  {{20.0f, variation_dm0_pt20to40},  {40.0f, variation_dm0_pt40toInf}}},
        {1,  {{20.0f, variation_dm1_pt20to40},  {40.0f, variation_dm1_pt40toInf}}},
        {10, {{20.0f, variation_dm10_pt20to40}, {40.0f, variation_dm10_pt40toInf}}},
        {11, {{20.0f, variation_dm11_pt20to40}, {40.0f, variation_dm11_pt40toInf}}},
    };

    // 2. Setup
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::VecOps::RVec<UChar_t>", decay_mode);

    auto evaluator = correction_manager.loadCorrection(es_file, correction_name);

    // 3. Calculation Lambda
    auto correction_lambda =
        [evaluator, id_algorithm, variations](const ROOT::RVec<float> &pts,
                         const ROOT::RVec<float> &etas,
                         const ROOT::RVec<UChar_t> &decay_modes_v12,
                         const ROOT::RVec<UChar_t> &gen_matches_char) {
            
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
            ROOT::RVec<float> corrected_pts(pts.size());
            
            for (size_t i = 0; i < pts.size(); i++) {
                float current_pt = pts.at(i);
                int current_dm = decay_modes.at(i);
                int current_gen = gen_matches.at(i);
                
                // Default: No correction
                corrected_pts[i] = current_pt;
                std::string variation_to_use = "nom"; 

                // Only correct if GenMatch is 5 (Genuine Tau)
                if (current_gen == 5) {
                    // Look for DM in our allowed map
                    auto dm_it = variations.find(current_dm);
                    
                    // CRITICAL FIX: Only call evaluate if DM is supported (0, 1, 10, 11)
                    if (dm_it != variations.end()) {
                        const auto &pt_map = dm_it->second;
                        
                        // Determine variation based on Pt
                        auto pt_it = pt_map.upper_bound(current_pt);
                        if (pt_it != pt_map.begin()) {
                            variation_to_use = std::prev(pt_it)->second;
                        }

                        // Evaluate now safe because current_dm is guaranteed valid
                        auto correction_factor = evaluator->evaluate(
                            {current_pt, std::abs(etas.at(i)), current_dm,
                             current_gen, id_algorithm,
                             variation_to_use});
                        
                        corrected_pts[i] = current_pt * correction_factor;
                    }
                }

                // Logging
                Logger::get("physicsobject::tau::PtCorrection_genuineTau")
                    ->debug("tau pt before {}, tau pt after {}, decaymode {}, var {}",
                            current_pt, corrected_pts.at(i), current_dm, variation_to_use);
            }
            return corrected_pts;
        };

    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}



/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of genuine hadronic taus. The energy scale correction for these
 * objects is measured for four tau decay modes (dm0, dm1, dm10 and dm11) and
 * depends on the transverse momentum of the hadronic tau.
 *
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param wp working point of the vsJet ID
 * @param vsele_wp working point of the vsEle ID
 * @param variation_dm0 variation for decay mode 0, options are "nom", "up",
 * "down"
 * @param variation_dm1 variation for decay mode 1, options are "nom", "up",
 * "down"
 * @param variation_dm10 variation for decay mode 10, options are "nom", "up",
 * "down"
 * @param variation_dm11 variation for decay mode 11, options are "nom", "up",
 * "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to genuine hadronic taus
 * (`gen_match=5`).
 */
ROOT::RDF::RNode PtCorrectionMC_genuineTau(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &pt,
    const std::string &eta, const std::string &decay_mode,
    const std::string &gen_match, const std::string &es_file,
    const std::string &correction_name, const std::string &id_algorithm,
    const std::string &wp, const std::string &vsele_wp,
    const std::string &variation_dm0, const std::string &variation_dm1,
    const std::string &variation_dm10, const std::string &variation_dm11) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::VecOps::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda =
        [evaluator, id_algorithm, wp, vsele_wp, variation_dm0, variation_dm1, variation_dm10,
         variation_dm11](const ROOT::RVec<float> &pts,
                         const ROOT::RVec<float> &etas,
                         const ROOT::RVec<UChar_t> &decay_modes_v12,
                         const ROOT::RVec<UChar_t> &gen_matches_char) {
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (gen_matches.at(i) == 5) {
                    if (decay_modes.at(i) == 0) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm, wp, vsele_wp,
                             variation_dm0});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 1) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm, wp, vsele_wp,
                             variation_dm1});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 10) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm, wp, vsele_wp,
                             variation_dm10});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 11) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm, wp, vsele_wp,
                             variation_dm11});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    }
                } else {
                    corrected_pts[i] = pts.at(i);
                }
                Logger::get("physicsobject::tau::PtCorrection_genuineTau")
                    ->debug("tau pt before {}, tau pt after {}, decaymode {}",
                            pts.at(i), corrected_pts.at(i), decay_modes.at(i));
            }
            return corrected_pts;
        };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}


/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of genuine hadronic taus. The energy scale correction for these
 * objects is measured for four tau decay modes (dm0, dm1, dm10 and dm11) and
 * depends on the transverse momentum of the hadronic tau.
 *
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param wp working point of the vsJet ID
 * @param vsele_wp working point of the vsEle ID
 * @param variation_dm0_20to40 variation for decay mode 0, options are "nom", "up",
 * "down"
 * @param variation_dm1_20to40 variation for decay mode 1, options are "nom", "up",
 * "down"
 * @param variation_dm10_20to40 variation for decay mode 10, options are "nom", "up",
 * "down"
 * @param variation_dm11_20to40 variation for decay mode 11, options are "nom", "up",
 * "down"
 * @param variation_dm0_40toInf variation for decay mode 0, options are "nom", "up",
 * "down"
 * @param variation_dm1_40toInf variation for decay mode 1, options are "nom", "up",
 * "down"
 * @param variation_dm10_40toInf variation for decay mode 10, options are "nom", "up",
 * "down"
 * @param variation_dm11_40toInf variation for decay mode 11, options are "nom", "up",
 * "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to genuine hadronic taus
 * (`gen_match=5`).
 */
ROOT::RDF::RNode PtCorrectionMC_genuineTau(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &pt,
    const std::string &eta, const std::string &decay_mode,
    const std::string &gen_match, const std::string &es_file,
    const std::string &correction_name, const std::string &id_algorithm,
    const std::string &wp, const std::string &vsele_wp,
    const std::string &variation_dm0_20to40, const std::string &variation_dm1_20to40,
    const std::string &variation_dm10_20to40, const std::string &variation_dm11_20to40,
    const std::string &variation_dm0_40toInf, const std::string &variation_dm1_40toInf,
    const std::string &variation_dm10_40toInf, const std::string &variation_dm11_40toInf) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    const std::map<float, std::string> variations_dm0 = {
        {20.0f,     variation_dm0_20to40},
        {40.0f,     variation_dm0_40toInf},
        {100000.0f, variation_dm0_40toInf},
    };
    const std::map<float, std::string> variations_dm1 = {
        {20.0f,     variation_dm1_20to40},
        {40.0f,     variation_dm1_40toInf},
        {100000.0f, variation_dm1_40toInf},
    };
    const std::map<float, std::string> variations_dm10 = {
        {20.0f,     variation_dm10_20to40},
        {40.0f,     variation_dm10_40toInf},
        {100000.0f, variation_dm10_40toInf},
    };
    const std::map<float, std::string> variations_dm11 = {
        {20.0f,     variation_dm11_20to40},
        {40.0f,     variation_dm11_40toInf},
        {100000.0f, variation_dm11_40toInf},
    };

    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::VecOps::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda =
        [evaluator, id_algorithm, wp, vsele_wp, variations_dm0, variations_dm1, variations_dm10, variations_dm11](const ROOT::RVec<float> &pts,
                         const ROOT::RVec<float> &etas,
                         const ROOT::RVec<UChar_t> &decay_modes_v12,
                         const ROOT::RVec<UChar_t> &gen_matches_char) {
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (gen_matches.at(i) == 5) {
                    if (decay_modes.at(i) == 0) {
                        auto it = variations_dm0.upper_bound(pts.at(i));
                        if (it != variations_dm0.begin()){
                            it = std::prev(it);
                            std::string variation_dm0 = it->second;
                            auto correction_factor = evaluator->evaluate(
                                {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                                gen_matches.at(i), id_algorithm, wp, vsele_wp,
                                variation_dm0});
                            corrected_pts[i] = pts.at(i) * correction_factor;
                        } else {
                            corrected_pts[i] = pts.at(i);
                        }
                    } else if (decay_modes.at(i) == 1) {
                        auto it = variations_dm1.upper_bound(pts.at(i));
                        if (it != variations_dm1.begin()){
                            it = std::prev(it);
                            std::string variation_dm1 = it->second;
                            auto correction_factor = evaluator->evaluate(
                                {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                                gen_matches.at(i), id_algorithm, wp, vsele_wp,
                                variation_dm1});
                            corrected_pts[i] = pts.at(i) * correction_factor;
                        } else {
                            corrected_pts[i] = pts.at(i);
                        }
                    } else if (decay_modes.at(i) == 10) {
                        auto it = variations_dm10.upper_bound(pts.at(i));
                        if (it != variations_dm10.begin()){
                            it = std::prev(it);
                            std::string variation_dm10 = it->second;
                            auto correction_factor = evaluator->evaluate(
                                {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                                gen_matches.at(i), id_algorithm, wp, vsele_wp,
                                variation_dm10});
                            corrected_pts[i] = pts.at(i) * correction_factor;
                        } else {
                            corrected_pts[i] = pts.at(i);
                        }
                    } else if (decay_modes.at(i) == 11) {
                        auto it = variations_dm11.upper_bound(pts.at(i));
                        if (it != variations_dm11.begin()){
                            it = std::prev(it);
                            std::string variation_dm11 = it->second;
                            auto correction_factor = evaluator->evaluate(
                                {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                                gen_matches.at(i), id_algorithm, wp, vsele_wp,
                                variation_dm11});
                            corrected_pts[i] = pts.at(i) * correction_factor;
                        } else {
                            corrected_pts[i] = pts.at(i);
                        }
                    } else {
                        corrected_pts[i] = pts.at(i);
                    }
                } else {
                    corrected_pts[i] = pts.at(i);
                }
                Logger::get("physicsobject::tau::PtCorrection_genuineTau")
                    ->debug("tau pt before {}, tau pt after {}, decaymode {}",
                            pts.at(i), corrected_pts.at(i), decay_modes.at(i));
            }
            return corrected_pts;
        };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}
namespace quantity {

/**
 * @brief This function writes out a flag that is true if a tau passes a specific 
 * tau ID working point (`wp`). The working point is defined by the bit value. The
 * bit values can be found e.g. in the description of the tau ID scale factors.
 *
 * @note This function should be used only for `nanoAODv9`. Starting with `nanoAODv12`
 * `physicsobject::tau::quantity::IDFlag_v12` should be used instead because the content
 * of the tau ID branches was changed.
 * 
 * @param df input dataframe
 * @param outputname name of the output column containing the flag
 * @param ID name of the column containing the tau ID values 
 * @param index_vector name of the column containing the vector with the relevant
 * tau pair indices
 * @param position position in the index vector of the relevant tau in the pair
 * @param wp bit value of the WP that has to be passed
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode IDFlag_v9(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &ID, const std::string &index_vector,
                        const int &position, const int &wp) {
    return df.Define(
        outputname,
        [position, wp](const ROOT::RVec<int> &index_vector,
                          const ROOT::RVec<UChar_t> &IDs) {
            Logger::get("physicsobject::tau::quantity::IDFlag_v9")
                ->debug(
                    "position tau in pair {}, tau pair {}, id wp bit {}, vsjet ids {}",
                    position, index_vector, wp, IDs);
            const int index = index_vector.at(position);
            const int id_value = IDs.at(index, default_int);
            if (id_value != default_int)
                return std::min(1, int(id_value & 1 << (wp - 1)));
            else
                return int(id_value);
        },
        {index_vector, ID});
}

/**
 * @brief This function writes out a flag that is true if a tau passes a specific 
 * tau ID working point (`wp`). The content of the tau ID branch changed in 
 * `nanoAODv12`. The working points are defied in integer steps (still saved as 
 * `UChar_t`).
 * 
 * For `vsJet` and `vsEle`: 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 
 * 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight
 *
 * For `vsMu`: 1 = VLoose, 2 = Loose, 3 = Medium, 4 = Tight
 * 
 * @param df input dataframe
 * @param outputname name of the output column containing the flag
 * @param ID name of the column containing the tau ID values 
 * @param index_vector name of the column containing the vector with the relevant
 * tau pair indices
 * @param position position in the index vector of the relevant tau in the pair
 * @param wp bit value of the WP that has to be passed
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode IDFlag_v12(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &ID, const std::string &index_vector,
                        const int &position, const int &wp) {
    return df.Define(
        outputname,
        [position, wp](const ROOT::RVec<int> &index_vector,
                          const ROOT::RVec<UChar_t> &IDs) {
            Logger::get("physicsobject::tau::quantity::IDFlag_v12")
                ->debug(
                    "position tau in pair {}, tau pair {}, id wp bit {}, vsjet ids {}",
                    position, index_vector, wp, IDs);
            const int index = index_vector.at(position);
            const int id_value = static_cast<int>(IDs.at(index, default_int));
            if (id_value != default_int)
                return int(id_value >= wp);
            else
                return id_value;
        },
        {index_vector, ID});
}
} // end namespace quantity

namespace scalefactor {

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against jets (`vsJet`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in transverse momenta (\f$p_T\f$) of hadronic taus 
 * for this function. This dependence is usually used in semi-leptonic channels 
 * (\f$e\tau_h\f$, \f$\mu\tau_h\f$).
 *
 * Description of the bit map used to define the tau ID against jets working points of the
 * DeepTau v2.1 tagger.
 * vsJets                              | Value | Bit (value used in the config)
 * ------------------------------------|-------|-------
 * no ID selection (takes every tau)   |  0    | -
 * VVVLoose                            |  1    | 1
 * VVLoose                             |  2    | 2
 * VLoose                              |  4    | 3
 * Loose                               |  8    | 4
 * Medium                              |  16   | 5
 * Tight                               |  32   | 6
 * VTight                              |  64   | 7
 * VVTight                             |  128  | 8
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the vsJets ID scale factor
 * @param pt name of the column containing the transverse momentum of a tau
 * @param decay_mode name of the column containing the decay mode of the tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsJet ID correction
 * @param selected_dms list of allowed decay modes for which a scale factor
 * should be calculated
 * @param wp working point of the vsJet ID
 * @param vsele_wp working point of the vsEle ID
 * @param sf_dependence variable dependence of the scale factor, opions are "pt" or "dm"
 * @param variation_pt30to35 name of the scale factor variation for \f$30 \leq p_T <35\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt35to40 name of the scale factor variation for \f$35 \leq p_T <40\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt40to500 name of the scale factor variation for \f$40 \leq p_T <500\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt500to1000 name of the scale factor variation for \f$500 \leq p_T <1000\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt1000toInf name of the scale factor variation for \f$1000 \leq p_T < \infty \f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Id_vsJet_lt(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correction_manager,
            const std::string &outputname,
            const std::string &pt, const std::string &decay_mode,
            const std::string &gen_match, 
            const std::string &sf_file,
            const std::string &sf_name,
            const std::vector<int> &selected_dms,
            const std::string &wp, const std::string &vsele_wp,
            const std::string &sf_dependence,
            const std::string &variation_pt30to35,
            const std::string &variation_pt35to40,
            const std::string &variation_pt40to500,
            const std::string &variation_pt500to1000,
            const std::string &variation_pt1000toInf) {
    
    const std::map<float, std::string> variations = {
        {30.0f,     variation_pt30to35},
        {35.0f,     variation_pt35to40},
        {40.0f,     variation_pt40to500},
        {500.0f,    variation_pt500to1000},
        {1000.0f,   variation_pt1000toInf},
        {100000.0f, variation_pt1000toInf},
    };
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, vsele_wp, variations,
                            sf_dependence, selected_dms,
                            sf_name](const float &pt, const int &decay_mode,
                                         const int &gen_match) {
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("ID - decayMode {}", decay_mode);
        // only calculate SFs for allowed tau decay modes (also excludes default
        // values due to tau energy correction shifts below good tau pt
        // selection)
        double sf = 1.;
        if (std::find(selected_dms.begin(), selected_dms.end(), decay_mode) !=
            selected_dms.end()) {
            auto it = variations.upper_bound(pt);
            if (it != variations.begin()){
                it = std::prev(it);
                std::string variation = it->second;
                Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation, sf_dependence);
                sf = evaluator->evaluate({pt, decay_mode, gen_match, wp, vsele_wp, variation, sf_dependence});
            } else {
                sf = 1.;
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, gen_match});
    return df1;
}



ROOT::RDF::RNode
Id_vsJet_lt(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &pt, const std::string &decay_mode,
        const std::string &gen_match, 
        const std::string &sf_file,
        const std::string &sf_name,
        const std::vector<int> &selected_dms,
        const std::string &wp, const std::string &vsele_wp,
        const std::string &sf_dependence,
        // DM + pT split
        const std::string &variation_dm0_pt20to40,
        const std::string &variation_dm0_pt40toInf,
        const std::string &variation_dm1_pt20to40,
        const std::string &variation_dm1_pt40toInf,
        const std::string &variation_dm10_pt20to40,
        const std::string &variation_dm10_pt40toInf,
        const std::string &variation_dm11_pt20to40,
        const std::string &variation_dm11_pt40toInf) {

    const std::unordered_map<int, std::map<float, std::string>> variations = {
        {0,  {{20.0f, variation_dm0_pt20to40},  {40.0f, variation_dm0_pt40toInf}}},
        {1,  {{20.0f, variation_dm1_pt20to40},  {40.0f, variation_dm1_pt40toInf}}},
        {10, {{20.0f, variation_dm10_pt20to40}, {40.0f, variation_dm10_pt40toInf}}},
        {11, {{20.0f, variation_dm11_pt20to40}, {40.0f, variation_dm11_pt40toInf}}},
    };

    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")
        ->debug("Setting up function for tau id vsJet sf (DM & Pt binned)");
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("ID - Name {}", sf_name);
    
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    
    auto sf_calculator = [evaluator, wp, vsele_wp, variations,
                            sf_dependence, selected_dms,
                            sf_name](const float &pt, const int &decay_mode,
                                         const int &gen_match) {
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("ID - decayMode {}", decay_mode);
        
        double sf = 1.;
        
        if (std::find(selected_dms.begin(), selected_dms.end(), decay_mode) != selected_dms.end()) {
            auto dm_it = variations.find(decay_mode);  // specific map for this Decay Mode
            if (dm_it != variations.end()) {
                const auto &pt_map = dm_it->second;
                
                auto pt_it = pt_map.upper_bound(pt);
                
                if (pt_it != pt_map.begin()){
                    pt_it = std::prev(pt_it); // Move back to the lower bound key
                    std::string variation = pt_it->second;
                    
                    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")
                    ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                            "vsele_wp {}, variation {}, sf_dependence {}",
                            sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                            variation, sf_dependence);
                    
                    sf = evaluator->evaluate({pt, decay_mode, gen_match, wp, vsele_wp, variation, sf_dependence});
                }
                // If pt_it == pt_map.begin(), it means pt < 20.0f (the first key), so sf remains 1.
            }
        }
        
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("Scale Factor {}", sf);
        return sf;
    };

    bool all_nominal = true;
    for (const auto& [dm, pt_map] : variations) {
        for (const auto& [pt_val, var_name] : pt_map) {
            if (var_name != "nom") {
                all_nominal = false;
                break;
            }
        }
        if (!all_nominal) break;
    }

    if (all_nominal) {
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")
            ->debug("All variations are 'nom'. Skipping definition of {} to avoid collision.", outputname);
        return df;
    }

    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, gen_match});
    return df1;
}




/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against jets (`vsJet`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in decay modes of hadronic taus for this function. This 
 * dependence is usually used in the fully hadronic channel (\f$\tau_h\tau_h\f$).
 *
 * Description of the bit map used to define the tau ID against jets working points of the
 * DeepTau v2.1 tagger.
 * vsJets                              | Value | Bit (value used in the config)
 * ------------------------------------|-------|-------
 * no ID selection (takes every tau)   |  0    | -
 * VVVLoose                            |  1    | 1
 * VVLoose                             |  2    | 2
 * VLoose                              |  4    | 3
 * Loose                               |  8    | 4
 * Medium                              |  16   | 5
 * Tight                               |  32   | 6
 * VTight                              |  64   | 7
 * VVTight                             |  128  | 8
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the vsJets ID scale factor
 * @param pt name of the column containing the transverse momentum of a tau
 * @param decay_mode name of the column containing the decay mode of the tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsJet ID correction
 * @param wp working point of the vsJet ID
 * @param vsele_wp working point of the vsEle ID
 * @param sf_dependence variable dependence of the scale factor, opions are "pt" or "dm"
 * @param variation_dm0 name of the scale factor variation for decay mode 0, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm1 name of the scale factor variation for decay mode 1, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm10 name of the scale factor variation for decay mode 10, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm11 name of the scale factor variation for decay mode 11, "nom" for nominal
 * and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Id_vsJet_tt(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, 
    const std::string &pt, const std::string &decay_mode,
    const std::string &gen_match, 
    const std::string &sf_file, const std::string &sf_name,
    const std::string &wp, const std::string &vsele_wp,
    const std::string &sf_dependence,
    const std::string &variation_dm0,
    const std::string &variation_dm1, 
    const std::string &variation_dm10,
    const std::string &variation_dm11) {
    
    const std::unordered_map<int, std::string> variations = {
        {0,  variation_dm0},
        {1,  variation_dm1},
        {10, variation_dm10},
        {11, variation_dm11},
    };
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, vsele_wp, variations,
                            sf_dependence, sf_name](const float &pt, const int &decay_mode,
                                         const int &gen_match) {
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("ID - decayMode {}", decay_mode);
        // only calculate SFs for allowed tau decay modes (also excludes default
        // values due to tau energy correction shifts below good tau pt
        // selection)
        double sf = 1.;
        if (auto it = variations.find(decay_mode); it != variations.end()) {
            std::string variation = it->second;
            Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation, sf_dependence);
            sf = evaluator->evaluate(
                {pt, decay_mode, gen_match, wp, vsele_wp, 
                variation, sf_dependence});
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against jets (`vsJet`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in decay modes of hadronic taus for this function. This 
 * dependence is usually used in the fully hadronic channel (\f$\tau_h\tau_h\f$).
 *
 * Description of the bit map used to define the tau ID against jets working points of the
 * DeepTau v2.1 tagger.
 * vsJets                              | Value | Bit (value used in the config)
 * ------------------------------------|-------|-------
 * no ID selection (takes every tau)   |  0    | -
 * VVVLoose                            |  1    | 1
 * VVLoose                             |  2    | 2
 * VLoose                              |  4    | 3
 * Loose                               |  8    | 4
 * Medium                              |  16   | 5
 * Tight                               |  32   | 6
 * VTight                              |  64   | 7
 * VVTight                             |  128  | 8
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the vsJets ID scale factor
 * @param pt name of the column containing the transverse momentum of a tau
 * @param decay_mode name of the column containing the decay mode of the tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsJet ID correction
 * @param wp working point of the vsJet ID
 * @param vsele_wp working point of the vsEle ID
 * @param sf_dependence variable dependence of the scale factor, opions are "pt" or "dm"
 * @param variation_dm0 name of the scale factor variation for decay mode 0, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm1 name of the scale factor variation for decay mode 1, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm10 name of the scale factor variation for decay mode 10, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm11 name of the scale factor variation for decay mode 11, "nom" for nominal
 * and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Id_vsJet_tt(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, 
    const std::string &pt, const std::string &decay_mode,
    const std::string &gen_match, 
    const std::string &sf_file, const std::string &sf_name,
    const std::string &wp, const std::string &vsele_wp,
    const std::string &sf_dependence,
    const std::string &variation_dm0_20to40, const std::string &variation_dm1_20to40,
    const std::string &variation_dm10_20to40, const std::string &variation_dm11_20to40,
    const std::string &variation_dm0_40toInf, const std::string &variation_dm1_40toInf,
    const std::string &variation_dm10_40toInf, const std::string &variation_dm11_40toInf) {
    
    // const std::unordered_map<int, std::string> variations = {
    const std::map<float, std::string> variations_dm0 = {
        {20.0f,     variation_dm0_20to40},
        {40.0f,     variation_dm0_40toInf},
        {100000.0f, variation_dm0_40toInf},
    };
    const std::map<float, std::string> variations_dm1 = {
        {20.0f,     variation_dm1_20to40},
        {40.0f,     variation_dm1_40toInf},
        {100000.0f, variation_dm1_40toInf},
    };
    const std::map<float, std::string> variations_dm10 = {
        {20.0f,     variation_dm10_20to40},
        {40.0f,     variation_dm10_40toInf},
        {100000.0f, variation_dm10_40toInf},
    };
    const std::map<float, std::string> variations_dm11 = {
        {20.0f,     variation_dm11_20to40},
        {40.0f,     variation_dm11_40toInf},
        {100000.0f, variation_dm11_40toInf},
    };
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, vsele_wp, variations_dm0, variations_dm1, variations_dm10, variations_dm11,
                            sf_dependence, sf_name](const float &pt, const int &decay_mode,
                                         const int &gen_match) {
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("ID - decayMode {}", decay_mode);
        // only calculate SFs for allowed tau decay modes (also excludes default
        // values due to tau energy correction shifts below good tau pt
        // selection)
        double sf = 1.;
        if (decay_mode == 0) {
            auto it = variations_dm0.upper_bound(pt);
            if (it != variations_dm0.begin()){
                it = std::prev(it);
                std::string variation_dm0 = it->second;
                Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation_dm0, sf_dependence);
            sf = evaluator->evaluate(
                {pt, decay_mode, gen_match, wp, vsele_wp, 
                variation_dm0, sf_dependence});
            }
        } else if (decay_mode == 1) {
            auto it = variations_dm1.upper_bound(pt);
            if (it != variations_dm1.begin()){
                it = std::prev(it);
                std::string variation_dm1 = it->second;
                Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation_dm1, sf_dependence);
            sf = evaluator->evaluate(
                {pt, decay_mode, gen_match, wp, vsele_wp, 
                variation_dm1, sf_dependence});
            }
        } else if (decay_mode == 10) {
            auto it = variations_dm10.upper_bound(pt);
            if (it != variations_dm10.begin()){
                it = std::prev(it);
                std::string variation_dm10 = it->second;
                Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation_dm10, sf_dependence);
            sf = evaluator->evaluate(
                {pt, decay_mode, gen_match, wp, vsele_wp, 
                variation_dm10, sf_dependence});
            }
        } else if (decay_mode == 11) {
            auto it = variations_dm11.upper_bound(pt);
            if (it != variations_dm11.begin()){
                it = std::prev(it);
                std::string variation_dm11 = it->second;
                Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation_dm11, sf_dependence);
            sf = evaluator->evaluate(
                {pt, decay_mode, gen_match, wp, vsele_wp, 
                variation_dm11, sf_dependence});
            }
        }
        // if (auto it = variations.find(decay_mode); it != variations.end()) {
        //     std::string variation = it->second;
        //     Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
        //         ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
        //                 "vsele_wp {}, variation {}, sf_dependence {}",
        //                 sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
        //                 variation, sf_dependence);
        //     sf = evaluator->evaluate(
        //         {pt, decay_mode, gen_match, wp, vsele_wp, 
        //         variation, sf_dependence});
        // }
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, gen_match});
    return df1;
}


/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against electrons (`vsEle`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in pseudorapidities of hadronic taus for this function.
 *
 * Description of the bit map used to define the tau ID against electrons working points of the
 * DeepTau v2.1 tagger.
 * vsElectrons                         | Value | Bit (value used in the config)
 * ------------------------------------|-------|-------
 * no ID selection (takes every tau)   |  0    | -
 * VVVLoose                            |  1    | 1
 * VVLoose                             |  2    | 2
 * VLoose                              |  4    | 3
 * Loose                               |  8    | 4
 * Medium                              |  16   | 5
 * Tight                               |  32   | 6
 * VTight                              |  64   | 7
 * VVTight                             |  128  | 8
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the vsEle ID scale factor
 * @param eta name of the column containing the pseudorapidity of a tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsEle ID correction
 * @param wp working point of the vsEle ID
 * @param variation_barrel name of the scale factor variation for the barrel region
 * (\f$|\eta| <1.46\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_endcap name of the scale factor variation for the endcap region
 * (\f$1.558 \leq |\eta| <2.3\f$), "nom" for nominal and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 * 
 * @note This function is intended to be used with Run 2 analyses. The scale factor additionally
 * depends on the decay mode of the tau in Run 3. For Run 3 analyses, use the overloaded version
 * of this function.
 */
ROOT::RDF::RNode
Id_vsEle(ROOT::RDF::RNode df,
         correctionManager::CorrectionManager &correction_manager,
         const std::string &outputname,
         const std::string &eta,
         const std::string &gen_match, 
         const std::string &sf_file, const std::string &sf_name,
         const std::string &wp, 
         const std::string &variation_barrel,
         const std::string &variation_endcap) {

    Logger::get("physicsobject::tau::scalefactor::Id_vsEle")
        ->debug("Setting up function for tau id vsEle sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsEle")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, variation_barrel, variation_endcap,
                            sf_name](const float &eta, const int &gen_match) {
        double sf = 1.;

        // set edges of barrel and endcap region
        double max_abs_eta_barrel = 1.46;
        double min_abs_eta_endcap = 1.558;
        double max_abs_eta_endcap = 2.3;

        // exclude default values due to tau energy correction shifts below good tau
        // pt selection
        if (eta > -5.0) {
            Logger::get("physicsobject::tau::scalefactor::Id_vsEle")
                ->debug("ID {} - eta {}, gen_match {}, wp {}, variation_barrel "
                        "{}, variation_endcap {}",
                        sf_name, eta, gen_match, wp, variation_barrel,
                        variation_endcap);
            // the eta cuts are taken from the correctionlib json file to define 
            // barrel and endcap
            if (std::abs(eta) < max_abs_eta_barrel) {
                sf = evaluator->evaluate(
                    {eta, gen_match, wp, variation_barrel});
            } else if (std::abs(eta) >= min_abs_eta_endcap && std::abs(eta) < max_abs_eta_endcap) {
                sf = evaluator->evaluate(
                    {eta, gen_match, wp, variation_endcap});
            } else {
                sf = 1.;
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsEle")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 =
        df.Define(outputname, sf_calculator, {eta, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against electrons (`vsEle`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in pseudorapidities of hadronic taus for this function. This
 * function is intended to be used with Run 3 analyses.
 *
 * Description of the bit map used to define the tau ID against electrons working points of the
 * DeepTau v2.5 tagger.
 * vsElectrons                         |  ID value (value used in the config)
 * ------------------------------------|--------
 * no ID selection (takes every tau)   |  -
 * VVVLoose                            |  1
 * VVLoose                             |  2
 * VLoose                              |  3
 * Loose                               |  4
 * Medium                              |  5
 * Tight                               |  6
 * VTight                              |  7
 * VVTight                             |  8
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the vsEle ID scale factor
 * @param eta name of the column containing the pseudorapidity of a tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param decay_mode name of the column containing the decay mode of the tau
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsEle ID correction
 * @param wp working point of the vsEle ID
 * @param variation_barrel name of the scale factor variation for the barrel region
 * (\f$|\eta| <1.46\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_endcap name of the scale factor variation for the endcap region
 * (\f$1.558 \leq |\eta| <2.3\f$), "nom" for nominal and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 * 
 * @note This function is intended to be used with Run 3 analyses. The scale factor additionally
 * depends on the decay mode of the tau. For Run 2 analyses, use the overloaded version of this
 * function.
 */
ROOT::RDF::RNode
Id_vsEle(ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname,
    const std::string &eta,
    const std::string &decay_mode,
    const std::string &gen_match, 
    const std::string &sf_file, const std::string &sf_name,
    const std::string &wp, 
    const std::string &variation_barrel,
    const std::string &variation_endcap
) {
    Logger::get("physicsobject::tau::scalefactor::Id_vsEle")
        ->debug("Setting up function for tau id vsEle sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsEle")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, variation_barrel, variation_endcap,
                            sf_name](const float &eta, const int &dm, const int &gen_match) {
        double sf = 1.;

        // set edges of barrel and endcap region
        double max_abs_eta_barrel = 1.46;
        double min_abs_eta_endcap = 1.558;
        double max_abs_eta_endcap = 2.5;

        // exclude default values due to tau energy correction shifts below good tau
        // pt selection
        if (eta > -5.0) {
            Logger::get("physicsobject::tau::scalefactor::Id_vsEle")
                ->debug("ID {} - eta {}, dm {}, gen_match {}, wp {}, variation_barrel "
                        "{}, variation_endcap {}",
                        sf_name, eta, dm, gen_match, wp, variation_barrel,
                        variation_endcap);
            // the eta cuts are taken from the correctionlib json file to define 
            // barrel and endcap
            if (std::abs(eta) < max_abs_eta_barrel) {
                sf = evaluator->evaluate(
                    {eta, dm, gen_match, wp, variation_barrel});
            } else if (std::abs(eta) >= min_abs_eta_endcap && std::abs(eta) < max_abs_eta_endcap) {
                sf = evaluator->evaluate(
                    {eta, dm, gen_match, wp, variation_endcap});
            } else {
                sf = 1.;
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsEle")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 =
        df.Define(outputname, sf_calculator, {eta, decay_mode, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against muons (`vsMu`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in pseudorapidities of hadronic taus for this function.
 *
 * Description of the bit map used to define the tau ID against muons working points of the
 * DeepTau v2.1 tagger.
 * vsMuons                             | Value | Bit (value used in the config)
 * ------------------------------------|-------|-------
 * no ID selection (takes every tau)   |  0    | -
 * VLoose                              |  1    | 1
 * Loose                               |  2    | 2
 * Medium                              |  4    | 3
 * Tight                               |  8    | 4
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the vsMu ID scale factor
 * @param eta name of the column containing the pseudorapidity of a tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsMu ID correction
 * @param wp working point of the vsMu ID
 * @param variation_wheel1 name of the scale factor variation for the muon wheel
 * (\f$|\eta| <0.4\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel2 name of the scale factor variation for the muon wheel
 * (\f$0.4 \leq |\eta| <0.8\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel3 name of the scale factor variation for the muon wheel
 * (\f$0.8 \leq |\eta| <1.2\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel4 name of the scale factor variation for the muon wheel
 * (\f$1.2 \leq |\eta| <1.7\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel5 name of the scale factor variation for the muon wheel
 * (\f$1.7 \leq |\eta| <2.3\f$), "nom" for nominal and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Id_vsMu(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &eta,
        const std::string &gen_match, 
        const std::string &sf_file,
        const std::string &sf_name,
        const std::string &wp,  
        const std::string &variation_wheel1,
        const std::string &variation_wheel2, 
        const std::string &variation_wheel3,
        const std::string &variation_wheel4, 
        const std::string &variation_wheel5) {
    
    const std::map<float, std::string> variations = {
        {0.0f, variation_wheel1},
        {0.4f, variation_wheel2},
        {0.8f, variation_wheel3},
        {1.2f, variation_wheel4},
        {1.7f, variation_wheel5},
        {2.4f, variation_wheel5},  // 2.4 to cover full muon system acceptance for Run 3 taus
                                   // should not affect Run 2 analyses, which cut on |eta| < 2.3
    };
    Logger::get("physicsobject::tau::scalefactor::Id_vsMu")->debug("Setting up function for tau id vsMu sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsMu")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, variations,
                            sf_name](const float &eta, const int &gen_match) {
        double sf = 1.;
        // exclude default values due to tau energy correction shifts below good tau
        // pt selection
        if (eta > -5.0) {
            auto it = variations.upper_bound(abs(eta));
            if (it != variations.begin()){
                it = std::prev(it);
                std::string variation = it->second;
                Logger::get("physicsobject::tau::scalefactor::Id_vsMu")
                    ->debug("ID {} - eta {}, gen_match {}, wp {}, variation {}",
                        sf_name, eta, gen_match, wp, variation);
                sf = evaluator->evaluate({eta, gen_match, wp, variation});
            } else {
                sf = 1.;
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsMu")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 =
        df.Define(outputname, sf_calculator, {eta, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau triggers. The scale factors 
 * are loaded from a correctionlib file using a specified scale factor name and variation. 
 * The scale factor is binned in \f$p_T\f$ and decay modes of hadronic taus for this function.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the trigger scale factor
 * @param pt name of the column containing the transverse momentum of a tau
 * @param decay_mode name of the column containing the decay mode of the tau
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the trigger correction
 * @param trigger_name name of the trigger, e.g. "ditau", "etau", "mutau"
 * @param wp working point of the vsJet ID
 * @param corr_type type of the value to be read out, options are "sf" (for 
 * scale factors), "eff_data", "eff_mc"
 * @param variation name the scale factor variation, "nom" for the nominal
 * scale factor and "up"/"down" for the up/down variation
 *
 * @return a new dataframe containing the new column
 *
 * @warning This function is deprecated. Use the overloaded function with the additional
 * parameter `trigger_flag` instead.
 */
ROOT::RDF::RNode
Trigger(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &pt, const std::string &decay_mode, 
        const std::string &sf_file,
        const std::string &sf_name,
        const std::string &trigger_name, const std::string &wp,
        const std::string &corr_type, const std::string &variation) {
    Logger::get("physicsobject::tau::scalefactor::Trigger")
        ->warn("Function is deprecated, use the overloaded version instead");
    Logger::get("physicsobject::tau::scalefactor::Trigger")
        ->debug("Setting up function for tau trigger sf");
    Logger::get("physicsobject::tau::scalefactor::Trigger")
        ->debug("ID - Name {}, file {}", sf_name, sf_file);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    Logger::get("physicsobject::tau::scalefactor::Trigger")->debug("WP {} - type {}", wp, corr_type);
    auto sf_calculator = [evaluator, trigger_name, wp, corr_type, variation, sf_name](
                                     const float &pt, const int &decay_mode) {
        float sf = 1.;
        try {
            Logger::get("physicsobject::tau::scalefactor::Trigger")
                ->debug("ID {} - decaymode {}, wp {} "
                    "pt {}, type {}, variation {}",
                    sf_name, decay_mode, wp, pt, corr_type, variation);
            if (pt >= 0.) {
                if (decay_mode == 0 || decay_mode == 1 || decay_mode == 10 ||
                    decay_mode == 11) {
                    sf = evaluator->evaluate(
                        {pt, decay_mode, trigger_name, wp, corr_type, variation});
                } else {
                    sf = evaluator->evaluate({pt, -1, trigger_name, wp, corr_type, variation});
                }
            }
        } catch (const std::runtime_error &e) {
            Logger::get("physicsobject::tau::scalefactor::Trigger")
                ->debug("SF evaluation for {} failed for pt {}", sf_name, pt);
        }
        Logger::get("physicsobject::tau::scalefactor::Trigger")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau triggers. The scale factors 
 * are loaded from a correctionlib file using a specified scale factor name and variation. 
 * The scale factor is binned in \f$p_T\f$ and decay modes of hadronic taus for this function.
 * This function only uses the scale factor from the correctionlib evaluation if the
 * corresponding trigger flag is set to `true`. Otherwise, it returns a scale factor of
 * 1.0.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the trigger scale factor
 * @param pt name of the column containing the transverse momentum of a tau
 * @param decay_mode name of the column containing the decay mode of the tau
 * @param trigger_flag name of the column containing the trigger flag
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the trigger correction
 * @param trigger_name name of the trigger, e.g. "ditau", "etau", "mutau"
 * @param wp working point of the vsJet ID
 * @param corr_type type of the value to be read out, options are "sf" (for 
 * scale factors), "eff_data", "eff_mc"
 * @param variation name the scale factor variation, "nom" for the nominal
 * scale factor and "up"/"down" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Trigger(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &pt, const std::string &decay_mode, 
        const std::string &trigger_flag,
        const std::string &sf_file,
        const std::string &sf_name,
        const std::string &trigger_name, const std::string &wp,
        const std::string &corr_type, const std::string &variation) {

    Logger::get("physicsobject::tau::scalefactor::Trigger")
        ->debug("Setting up function for tau trigger sf");
    Logger::get("physicsobject::tau::scalefactor::Trigger")
        ->debug("ID - Name {}, file {}", sf_name, sf_file);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    Logger::get("physicsobject::tau::scalefactor::Trigger")->debug("WP {} - type {}", wp, corr_type);
    auto sf_calculator = [evaluator, trigger_name, wp, corr_type, variation, sf_name](
                                     const float &pt, const int &decay_mode,
                                     const bool &trigger_flag) {
        float sf = 1.;
        try {
            Logger::get("physicsobject::tau::scalefactor::Trigger")
            ->debug("ID {} - decaymode {}, wp {} "
                "pt {}, trigger_flag {}, type {}, variation {}",
                sf_name, decay_mode, wp, pt, trigger_flag, corr_type, variation);
            if (pt >= 0. && trigger_flag) {
                if (decay_mode == 0 || decay_mode == 1 || decay_mode == 10 ||
                    decay_mode == 11) {
                    sf = evaluator->evaluate(
                        {pt, decay_mode, trigger_name, wp, corr_type, variation});
                } else {
                    sf = evaluator->evaluate({pt, -1, trigger_name, wp, corr_type, variation});
                }
            }
        } catch (const std::runtime_error &e) {
            Logger::get("physicsobject::tau::scalefactor::Trigger")
                ->debug("SF evaluation for {} failed for pt {}", sf_name, pt);
        }
        Logger::get("physicsobject::tau::scalefactor::Trigger")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, trigger_flag});
    return df1;
}
} // end namespace scalefactor
} // end namespace tau
} // end namespace physicsobject
#endif /* GUARD_TAUS_H */
