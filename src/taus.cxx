#ifndef GUARD_TAUS_H
#define GUARD_TAUS_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "correction.h"

namespace physicsobject {
namespace tau {

/**
 * @brief This function corrects the transverse momentum (pT) in MC simulations of 
 * hadronic taus that originate from electrons that are misidentified. The energy 
 * scale correction for these objects is measured for two tau decay modes (dm0 and 
 * dm1) and depends on the barrel and endcap region of the detector. 
 *
 * The correction procedure is taken from the officially recommendation of the 
 * TauPOG:
 *
 * Run2 (UL): https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3 (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the 
 * correction file
 * @param corrected_pt name of the output column storing the corrected hadronic 
 * tau pT values
 * @param pt name of the input column containing hadronic tau pT values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the hadronic 
 * tau to generator-level particles
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
 * @note This correction is only applied to misidentified hadronic taus originating 
 * from prompt electrons (`gen_match=1`) and electrons that decayed from a tau lepton 
 * (`gem_match=3`).
 */
ROOT::RDF::RNode
PtCorrectionMC_eleFake(ROOT::RDF::RNode df,
                       correctionManager::CorrectionManager &correction_manager,
                       const std::string &corrected_pt, 
                       const std::string &pt,
                       const std::string &eta, 
                       const std::string &decay_mode,
                       const std::string &gen_match, 
                       const std::string &es_file,
                       const std::string &correction_name,
                       const std::string &id_algorithm,
                       const std::string &variation_dm0_barrel, 
                       const std::string &variation_dm1_barrel,
                       const std::string &variation_dm0_endcap, 
                       const std::string &variation_dm1_endcap) {
    auto evaluator = correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda = [evaluator, id_algorithm, variation_dm0_barrel, 
                              variation_dm1_barrel, variation_dm0_endcap, 
                              variation_dm1_endcap](
                                const ROOT::RVec<float> &pts,
                                const ROOT::RVec<float> &etas,
                                const ROOT::RVec<int> &decay_modes,
                                const ROOT::RVec<UChar_t> &gen_matches) {
        ROOT::RVec<float> corrected_pts(pts.size());
        for (int i = 0; i < pts.size(); i++) {
            if (gen_matches.at(i) == 1 || gen_matches.at(i) == 3) {
                if (decay_modes.at(i) == 0 && 
                    std::abs(etas.at(i)) <= 1.5) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm0_barrel});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                } else if (decay_modes.at(i) == 0 &&
                           std::abs(etas.at(i)) > 1.5 &&
                           std::abs(etas.at(i)) <= 2.5) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm0_endcap});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(etas.at(i)) <= 1.5) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm1_barrel});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(etas.at(i)) > 1.5 &&
                           std::abs(etas.at(i)) <= 2.5) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm1_endcap});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                }
            } else {
                corrected_pts[i] = pts.at(i);
            }
            Logger::get("PtCorrectionMC_eleFake")
                ->debug("tau pt before {}, tau pt after {}", pts.at(i),
                        corrected_pts.at(i));
        }
        return corrected_pts;
    };
    auto df1 = df.Define(corrected_pt, correction_lambda,
                         {pt, eta, decay_mode, gen_match});
    return df1;
}

/**
 * @brief This function corrects the transverse momentum (pT) in MC simulations of 
 * hadronic taus that originate from muons that are misidentified. The energy 
 * scale correction is always set to `1` and for the uncertainty varied by `1%`.
 *
 * The correction procedure is taken from the officially recommendation of the 
 * TauPOG:
 *
 * Run2 (UL): https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3 (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the 
 * correction file
 * @param corrected_pt name of the output column storing the corrected hadronic 
 * tau pT values
 * @param pt name of the input column containing hadronic tau pT values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the hadronic 
 * tau to generator-level particles
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param variation variation of the correction, options are "nom", "up", "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to misidentified hadronic taus originating 
 * from prompt muons (`gen_match=2`) and muons that decayed from a tau lepton 
 * (`gem_match=4`).
 */
ROOT::RDF::RNode
PtCorrectionMC_muFake(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correction_manager,
                      const std::string &corrected_pt, 
                      const std::string &pt,
                      const std::string &eta, 
                      const std::string &decay_mode,
                      const std::string &gen_match, 
                      const std::string &es_file,
                      const std::string &correction_name,
                      const std::string &id_algorithm, 
                      const std::string &variation) {
    auto evaluator = correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda =
        [evaluator, id_algorithm, variation](const ROOT::RVec<float> &pts,
                                        const ROOT::RVec<float> &etas,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &gen_matches) {
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (gen_matches.at(i) == 2 || gen_matches.at(i) == 4) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                } else {
                    corrected_pts[i] = pts.at(i);
                }
                Logger::get("PtCorrectionMC_muFake")->debug(
                    "tau pt before {}, tau pt after {}", pts.at(i),
                    corrected_pts.at(i));
            }
            return corrected_pts;
        };
    auto df1 = df.Define(corrected_pt, correction_lambda,
                         {pt, eta, decay_mode, gen_match});
    return df1;
}

/**
 * @brief This function corrects the transverse momentum (pT) in MC simulations of 
 * genuine hadronic taus. The energy scale correction for these objects is measured 
 * for four tau decay modes (dm0, dm1, dm10 and dm11) and depends on the transverse 
 * momentum of the hadronic tau.
 *
 * The correction procedure is taken from the officially recommendation of the 
 * TauPOG:
 *
 * Run2 (UL): https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * - https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3 (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the 
 * correction file
 * @param corrected_pt name of the output column storing the corrected hadronic 
 * tau pT values
 * @param pt name of the input column containing hadronic tau pT values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the hadronic 
 * tau to generator-level particles
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param variation_dm0 variation for decay mode 0, options are "nom", "up", "down"
 * @param variation_dm1 variation for decay mode 1, options are "nom", "up", "down"
 * @param variation_dm10 variation for decay mode 10, options are "nom", "up", "down"
 * @param variation_dm11 variation for decay mode 11, options are "nom", "up", "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to genuine hadronic taus (`gen_match=5`).
 */
ROOT::RDF::RNode
PtCorrectionMC_genuineTau(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correction_manager,
                    const std::string &corrected_pt, 
                    const std::string &pt,
                    const std::string &eta, 
                    const std::string &decay_mode,
                    const std::string &gen_match, 
                    const std::string &es_file,
                    const std::string &correction_name,
                    const std::string &id_algorithm, 
                    const std::string &variation_dm0,
                    const std::string &variation_dm1, 
                    const std::string &variation_dm10,
                    const std::string &variation_dm11) {
    auto evaluator = correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda = [evaluator, id_algorithm, variation_dm0, variation_dm1, 
                              variation_dm10, variation_dm11](
                                const ROOT::RVec<float> &pts,
                                const ROOT::RVec<float> &etas,
                                const ROOT::RVec<int> &decay_modes,
                                const ROOT::RVec<UChar_t> &gen_matches) {
        ROOT::RVec<float> corrected_pts(pts.size());
        for (int i = 0; i < pts.size(); i++) {
            if (gen_matches.at(i) == 5) {
                if (decay_modes.at(i) == 0) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm0});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                } else if (decay_modes.at(i) == 1) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm1});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                } else if (decay_modes.at(i) == 10) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm10});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                } else if (decay_modes.at(i) == 11) {
                    auto correction_factor = evaluator->evaluate(
                        {pts.at(i), std::abs(etas.at(i)),
                         decay_modes.at(i), static_cast<int>(gen_matches.at(i)),
                         id_algorithm, variation_dm11});
                    corrected_pts[i] = pts.at(i) * correction_factor;
                }
            } else {
                corrected_pts[i] = pts.at(i);
            }
            Logger::get("PtCorrection_genuineTau")
                ->debug("tau pt before {}, tau pt after {}, decaymode {}",
                        pts.at(i), corrected_pts.at(i),
                        decay_modes.at(i));
        }
        return corrected_pts;
    };
    auto df1 = df.Define(corrected_pt, correction_lambda,
                         {pt, eta, decay_mode, gen_match});
    return df1;
}
} // end namespace tau
} // namespace physicsobject
#endif /* GUARD_TAUS_H */