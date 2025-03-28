#ifndef GUARD_TAUS_H
#define GUARD_TAUS_H

#include "../include/RoccoR.hxx"
#include "../include/basefilters.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace physicsobject {
namespace tau {
/// Function to cut on taus based on the tau decay mode
///
/// \param[in] df the input dataframe
/// \param[in] tau_dms name of the column with tau decay modes
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe
/// \param[in] SelectedDecayModes a `std::vector<int>` containing the
/// decay modes, that should pass the cut
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutDecayModes(ROOT::RDF::RNode df, const std::string &maskname,
                               const std::string &tau_dms,
                               const std::vector<int> &SelectedDecayModes) {
    auto df1 = df.Define(
        maskname,
        [SelectedDecayModes](const ROOT::RVec<Int_t> &decaymodes) {
            ROOT::RVec<int> mask;
            for (auto n : decaymodes) {
                mask.push_back(int(std::find(SelectedDecayModes.begin(),
                                             SelectedDecayModes.end(),
                                             n) != SelectedDecayModes.end()));
            }
            return mask;
        },
        {tau_dms});
    return df1;
}
/// Function to cut taus based on the tau ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the dataframe 
/// \param[in] nameID name of the ID column in the NanoAOD 
/// \param[in] idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutTauID(ROOT::RDF::RNode df, const std::string &maskname,
                          const std::string &nameID, const int &idxID) {
    auto df1 = df.Define(maskname, basefilters::Bitmask<UChar_t>(idxID), {nameID});
    return df1;
}
/// Function to correct e to tau fake pt
///
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] df the input dataframe
/// \param[in] correctionManager the correction manager instance
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau)
/// \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_dm0_b scale factor to be applied to taus with decay mode 0 and
/// eta region barrel
/// \param[in] sf_dm1_b scale factor to be applied to taus
/// with decay mode 1 and eta region barrel
/// \param[in] sf_dm0_e scale factor to
/// be applied to taus with decay mode 0 and eta region endcap
/// \param[in] sf_dm1_e scale factor to be applied to taus with decay mode 1 and
/// eta region endcap name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_eleFake(ROOT::RDF::RNode df,
                     correctionManager::CorrectionManager &correctionManager,
                     const std::string &corrected_pt, const std::string &pt,
                     const std::string &eta, const std::string &decayMode,
                     const std::string &genMatch, const std::string &sf_file,
                     const std::string &jsonESname,
                     const std::string &idAlgorithm,
                     const std::string &sf_dm0_b, const std::string &sf_dm1_b,
                     const std::string &sf_dm0_e, const std::string &sf_dm1_e) {
    auto evaluator = correctionManager.loadCorrection(sf_file, jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, sf_dm0_b, sf_dm1_b,
                                     sf_dm0_e, sf_dm1_e](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 1 || genmatch.at(i) == 3) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0 &&
                    std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 0 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("ptcorrection ele fake")
                ->debug("tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct e to tau fake pt
/// WARNING: The function without the CorrectionManager is deprecated and will
/// be removed in the future \param[out] corrected_pt name of the corrected tau
/// pt to be calculated \param[in] df the input dataframe \param[in] pt name of
/// the raw tau pt \param[in] eta name of raw tau eta \param[in] decayMode decay
/// mode of the tau \param[in] genMatch column with genmatch values (from prompt
/// e, prompt mu, tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_dm0_b scale factor to be applied to taus with decay mode 0 and
/// eta region barrel
/// \param[in] sf_dm1_b scale factor to be applied to taus
/// with decay mode 1 and eta region barrel
/// \param[in] sf_dm0_e scale factor to
/// be applied to taus with decay mode 0 and eta region endcap
/// \param[in] sf_dm1_e scale factor to be applied to taus with decay mode 1 and
/// eta region endcap name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_eleFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const std::string &decayMode, const std::string &genMatch,
                     const std::string &sf_file, const std::string &jsonESname,
                     const std::string &idAlgorithm,
                     const std::string &sf_dm0_b, const std::string &sf_dm1_b,
                     const std::string &sf_dm0_e, const std::string &sf_dm1_e) {
    Logger::get("ptcorrection ele fake")
        ->warn("The function  without CorrectionManager is deprecated and will "
               "be removed in the future. Please use the function with "
               "CorrectionManager instead.");
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, sf_dm0_b, sf_dm1_b,
                                     sf_dm0_e, sf_dm1_e](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 1 || genmatch.at(i) == 3) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0 &&
                    std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 0 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("ptcorrection ele fake")
                ->debug("tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct mu to tau fake pt
///
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] correctionManager the correction manager instance
/// \param[in] df the input dataframe
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_es scale factor to be applied to taus faked by muons
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_muFake(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &corrected_pt, const std::string &pt,
                    const std::string &eta, const std::string &decayMode,
                    const std::string &genMatch, const std::string &sf_file,
                    const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &sf_es) {
    auto evaluator = correctionManager.loadCorrection(sf_file, jsonESname);
    auto tau_pt_correction_lambda =
        [evaluator, idAlgorithm, sf_es](const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    // only considering wanted tau decay modes
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_es});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    Logger::get("mu fake")->debug(
                        "tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
                }
            }
            return corrected_pt_values;
        };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct mu to tau fake pt
/// WARNING: The function without the CorrectionManager is deprecated and will
/// be removed in the future \param[out] corrected_pt name of the corrected tau
/// pt to be calculated \param[in] df the input dataframe \param[in] pt name of
/// the raw tau pt \param[in] eta name of raw tau eta \param[in] decayMode decay
/// mode of the tau \param[in] genMatch column with genmatch values (from prompt
/// e, prompt mu, tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_es scale factor to be applied to taus faked by muons
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_muFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &sf_es) {
    Logger::get("ptcorrection mu fake")
        ->warn("The function without CorrectionManager is deprecated and will "
               "be removed in the future. Please use the function with "
               "CorrectionManager instead.");
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda =
        [evaluator, idAlgorithm, sf_es](const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    // only considering wanted tau decay modes
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_es});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    Logger::get("mu fake")->debug(
                        "tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
                }
            }
            return corrected_pt_values;
        };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct tau pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt
/// \param[in] decayMode decay mode of the tau
/// \param[in] sf_dm0 scale factor to be applied to taus with decay mode 0
/// \param[in] sf_dm1 scale factor to be applied to other 1 prong taus
/// \param[in] sf_dm10 scale factor to be applied to taus with decay mode 10
/// \param[in] sf_dm11 scale factor to be applied to other 3 prong taus
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &decayMode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11) {
    auto tau_pt_correction_lambda =
        [sf_dm0, sf_dm1, sf_dm10, sf_dm11](const ROOT::RVec<float> &pt_values,
                                           const ROOT::RVec<int> &decay_modes) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (decay_modes.at(i) == 0)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm0;
                else if (decay_modes.at(i) > 0 && decay_modes.at(i) < 5)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm1;
                else if (decay_modes.at(i) == 10)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm10;
                else if (decay_modes.at(i) > 10 && decay_modes.at(i) < 15)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm11;
                else
                    corrected_pt_values[i] = pt_values.at(i);
            }
            return corrected_pt_values;
        };
    auto df1 =
        df.Define(corrected_pt, tau_pt_correction_lambda, {pt, decayMode});
    return df1;
}
/// Function to correct tau pt with correctionlib
///
/// \param[in] df the input dataframe
/// \param[in] correctionManager the correction manager instance
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] DM0 variation decay mode 0
/// \param[in] DM1 variation decay mode 1
/// \param[in] DM10 variation decay mode 10
/// \param[in] DM11 variation decay mode 11
/// DM values: "nom","up","down"
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_genTau(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &corrected_pt, const std::string &pt,
                    const std::string &eta, const std::string &decayMode,
                    const std::string &genMatch, const std::string &sf_file,
                    const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &DM0,
                    const std::string &DM1, const std::string &DM10,
                    const std::string &DM11) {
    auto evaluator = correctionManager.loadCorrection(sf_file, jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, DM0, DM1, DM10,
                                     DM11](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 5) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM0});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM1});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 10) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM10});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 11) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM11});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("tauEnergyCorrection")
                ->debug("tau pt before {}, tau pt after {}, decaymode {}",
                        pt_values.at(i), corrected_pt_values.at(i),
                        decay_modes.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct tau pt with correctionlib
/// WARNING: The function without the CorrectionManager is deprecated and will
/// be removed in the future \param[in] df the input dataframe \param[out]
/// corrected_pt name of the corrected tau pt to be calculated \param[in] pt
/// name of the raw tau pt \param[in] eta name of raw tau eta \param[in]
/// decayMode decay mode of the tau \param[in] genMatch column with genmatch
/// values (from prompt e, prompt mu, tau->e, tau->mu, had. tau) \param[in]
/// sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] DM0 variation decay mode 0
/// \param[in] DM1 variation decay mode 1
/// \param[in] DM10 variation decay mode 10
/// \param[in] DM11 variation decay mode 11
/// DM values: "nom","up","down"
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_genTau(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &DM0,
                    const std::string &DM1, const std::string &DM10,
                    const std::string &DM11) {
    Logger::get("ptcorrection gen tau")
        ->warn("The function without CorrectionManager is deprecated and will "
               "be removed in the future. Please use the function with "
               "CorrectionManager instead.");
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, DM0, DM1, DM10,
                                     DM11](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 5) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM0});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM1});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 10) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM10});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 11) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM11});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("tauEnergyCorrection")
                ->debug("tau pt before {}, tau pt after {}, decaymode {}",
                        pt_values.at(i), corrected_pt_values.at(i),
                        decay_modes.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
} // end namespace tau
} // namespace physicsobject
#endif /* GUARD_TAUS_H */