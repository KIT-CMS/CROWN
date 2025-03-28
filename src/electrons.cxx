#ifndef GUARD_ELECTRONS_H
#define GUARD_ELECTRONS_H

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
namespace electron {
/// Function to correct electron pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected electron pt to be calculated
/// \param[in] pt name of the raw electron pt
/// \param[in] eta the name of the raw electron eta
/// \param[in] sf_barrel scale factor to be applied to electrons in the barrel
/// \param[in] sf_endcap scale factor to be applied to electrons in the endcap
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap) {
    auto electron_pt_correction_lambda =
        [sf_barrel, sf_endcap](const ROOT::RVec<float> &pt_values,
                               const ROOT::RVec<float> &eta) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (abs(eta.at(i)) <= 1.479) {
                    corrected_pt_values[i] = pt_values.at(i) * sf_barrel;
                } else if (abs(eta.at(i)) > 1.479) {
                    corrected_pt_values[i] = pt_values.at(i) * sf_endcap;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
            }
            return corrected_pt_values;
        };
    auto df1 =
        df.Define(corrected_pt, electron_pt_correction_lambda, {pt, eta});
    return df1;
}
/// Function to correct electron pt, based on correctionlib file
///
/// \param[in] df the input dataframe
/// \param[in] correctionManager the correction manager instance
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] sf_barrel scale factor to be applied to electrons in the barrel
/// \param[in] sf_endcap scale factor to be applied to electrons in the endcap
/// \param[in] sf_file:
/// \param[in] jsonESname name of the tau energy correction in the json file
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correctionManager,
             const std::string &corrected_pt, const std::string &pt,
             const std::string &eta, const std::string &sf_barrel,
             const std::string &sf_endcap, const std::string &sf_file,
             const std::string &jsonESname) {
    auto evaluator = correctionManager.loadCorrection(sf_file, jsonESname);
    auto electron_pt_correction_lambda = [evaluator, sf_barrel, sf_endcap](
                                             const ROOT::RVec<float> &pt_values,
                                             const ROOT::RVec<float> &eta) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (abs(eta.at(i)) <= 1.479) {
                auto sf = evaluator->evaluate({"barrel", sf_barrel});
                corrected_pt_values[i] = pt_values.at(i) * sf;
                Logger::get("eleEnergyCorrection")
                    ->debug("barrel: ele pt before {}, ele pt after {}, sf {}",
                            pt_values.at(i), corrected_pt_values.at(i), sf);
            } else if (abs(eta.at(i)) > 1.479) {
                auto sf = evaluator->evaluate({"endcap", sf_endcap});
                corrected_pt_values[i] = pt_values.at(i) * sf;
                Logger::get("eleEnergyCorrection")
                    ->debug("endcap: ele pt before {}, ele pt after {}, sf {}",
                            pt_values.at(i), corrected_pt_values.at(i), sf);
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
        }
        return corrected_pt_values;
    };
    auto df1 =
        df.Define(corrected_pt, electron_pt_correction_lambda, {pt, eta});
    return df1;
}
/// Function to correct electron pt, based on correctionlib file
/// WARNING: The function without the CorrectionManager is deprecated and will
/// be removed in the future \param[in] df the input dataframe \param[out]
/// corrected_pt name of the corrected tau pt to be calculated \param[in] pt
/// name of the raw tau pt \param[in] eta name of raw tau eta \param[in]
/// sf_barrel scale factor to be applied to electrons in the barrel \param[in]
/// sf_endcap scale factor to be applied to electrons in the endcap \param[in]
/// sf_file: \param[in] jsonESname name of the tau energy correction in the json
/// file
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df, const std::string &corrected_pt,
             const std::string &pt, const std::string &eta,
             const std::string &sf_barrel, const std::string &sf_endcap,
             const std::string &sf_file, const std::string &jsonESname) {
    Logger::get("eleEnergyCorrection")
        ->warn("The function without CorrectionManager is deprecated and will "
               "be removed in the future. Please use the function with "
               "CorrectionManager instead.");
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto electron_pt_correction_lambda = [evaluator, sf_barrel, sf_endcap](
                                             const ROOT::RVec<float> &pt_values,
                                             const ROOT::RVec<float> &eta) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (abs(eta.at(i)) <= 1.479) {
                auto sf = evaluator->evaluate({"barrel", sf_barrel});
                corrected_pt_values[i] = pt_values.at(i) * sf;
                Logger::get("eleEnergyCorrection")
                    ->debug("barrel: ele pt before {}, ele pt after {}, sf {}",
                            pt_values.at(i), corrected_pt_values.at(i), sf);
            } else if (abs(eta.at(i)) > 1.479) {
                auto sf = evaluator->evaluate({"endcap", sf_endcap});
                corrected_pt_values[i] = pt_values.at(i) * sf;
                Logger::get("eleEnergyCorrection")
                    ->debug("endcap: ele pt before {}, ele pt after {}, sf {}",
                            pt_values.at(i), corrected_pt_values.at(i), sf);
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
        }
        return corrected_pt_values;
    };
    auto df1 =
        df.Define(corrected_pt, electron_pt_correction_lambda, {pt, eta});
    return df1;
}
/// Function to calculate uncertainties for electron pt correction. The electron
/// energy correction is already applied in nanoAOD and in general there are
/// branches in nanoAOD with the energy shifts for scale and resolution, but due
/// to a bug the scale shifts are all 0 and have to be calculated from a json
/// file. Information taken from
/// https://cms-talk.web.cern.ch/t/electron-scale-smear-variables-in-nanoaod/20210
/// and https://twiki.cern.ch/twiki/bin/view/CMS/EgammaSFJSON
///
/// \param[in] df the input dataframe
/// \param[in] correctionManager the correction manager instance
/// \param[out] corrected_pt name of the corrected electron pt to be calculated
/// \param[in] pt name of the electron pt
/// \param[in] eta name of electron eta
/// \param[in] gain name of electron seedGain
/// \param[in] ES_sigma_up name of electron energy smearing value 1 sigma up
/// shifted \param[in] ES_sigma_down name of electron energy smearing value 1
/// sigma down shifted \param[in] era era of the electron measurement e.g.
/// "2018" \param[in] variation name of the variation to be calculated (nominal
/// correction is already applied) \param[in] ES_file name of the json file with
/// the energy scale uncertainties
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correctionManager,
               const std::string &corrected_pt, const std::string &pt,
               const std::string &eta, const std::string &gain,
               const std::string &ES_sigma_up, const std::string &ES_sigma_down,
               const std::string &era, const std::string &variation,
               const std::string &ES_file) {
    auto evaluator =
        correctionManager.loadCorrection(ES_file, "UL-EGM_ScaleUnc");
    auto electron_pt_correction_lambda =
        [evaluator, era, variation](const ROOT::RVec<float> &pt_values,
                                    const ROOT::RVec<float> &eta,
                                    const ROOT::RVec<UChar_t> &gain,
                                    const ROOT::RVec<float> &ES_sigma_up,
                                    const ROOT::RVec<float> &ES_sigma_down) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (variation == "resolutionUp") {
                    auto dpt = ES_sigma_up.at(i) / std::cosh(eta.at(i));
                    corrected_pt_values[i] = pt_values.at(i) + dpt;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, dpt {}",
                                pt_values.at(i), corrected_pt_values.at(i),
                                dpt);
                } else if (variation == "resolutionDown") {
                    auto dpt = ES_sigma_down.at(i) / std::cosh(eta.at(i));
                    corrected_pt_values[i] = pt_values.at(i) + dpt;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, dpt {}",
                                pt_values.at(i), corrected_pt_values.at(i),
                                dpt);
                } else if (variation == "scaleUp") {
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("inputs: era {}, eta {}, gain {}", era,
                                eta.at(i), static_cast<int>(gain.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaleup", eta.at(i),
                                             static_cast<int>(gain.at(i))});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, sf {}",
                                pt_values.at(i), corrected_pt_values.at(i), sf);
                } else if (variation == "scaleDown") {
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("inputs: era {}, eta {}, gain {}", era,
                                eta.at(i), static_cast<int>(gain.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaledown", eta.at(i),
                                             static_cast<int>(gain.at(i))});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, sf {}",
                                pt_values.at(i), corrected_pt_values.at(i), sf);
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}",
                                pt_values.at(i), corrected_pt_values.at(i));
                }
            }
            return corrected_pt_values;
        };
    auto df1 = df.Define(corrected_pt, electron_pt_correction_lambda,
                         {pt, eta, gain, ES_sigma_up, ES_sigma_down});
    return df1;
}
/// Function to calculate uncertainties for electron pt correction. The electron
/// energy correction is already applied in nanoAOD and in general there are
/// branches in nanoAOD with the energy shifts for scale and resolution, but due
/// to a bug the scale shifts are all 0 and have to be calculated from a json
/// file. Information taken from
/// https://cms-talk.web.cern.ch/t/electron-scale-smear-variables-in-nanoaod/20210
/// and https://twiki.cern.ch/twiki/bin/view/CMS/EgammaSFJSON
/// WARNING: The function without the CorrectionManager is deprecated and will
/// be removed in the future
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected electron pt to be calculated
/// \param[in] pt name of the electron pt
/// \param[in] eta name of electron eta
/// \param[in] gain name of electron seedGain
/// \param[in] ES_sigma_up name of electron energy smearing value 1 sigma up
/// shifted \param[in] ES_sigma_down name of electron energy smearing value 1
/// sigma down shifted \param[in] era era of the electron measurement e.g.
/// "2018" \param[in] variation name of the variation to be calculated (nominal
/// correction is already applied) \param[in] ES_file name of the json file with
/// the energy scale uncertainties
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df, const std::string &corrected_pt,
               const std::string &pt, const std::string &eta,
               const std::string &gain, const std::string &ES_sigma_up,
               const std::string &ES_sigma_down, const std::string &era,
               const std::string &variation, const std::string &ES_file) {
    Logger::get("ElectronPtCorrectionMC")
        ->warn("The function without CorrectionManager is deprecated and will "
               "be removed in the future. Please use the function with "
               "CorrectionManager instead.");
    auto evaluator =
        correction::CorrectionSet::from_file(ES_file)->at("UL-EGM_ScaleUnc");
    auto electron_pt_correction_lambda =
        [evaluator, era, variation](const ROOT::RVec<float> &pt_values,
                                    const ROOT::RVec<float> &eta,
                                    const ROOT::RVec<UChar_t> &gain,
                                    const ROOT::RVec<float> &ES_sigma_up,
                                    const ROOT::RVec<float> &ES_sigma_down) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (variation == "resolutionUp") {
                    auto dpt = ES_sigma_up.at(i) / std::cosh(eta.at(i));
                    corrected_pt_values[i] = pt_values.at(i) + dpt;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, dpt {}",
                                pt_values.at(i), corrected_pt_values.at(i),
                                dpt);
                } else if (variation == "resolutionDown") {
                    auto dpt = ES_sigma_down.at(i) / std::cosh(eta.at(i));
                    corrected_pt_values[i] = pt_values.at(i) + dpt;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, dpt {}",
                                pt_values.at(i), corrected_pt_values.at(i),
                                dpt);
                } else if (variation == "scaleUp") {
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("inputs: era {}, eta {}, gain {}", era,
                                eta.at(i), static_cast<int>(gain.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaleup", eta.at(i),
                                             static_cast<int>(gain.at(i))});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, sf {}",
                                pt_values.at(i), corrected_pt_values.at(i), sf);
                } else if (variation == "scaleDown") {
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("inputs: era {}, eta {}, gain {}", era,
                                eta.at(i), static_cast<int>(gain.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaledown", eta.at(i),
                                             static_cast<int>(gain.at(i))});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, sf {}",
                                pt_values.at(i), corrected_pt_values.at(i), sf);
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                    Logger::get("ElectronPtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}",
                                pt_values.at(i), corrected_pt_values.at(i));
                }
            }
            return corrected_pt_values;
        };
    auto df1 = df.Define(corrected_pt, electron_pt_correction_lambda,
                         {pt, eta, gain, ES_sigma_up, ES_sigma_down});
    return df1;
}
/// Function to cut electrons based on the electron MVA ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
/// Function to cut electrons based on the cut based electron ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe
/// \param[in] nameID name of the ID column in the NanoAOD
/// \param[in] IDvalue value of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const int &IDvalue) {
    auto df1 =
        df.Define(maskname, basefilters::Min<int>(IDvalue), {nameID});
    return df1;
}
/// Function to cut electrons based on failing the cut based electron ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe
/// \param[in] nameID name of the ID column in the NanoAOD
/// \param[in] IDvalue value of the WP the has to be failed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode AntiCutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                             const std::string &nameID, const int &IDvalue) {
    auto df1 =
        df.Define(maskname, basefilters::Max<int>(IDvalue), {nameID});
    return df1;
}

/// Function to cut electrons based on the electron isolation using
/// basefilters::Max
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to
/// the dataframe
/// \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold) {
    auto df1 = df.Define(maskname, basefilters::Max<float>(Threshold),
                         {isolationName});
    return df1;
}
/// Function to select electrons below an Dxy and Dz threshold, based on the
/// electrons supercluster
///
/// \param[in] df the input dataframe
/// \param[in] eta quantity name of the electron eta column in the NanoAOD
/// \param[in] detasc quantity name of the electron deltaEtaSC column in the
/// NanoAOD
/// \param[in] dxy quantity name of the Dxy column in the NanoAOD
/// \param[in] dz quantity name of the Dz column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe
/// \param[in] abseta_eb_ee abs(eta) of the EB-EE transition
/// \param[in] max_dxy_eb Threshold maximal Dxy value in the barrel
/// \param[in] max_dz_eb Threshold maximal Dz value in the barrel
/// \param[in] max_dxy_ee hreshold maximal Dxy value in the endcap
/// \param[in] max_dz_ee Threshold maximal Dz value in the endcap
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIP(ROOT::RDF::RNode df, const std::string &eta,
                       const std::string &detasc, const std::string &dxy,
                       const std::string &dz, const std::string &maskname,
                       const float &abseta_eb_ee, const float &max_dxy_eb,
                       const float &max_dz_eb, const float &max_dxy_ee,
                       const float &max_dz_ee) {
    auto lambda = [abseta_eb_ee, max_dxy_eb, max_dz_eb, max_dxy_ee,
                   max_dz_ee](const ROOT::RVec<float> &eta,
                              const ROOT::RVec<float> &detasc,
                              const ROOT::RVec<float> &dxy,
                              const ROOT::RVec<float> &dz) {
        ROOT::RVec<int> mask = (((abs(eta + detasc) < abseta_eb_ee) &&
                                 (dxy < max_dxy_eb) && (dz < max_dz_eb)) ||
                                ((abs(eta + detasc) >= abseta_eb_ee) &&
                                 (dxy < max_dxy_ee) && (dz < max_dz_ee)));
        return mask;
    };

    auto df1 = df.Define(maskname, lambda, {eta, detasc, dxy, dz});
    return df1;
}

/// Function to veto electrons in the transition region of EB and EE, based on
/// the electrons supercluster
///
/// \param[in] df the input dataframe
/// \param[in] eta quantity name of the electron eta column in the NanoAOD
/// \param[in] detasc quantity name of the electron deltaEtaSC column in the
/// NanoAOD
/// \param[out] maskname the name of the mask to be added as column to
/// the dataframe
/// \param[in] end_eb abs(eta) of the beginning of the transition
/// region
///\param[in] start_ee abs(eta) of the end of the transition region
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutGap(ROOT::RDF::RNode df, const std::string &eta,
                        const std::string &detasc, const std::string &maskname,
                        const float &end_eb, const float &start_ee) {
    auto lambda = [end_eb, start_ee](const ROOT::RVec<float> &eta,
                                     const ROOT::RVec<float> &detasc) {
        ROOT::RVec<int> mask =
            (abs(eta + detasc) < end_eb) || (abs(eta + detasc) > start_ee);
        return mask;
    };

    auto df1 = df.Define(maskname, lambda, {eta, detasc});
    return df1;
}

} // namespace electron
} // namespace physicsobject
#endif /* GUARD_ELECTRONS_H */