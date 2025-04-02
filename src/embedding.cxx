#ifndef GUARD_EMBEDDING_H
#define GUARD_EMBEDDING_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include <nlohmann/json.hpp>

namespace embedding {
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
} // namespace electron
} // namespace embedding

#endif /* GUARD_EMBEDDING_H */