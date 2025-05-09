#ifndef GUARD_EMBEDDING_H
#define GUARD_EMBEDDING_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include <nlohmann/json.hpp>

namespace embedding {
namespace electron {

/**
 * @brief This function scales the \f$p_T\f$ values of electrons depending on 
 * whether they fall within the barrel or endcap region of the detector. The 
 * correction factors are provided as input parameters.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the corrected \f$p_T\f$ 
 * values
 * @param pt name of the input column containing \f$p_T\f$ values
 * @param eta name of the column containing eta values
 * @param sf_barrel scale factor applied in the barrel region (|eta| <= 1.479)
 * @param sf_endcap scale factor applied in the endcap region (|eta| > 1.479)
 *
 * @return a dataframe with the corrected \f$p_T\f$ values
 */
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap) {
    auto correction_lambda =
        [sf_barrel, sf_endcap](const ROOT::RVec<float> &pts,
                               const ROOT::RVec<float> &etas) {
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (abs(etas.at(i)) <= 1.479)
                    corrected_pts[i] = pts.at(i) * sf_barrel;
                else if (abs(etas.at(i)) > 1.479)
                    corrected_pts[i] = pts.at(i) * sf_endcap;
                else
                    corrected_pts[i] = pts.at(i);
                
                Logger::get("embedding::electron::PtCorrection_byValue")
                ->debug("ele pt before {}, ele pt after {}",
                        pts.at(i), corrected_pts.at(i));
            }
            return corrected_pts;
        };
    auto df1 =
        df.Define(outputname, correction_lambda, {pt, eta});
    return df1;
}

/**
 * @brief This function scales the \f$p_T\f$ values of electrons based on their 
 * pseudorapidity (`eta`), applying different correction factors for 
 * the barrel and endcap regions of the detector.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the 
 * correction file
 * @param outputname name of the output column storing the corrected 
 * electron \f$p_T\f$pT values
 * @param pt name of the input column containing electron \f$p_T\f$ values
 * @param eta name of the column containing electron eta values
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param variation_barrel variation for the barrel region, options are 
 * "nom", "up", "down"
 * @param variation_endcap variation for the endcap region, options are 
 * "nom", "up", "down"
 *
 * @return a dataframe with the corrected \f$p_T\f$ values
 */
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correction_manager,
             const std::string &outputname, const std::string &pt,
             const std::string &eta, const std::string &es_file,
             const std::string &correction_name, 
             const std::string &variation_barrel,
             const std::string &variation_endcap) {
    auto evaluator = correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda = [evaluator, variation_barrel, variation_endcap](
                                const ROOT::RVec<float> &pts,
                                const ROOT::RVec<float> &etas) {
        ROOT::RVec<float> corrected_pts(pts.size());
        for (int i = 0; i < pts.size(); i++) {
            if (abs(etas.at(i)) <= 1.479) {
                auto correction_factor = evaluator->evaluate({"barrel", variation_barrel});
                corrected_pts[i] = pts.at(i) * correction_factor;
                Logger::get("embedding::electron::PtCorrection")
                    ->debug("barrel: ele pt before {}, ele pt after {}, factor {}",
                            pts.at(i), corrected_pts.at(i), correction_factor);
            } else if (abs(etas.at(i)) > 1.479) {
                auto correction_factor = evaluator->evaluate({"endcap", variation_endcap});
                corrected_pts[i] = pts.at(i) * correction_factor;
                Logger::get("embedding::electron::PtCorrection")
                    ->debug("endcap: ele pt before {}, ele pt after {}, factor {}",
                            pts.at(i), corrected_pts.at(i), correction_factor);
            } else {
                corrected_pts[i] = pts.at(i);
            }
        }
        return corrected_pts;
    };
    auto df1 =
        df.Define(outputname, correction_lambda, {pt, eta});
    return df1;
}
} // namespace electron

namespace tau {

/**
 * @brief This function scales the \f$p_T\f$ values of hadronic taus based on their 
 * decay mode. The correction factors are provided as input parameters.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the corrected \f$p_T\f$ values
 * @param pt name of the input column containing \f$p_T\f$ values
 * @param decay_mode name of the column containing the tau decay modes
 * @param sf_dm0 scale factor applied for decay mode 0
 * @param sf_dm1 scale factor applied for decay mode 1
 * @param sf_dm10 scale factor applied for decay mode 10
 * @param sf_dm11 scale factor applied for decay mode 11
 *
 * @return a dataframe with the corrected \f$p_T\f$ values
 */
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pt, const std::string &decay_mode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11) {
    auto correction_lambda =
        [sf_dm0, sf_dm1, sf_dm10, sf_dm11](const ROOT::RVec<float> &pts,
                                           const ROOT::RVec<int> &decay_modes) {
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (decay_modes.at(i) == 0)
                    corrected_pts[i] = pts.at(i) * sf_dm0;
                else if (decay_modes.at(i) == 1)
                    corrected_pts[i] = pts.at(i) * sf_dm1;
                else if (decay_modes.at(i) == 10)
                    corrected_pts[i] = pts.at(i) * sf_dm10;
                else if (decay_modes.at(i) == 11)
                    corrected_pts[i] = pts.at(i) * sf_dm11;
                else
                    corrected_pts[i] = pts.at(i);
                
                Logger::get("embedding::tau::PtCorrection_byValue")
                ->debug("tau pt before {}, tau pt after {}",
                        pts.at(i), corrected_pts.at(i));
            }
            return corrected_pts;
        };
    auto df1 =
        df.Define(outputname, correction_lambda, {pt, decay_mode});
    return df1;
}
} // namespace tau
} // namespace embedding

#endif /* GUARD_EMBEDDING_H */