#ifndef GUARD_ELECTRONS_H
#define GUARD_ELECTRONS_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "correction.h"

namespace physicsobject {
namespace electron {

/**
 * @brief This function calculates uncertainties for the electron energy scale
 * and resolution corrections that are already applied in nanoAOD (Run2 UL). The
 * resolution uncertainty is taken for dedicated variation branches in the
 * nanoAODs. The same cannot be done for the scale uncertainty due to a bug in
 * the nanoAOD (Run2 UL) production. Therefore, a patch is used based on a
 * correctionlib json file. This procedure is recommended by EGM POG and
 * described in
 *
 * https://cms-talk.web.cern.ch/t/electron-scale-smear-variables-in-nanoaod/20210
 *
 * and https://twiki.cern.ch/twiki/bin/view/CMS/EgammaSFJSON
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction scale uncertainty patch file
 * @param outputname name of the output column for corrected \f$p_T\f$ values
 * @param pt name of the column containing electron \f$p_T\f$ values
 * @param eta name of the column containing electron pseudorapidities
 * @param gain name of the column containing electron gain values
 * @param es_resolution_up name of the column containing the one sigma upward
 * energy smearing uncertainties
 * @param es_resolution_down name of the column containing the one sigma
 * downward energy smearing uncertainties
 * @param es_file path to the correction file for the energy scale uncertainties
 * @param era data-taking period of Run2, possible options are "2018", "2017",
 * "2016postVFP", "2016preVFP"
 * @param variation name of the energy correction variation that should be
 * calculated (e.g., "resolutionUp", "resolutionDown", "scaleUp", "scaleDown"),
 * for "nominal" nothing is done because energy correction is already applied
 *
 * @return a dataframe containing the varied electron transverse momenta
 *
 * @note TODO: This function should be extended to support Run3
 */
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &pt,
               const std::string &eta, const std::string &gain,
               const std::string &es_resolution_up,
               const std::string &es_resolution_down,
               const std::string &es_file, const std::string &era,
               const std::string &variation) {
    auto evaluator =
        correction_manager.loadCorrection(es_file, "UL-EGM_ScaleUnc");
    auto electron_pt_correction_lambda =
        [evaluator, era, variation](const ROOT::RVec<float> &pts,
                                    const ROOT::RVec<float> &etas,
                                    const ROOT::RVec<UChar_t> &gains,
                                    const ROOT::RVec<float> &es_reso_up,
                                    const ROOT::RVec<float> &es_reso_down) {
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (variation == "resolutionUp") {
                    auto dpt = es_reso_up.at(i) / std::cosh(etas.at(i));
                    corrected_pts[i] = pts.at(i) + dpt;
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, dpt {}",
                                pts.at(i), corrected_pts.at(i), dpt);
                } else if (variation == "resolutionDown") {
                    auto dpt = es_reso_down.at(i) / std::cosh(etas.at(i));
                    corrected_pts[i] = pts.at(i) + dpt;
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, dpt {}",
                                pts.at(i), corrected_pts.at(i), dpt);
                } else if (variation == "scaleUp") {
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("inputs: era {}, eta {}, gain {}", era,
                                etas.at(i), static_cast<int>(gains.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaleup", etas.at(i),
                                             static_cast<int>(gains.at(i))});
                    corrected_pts[i] = pts.at(i) * sf;
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, sf {}",
                                pts.at(i), corrected_pts.at(i), sf);
                } else if (variation == "scaleDown") {
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("inputs: era {}, eta {}, gain {}", era,
                                etas.at(i), static_cast<int>(gains.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaledown", etas.at(i),
                                             static_cast<int>(gains.at(i))});
                    corrected_pts[i] = pts.at(i) * sf;
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, sf {}",
                                pts.at(i), corrected_pts.at(i), sf);
                } else {
                    corrected_pts[i] = pts.at(i);
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}", pts.at(i),
                                corrected_pts.at(i));
                }
            }
            return corrected_pts;
        };
    auto df1 = df.Define(outputname, electron_pt_correction_lambda,
                         {pt, eta, gain, es_resolution_up, es_resolution_down});
    return df1;
}

/**
 * @brief This function defines a boolean mask that identifies electrons
 * falling within the ECAL barrel-endcap transition gap. The mask is `true`
 * for electrons outside the gap and `false` for those inside.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the veto mask
 * @param eta name of the column containing electron pseudorapidities
 * @param delta_eta_sc name of the column containing the distance in
 * pseudorapidity between supercluster and electron
 * @param end_ecal_barrel end of the ECAL barrel region (in pseudorapidity)
 * @param start_ecal_endcap begin of the ECAL endcap region (in pseudorapidity)
 *
 * @return a dataframe containing the new mask as a column
 */
ROOT::RDF::RNode VetoECALGap(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &eta,
                             const std::string &delta_eta_sc,
                             const float &end_ecal_barrel,
                             const float &start_ecal_endcap) {
    auto lambda = [end_ecal_barrel,
                   start_ecal_endcap](const ROOT::RVec<float> &eta,
                                      const ROOT::RVec<float> &delta_eta_sc) {
        ROOT::RVec<int> mask = (abs(eta + delta_eta_sc) < end_ecal_barrel) ||
                               (abs(eta + delta_eta_sc) >= start_ecal_endcap);
        return mask;
    };

    auto df1 = df.Define(outputname, lambda, {eta, delta_eta_sc});
    return df1;
}

/**
 * @brief This function creates a selection mask based on the transverse
 * impact parameter `dxy` and longitudinal impact parameter `dz` of electrons.
 * The selection criteria differ between the barrel and endcap regions of the
 * detector.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the selection mask
 * @param eta name of the column containing electron pseudorapidities
 * @param delta_eta_sc name of the column containing the distance in
 * pseudorapidity between supercluster and electron
 * @param dxy name of the column containing transverse impact parameter
 * @param dz name of the column containing longitudinal impact parameter
 * @param ecal_barrel_endcap_boundary absolute eta boundary separating the
 * barrel and endcap regions
 * @param max_dxy_barrel maximal dxy in the ECAL barrel region
 * @param max_dz_barrel maximal dz in the ECAL barrel region
 * @param max_dxy_endcap maximal dxy in the ECAL endcap region
 * @param max_dz_endcap maximal dz in the ECAL endcap region
 *
 * @return a dataframe containing the new mask as a column
 */
ROOT::RDF::RNode
CutInteractionPoint(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &eta, const std::string &delta_eta_sc,
                    const std::string &dxy, const std::string &dz,
                    const float &ecal_barrel_endcap_boundary,
                    const float &max_dxy_barrel, const float &max_dz_barrel,
                    const float &max_dxy_endcap, const float &max_dz_endcap) {
    auto lambda = [ecal_barrel_endcap_boundary, max_dxy_barrel, max_dz_barrel,
                   max_dxy_endcap,
                   max_dz_endcap](const ROOT::RVec<float> &eta,
                                  const ROOT::RVec<float> &delta_eta_sc,
                                  const ROOT::RVec<float> &dxy,
                                  const ROOT::RVec<float> &dz) {
        ROOT::RVec<int> mask =
            (((abs(eta + delta_eta_sc) < ecal_barrel_endcap_boundary) &&
              (dxy < max_dxy_barrel) && (dz < max_dz_barrel)) ||
             ((abs(eta + delta_eta_sc) >= ecal_barrel_endcap_boundary) &&
              (dxy < max_dxy_endcap) && (dz < max_dz_endcap)));
        return mask;
    };

    auto df1 = df.Define(outputname, lambda, {eta, delta_eta_sc, dxy, dz});
    return df1;
}

namespace scalefactor {

/**
 * @brief This function calculates electron ID scale factors (SFs) for a single
 * electron dependening on its pseudorapidity (\f$\eta\f$) and transverse
 * momentum
 * (\f$p_T\f$). The scale factors are loaded from a correctionlib file using a
 * specified scale factor name and variation.
 *
 * Recommendations by EgammaPOG:
 * - [Run2](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018)
 * - [Run3](https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * electron scale factor file
 * @param outputname name of the output column containing the ID scale factor
 * @param pt name of the column containing the transverse momentum of an
 * electron
 * @param eta name of the column containing the pseudorapidity of an electron
 * @param era string with the era name of a data taking period, e.g.
 * "2016preVFP"
 * @param wp working point of the electron id that should be used, e.g.
 * "Medium", "wp90noiso", ...
 * @param sf_file path to the file with the electron scale factors
 * @param sf_name name of the electron scale factor for the ID correction,
 * e.g. "UL-Electron-ID-SF"
 * @param variation name the scale factor variation, "sf" for the nominal
 * scale factor and "sfup"/"sfdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Id(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &outputname, const std::string &pt,
                    const std::string &eta, const std::string &era,
                    const std::string &wp, const std::string &sf_file,
                    const std::string &sf_name, const std::string &variation) {
    Logger::get("physicsobject::electron::scalefactor::Id")
        ->debug("Setting up functions for electron id sf with correctionlib");
    Logger::get("physicsobject::electron::scalefactor::Id")
        ->debug("ID - Name {}", sf_name);
    auto evaluator = correctionManager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, era, sf_name, wp, variation](const float &pt,
                                                 const float &eta) {
            Logger::get("physicsobject::electron::scalefactor::Id")
                ->debug("Era {}, Variation {}, WP {}", era, variation, wp);
            Logger::get("physicsobject::electron::scalefactor::Id")
                ->debug("ID - pt {}, eta {}", pt, eta);
            double sf = 1.;
            if (pt >= 0.0) {
                sf = evaluator->evaluate({era, variation, wp, eta, pt});
            }
            Logger::get("physicsobject::electron::scalefactor::Id")
                ->debug("Scale Factor {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
} // namespace scalefactor
} // namespace electron
} // namespace physicsobject
#endif /* GUARD_ELECTRONS_H */