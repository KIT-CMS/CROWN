#ifndef GUARD_ELECTRONS_H
#define GUARD_ELECTRONS_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TRandom3.h"
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
 * @param seed_gain name of the column containing electron gain values
 * @param es_resolution_up name of the column containing the one sigma upward
 * energy smearing uncertainties
 * @param es_resolution_down name of the column containing the one sigma
 * downward energy smearing uncertainties
 * @param es_file path to the correction file for the energy scale uncertainties
 * @param es_name name of the correction in the `es_file`
 * @param era data-taking period of Run2, possible options are "2018", "2017",
 * "2016postVFP", "2016preVFP"
 * @param variation name of the energy correction variation that should be
 * calculated (e.g., "resolutionUp", "resolutionDown", "scaleUp", "scaleDown"),
 * for "nominal" nothing is done because energy correction is already applied
 *
 * @return a dataframe containing the varied electron transverse momenta
 *
 * @note This function is intended for analyses working with Run 2 NanoAODv9
 * samples. For the corresponding function that can be used with
 * Run 3 NanoAODv12, look at the overloaded version of this function.
 */
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &pt,
               const std::string &eta, const std::string &seed_gain,
               const std::string &es_resolution_up,
               const std::string &es_resolution_down,
               const std::string &es_file,
               const std::string &es_name,
               const std::string &era,
               const std::string &variation) {
    auto evaluator =
        correction_manager.loadCorrection(es_file, es_name);
    auto electron_pt_correction_lambda =
        [evaluator, era, variation](const ROOT::RVec<float> &pts,
                                    const ROOT::RVec<float> &etas,
                                    const ROOT::RVec<UChar_t> &seed_gains,
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
                        ->debug("inputs: era {}, eta {}, seed gain {}", era,
                                etas.at(i), static_cast<int>(seed_gains.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaleup", etas.at(i),
                                             static_cast<int>(seed_gains.at(i))});
                    corrected_pts[i] = pts.at(i) * sf;
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("ele pt before {}, ele pt after {}, sf {}",
                                pts.at(i), corrected_pts.at(i), sf);
                } else if (variation == "scaleDown") {
                    Logger::get("physicsobject::electron::PtCorrectionMC")
                        ->debug("inputs: era {}, eta {}, seed gain {}", era,
                                etas.at(i), static_cast<int>(seed_gains.at(i)));
                    auto sf =
                        evaluator->evaluate({era, "scaledown", etas.at(i),
                                             static_cast<int>(seed_gains.at(i))});
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
                         {pt, eta, seed_gain, es_resolution_up, es_resolution_down});
    return df1;
}

/**
 * @brief This function applies energy scale and resolution corrections to MC.
 * The corrections are obtained from a dedicated correctionlib file.
 * 
 * For Run 3 samples, the electron energy scale correction has to be evaluated
 * using a centrally provided correctionlib file. The documentation of the file
 * content can be found here:
 * 
 * - [2022preEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2022_Summer22_electronSS_EtDependent.html)
 * - [2022postEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2022_Summer22EE_electronSS_EtDependent.html)
 * - [2023preBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2023_Summer23_electronSS_EtDependent.html)
 * - [2023postBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2023_Summer23BPix_electronSS_EtDependent.html)
 *
 * An implementation recipe is provided here:
 * [egmScaleAndSmearingExample.py](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/egmScaleAndSmearingExample.py).
 * 
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction scale uncertainty patch file
 * @param outputname name of the output column for corrected \f$p_T\f$ values
 * @param pt name of the column containing electron \f$p_T\f$ values
 * @param eta name of the column containing electron pseudorapidities
 * @param delta_eta_sc name of the column containing the distance in
 * pseudorapidity between supercluster and electron
 * @param r9 name of the column containing the R9 value of the electron's
 * supercluster
 * @param event_seed name of the column containing the event seed for the
 * smearing
 * @param sf_file path to the correction file for the energy scale corrections
 * and variations
 * @param sf_name path to the correction and the uncertainty shifts to be
 * accessed.
 * @param variation name of the energy correction variation that should be
 * calculated (e.g., "resolutionUp", "resolutionDown", "scaleUp", "scaleDown"),
 * for "nominal" nothing is done because energy correction is already applied
 *
 * @return a dataframe containing the varied electron transverse momenta
 * 
 * @note This function is intended for analyses working with Run 3 NanoAODv12
 * or higher. In the Run 2 NanoAODv12 samples, the scale correction in data
 * is already applied in the NanoAOD files. Look at the overloaded version of
 * this function for Run 2 analyses.
 */
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &pt,
               const std::string &eta,
               const std::string &delta_eta_sc,
               const std::string &r9,
               const std::string &event_seed,
               const std::string &sf_file,
               const std::string &sf_name,
               const std::string &variation) {

    // load corrections
    auto evaluator =
        correction_manager.loadCorrection(sf_file, sf_name);

    // lambda function to apply the scale and smearing corrections
    auto correction = [evaluator, variation] (
        const ROOT::RVec<float> &pt,
        const ROOT::RVec<float> &eta,
        const ROOT::RVec<float> &delta_eta_sc,
        const ROOT::RVec<float> &r9,
        const unsigned int &event_seed
    ) {
        // initialize the random number generator with the event seed
        TRandom3 rand_gen = TRandom3(event_seed);

        // container for corrected pt values
        auto pt_corrected = ROOT::RVec<float>(pt.size());
        for (int i = 0; i < pt.size(); ++i) {
            // calculate supercluster eta
            float eta_sc = eta.at(i) + delta_eta_sc.at(i);

            // calculate the nominal corrected pt by smearing the original pt
            auto random_number = rand_gen.Gaus(0.0, 1.0);
            auto smear_nom = evaluator->evaluate({
                "smear",
                pt.at(i),
                r9.at(i),
                eta_sc
            });

            // get the scale uncertainty
            auto smear_unc = evaluator->evaluate({
                "escale",
                pt.at(i),
                r9.at(i),
                eta_sc
            });

            // get the resolution uncertainty
            auto scale_unc = evaluator->evaluate({
                "esmear",
                pt.at(i),
                r9.at(i),
                eta_sc
            });

            // set the corrected pt based on the considered variation
            float sf = 1.0;
            if (variation == "nom") {
                sf = std::max(0.0, 1.0 + smear_nom * random_number);
                pt_corrected[i] = pt.at(i) * sf;
            } else if (variation == "resolutionUp") {
                sf = std::max(0.0, 1.0 + (smear_nom + smear_unc) * random_number);
                pt_corrected[i] = pt.at(i) * sf;
            } else if (variation == "resolutionDown") {
                sf = std::max(0.0, 1.0 + (smear_nom - smear_unc) * random_number);
                pt_corrected[i] = pt.at(i) * sf;
            } else if (variation == "scaleUp") {
                sf = 1.0 + scale_unc;
                pt_corrected[i] = pt.at(i) * sf;
            } else if (variation == "scaleDown") {
                sf = 1.0 - scale_unc;
                pt_corrected[i] = pt.at(i) * sf;
            } else {
                Logger::get("physicsobject::electron::PtCorrectionMC")
                    ->debug("unknown variation {}", variation);
                throw std::runtime_error("unknown variation");
            }

            // logging output
            Logger::get("physicsobject::electron::PtCorrectionMC")
                ->debug("ele pt before {}, ele pt after {}, sf {}, variation {}",
                        pt.at(i), pt_corrected.at(i), sf, variation);

        }

        return pt_corrected;
    };

    return df.Define(
        outputname,
        correction,
        {pt, eta, delta_eta_sc, r9, event_seed}
    );
}

/**
 * @brief This function applies energy scale corrections to data. The corrections are
 * obtained from a dedicated correctionlib file.
 * 
 * For Run 3 samples, the electron scale correction is not available in the NanoAOD files
 * and the corrections have to be evaluated using a centrally provided correctionlib file.
 * This function should only be used  The documentation of the file content
 * can be found here:
 * 
 * - [2022preEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2022_Summer22_electronSS_EtDependent.html)
 * - [2022postEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2022_Summer22EE_electronSS_EtDependent.html)
 * - [2023preBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2023_Summer23_electronSS_EtDependent.html)
 * - [2023postBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2023_Summer23BPix_electronSS_EtDependent.html)
 *
 * An implementation recipe is provided here:
 * [egmScaleAndSmearingExample.py](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/egmScaleAndSmearingExample.py).
 * 
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction scale uncertainty patch file
 * @param outputname name of the output column for corrected \f$p_T\f$ values
 * @param pt name of the column containing electron \f$p_T\f$ values
 * @param eta name of the column containing electron pseudorapidities
 * @param delta_eta_sc name of the column containing the distance in
 * pseudorapidity between supercluster and electron
 * @param seed_gain name of the column containing electron gain values
 * @param r9 name of the column containing the R9 value of the electron's
 * supercluster
 * @param run name of the column containing the run number
 * @param sf_file path to the correction file for the energy scale corrections
 * @param sf_name name of the correction to be accessed.
 * calculated (e.g., "resolutionUp", "resolutionDown", "scaleUp", "scaleDown"),
 * for "nominal" nothing is done because energy correction is already applied
 *
 * @return a dataframe containing the varied electron transverse momenta
 * 
 * @note This function is intended for analyses working with Run 3 NanoAODv12
 * or higher.
 */
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &pt,
               const std::string &eta, const std::string &delta_eta_sc,
               const std::string &seed_gain,
               const std::string &r9, const std::string &run,
               const std::string &sf_file,
               const std::string &sf_name) {
    // load the correctionlib evaluator
    auto evaluator =
        correction_manager.loadCompoundCorrection(sf_file, sf_name);

    // lambda function to apply the scale correction
    auto correction = [evaluator] (
        const ROOT::RVec<float> &pt,
        const ROOT::RVec<float> &eta,
        const ROOT::RVec<float> &delta_eta_sc,
        const ROOT::RVec<UChar_t> &seed_gain,
        const ROOT::RVec<float> &r9,
        const unsigned int &run
    ) {
        // for the data correction, we just need to multiply the original pt with the scale factor
        ROOT::RVec<float> corrected_pt(pt.size());
        for (int i = 0; i < pt.size(); ++i) {
            // calculate supercluster eta
            float eta_sc = eta.at(i) + delta_eta_sc.at(i);

            // evaluate the nominal correction scale factor from correctionlib
            auto sf = evaluator->evaluate({
                "scale",
                static_cast<float>(run),
                eta_sc,
                r9.at(i),
                pt.at(i),
                static_cast<float>(seed_gain.at(i))
            });
            corrected_pt[i] = sf * pt.at(i);
            Logger::get("physicsobject::electron::PtCorrectionData")
                ->debug("ele pt before {}, ele pt after {}, sf {}",
                        pt.at(i), corrected_pt.at(i), sf);
        }
        return corrected_pt;
    };

    return df.Define(
        outputname,
        correction,
        {pt, eta, delta_eta_sc, seed_gain, r9, run}
    );
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
 * momentum (\f$p_T\f$). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation.
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
 * @param phi name of the column containing the azimuthal angle of an electron
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
 *
 * @note This function needs the dependence on phi only in case of 2023 data
 * because for whatever reason EGM POG intoduced it only in that era. 
 */
ROOT::RDF::RNode Id(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correction_manager,
                    const std::string &outputname, const std::string &pt,
                    const std::string &eta, const std::string &phi,
                    const std::string &era,
                    const std::string &wp, const std::string &sf_file,
                    const std::string &sf_name, const std::string &variation) {
    Logger::get("physicsobject::electron::scalefactor::Id")
        ->debug("Setting up functions for electron id sf with correctionlib");
    Logger::get("physicsobject::electron::scalefactor::Id")
        ->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, era, sf_name, wp, variation](const float &pt,
                                                 const float &eta,
                                                 const float &phi) {
            Logger::get("physicsobject::electron::scalefactor::Id")
                ->debug("Era {}, Variation {}, WP {}", era, variation, wp);
            Logger::get("physicsobject::electron::scalefactor::Id")
                ->debug("ID - pt {}, eta {}, phi {}", pt, eta, phi);
            double sf = 1.;
            if (pt >= 0.0) {
                if (era.find("2023") != std::string::npos) {
                    // for 2023, phi is needed as input
                    sf = evaluator->evaluate({era, variation, wp, eta, pt, phi});
                } 
                else sf = evaluator->evaluate({era, variation, wp, eta, pt});
            }
            Logger::get("physicsobject::electron::scalefactor::Id")
                ->debug("Scale Factor {}", sf);
            return sf;
        },
        {pt, eta, phi});
    return df1;
}

/**
 * @brief This function calculates single electron trigger scale factors (SFs)
 * for a single electron dependening on its pseudorapidity (\f$\eta\f$), its transverse
 * momentum (\f$p_T\f$), and the electron identification working point. The scale factors
 * are loaded from a correctionlib file using a specified scale factor name and variation.
 * This function only uses the scale factor from the correctionlib evaluation if the
 * corresponding trigger flag is set to `true`. Otherwise, it returns a scale factor of
 * 1.0.
 *
 * Recommendations by EgammaPOG:
 * - [Run3 triggers](https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary)
 * - [Run3 scale factors](https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3)
 * 
 * The documentation of the corresponding jsonPOG files can be found here:
 * - [2022preEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2022_Summer22_electronHlt.html)
 * - [2022postEE](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2022_Summer22EE_electronHlt.html)
 * - [2023preBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2023_Summer23_electronHlt.html)
 * - [2023postBPix](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2023_Summer23BPix_electronHlt.html)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * electron scale factor file
 * @param outputname name of the output column containing the ID scale factor
 * @param pt name of the column containing the transverse momentum of an
 * electron
 * @param eta name of the column containing the pseudorapidity of an electron
 * @param trigger_flag name of the column containing the trigger flag
 * @param era string with the era name of a data taking period, e.g.
 * "2016preVFP"
 * @param path_id_name string that serves as an identifier for the used
 * combination of trigger path and electron ID, e.g. "HLT_SF_Ele30_TightID",
 * "HLT_SF_Ele30_MVAiso90ID", ...
 * @param sf_file path to the file with the electron scale factors
 * @param sf_name name of the electron scale factor for the ID correction,
 * e.g. "UL-Electron-ID-SF"
 * @param variation name the scale factor variation, "sf" for the nominal
 * scale factor and "sfup"/"sfdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Trigger(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correction_manager,
                    const std::string &outputname, const std::string &pt,
                    const std::string &eta, const std::string &trigger_flag,
                    const std::string &era, const std::string &path_id_name,
                    const std::string &sf_file, const std::string &sf_name,
                    const std::string &variation) {
    Logger::get("physicsobject::electron::scalefactor::Trigger")
        ->debug("Setting up functions for electron trigger sf with correctionlib");
    Logger::get("physicsobject::electron::scalefactor::Trigger")
        ->debug("Scale factor - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, era, sf_name, path_id_name, variation](
            const float &pt, const float &eta, const bool &trigger_flag) {
            Logger::get("physicsobject::electron::scalefactor::Trigger")
                ->debug("Era {}, Variation {}, Trigger at electron ID {}", era, variation, path_id_name);
            Logger::get("physicsobject::electron::scalefactor::Trigger")
                ->debug("ID - pt {}, eta {}, trigger flag {}", pt, eta, trigger_flag);
            double sf = 1.;
            // check to prevent electrons with default values due to tau energy
            // correction shifts below good tau pt selection
            try {
                if (pt >= 0.0 && std::abs(eta) >= 0.0 && trigger_flag) {
                    sf = evaluator->evaluate({era, variation, path_id_name, eta, pt});
                    Logger::get("physicsobject::electron::scalefactor::Trigger")
                        ->debug("Scale Factor {}", sf);
                }
            } catch (const std::runtime_error &e) {
                // this error can occur because the pt range starts at different
                // values for different triggers
                Logger::get("physicsobject::electron::scalefactor::Trigger")
                    ->debug("SF evaluation for {} failed for pt {}", sf_name,
                            pt);
            }
            return sf;
        },
        {pt, eta, trigger_flag});
    return df1;
}

} // end namespace scalefactor
} // end namespace electron
} // end namespace physicsobject
#endif /* GUARD_ELECTRONS_H */
