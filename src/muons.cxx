#ifndef GUARD_MUONS_H
#define GUARD_MUONS_H

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/RoccoR.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "correction.h"

namespace physicsobject {
namespace muon {

/**
 * @brief This function defines a new column in the dataframe that contains the
 * corrected transverse momentum (\f$p_T\f$) of a muon, using the Muon-POG
 * `muonscarekit` scale-and-resolution ("ScaRe") corrections for Run3.
 *
 * The correction is taken from
 * [MuonScaRe.py](https://gitlab.cern.ch/cms-muonPOG/muonscarekit/-/blob/master/scripts/MuonScaRe.py)
 * and its application is taken from
 * [apply_corrections_rdf.py](https://gitlab.cern.ch/cms-muonPOG/muonscarekit/-/blob/master/scripts/apply_corrections_rdf.py).
 *
 * A charge- and \f$(\eta,\phi)\f$-dependent scale correction (`pt_scale`) is
 * applied to both data and MC. For MC, a resolution smearing (`pt_resol`) is
 * additionally applied on top of the scale-corrected \f$p_T\f$, using a random
 * number drawn from an asymmetric Crystal Ball distribution (`CrystalBall`,
 * see `RoccoR.hxx`) and a residual data/MC resolution difference `k`. If
 * `shift` requests a systematic variation ("ScaleUp"/"ScaleDown" for the
 * scale uncertainty, "ResoUp"/"ResoDown" for the resolution uncertainty), the
 * corresponding uncertainty (`pt_scale_var`/`pt_resol_var`) is propagated
 * instead of the nominal correction. Muons corrected outside of
 * `[min_muon_pt, 200]` GeV, or for which a correction evaluates to NaN or an
 * unreasonable ratio to the input \f$p_T\f$, are reset to their pre-step
 * value.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * muon scale-and-resolution correction file
 * @param outputname name of the new column containing the corrected \f$p_T\f$
 * values
 * @param pt name of the column containing the muon transverse momenta
 * @param eta name of the column containing the muon eta values
 * @param phi name of the column containing the muon phi values
 * @param charge name of the column containing the muon charges
 * @param n_tracker_layers name of the column containing the number of tracker
 * layers
 * @param lumi name of the column containing the luminosity block number
 * @param event name of the column containing the event number
 * @param min_muon_pt lower \f$p_T\f$ boundary below which a corrected value is
 * considered invalid and reset to the input value
 * @param scale_file path to the muonscarekit scale-and-resolution correction
 * file
 * @param shift name of the requested variation, "nom" for the nominal
 * correction, "ScaleUp"/"ScaleDown" for the scale uncertainty, or
 * "ResoUp"/"ResoDown" for the resolution uncertainty
 * @param is_data flag that is `true` if the input is data and `false` if it
 * is MC; resolution smearing is only applied to MC
 *
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df,
                correctionManager::CorrectionManager &correction_manager,
                const std::string &outputname,
                const std::string &pt, const std::string &eta,
                const std::string &phi, const std::string &charge,
                const std::string &n_tracker_layers,
                const std::string &lumi,
                const std::string &event,
                const float &min_muon_pt,
                const std::string &scale_file,
                const std::string &shift, bool is_data) {

    std::string name = is_data ? "data" : "mc";

    auto scale_evaluator_a = correction_manager.loadCorrection(
        scale_file, "a_" + name);

    auto scale_evaluator_m = correction_manager.loadCorrection(
        scale_file, "m_" + name);

    auto reso_evaluator_kmc = correction_manager.loadCorrection(
        scale_file, "k_mc");

    auto reso_evaluator_kdata = correction_manager.loadCorrection(
        scale_file, "k_data");

    auto reso_evaluator_cb = correction_manager.loadCorrection(
        scale_file, "cb_params");

    auto reso_evaluator_poly = correction_manager.loadCorrection(
        scale_file, "poly_params");

    auto smear_evaluator = correction_manager.loadCorrection(
        scale_file, "RandomSmearing");

    auto lambda = [scale_evaluator_a, scale_evaluator_m, reso_evaluator_kmc,
                   reso_evaluator_kdata, reso_evaluator_cb,
                   reso_evaluator_poly, smear_evaluator, min_muon_pt, shift,
                   is_data](
                      const ROOT::RVec<float> &pt,
                      const ROOT::RVec<float> &eta,
                      const ROOT::RVec<float> &phi,
                      const ROOT::RVec<int> &charge,
                      const ROOT::RVec<UChar_t> &n_tracker_layers,
                      const unsigned int &lumi,
                      const unsigned long long &event) {

        ROOT::RVec<float> corrected_pts(pt.size());

        for (std::size_t i = 0; i < pt.size(); ++i) {

            float corr_pt = pt.at(i);
            float scale_pt = pt.at(i);

            float a_f = scale_evaluator_a->evaluate({eta.at(i), phi.at(i), "nom"});
            float m_f = scale_evaluator_m->evaluate({eta.at(i), phi.at(i), "nom"});

            // apply scale correction to both data and mc
            scale_pt = 1. / (m_f / pt.at(i) + charge.at(i) * a_f);

            if (scale_pt < min_muon_pt || scale_pt > 200.0 || std::isnan(scale_pt)) {
                // set pt to initial value if correction is outside of valid boundaries
                scale_pt = pt.at(i);
            }

            corr_pt = scale_pt;

            if (is_data == false) {
                // apply resolution smearing
                float mean_f = reso_evaluator_cb->evaluate({std::abs(eta.at(i)), float(n_tracker_layers.at(i)), 0});
                float sigma_f = reso_evaluator_cb->evaluate({std::abs(eta.at(i)), float(n_tracker_layers.at(i)), 1});
                float n_f = reso_evaluator_cb->evaluate({std::abs(eta.at(i)), float(n_tracker_layers.at(i)), 2});
                float alpha_f = reso_evaluator_cb->evaluate({std::abs(eta.at(i)), float(n_tracker_layers.at(i)), 3});

                float rndm_f = smear_evaluator->evaluate(
                    {static_cast<int>(event), static_cast<int>(lumi), phi.at(i)});

                CrystalBall cb;
                cb.m = mean_f;
                cb.s = sigma_f;
                cb.a = alpha_f;
                cb.n = n_f;
                cb.init();
                float rndm = cb.invcdf(rndm_f);

                float param0_f = reso_evaluator_poly->evaluate({std::abs(eta.at(i)), float(n_tracker_layers.at(i)), 0});
                float param1_f = reso_evaluator_poly->evaluate({std::abs(eta.at(i)), float(n_tracker_layers.at(i)), 1});
                float param2_f = reso_evaluator_poly->evaluate({std::abs(eta.at(i)), float(n_tracker_layers.at(i)), 2});

                float sigma_std = param0_f + param1_f * scale_pt + param2_f * scale_pt * scale_pt;
                float std_dev = std::max(0.f, sigma_std);

                float k_data_f = reso_evaluator_kdata->evaluate({std::abs(eta.at(i)), "nom"});
                float k_mc_f = reso_evaluator_kmc->evaluate({std::abs(eta.at(i)), "nom"});
                float k = 0.0;
                if (k_mc_f < k_data_f) {
                    k = std::sqrt(k_data_f * k_data_f - k_mc_f * k_mc_f);
                }

                float reso_pt = scale_pt * (1 + k * std_dev * rndm);

                if (reso_pt < min_muon_pt || reso_pt > 200.0 || std::isnan(reso_pt)) {
                    // set pt to initial value if correction is outside of valid boundaries
                    reso_pt = scale_pt;
                }

                if (reso_pt / scale_pt > 2 || reso_pt / scale_pt < 0.1 || reso_pt < 0) {
                    // set pt to initial value if correction is outside of valid boundaries,
                    // not to be merged with the step before
                    reso_pt = scale_pt;
                }

                corr_pt = reso_pt;

                if (shift == "ScaleStatUp" || shift == "ScaleStatDown") {
                    // apply scale uncertainty on top of the scale+reso corrected pt
                    float stat_a_f = scale_evaluator_a->evaluate({eta.at(i), phi.at(i), "stat"});
                    float stat_m_f = scale_evaluator_m->evaluate({eta.at(i), phi.at(i), "stat"});
                    float stat_rho_f = scale_evaluator_m->evaluate({eta.at(i), phi.at(i), "rho_stat"});

                    float unc = corr_pt * corr_pt * std::sqrt(
                        stat_m_f * stat_m_f / (corr_pt * corr_pt) +
                        stat_a_f * stat_a_f +
                        2 * charge.at(i) * stat_rho_f * stat_m_f * stat_a_f / corr_pt);

                    if (shift == "ScaleStatUp") corr_pt = corr_pt + unc;
                    else corr_pt = corr_pt - unc;
                } else if (shift == "ScaleSystUp" || shift == "ScaleSystDown") {
                    // apply scale uncertainty on top of the scale+reso corrected pt
                    float syst_a_f = scale_evaluator_a->evaluate({eta.at(i), phi.at(i), "syst"});
                    float syst_m_f = scale_evaluator_m->evaluate({eta.at(i), phi.at(i), "syst"});
                    float syst_rho_f = scale_evaluator_m->evaluate({eta.at(i), phi.at(i), "rho_syst"});

                    float unc = corr_pt * corr_pt * std::sqrt(
                        syst_m_f * syst_m_f / (corr_pt * corr_pt) +
                        syst_a_f * syst_a_f +
                        2 * charge.at(i) * syst_rho_f * syst_m_f * syst_a_f / corr_pt);

                    if (shift == "ScaleSystUp") corr_pt = corr_pt + unc;
                    else corr_pt = corr_pt - unc;
                } else if (shift == "ResoStatUp" || shift == "ResoStatDown") {
                    // apply resolution uncertainty on top of the scale corrected pt
                    float k_unc_f = reso_evaluator_kmc->evaluate({std::abs(eta.at(i)), "stat"});

                    if (k_mc_f > 0) {
                        float std_x_cb = (reso_pt / scale_pt - 1) / k_mc_f;
                        if (shift == "ResoStatUp") corr_pt = scale_pt * (1 + (k_mc_f + k_unc_f) * std_x_cb);
                        else corr_pt = scale_pt * (1 + (k_mc_f - k_unc_f) * std_x_cb);
                    }

                    if (corr_pt / scale_pt > 2 || corr_pt / scale_pt < 0.1 || corr_pt < 0) corr_pt = scale_pt;
                } else if (shift == "ResoSystUp" || shift == "ResoSystDown") {
                    // apply resolution uncertainty on top of the scale corrected pt
                    float k_unc_f = reso_evaluator_kmc->evaluate({std::abs(eta.at(i)), "syst"});

                    if (k_mc_f > 0) {
                        float std_x_cb = (reso_pt / scale_pt - 1) / k_mc_f;
                        if (shift == "ResoSystUp") corr_pt = scale_pt * (1 + (k_mc_f + k_unc_f) * std_x_cb);
                        else corr_pt = scale_pt * (1 + (k_mc_f - k_unc_f) * std_x_cb);
                    }

                    if (corr_pt / scale_pt > 2 || corr_pt / scale_pt < 0.1 || corr_pt < 0) corr_pt = scale_pt;
                }

            }

            Logger::get("physicsobject::muon::PtCorrectionMC")
                ->debug("muon pt before {}, muon pt after {} shift {}", pt.at(i),
                        corr_pt, shift);

            corrected_pts.at(i) = corr_pt;
        }

        return corrected_pts;
    };

    return df.Define(outputname, lambda,
                     {pt, eta, phi, charge, n_tracker_layers,
                      lumi, event});
}

/**
 * @brief This function defines a new column in the dataframe that contains the
 * corrected transverse momentum (\f$p_T\f$) of a muon, calculated using the
 * Rochester correction for Monte Carlo (MC) simulation.
 *
 * The correction is taken from https://gitlab.cern.ch/akhukhun/roccor (Run2)
 * and is evaluated using the `RoccoR` correction function which is also taken
 * from this repository.
 *
 * The correction is applied depending on whether the muon could be matched to
 * a generator level muon or not. If the muon is matched, the correction is done
 * using `kSpreadMC`, otherwise `kSmearMC` is used.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected \f$p_T\f$
 * values
 * @param charge name of the column containing the muon charges
 * @param pt name of the column containing the muon transverse momenta
 * @param eta name of the column containing the muon eta values
 * @param phi name of the column containing the muon phi values
 * @param gen_pt name of the column containing the generator level transverse
 * momentum of the matched muon, if no muon can be match a negative default
 * value has to be provide
 * @param n_tracker_layers name of the column containing the number of tracker
 * layers
 * @param rndm_vector name of the column containing random values for the
 * smearing
 * @param index_vector name of the column containing index values
 * @param position position within the index vector used to retrieve the index
 * of the wanted muon
 * @param filename file path to the Rochester correction
 * @param error_set error set number that should be used
 * @param error_member error member number that should be used
 *
 * @return a dataframe with the new column
 *
 * @note TODO: Corrections for Run3 are not yet implemented
 */
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df, const std::string &outputname,
               const std::string &charge, const std::string &pt,
               const std::string &eta, const std::string &phi,
               const std::string &gen_pt, const std::string &n_tracker_layers,
               const std::string &rndm_vector, const std::string &index_vector,
               const int &position, const std::string &filename,
               const int error_set, const int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set, error_member](
                      const ROOT::RVec<int> &charges,
                      const ROOT::RVec<float> &pts,
                      const ROOT::RVec<float> &etas,
                      const ROOT::RVec<float> &phis, const float &gen_pt,
                      const ROOT::RVec<int> &n_tracker_layers,
                      const ROOT::RVec<float> &rndm_values,
                      const ROOT::RVec<int> &index_vector) {
        double corrected_pt = default_float;
        const int index = index_vector.at(position);
        if (gen_pt > 0.) {
            corrected_pt =
                pts.at(index) * rc.kSpreadMC(charges.at(index), pts.at(index),
                                             etas.at(index), phis.at(index),
                                             gen_pt, error_set, error_member);
        } else {
            corrected_pt =
                pts.at(index) *
                rc.kSmearMC(charges.at(index), pts.at(index), etas.at(index),
                            phis.at(index), n_tracker_layers.at(index),
                            rndm_values.at(position), error_set, error_member);
        }
        Logger::get("physicsobject::muon::PtCorrectionMC")
            ->debug("muon pt before {}, muon pt after {}", pts.at(index),
                    corrected_pt);
        return corrected_pt;
    };

    return df.Define(outputname, lambda,
                     {charge, pt, eta, phi, gen_pt, n_tracker_layers,
                      rndm_vector, index_vector});
}

/**
 * @brief This function defines a new column in the dataframe that contains the
 * corrected transverse momentum (\f$p_T\f$) of a muon, calculated using the
 * Rochester correction for data.
 *
 * The correction is taken from https://gitlab.cern.ch/akhukhun/roccor (Run2)
 * and is evaluated using the `RoccoR` correction function which is also taken
 * from this repository.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected \f$p_T\f$
 * values
 * @param charge name of the column containing the muon charges
 * @param pt name of the column containing the muon transverse momenta
 * @param eta name of the column containing the muon eta values
 * @param phi name of the column containing the muon phi values
 * @param index_vector name of the column containing index values
 * @param position position within the index vector used to retrieve the index
 * of the wanted muon
 * @param filename file path to the Rochester correction
 * @param error_set error set number that should be used
 * @param error_member error member number that should be used
 *
 * @return a dataframe with the new column
 *
 * @note TODO: Corrections for Run3 are not yet implemented
 */
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df, const std::string &outputname,
                 const std::string &charge, const std::string &pt,
                 const std::string &eta, const std::string &phi,
                 const std::string &index_vector, const int &position,
                 const std::string &filename, int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set,
                   error_member](const ROOT::RVec<int> &charges,
                                 const ROOT::RVec<float> &pts,
                                 const ROOT::RVec<float> &etas,
                                 const ROOT::RVec<float> &phis,
                                 const ROOT::RVec<int> &index_vector) {
        const int index = index_vector.at(position);
        double corrected_pt =
            pts.at(index) * rc.kScaleDT(charges.at(index), pts.at(index),
                                        etas.at(index), phis.at(index),
                                        error_set, error_member);
        Logger::get("physicsobject::muon::PtCorrectionData")
            ->debug("muon pt before {}, muon pt after {}", pts.at(index),
                    corrected_pt);
        return corrected_pt;
    };

    return df.Define(outputname, lambda, {charge, pt, eta, phi, index_vector});
}

namespace scalefactor {

/**
 * @brief This function calculates muon reco scale factors (SFs) for a single
 * muon dependening on its pseudorapidity (\f$\eta\f$) and transverse momentum
 * (\f$p_T\f$). The scale factors are loaded from a correctionlib file using a
 * specified scale factor name and variation.
 *
 * Recommendations by MuonPOG:
 * - [Run2](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_5_1)
 * - [Run3](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_5_2)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * muon scale factor file
 * @param outputname name of the output column containing the reco scale factor
 * @param pt name of the column containing the transverse momentum of a muon
 * @param eta name of the column containing the pseudorapidity of a muon
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon scale factor for the reco correction,
 * e.g. "NUM_TrackerMuons_DEN_genTracks"
 * @param variation name the scale factor variation, "nominal" for the nominal
 * scale factor and "systup"/"systdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 *
 * @note The \f$p_T\f$ dependence of the reco scale factor is only for
 * consistency with the other scale factors. It was only derived in one
 * \f$p_T\f$ bin and should be applied for all muons \f$p_T\f$'s
 * (recommendation: [10-200] GeV).
 */
ROOT::RDF::RNode Reco(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correction_manager,
                      const std::string &outputname, const std::string &pt,
                      const std::string &eta, const std::string &sf_file,
                      const std::string &sf_name,
                      const std::string &variation) {
    Logger::get("physicsobject::muon::scalefactor::Reco")
        ->debug("Setting up functions for muon reco sf");
    Logger::get("physicsobject::muon::scalefactor::Reco")
        ->debug("Reco - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, variation](const float &pt, const float &eta) {
            Logger::get("physicsobject::muon::scalefactor::Reco")
                ->debug("Reco - pt {}, eta {}", pt, eta);
            double sf = 1.;
            // check to prevent muons with default values due to tau energy
            // correction shifts below good tau pt selection
            if (pt >= 40.0 && std::abs(eta) >= 0.0) {
                sf = evaluator->evaluate({std::abs(eta), pt, variation});
            }
            // the reco scale factor in the json file is only defined above 40
            // GeV
            else if (pt >= 0.0 && pt < 40.0 && std::abs(eta) >= 0.0) {
                sf = evaluator->evaluate({std::abs(eta), 40.0, variation});
            }
            return sf;
        },
        {pt, eta});
    return df1;
}

/**
 * @brief This function calculates muon iso and ID scale factors (SFs) for a single
 * muon dependening on its pseudorapidity (\f$\eta\f$) and transverse momentum
 * (\f$p_T\f$). The scale factors are loaded from a correctionlib file using a
 * specified scale factor name and variation.
 *
 * Recommendations by MuonPOG:
 * - [Run2](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_7_1)
 * - [Run3](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_7_2)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * muon scale factor file
 * @param outputname name of the output column containing the scale factor
 * @param pt name of the column containing the transverse momentum of a muon
 * @param eta name of the column containing the pseudorapidity of a muon
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon scale factor for the correction,
 * e.g. "NUM_TightRelIso_DEN_MediumID"
 * @param variation name the scale factor variation, "nominal" for the nominal
 * scale factor and "systup"/"systdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode IsoAndID(ROOT::RDF::RNode df,
                     correctionManager::CorrectionManager &correction_manager,
                     const std::string &outputname, const std::string &pt,
                     const std::string &eta, const std::string &sf_file,
                     const std::string &sf_name, const std::string &variation) {
    Logger::get("physicsobject::muon::scalefactor::IsoAndID")
        ->debug("Setting up functions for muon iso sf");
    Logger::get("physicsobject::muon::scalefactor::IsoAndID")
        ->debug("SF - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, variation](const float &pt, const float &eta) {
            Logger::get("physicsobject::muon::scalefactor::IsoAndID")
                ->debug("SF - pt {}, eta {}", pt, eta);
            double sf = 1.;
            // check to prevent muons with default values due to tau energy
            // correction shifts below good tau pt selection
            if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                if (variation=="nominal") {
                    sf = evaluator->evaluate({std::abs(eta), pt, "nominal"});
                }
                else if (variation=="systup") {
                    sf = sf + evaluator->evaluate({std::abs(eta), pt, "syst"});
                }
                else if (variation=="systdown") {
                    sf = sf - evaluator->evaluate({std::abs(eta), pt, "syst"});
                }
                else if (variation=="statup") {
                    sf = sf + evaluator->evaluate({std::abs(eta), pt, "stat"});
                }
                else if (variation=="statdown") {
                    sf = sf - evaluator->evaluate({std::abs(eta), pt, "stat"});
                }
                else { 
                    Logger::get("physicsobject::muon::scalefactor::IsoAndID")
                         ->debug("variation {} not implemented, check your code");
                }
            }
            return sf;
        },
        {pt, eta});
    return df1;
}

/**
 * @brief This function calculates muon trigger scale factors (SFs) for a single
 * muon dependening on its pseudorapidity (\f$\eta\f$) and transverse momentum
 * (\f$p_T\f$). The scale factors are loaded from a correctionlib file using a
 * specified scale factor name and variation.
 *
 * Recommendations by MuonPOG:
 * - [Run2](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_8_1)
 * - [Run3](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_8_2)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * muon scale factor file
 * @param outputname name of the output column containing the trigger scale
 * factor
 * @param pt name of the column containing the transverse momentum of a muon
 * @param eta name of the column containing the pseudorapidity of a muon
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon scale factor for the trigger correction,
 * e.g. "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"
 * @param variation name the scale factor variation, "nominal" for the nominal
 * scale factor and "systup"/"systdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 *
 * @warning This function is deprecated. Use the overloaded function with the
 * additional parameter `trigger_flag` instead.
 */
ROOT::RDF::RNode
Trigger(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname, const std::string &pt,
        const std::string &eta, const std::string &sf_file,
        const std::string &sf_name, const std::string &variation) {
    Logger::get("physicsobject::muon::scalefactor::Trigger")
        ->warn("Function is deprecated, use the overloaded version instead");
    Logger::get("physicsobject::muon::scalefactor::Trigger")
        ->debug("Setting up functions for muon trigger sf");
    Logger::get("physicsobject::muon::scalefactor::Trigger")
        ->debug("Trigger - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, variation, sf_name](const float &pt, const float &eta) {
            Logger::get("physicsobject::muon::scalefactor::Trigger")
                ->debug("Trigger - pt {}, eta {}", pt, eta);
            double sf = 1.;
            // check to prevent muons with default values due to tau energy
            // correction shifts below good tau pt selection
            try {
                if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                    sf = evaluator->evaluate({std::abs(eta), pt, variation});
                }
            } catch (const std::runtime_error &e) {
                // this error can occur because the pt range starts at different
                // values for different triggers
                Logger::get("physicsobject::muon::scalefactor::Trigger")
                    ->debug("SF evaluation for {} failed for pt {}", sf_name,
                            pt);
            }
            return sf;
        },
        {pt, eta});
    return df1;
}

/**
 * @brief This function calculates muon trigger scale factors (SFs) for a single
 * muon depending on its pseudorapidity (\f$\eta\f$) and transverse momentum
 * (\f$p_T\f$). The scale factors are loaded from a correctionlib file using a
 * specified scale factor name and variation. This function only uses the scale
 * factor from the correctionlib evaluation if the corresponding trigger flag
 * is set to `true`. Otherwise, it returns a scale factor of 1.0.
 *
 * Recommendations by MuonPOG:
 * - [Run2](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_8_1)
 * - [Run3](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_8_2)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * muon scale factor file
 * @param outputname name of the output column containing the trigger scale
 * factor
 * @param pt name of the column containing the transverse momentum of a muon
 * @param eta name of the column containing the pseudorapidity of a muon
 * @param trigger_flag name of the column containing the trigger flag
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon scale factor for the trigger correction,
 * e.g. "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"
 * @param variation name the scale factor variation, "nominal" for the nominal
 * scale factor and "systup"/"systdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Trigger(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname, const std::string &pt,
        const std::string &eta, const std::string &trigger_flag,
        const std::string &sf_file, const std::string &sf_name,
        const std::string &variation) {
    Logger::get("physicsobject::muon::scalefactor::Trigger")
        ->debug("Setting up functions for muon trigger sf");
    Logger::get("physicsobject::muon::scalefactor::Trigger")
        ->debug("Trigger - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, variation, sf_name](const float &pt, const float &eta,
                                        const bool &trigger_flag) {
            Logger::get("physicsobject::muon::scalefactor::Trigger")
                ->debug("Trigger - pt {}, eta {}, trigger flag {}", pt, eta,
                        trigger_flag);
            double sf = 1.;
            // check to prevent muons with default values due to tau energy
            // correction shifts below good tau pt selection
            try {
                if (trigger_flag) {
                    sf = evaluator->evaluate({std::abs(eta), pt, variation});
                }
            } catch (const std::runtime_error &e) {
                // this error can occur because the pt range starts at different
                // values for different triggers
                Logger::get("physicsobject::muon::scalefactor::Trigger")
                    ->debug("SF evaluation for {} failed for pt {}", sf_name,
                            pt);
            }
            return sf;
        },
        {pt, eta, trigger_flag});
    return df1;
}
} // end namespace scalefactor
} // end namespace muon
} // end namespace physicsobject
#endif /* GUARD_MUONS_H */
