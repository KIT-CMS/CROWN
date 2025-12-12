#ifndef GUARD_REWEIGHTING_H
#define GUARD_REWEIGHTING_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "correction.h"
#include <Math/Vector4D.h>


namespace event {
namespace reweighting {

/**
 * @brief This function is used to correct Monte Carlo (MC) simulations for 
 * differences in the pileup distribution compared to the one measured in data. 
 * It retrieves a per-event weight from a correction file based on the true number 
 * of pileup interactions in an event.
 *
 * The correction files are provided by the Luminosity POG and more information 
 * about the pileup reweighting can be found here: 
 * https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * pileup weights file
 * @param outputname name of the output column containing the pileup event weight
 * @param true_pileup_number name of the column containing the true mean number 
 * of the poisson distribution for an event from which the number of interactions 
 * each bunch crossing has been sampled
 * @param corr_file path to the file with the pileup weights
 * @param corr_name name of the pileup correction in the file, 
 * e.g. "Collisions18_UltraLegacy_goldenJSON"
 * @param variation name of the pileup weight variation, options are "nominal", 
 * "up" and "down"
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Pileup(ROOT::RDF::RNode df,
          correctionManager::CorrectionManager &correction_manager,
          const std::string &outputname, const std::string &true_pileup_number,
          const std::string &corr_file, const std::string &corr_name,
          const std::string &variation) {
    auto evaluator = correction_manager.loadCorrection(corr_file, corr_name);
    auto df1 = df.Define(outputname,
                         [evaluator, variation](const float &pu) {
                             double weight =
                                 evaluator->evaluate({pu, variation});
                             return weight;
                         },
                         {true_pileup_number});
    return df1;
}

/**
 * @brief This function is used to evaluate the parton shower (PS) weight of an event. 
 * The weights are stored in the nanoAOD files and defined as 
 * \f$w_{variation}\f$ / \f$w_{nominal}\f$. The nominal weight is already applied, 
 * therefore, the main use of this function is to get the initial state radiation (ISR) 
 * and final state radiation (FSR) variations to the nominal PS weight.
 *
 * Depending on the selected ISR and FSR value, a specific index has to be identified. 
 * The mapping between the index and the ISR and FSR values is:
 *  ISR       | FSR        | index
 * -----------|------------|---------
 *  2.0       | 1.0        | 0
 *  1.0       | 2.0        | 1
 *  0.5       | 1.0        | 2
 *  1.0       | 0.5        | 3
 *
 * @note For some simulated samples this mapping might be defined differently, 
 * therefore, it is advisable to check the documentation of the `PSWeight` 
 * branch in the nanoAOD files of the samples if issues occur.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the ISR/FSR event weight
 * @param ps_weights name of the column containing the parton shower (ISR/FSR) weights
 * @param isr value of the ISR variation, possible values are 0.5, 1.0, 2.0
 * @param fsr value of the FSR variation, possible values are 0.5, 1.0, 2.0
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode PartonShower(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &ps_weights,
                            const float isr, const float fsr) {
    // find the index we have to use, first check if the isr and fsr values are
    // valid, only 0.5, 1.0, 2.0 are allowed
    std::vector<float> allowed_values = {0.5, 1.0, 2.0};
    if (std::find(allowed_values.begin(), allowed_values.end(), isr) ==
        allowed_values.end()) {
        Logger::get("event::reweighting::PartonShower")
            ->error("Invalid value for isr: {}", isr);
        throw std::runtime_error("Invalid value for isr");
    }
    if (std::find(allowed_values.begin(), allowed_values.end(), fsr) ==
        allowed_values.end()) {
        Logger::get("event::reweighting::PartonShower")
            ->error("Invalid value for fsr: {}", fsr);
        throw std::runtime_error("Invalid value for fsr");
    }
    
    auto ps_weights_lambda =
        [isr, fsr](const ROOT::RVec<float> ps_weights) {
            if (isr == 1.0 && fsr == 1.0) {
                // if the ISR and FSR are both 1.0, we return the nominal weight
                return (float)1.0;
            }
            // now find the index
            std::map<std::pair<const float, const float>, int> index_map;
            if (ps_weights.size() == 4) {
                index_map = {
                    {{2.0, 1.0}, 0}, {{1.0, 2.0}, 1}, 
                    {{0.5, 1.0}, 2}, {{1.0, 0.5}, 3}
                };
            } else {
                Logger::get("event::reweighting::PartonShower")
                    ->error("Invalid number of PS weights: {}",
                            ps_weights.size());
                throw std::runtime_error("Invalid number of PS weights");
            }
            std::pair<const float, const float> variations = {isr, fsr};
            int index = index_map[variations];
            return ps_weights.at(index);
        };
    auto df1 =
        df.Define(outputname, ps_weights_lambda, {ps_weights});
    return df1;
}

/**
 * @brief This function is used to evaluate the LHE scale weight of an event. The weights
 * are stored in the nanoAOD files and defined as \f$w_{variation}\f$ / \f$w_{nominal}\f$. 
 * The nominal weight is already applied, therefore, the main use of this function is 
 * to get the factorization and renormalization scale variations to the nominal scale 
 * weight.
 *
 * Depending on the selected \f$\mu_R\f$ and \f$\mu_F\f$ value, a specific index has 
 * to be identified. The mapping between the index and the \f$\mu_R\f$ and \f$\mu_F\f$ 
 * values is:
 *  mu_f       | mu_r        | index
 * ------------|-------------|---------
 *  0.5        | 0.5         | 0
 *  1.0        | 0.5         | 1
 *  2.0        | 0.5         | 2
 *  0.5        | 1.0         | 3
 *  1.0        | 1.0         | 4 (not always included)
 *  2.0        | 1.0         | 5 (4)
 *  0.5        | 2.0         | 6 (5)
 *  1.0        | 2.0         | 7 (6)
 *  2.0        | 2.0         | 8 (7)
 *
 * @note For some simulated samples this mapping might be defined differently, 
 * therefore, it is advisable to check the documentation of the `LHEScaleWeight` 
 * branch in the nanoAOD files of the samples if issues occur.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the LHE scale event weight
 * @param lhe_scale_weights name of the column containing the LHE scale weights
 * @param mu_r value of \f$\mu_R\f$ variation, possible values are 0.5, 1.0, 2.0
 * @param mu_f value of \f$\mu_F\f$ variation, possible values are 0.5, 1.0, 2.0
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode LHEscale(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &lhe_scale_weights,
                            const float mu_r, const float mu_f) {
    // find the index we have to use, first check if the mu_r and mu_f values are
    // valid, only 0.5, 1.0, 2.0 are allowed
    std::vector<float> allowed_values = {0.5, 1.0, 2.0};
    if (std::find(allowed_values.begin(), allowed_values.end(), mu_r) ==
        allowed_values.end()) {
        Logger::get("event::reweighting::LHEscale")
            ->error("Invalid value for mu_r: {}", mu_r);
        throw std::runtime_error("Invalid value for mu_r");
    }
    if (std::find(allowed_values.begin(), allowed_values.end(), mu_f) ==
        allowed_values.end()) {
        Logger::get("event::reweighting::LHEscale")
            ->error("Invalid value for mu_f: {}", mu_f);
        throw std::runtime_error("Invalid value for mu_f");
    }
    
    auto lhe_scale_weights_lambda =
        [mu_r, mu_f](const ROOT::RVec<float> scale_weights) {
            // now find the index
            std::map<std::pair<const float, const float>, int> index_map;
            if (scale_weights.size() == 9) {
                index_map = {
                    {{0.5, 0.5}, 0}, {{1.0, 0.5}, 1}, {{2.0, 0.5}, 2},
                    {{0.5, 1.0}, 3}, {{1.0, 1.0}, 4}, {{2.0, 1.0}, 5},
                    {{0.5, 2.0}, 6}, {{1.0, 2.0}, 7}, {{2.0, 2.0}, 8}
                };
            }
            else if (scale_weights.size() == 8) {
                index_map = {
                    {{0.5, 0.5}, 0}, {{1.0, 0.5}, 1}, {{2.0, 0.5}, 2},
                    {{0.5, 1.0}, 3}, {{2.0, 1.0}, 4}, {{0.5, 2.0}, 5}, 
                    {{1.0, 2.0}, 6}, {{2.0, 2.0}, 7}
                };
            }
            else {
                Logger::get("event::reweighting::LHEscale")
                    ->error("Invalid number of LHE scale weights: {}",
                            scale_weights.size());
                throw std::runtime_error("Invalid number of LHE scale weights");
            }
            std::pair<const float, const float> variations = {mu_f, mu_r};
            int index = index_map[variations];
            return scale_weights.at(index);
        };
    auto df1 =
        df.Define(outputname, lhe_scale_weights_lambda, {lhe_scale_weights});
    return df1;
}

/**
 * @brief This function is used to evaluate the LHE PDF weight of an event. The weights
 * are stored in the nanoAOD files and defined as \f$w_{variation}\f$ / \f$w_{nominal}\f$. 
 * The nominal weight is already applied, therefore, the main use of this function is 
 * to get the variation of the PDF weights to the nominal PDF weight.
 *
 * The PDF weights consist of 101 weights, where the first weight is the nominal weight 
 * and the remaining 100 weights correspond to alternative PDF sets. 
 *
 * @note The proper procedure is to use each alternative PDF set as an independent 
 * systematic vatiation. However, in case of this function, a simplified approach is used 
 * to calculate a single PDF weight variation. The standard deviation of the 100 
 * alternative PDF weights is calculated and used to define the up and down variations as 
 * follows: \f$w_{up/down} = 1 \pm \sqrt{\sum_{i=1}^{100} (w_i - 1)^2}\f$
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the LHE PDF event weight
 * @param lhe_pdf_weights name of the column containing the LHE PDF weights
 * @param variation name of the variation that should be evaluated, possible values 
 * are "nominal", "up", "down"
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode LHEpdf(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &lhe_pdf_weights,
                        const std::string &variation) {
    auto lhe_pdf_weights_lambda =
        [variation](const ROOT::RVec<float> pdf_weights) {
            // the nominal weight is already applied, so we can return 1.0
            if (variation == "nominal") {
                return (float)1.0;
            }
            const int n_pdfs = pdf_weights.size();
            if (n_pdfs == 101 || n_pdfs == 103) {
                float sum = 0.0;
                for (size_t i = 1; i < 101; i++) {
                    float diff = pdf_weights[i] - 1;
                    sum += diff * diff;        
                }
                if (variation == "up") {
                    return (float)(1.0 + std::sqrt(sum));
                } else if (variation == "down") {
                    return (float)(1.0 - std::sqrt(sum));
                } else {
                    Logger::get("event::reweighting::LHEpdf")
                        ->error("Invalid variation: {}", variation);
                    throw std::runtime_error("Invalid variation for LHE PDF weights");
                }
            } else {
                Logger::get("event::reweighting::LHEpdf")
                    ->error("Invalid number of LHE PDF weights: {}",
                            n_pdfs);
                throw std::runtime_error("Invalid number of LHE PDF weights");
            }
        };

    auto df1 =
        df.Define(outputname, lhe_pdf_weights_lambda, {lhe_pdf_weights});
    return df1;
}

/**
 * @brief This function is used to evaluate the LHE \f$\alpha_S\f$ weight of an event. 
 * The weights are stored in the nanoAOD files and defined as 
 * \f$w_{variation}\f$ / \f$w_{nominal}\f$. The nominal weight is already applied, 
 * therefore, the main use of this function is to get the variation of the \f$\alpha_S\f$ 
 * weight to the nominal weight.
 *
 * For some samples the \f$\alpha_S\f$ weight is included in the PDF weights vector. In 
 * that case the full PDF weights vector is expected to contains 103 entries, where the
 * first 101 entries are PDF weights and the last two entries correspond to the up and 
 * down varied \f$\alpha_S\f$ weight.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the LHE \f$\alpha_S\f$ event weight
 * @param lhe_pdf_weights name of the column containing the LHE \f$\alpha_S\f$ weights (it is 
 * part of the LHE PDF weights)
 * @param variation name of the variation that should be evaluated, possible values 
 * are "nominal", "up", "down"
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode LHEalphaS(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &lhe_pdf_weights,
                        const std::string &variation) {
    auto lhe_alphaS_weights_lambda =
        [variation](const ROOT::RVec<float> pdf_weights) {
            // the nominal weight is already applied, so we can return 1.0
            if (variation == "nominal") {
                return (float)1.0;
            }
            const int n_pdfs = pdf_weights.size();
            if (n_pdfs == 103) {
                if (variation == "down") {
                    return pdf_weights[101];
                } else if (variation == "up") {
                    return pdf_weights[102];
                } else {
                    Logger::get("event::reweighting::LHEalphaS")
                        ->error("Invalid variation: {}", variation);
                    throw std::runtime_error("Invalid variation for LHE PDF weights");
                }
            } else {
                Logger::get("event::reweighting::LHEalphaS")
                    ->debug("LHE PDF weights do not include alphaS uncertainty: {}",
                            n_pdfs);
                return (float)1.0;            
            }
        };

    auto df1 =
        df.Define(outputname, lhe_alphaS_weights_lambda, {lhe_pdf_weights});
    return df1;
}

/**
 * @brief This function is used to calculate an event weight to correct the top
 * quark \f$p_T\f$ mismodeling in simulated \f$t\bar{t}\f$ events. The correction
 * is provided by the Top POG and in case of this function the calculated weight 
 * corrects NLO simulation (POWHEG+Pythia8) to data.
 *
 * For reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
 *
 * The weight is calculated as \f$w=\sqrt{SF(t)\cdot SF(\bar{t})}\f$ 
 *
 * with \f$SF= \exp(0.0615-0.0005\cdot p_T)\f$
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the derived event weight
 * @param genparticles_pdg_id name of the column containing the PDG IDs of the generator
 * particles
 * @param genparticles_status_flags name of the column containing the status flags of the
 * generator particles, where bit 13 contains the isLastCopy flag
 * @param genparticles_pt name of the column containing the pt of the generator particles
 *
 * @return a new dataframe containing the new column
 *
 * @note The Top POG also provides other reweighting functions, e.g. for NNLO to 
 * data or NLO to NNLO which could be preferred depending on the use case.
 */
ROOT::RDF::RNode TopPt(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &genparticles_pdg_id,
                        const std::string &genparticles_status_flags,
                        const std::string &genparticles_pt) {
    // In nanoAODv12 the type of genparticle status flags was changed to UShort_t
    // For v9 compatibility a type casting is applied
    auto [df1, genparticles_status_flags_column] = utility::Cast<ROOT::RVec<UShort_t>, ROOT::RVec<Int_t>>(
            df, genparticles_status_flags+"_v12", "ROOT::VecOps::RVec<UShort_t>", genparticles_status_flags);

    auto ttbarreweightlambda = [](const ROOT::RVec<int> pdg_ids,
                                  const ROOT::RVec<UShort_t> status_flags_v12,
                                  const ROOT::RVec<float> pts) {
        auto status_flags = static_cast<ROOT::RVec<int>>(status_flags_v12);
        std::vector<float> top_pts;
        for (size_t i = 0; i < pdg_ids.size(); i++) {
            if (std::abs(pdg_ids[i]) == 6 && ((status_flags[i] >> 13) & 1) == 1)
                top_pts.push_back(pts[i]);
        }
        if (top_pts.size() != 2) {
            std::cout << top_pts.size();
            Logger::get("event::reweighting::TopPt")
                ->error("TTbar reweighting applied to event with not exactly "
                        "two top quarks. Probably due to wrong sample type. "
                        "n_top: {}",
                        top_pts.size());
            throw std::runtime_error("Bad number of top quarks.");
        }

        if (top_pts[0] > 500.0)
            top_pts[0] = 500.0;
        if (top_pts[1] > 500.0)
            top_pts[1] = 500.0;
        const float parameter_a = 0.0615;
        const float parameter_b = -0.0005;
        return sqrt(exp(parameter_a + parameter_b * top_pts[0]) *
                    exp(parameter_a + parameter_b * top_pts[1]));
    };
    auto df2 = df1.Define(outputname, ttbarreweightlambda,
                {genparticles_pdg_id, genparticles_status_flags_column, genparticles_pt});
    return df2;
}

/**
 * @brief This function is used to calculate an event weight to correct the Z boson 
 * \f$p_T\f$. These corrections are recommended especially
 * for LO Drell-Yan samples, where the \f$p_T\f$ and mass of the Z boson are
 * mismodeled compared to data. This function is defined for the corrections provided
 * by the CMS HLepRare group. More details can be found here:
 * https://cms-higgs-leprare.docs.cern.ch/htt-common/DY_reweight/
 *
 * @note HLepRare only provides corrections for Run3. For Run2 see
 *`event::reweighting::ZPtMass`.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column containing the derived event weight
 * @param gen_boson name of the column containing the Lorentz vector of the
 * generator-level boson
 * @param corr_file path to the correction file containing the Z boson \f$p_T\f$
 * corrections
 * @param corr_name name of the correction in the json file
 * @param order order of the used DY samples: "LO" for madgraph, "NLO" for amcatnlo,
 * "NNLO" for powheg
 * @param variation name of the variation that should be evaluated, options are
 * "nom", "up", "down" or "upX", "downX". For "up" and "down" the uncertainty is
 * defined by the envelope of all provided uncertainty sources in the correction
 * file. Otherwise the specific uncertainty source "X" is used (where X is a number
 * e.g. 1,2,3,...).
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
ZBosonPt(ROOT::RDF::RNode df,
         correctionManager::CorrectionManager &correction_manager,
         const std::string &outputname,
         const std::string &gen_boson,
         const std::string &corr_file,
         const std::string &corr_name,
         const std::string &order,
         const std::string &variation) {

    auto corrDY = correction_manager.loadCorrection(corr_file, corr_name);
    auto corrDYunc = correction_manager.loadCorrection(corr_file, corr_name + "_N_uncertainty");

    // Get number of uncertainties
    auto n_unc = corrDYunc->evaluate({order});

    // Get weight depending on variation
    auto corr = [corrDY, order, variation, n_unc](const ROOT::Math::PtEtaPhiMVector &gen_boson) {
        float weight = 1.0;
        if (variation == "up" || variation == "down") {
            float weight_nom = corrDY->evaluate({order, (float)gen_boson.Pt(), "nom"});
            weight = weight_nom;

            for (int i = 0; i < n_unc; i++) {
                float unc_up = corrDY->evaluate({order, (float)gen_boson.Pt(), "up" + std::to_string(i)});
                float unc_down = corrDY->evaluate({order, (float)gen_boson.Pt(), "down" + std::to_string(i)});

                if (variation == "up" && unc_up > weight_nom && unc_up > weight) {
                    weight = unc_up;
                } else if (variation == "up" && unc_down > weight_nom && unc_down > weight) {
                    weight = unc_down;
                } else if (variation == "down" && unc_down < weight_nom && unc_down < weight) {
                    weight = unc_down;
                } else if (variation == "down" && unc_up < weight_nom && unc_up < weight) {
                    weight = unc_up;
                }
            }
        } 
        else if (variation == "nom" || variation.rfind("up") != std::string::npos || 
                 variation.rfind("down") != std::string::npos) {
            weight = corrDY->evaluate({order, (float)gen_boson.Pt(), variation});
        }
        else {
            Logger::get("event::reweighting::ZBosonPt")
                ->error("Invalid variation: {}", variation);
            throw std::runtime_error("Invalid variation for Z boson pT weight");
        }
        return weight;
    };
    auto df1 = df.Define(outputname, corr, {gen_boson});
    return df1;
}

/**
 * @brief This function is used to calculate an event weight based on Z boson 
 * \f$p_T\f$ and mass corrections. These corrections are recommended especially 
 * for LO Drell-Yan samples, where the \f$p_T\f$ and mass of the Z boson are
 * mismodeled compared to data. 
 *
 * @warning This function is based on workspaces and functions that were derived
 * for the legacy \f$H(\tau\tau)\f$ analysis and therefore not up-to-date anymore 
 * for UL or Run3.
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the derived event weight
 * @param gen_boson name of the column containing the Lorentz vector of the 
 * generator-level boson
 * @param workspace_file path to the file which contains the workspace that should be
 * used
 * @param functor_name name of the function in the workspace that should be used
 * @param argset additional arguments that are needed for the function
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode ZPtMass(ROOT::RDF::RNode df,
                         const std::string &outputname,
                         const std::string &gen_boson,
                         const std::string &workspace_file,
                         const std::string &functor_name,
                         const std::string &argset) {
    // retrieve pt and mass of gen boson reconstructed with the method used by
    // recoil corrections; resulting quantities are only for the purpose of this
    // method
    auto df1 = df.Define(gen_boson + "_pt",
                         [](const ROOT::Math::PtEtaPhiMVector &p4) {
                             return (float)p4.pt();
                         },
                         {gen_boson});
    auto df2 = df1.Define(gen_boson + "_mass",
                          [](const ROOT::Math::PtEtaPhiMVector &p4) {
                              return (float)p4.mass();
                          },
                          {gen_boson});

    // set up workspace
    Logger::get("event::reweighting::ZPtMass")
        ->debug("Setting up functions for zPtMassReweighting");
    Logger::get("event::reweighting::ZPtMass")
        ->debug("zPtMassReweighting - Function {} // argset {}", functor_name,
                argset);

    const std::shared_ptr<RooFunctorThreadsafe> weight_function =
        loadFunctor(workspace_file, functor_name, argset);
    auto df3 = utility::EvaluateWorkspaceFunction(
        df2, outputname, weight_function, gen_boson + "_mass",
        gen_boson + "_pt");
    return df3;
}
} // end namespace reweighting
} // end namespace event
#endif /* GUARD_REWEIGHTING_H */
