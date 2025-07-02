#ifndef GUARD_EMBEDDING_H
#define GUARD_EMBEDDING_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include <nlohmann/json.hpp>

namespace embedding {

namespace scalefactor {

/**
 * @brief This function calculates a scale factor (SF) that corrects the bias 
 * in the trigger selection of events for the embedding method. The scale factors 
 * are measured by the tau embedding group and are loaded from a correctionlib 
 * file using a specified scale factor name.
 *
 * For the calculation the two initially selected muons in data need to be 
 * provided.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column containing the embedding 
 * selection trigger scale factor
 * @param pt_1 \f$p_T\f$ of the leading generator particle which corresponds 
 * to the leading muon selected by the embedding selection
 * @param eta_1 \f$\eta\f$ of the leading generator particle which corresponds 
 * to the leading muon selected by the embedding selection
 * @param pt_2 \f$p_T\f$ of the subleading generator particle which corresponds 
 * to the subleading muon selected by the embedding selection
 * @param eta_2 \f$\eta\f$ of the subleading generator particle which corresponds 
 * to the subleading muon selected by the embedding selection
 * @param sf_file path to the file with the embedding selection trigger scale 
 * factors
 * @param sf_name name of the embedding selection trigger scale factor correction
 * e.g. "m_sel_trg_kit_ratio"
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
SelectionTrigger(ROOT::RDF::RNode df,
                  correctionManager::CorrectionManager &correction_manager,
                  const std::string &outputname,
                  const std::string &pt_1, const std::string &eta_1,
                  const std::string &pt_2, const std::string &eta_2,
                  const std::string &sf_file,
                  const std::string &sf_name) {
    Logger::get("embedding::scalefactor::SelectionTrigger")
        ->debug("Correction - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator](const float &pt_1, const float &eta_1, const float &pt_2,
                    const float &eta_2) {
            Logger::get("embedding::scalefactor::SelectionTrigger")
                ->debug(" pt_1 {}, eta_1 {}, pt_2 {}, eta_2 {}", pt_1, eta_1,
                        pt_2, eta_2);
            double sf = 1.;
            sf = evaluator->evaluate(
                {pt_1, std::abs(eta_1), pt_2, std::abs(eta_2)});
            Logger::get("embedding::scalefactor::SelectionTrigger")->debug("sf {}", sf);
            return sf;
        },
        {pt_1, eta_1, pt_2, eta_2});
    return df1;
}

/**
 * @brief This function calculates a scale factor (SF) that corrects the bias 
 * in the ID selection of muons in events for the embedding method. The scale 
 * factors are measured by the tau embedding group and are loaded from a 
 * correctionlib file using a specified scale factor name.
 *
 * For the calculation the one of the two initially selected muons in data 
 * need to be provided.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column containing the embedding 
 * selection ID scale factor
 * @param pt \f$p_T\f$ of a generator particle which corresponds to
 * one of the muons selected by the embedding selection
 * @param eta \f$\eta\f$ of a generator particle which corresponds to
 * one of the muons selected by the embedding selection
 * @param sf_file tpath to the file with the embedding selection ID scale 
 * factors
 * @param sf_name name of the embedding selection trigger scale factor correction
 * e.g. "EmbID_pt_eta_bins"
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
SelectionId(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correction_manager,
             const std::string &outputname,
             const std::string &pt, const std::string &eta,
             const std::string &sf_file,
             const std::string &sf_name) {
    Logger::get("embedding::scalefactor::SelectionId")
        ->debug("Correction - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 =
        df.Define(outputname,
                  [evaluator](const float &pt, const float &eta) {
                      Logger::get("embedding::scalefactor::SelectionId")
                          ->debug(" pt {}, eta {},", pt, eta);
                      double sf = 1.;
                      sf = evaluator->evaluate({pt, std::abs(eta)});
                      Logger::get("embedding::scalefactor::SelectionId")->debug("sf {}", sf);
                      return sf;
                  },
                  {pt, eta});
    return df1;
}
} // end namespace scalefactor

namespace muon {

/**
 * @brief This function calculates muon scale factors (SFs) that were 
 * measured by the tau embedding group for embedding samples. The scale factors 
 * are loaded from a correctionlib file using a specified scale factor name 
 * which defines the correction (e.g. ID, Iso, Trigger).
 *
 * The tau embedding group also measured muon scale factors for simulated
 * samples with the same T&P method as for embedding. Which scale factor is 
 * used can be specified by the `correction_type` parameter.
 *
 * The variation is introduced via an extrapolation factor. This factor is 
 * applied on top of the scale factor, thereby varying the scale factor by a
 * fixed percentage.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * muon scale factor file
 * @param outputname name of the output column containing the muon scale factor
 * @param pt name of the column containing the transverse momentum of a muon
 * @param eta name of the column containing the pseudorapidity of a muon
 * @param sf_file path to the file containing the muon scale factors
 * @param sf_name name of the muon scale factor correction, this parameter
 * defines what is corrected ID, Iso or Trigger, e.g. "Trg_IsoMu27_pt_eta_bins"
 * @param correction_type type of the correction, there are two options: `emb` 
 * for embedding and `mc` for simulation
 * @param extrapolation_factor factor that is used to introduce percent type 
 * variations of the scale factor, default value is `1.0`
 *
 * @return a new dataframe containing the new column
 *
 * @note If possible, the scale factors from the CMS POGs should be used for
 * simulated samples because they have a more sofisticated treatment of
 * systematic uncertainties.
 */
ROOT::RDF::RNode
Scalefactor(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &pt, const std::string &eta,
        const std::string &sf_file, const std::string &sf_name,
        const std::string correction_type,
        const float &extrapolation_factor = 1.0) {

    Logger::get("embedding::muon::Scalefactor")->debug("Correction - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, correction_type, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("embedding::muon::Scalefactor")
                ->debug(" pt {}, eta {}, correction_type {}, extrapolation "
                        "factor {}",
                        pt, eta, correction_type, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, std::abs(eta), correction_type});
            Logger::get("embedding::muon::Scalefactor")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
} // end namespace muon

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
    auto correction_lambda = [sf_barrel,
                              sf_endcap](const ROOT::RVec<float> &pts,
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
                ->debug("ele pt before {}, ele pt after {}", pts.at(i),
                        corrected_pts.at(i));
        }
        return corrected_pts;
    };
    auto df1 = df.Define(outputname, correction_lambda, {pt, eta});
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
    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda = [evaluator, variation_barrel,
                              variation_endcap](const ROOT::RVec<float> &pts,
                                                const ROOT::RVec<float> &etas) {
        ROOT::RVec<float> corrected_pts(pts.size());
        for (int i = 0; i < pts.size(); i++) {
            if (abs(etas.at(i)) <= 1.479) {
                auto correction_factor =
                    evaluator->evaluate({"barrel", variation_barrel});
                corrected_pts[i] = pts.at(i) * correction_factor;
                Logger::get("embedding::electron::PtCorrection")
                    ->debug(
                        "barrel: ele pt before {}, ele pt after {}, factor {}",
                        pts.at(i), corrected_pts.at(i), correction_factor);
            } else if (abs(etas.at(i)) > 1.479) {
                auto correction_factor =
                    evaluator->evaluate({"endcap", variation_endcap});
                corrected_pts[i] = pts.at(i) * correction_factor;
                Logger::get("embedding::electron::PtCorrection")
                    ->debug(
                        "endcap: ele pt before {}, ele pt after {}, factor {}",
                        pts.at(i), corrected_pts.at(i), correction_factor);
            } else {
                corrected_pts[i] = pts.at(i);
            }
        }
        return corrected_pts;
    };
    auto df1 = df.Define(outputname, correction_lambda, {pt, eta});
    return df1;
}

/**
 * @brief This function calculates electron scale factors (SFs) that were 
 * measured by the tau embedding group for embedding samples. The scale factors 
 * are loaded from a correctionlib file using a specified scale factor name 
 * which defines the correction (e.g. ID, Iso, Trigger).
 *
 * The tau embedding group also measured electron scale factors for simulated
 * samples with the same T&P method as for embedding. Which scale factor is 
 * used can be specified by the `correction_type` parameter.
 *
 * The variation is introduced via an extrapolation factor. This factor is 
 * applied on top of the scale factor, thereby varying the scale factor by a
 * fixed percentage.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * electron scale factor file
 * @param outputname name of the output column containing the electron scale factor
 * @param pt name of the column containing the transverse momentum of an electron
 * @param eta name of the column containing the pseudorapidity of an electron
 * @param sf_file path to the file containing the electron scale factors
 * @param sf_name name of the electron scale factor correction, this parameter
 * defines what is corrected ID, Iso or Trigger, e.g. "ID90_pt_bins_inc_eta"
 * @param correction_type type of the correction, there are two options: `emb` 
 * for embedding and `mc` for simulation
 * @param extrapolation_factor factor that is used to introduce percent type 
 * variations of the scale factor, default value is `1.0`
 *
 * @return a new dataframe containing the new column
 *
 * @note If possible, the scale factors from the CMS POGs should be used for
 * simulated samples because they have a more sofisticated treatment of
 * systematic uncertainties.
 */
ROOT::RDF::RNode
Scalefactor(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correction_manager,
            const std::string &outputname,
            const std::string &pt, const std::string &eta,
            const std::string &sf_file, const std::string &sf_name,
            const std::string correction_type, 
            const float &extrapolation_factor = 1.0) {

    Logger::get("embedding::electron::Scalefactor")
        ->debug("Correction - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, correction_type, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("embedding::electron::Scalefactor")
                ->debug(" pt {}, eta {}, correction_type {}, extrapolation "
                        "factor {}",
                        pt, eta, correction_type, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, eta, correction_type});
            Logger::get("embedding::electron::Scalefactor")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
} // end namespace electron

namespace tau {

/**
 * @brief This function scales the \f$p_T\f$ values of hadronic taus based on
 * their decay mode. The correction factors are provided as input parameters.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the corrected \f$p_T\f$
 * values
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
    auto correction_lambda = [sf_dm0, sf_dm1, sf_dm10,
                              sf_dm11](const ROOT::RVec<float> &pts,
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
                ->debug("tau pt before {}, tau pt after {}", pts.at(i),
                        corrected_pts.at(i));
        }
        return corrected_pts;
    };
    auto df1 = df.Define(outputname, correction_lambda, {pt, decay_mode});
    return df1;
}

namespace scalefactor {

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against jets (`vsJet`) for embedding samples. The scale factors are loaded from 
 * a correctionlib file using a specified scale factor name and variation. The 
 * variation and the scale factor itself is binned in transverse momenta (\f$p_T\f$) 
 * of hadronic taus for this function. This dependence is usually used in semi-leptonic 
 * channels (\f$e\tau_h\f$, \f$\mu\tau_h\f$).
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
 * @param variation_pt20to25 name of the scale factor variation for \f$20 \leq p_T <25\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt25to30 name of the scale factor variation for \f$25 \leq p_T <30\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt30to35 name of the scale factor variation for \f$30 \leq p_T <35\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt35to40 name of the scale factor variation for \f$35 \leq p_T <40\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt40toInf name of the scale factor variation for \f$40 \leq p_T < \infty \f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 *
 * @note The only differnce to the `physicsobject::tau::scalefactor::Id_vsJet_lt` function,
 * which is used for simulated samples, is the \f$p_T\f$ binning of the variations.
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
            const std::string &variation_pt20to25,
            const std::string &variation_pt25to30,
            const std::string &variation_pt30to35,
            const std::string &variation_pt35to40,
            const std::string &variation_pt40toInf) {
    
    const std::map<float, std::string> variations = {
        {20.0f,     variation_pt20to25},
        {25.0f,     variation_pt25to30},
        {30.0f,     variation_pt30to35},
        {35.0f,     variation_pt35to40},
        {40.0f,     variation_pt40toInf},
        {100000.0f, variation_pt40toInf},
    };
    Logger::get("embedding::tau::scalefactor::Id_vsJet_lt")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("embedding::tau::scalefactor::Id_vsJet_lt")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, vsele_wp, variations,
                            sf_dependence, selected_dms,
                            sf_name](const float &pt, const int &decay_mode,
                                         const int &gen_match) {
        Logger::get("embedding::tau::scalefactor::Id_vsJet_lt")->debug("ID - decayMode {}", decay_mode);
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
                Logger::get("embedding::tau::scalefactor::Id_vsJet_lt")
                    ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation, sf_dependence);
                sf = evaluator->evaluate({pt, decay_mode, gen_match, wp, vsele_wp, variation, sf_dependence});
            } else {
                sf = 1.;
            }
        }
        Logger::get("embedding::tau::scalefactor::Id_vsJet_lt")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, gen_match});
    return df1;
}
} // end namespace scalefactor
} // end namespace tau
} // end namespace embedding

#endif /* GUARD_EMBEDDING_H */