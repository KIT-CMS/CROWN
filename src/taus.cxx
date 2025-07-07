#ifndef GUARD_TAUS_H
#define GUARD_TAUS_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "../include/defaults.hxx"
#include "ROOT/RDataFrame.hxx"
#include "correction.h"
#include <unordered_map>

namespace physicsobject {
namespace tau {

/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of hadronic taus that originate from electrons that are
 * misidentified. The energy scale correction for these objects is measured for
 * two tau decay modes (dm0 and dm1) and depends on the barrel and endcap region
 * of the detector.
 *
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
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
 * @note This correction is only applied to misidentified hadronic taus
 * originating from prompt electrons (`gen_match=1`) and electrons that decayed
 * from a tau lepton
 * (`gem_match=3`).
 */
ROOT::RDF::RNode
PtCorrectionMC_eleFake(ROOT::RDF::RNode df,
                       correctionManager::CorrectionManager &correction_manager,
                       const std::string &outputname, const std::string &pt,
                       const std::string &eta, const std::string &decay_mode,
                       const std::string &gen_match, const std::string &es_file,
                       const std::string &correction_name,
                       const std::string &id_algorithm,
                       const std::string &variation_dm0_barrel,
                       const std::string &variation_dm1_barrel,
                       const std::string &variation_dm0_endcap,
                       const std::string &variation_dm1_endcap) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda =
        [evaluator, id_algorithm, variation_dm0_barrel, variation_dm1_barrel,
         variation_dm0_endcap, variation_dm1_endcap](
            const ROOT::RVec<float> &pts, const ROOT::RVec<float> &etas,
            const ROOT::RVec<UChar_t> &decay_modes_v12,
            const ROOT::RVec<UChar_t> &gen_matches_char) {
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
            ROOT::RVec<float> corrected_pts(pts.size());
            const float barrel_end_cut = 1.5;
            const float endcap_end_cut = 2.5;
            const int dm0 = 0;
            const int dm1 = 1;

            for (int i = 0; i < pts.size(); i++) {
                if (gen_matches.at(i) == 1 || gen_matches.at(i) == 3) {
                    if (decay_modes.at(i) == dm0 &&
                        std::abs(etas.at(i)) <= barrel_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm0_barrel});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == dm0 &&
                               std::abs(etas.at(i)) > barrel_end_cut &&
                               std::abs(etas.at(i)) <= endcap_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm0_endcap});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == dm1 &&
                               std::abs(etas.at(i)) <= barrel_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm1_barrel});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == dm1 &&
                               std::abs(etas.at(i)) > barrel_end_cut &&
                               std::abs(etas.at(i)) <= endcap_end_cut) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm1_endcap});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    }
                } else {
                    corrected_pts[i] = pts.at(i);
                }
                Logger::get("physicsobject::tau::PtCorrectionMC_eleFake")
                    ->debug("tau pt before {}, tau pt after {}", pts.at(i),
                            corrected_pts.at(i));
            }
            return corrected_pts;
        };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}

/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of hadronic taus that originate from muons that are
 * misidentified. The energy scale correction is always set to `1` and for the
 * uncertainty varied by `1%`.
 *
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param variation variation of the correction, options are "nom", "up", "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to misidentified hadronic taus
 * originating from prompt muons (`gen_match=2`) and muons that decayed from a
 * tau lepton
 * (`gem_match=4`).
 */
ROOT::RDF::RNode
PtCorrectionMC_muFake(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correction_manager,
                      const std::string &outputname, const std::string &pt,
                      const std::string &eta, const std::string &decay_mode,
                      const std::string &gen_match, const std::string &es_file,
                      const std::string &correction_name,
                      const std::string &id_algorithm,
                      const std::string &variation) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda = [evaluator, id_algorithm, variation](
                                 const ROOT::RVec<float> &pts,
                                 const ROOT::RVec<float> &etas,
                                 const ROOT::RVec<UChar_t> &decay_modes_v12,
                                 const ROOT::RVec<UChar_t> &gen_matches_char) {
        auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
        auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
        ROOT::RVec<float> corrected_pts(pts.size());
        for (int i = 0; i < pts.size(); i++) {
            if (gen_matches.at(i) == 2 || gen_matches.at(i) == 4) {
                auto correction_factor = evaluator->evaluate(
                    {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                     gen_matches.at(i), id_algorithm,
                     variation});
                corrected_pts[i] = pts.at(i) * correction_factor;
            } else {
                corrected_pts[i] = pts.at(i);
            }
            Logger::get("physicsobject::tau::PtCorrectionMC_muFake")
                ->debug("tau pt before {}, tau pt after {}", pts.at(i),
                        corrected_pts.at(i));
        }
        return corrected_pts;
    };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}

/**
 * @brief This function corrects the transverse momentum (\f$p_T\f$) in MC
 * simulations of genuine hadronic taus. The energy scale correction for these
 * objects is measured for four tau decay modes (dm0, dm1, dm10 and dm11) and
 * depends on the transverse momentum of the hadronic tau.
 *
 * The correction procedure is taken from the officially recommendation of the
 * TauPOG:
 *
 * Run2 (UL):
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2018_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2017_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016postVFP_UL_tau.html
 * -
 * https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/TAU_2016preVFP_UL_tau.html
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
 * (not added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * correction file
 * @param outputname name of the output column storing the corrected hadronic
 * tau \f$p_T\f$ values
 * @param pt name of the input column containing hadronic tau \f$p_T\f$ values
 * @param eta name of the column containing hadronic tau eta values
 * @param decay_mode name of the column containing hadronic tau decay modes
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param es_file path to the correction file for the energy scale correction
 * @param correction_name name of the correction in `es_file`
 * @param id_algorithm identification algorithm used for hadronic tau ID
 * @param variation_dm0 variation for decay mode 0, options are "nom", "up",
 * "down"
 * @param variation_dm1 variation for decay mode 1, options are "nom", "up",
 * "down"
 * @param variation_dm10 variation for decay mode 10, options are "nom", "up",
 * "down"
 * @param variation_dm11 variation for decay mode 11, options are "nom", "up",
 * "down"
 *
 * @return a dataframe containing the corrected transverse momenta
 *
 * @note This correction is only applied to genuine hadronic taus
 * (`gen_match=5`).
 */
ROOT::RDF::RNode PtCorrectionMC_genuineTau(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &pt,
    const std::string &eta, const std::string &decay_mode,
    const std::string &gen_match, const std::string &es_file,
    const std::string &correction_name, const std::string &id_algorithm,
    const std::string &variation_dm0, const std::string &variation_dm1,
    const std::string &variation_dm10, const std::string &variation_dm11) {
    // In nanoAODv12 the type of tau decay mode was changed to UChar_t
    // For v9 compatibility a type casting is applied
    auto [df1, decay_mode_column] = utility::Cast<ROOT::RVec<UChar_t>, ROOT::RVec<Int_t>>(
            df, decay_mode+"_v12", "ROOT::RVec<UChar_t>", decay_mode);

    auto evaluator =
        correction_manager.loadCorrection(es_file, correction_name);
    auto correction_lambda =
        [evaluator, id_algorithm, variation_dm0, variation_dm1, variation_dm10,
         variation_dm11](const ROOT::RVec<float> &pts,
                         const ROOT::RVec<float> &etas,
                         const ROOT::RVec<UChar_t> &decay_modes_v12,
                         const ROOT::RVec<UChar_t> &gen_matches_char) {
            auto decay_modes = static_cast<ROOT::RVec<int>>(decay_modes_v12);
            auto gen_matches = static_cast<ROOT::RVec<int>>(gen_matches_char);
            ROOT::RVec<float> corrected_pts(pts.size());
            for (int i = 0; i < pts.size(); i++) {
                if (gen_matches.at(i) == 5) {
                    if (decay_modes.at(i) == 0) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm0});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 1) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm1});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 10) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm10});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    } else if (decay_modes.at(i) == 11) {
                        auto correction_factor = evaluator->evaluate(
                            {pts.at(i), std::abs(etas.at(i)), decay_modes.at(i),
                             gen_matches.at(i), id_algorithm,
                             variation_dm11});
                        corrected_pts[i] = pts.at(i) * correction_factor;
                    }
                } else {
                    corrected_pts[i] = pts.at(i);
                }
                Logger::get("physicsobject::tau::PtCorrection_genuineTau")
                    ->debug("tau pt before {}, tau pt after {}, decaymode {}",
                            pts.at(i), corrected_pts.at(i), decay_modes.at(i));
            }
            return corrected_pts;
        };
    auto df2 = df1.Define(outputname, correction_lambda,
                         {pt, eta, decay_mode_column, gen_match});
    return df2;
}

namespace quantity {

/**
 * @brief This function writes out a flag that is true if a tau passes a specific 
 * tau ID working point (`wp`). The working point is defined by the bit value. The
 * bit values can be found e.g. in the description of the tau ID scale factors. 
 * 
 * @param df input dataframe
 * @param outputname name of the output column containing the flag
 * @param ID name of the column containing the tau ID values 
 * @param index_vector name of the column containing the vector with the relevant
 * tau pair indices
 * @param position position in the index vector of the relevant tau in the pair
 * @param wp bit value of the WP that has to be passed
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode IDFlag(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &ID, const std::string &index_vector,
                        const int &position, const int &wp) {
    return df.Define(
        outputname,
        [position, wp](const ROOT::RVec<int> &index_vector,
                          const ROOT::RVec<UChar_t> &IDs) {
            Logger::get("physicsobject::tau::quantity::IDFlag")
                ->debug(
                    "position tau in pair {}, tau pair {}, id wp bit {}, vsjet ids {}",
                    position, index_vector, wp, IDs);
            const int index = index_vector.at(position);
            const int id_value = IDs.at(index, default_int);
            if (id_value != default_int)
                return std::min(1, int(id_value & 1 << (wp - 1)));
            else
                return int(id_value);
        },
        {index_vector, ID});
}
} // end namespace quantity

namespace scalefactor {

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against jets (`vsJet`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in transverse momenta (\f$p_T\f$) of hadronic taus 
 * for this function. This dependence is usually used in semi-leptonic channels 
 * (\f$e\tau_h\f$, \f$\mu\tau_h\f$).
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
 * @param variation_pt30to35 name of the scale factor variation for \f$30 \leq p_T <35\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt35to40 name of the scale factor variation for \f$35 \leq p_T <40\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt40to500 name of the scale factor variation for \f$40 \leq p_T <500\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt500to1000 name of the scale factor variation for \f$500 \leq p_T <1000\f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_pt1000toInf name of the scale factor variation for \f$1000 \leq p_T < \infty \f$ GeV, "nom" for nominal
 * and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
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
            const std::string &variation_pt30to35,
            const std::string &variation_pt35to40,
            const std::string &variation_pt40to500,
            const std::string &variation_pt500to1000,
            const std::string &variation_pt1000toInf) {
    
    const std::map<float, std::string> variations = {
        {30.0f,     variation_pt30to35},
        {35.0f,     variation_pt35to40},
        {40.0f,     variation_pt40to500},
        {500.0f,    variation_pt500to1000},
        {1000.0f,   variation_pt1000toInf},
        {100000.0f, variation_pt1000toInf},
    };
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, vsele_wp, variations,
                            sf_dependence, selected_dms,
                            sf_name](const float &pt, const int &decay_mode,
                                         const int &gen_match) {
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("ID - decayMode {}", decay_mode);
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
                Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation, sf_dependence);
                sf = evaluator->evaluate({pt, decay_mode, gen_match, wp, vsele_wp, variation, sf_dependence});
            } else {
                sf = 1.;
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_lt")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against jets (`vsJet`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in decay modes of hadronic taus for this function. This 
 * dependence is usually used in the fully hadronic channel (\f$\tau_h\tau_h\f$).
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
 * @param wp working point of the vsJet ID
 * @param vsele_wp working point of the vsEle ID
 * @param sf_dependence variable dependence of the scale factor, opions are "pt" or "dm"
 * @param variation_dm0 name of the scale factor variation for decay mode 0, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm1 name of the scale factor variation for decay mode 1, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm10 name of the scale factor variation for decay mode 10, "nom" for nominal
 * and "up"/"down" the up/down variation
 * @param variation_dm11 name of the scale factor variation for decay mode 11, "nom" for nominal
 * and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Id_vsJet_tt(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, 
    const std::string &pt, const std::string &decay_mode,
    const std::string &gen_match, 
    const std::string &sf_file, const std::string &sf_name,
    const std::string &wp, const std::string &vsele_wp,
    const std::string &sf_dependence,
    const std::string &variation_dm0,
    const std::string &variation_dm1, 
    const std::string &variation_dm10,
    const std::string &variation_dm11) {
    
    const std::unordered_map<int, std::string> variations = {
        {0,  variation_dm0},
        {1,  variation_dm1},
        {10, variation_dm10},
        {11, variation_dm11},
    };
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, vsele_wp, variations,
                            sf_dependence, sf_name](const float &pt, const int &decay_mode,
                                         const int &gen_match) {
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("ID - decayMode {}", decay_mode);
        // only calculate SFs for allowed tau decay modes (also excludes default
        // values due to tau energy correction shifts below good tau pt
        // selection)
        double sf = 1.;
        if (auto it = variations.find(decay_mode); it != variations.end()) {
            std::string variation = it->second;
            Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")
                ->debug("ID {} - pt {}, decay_mode {}, gen_match {}, wp {}, "
                        "vsele_wp {}, variation {}, sf_dependence {}",
                        sf_name, pt, decay_mode, gen_match, wp, vsele_wp,
                        variation, sf_dependence);
            sf = evaluator->evaluate(
                {pt, decay_mode, gen_match, wp, vsele_wp, 
                variation, sf_dependence});
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsJet_tt")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against electrons (`vsEle`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in pseudorapidities of hadronic taus for this function.
 *
 * Description of the bit map used to define the tau ID against electrons working points of the
 * DeepTau v2.1 tagger.
 * vsElectrons                         | Value | Bit (value used in the config)
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
 * @param outputname name of the output column containing the vsEle ID scale factor
 * @param eta name of the column containing the pseudorapidity of a tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsEle ID correction
 * @param wp working point of the vsEle ID
 * @param variation_barrel name of the scale factor variation for the barrel region
 * (\f$|\eta| <1.46\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_endcap name of the scale factor variation for the endcap region
 * (\f$1.558 \leq |\eta| <2.3\f$), "nom" for nominal and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Id_vsEle(ROOT::RDF::RNode df,
         correctionManager::CorrectionManager &correction_manager,
         const std::string &outputname,
         const std::string &eta,
         const std::string &gen_match, 
         const std::string &sf_file, const std::string &sf_name,
         const std::string &wp, 
         const std::string &variation_barrel,
         const std::string &variation_endcap) {

    Logger::get("physicsobject::tau::scalefactor::Id_vsEle")
        ->debug("Setting up function for tau id vsEle sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsEle")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, variation_barrel, variation_endcap,
                            sf_name](const float &eta, const int &gen_match) {
        double sf = 1.;
        // exclude default values due to tau energy correction shifts below good tau
        // pt selection
        if (eta > -5.0) {
            Logger::get("physicsobject::tau::scalefactor::Id_vsEle")
                ->debug("ID {} - eta {}, gen_match {}, wp {}, variation_barrel "
                        "{}, variation_endcap {}",
                        sf_name, eta, gen_match, wp, variation_barrel,
                        variation_endcap);
            // the eta cuts are taken from the correctionlib json file to define 
            // barrel and endcap
            if (std::abs(eta) < 1.46) {
                sf = evaluator->evaluate(
                    {eta, gen_match, wp, variation_barrel});
            } else if (std::abs(eta) >= 1.558 && std::abs(eta) < 2.3) {
                sf = evaluator->evaluate(
                    {eta, gen_match, wp, variation_endcap});
            } else {
                sf = 1.;
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsEle")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 =
        df.Define(outputname, sf_calculator, {eta, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau identification (ID) 
 * against muons (`vsMu`). The scale factors are loaded from a correctionlib file 
 * using a specified scale factor name and variation. The variation and the scale 
 * factor itself is binned in pseudorapidities of hadronic taus for this function.
 *
 * Description of the bit map used to define the tau ID against muons working points of the
 * DeepTau v2.1 tagger.
 * vsMuons                             | Value | Bit (value used in the config)
 * ------------------------------------|-------|-------
 * no ID selection (takes every tau)   |  0    | -
 * VLoose                              |  1    | 1
 * Loose                               |  2    | 2
 * Medium                              |  4    | 3
 * Tight                               |  8    | 4
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the vsMu ID scale factor
 * @param eta name of the column containing the pseudorapidity of a tau
 * @param gen_match name of the column with the matching information of the
 * hadronic tau to generator-level particles (matches are: 1=prompt e, 2=prompt mu,
 * 3=tau->e, 4=tau->mu, 5=had. tau, 0=unmatched)
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the vsMu ID correction
 * @param wp working point of the vsMu ID
 * @param variation_wheel1 name of the scale factor variation for the muon wheel
 * (\f$|\eta| <0.4\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel2 name of the scale factor variation for the muon wheel
 * (\f$0.4 \leq |\eta| <0.8\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel3 name of the scale factor variation for the muon wheel
 * (\f$0.8 \leq |\eta| <1.2\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel4 name of the scale factor variation for the muon wheel
 * (\f$1.2 \leq |\eta| <1.7\f$), "nom" for nominal and "up"/"down" the up/down variation
 * @param variation_wheel5 name of the scale factor variation for the muon wheel
 * (\f$1.7 \leq |\eta| <2.3\f$), "nom" for nominal and "up"/"down" the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Id_vsMu(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &eta,
        const std::string &gen_match, 
        const std::string &sf_file,
        const std::string &sf_name,
        const std::string &wp,  
        const std::string &variation_wheel1,
        const std::string &variation_wheel2, 
        const std::string &variation_wheel3,
        const std::string &variation_wheel4, 
        const std::string &variation_wheel5) {
    
    const std::map<float, std::string> variations = {
        {0.0f, variation_wheel1},
        {0.4f, variation_wheel2},
        {0.8f, variation_wheel3},
        {1.2f, variation_wheel4},
        {1.7f, variation_wheel5},
        {2.3f, variation_wheel5},
    };
    Logger::get("physicsobject::tau::scalefactor::Id_vsMu")->debug("Setting up function for tau id vsMu sf");
    Logger::get("physicsobject::tau::scalefactor::Id_vsMu")->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto sf_calculator = [evaluator, wp, variations,
                            sf_name](const float &eta, const int &gen_match) {
        double sf = 1.;
        // exclude default values due to tau energy correction shifts below good tau
        // pt selection
        if (eta > -5.0) {
            auto it = variations.upper_bound(abs(eta));
            if (it != variations.begin()){
                it = std::prev(it);
                std::string variation = it->second;
                Logger::get("physicsobject::tau::scalefactor::Id_vsMu")
                    ->debug("ID {} - eta {}, gen_match {}, wp {}, variation {}",
                        sf_name, eta, gen_match, wp, variation);
                sf = evaluator->evaluate({eta, gen_match, wp, variation});
            } else {
                sf = 1.;
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Id_vsMu")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 =
        df.Define(outputname, sf_calculator, {eta, gen_match});
    return df1;
}

/**
 * @brief This function calculates scale factors (SFs) for tau triggers. The scale factors 
 * are loaded from a correctionlib file using a specified scale factor name and variation. 
 * The scale factor is binned in \f$p_T\f$ and decay modes of hadronic taus for this function.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * tau scale factor file
 * @param outputname name of the output column containing the trigger scale factor
 * @param pt name of the column containing the transverse momentum of a tau
 * @param decay_mode name of the column containing the decay mode of the tau
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau scale factor for the trigger correction
 * @param trigger_name name of the trigger, e.g. "ditau", "etau", "mutau"
 * @param wp working point of the vsJet ID
 * @param corr_type type of the value to be read out, options are "sf" (for 
 * scale factors), "eff_data", "eff_mc"
 * @param variation name the scale factor variation, "nom" for the nominal
 * scale factor and "up"/"down" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
Trigger(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &pt, const std::string &decay_mode, 
        const std::string &sf_file,
        const std::string &sf_name,
        const std::string &trigger_name, const std::string &wp,
        const std::string &corr_type, const std::string &variation) {

    Logger::get("physicsobject::tau::scalefactor::Trigger")
        ->debug("Setting up function for tau trigger sf");
    Logger::get("physicsobject::tau::scalefactor::Trigger")
        ->debug("ID - Name {}, file {}", sf_name, sf_file);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    Logger::get("physicsobject::tau::scalefactor::Trigger")->debug("WP {} - type {}", wp, corr_type);
    auto sf_calculator = [evaluator, trigger_name, wp, corr_type, variation, sf_name](
                                     const float &pt, const int &decay_mode) {
        float sf = 1.;
        if (pt > 0.) {
            Logger::get("physicsobject::tau::scalefactor::Trigger")
                ->debug("ID {} - decaymode {}, wp {} "
                    "pt {}, type {}, variation {}",
                    sf_name, decay_mode, wp, pt, corr_type, variation);
            if (decay_mode == 0 || decay_mode == 1 || decay_mode == 10 ||
                decay_mode == 11) {
                sf = evaluator->evaluate(
                    {pt, decay_mode, trigger_name, wp, corr_type, variation});
            } else {
                sf = evaluator->evaluate({pt, -1, trigger_name, wp, corr_type, variation});
            }
        }
        Logger::get("physicsobject::tau::scalefactor::Trigger")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(outputname, sf_calculator, {pt, decay_mode});
    return df1;
}
} // end namespace scalefactor
} // end namespace tau
} // end namespace physicsobject
#endif /* GUARD_TAUS_H */