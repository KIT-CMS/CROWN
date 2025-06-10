#ifndef GUARD_SCALEFACTORS_H
#define GUARD_SCALEFACTORS_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "RooFunctor.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "correction.h"
/// namespace used for scale factor related functions
namespace scalefactor {

namespace tau {
/**
 * @brief Function used to evaluate vsJets tau id scale factors in the lt
channel with the correctionlib for tauembedded samples


 * @param df The input dataframe
 * @param correctionManager The CorrectionManager object
 * @param pt tau pt
 * @param wp working point of the ID cut
 * @param sf_vsjet_tau20to25 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau25to30 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau30to35 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau35to40 id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau40toInf id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param correctionset name of the correction set containing the tau id scale
factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
id_vsJet_lt_embedding(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correctionManager,
                      const std::string &pt, const std::string &wp,
                      const std::string &sf_vsjet_tau20to25,
                      const std::string &sf_vsjet_tau25to30,
                      const std::string &sf_vsjet_tau30to35,
                      const std::string &sf_vsjet_tau35to40,
                      const std::string &sf_vsjet_tau40toInf,
                      const std::string &id_output, const std::string &sf_file,
                      const std::string &correctionset) {

    Logger::get("TauIDvsJet_lt_SF_embedding")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("TauIDvsJet_lt_SF_embedding")
        ->debug("ID - Name {}", correctionset);
    auto evaluator = correctionManager.loadCorrection(sf_file, correctionset);
    auto idSF_calculator = [evaluator, wp, sf_vsjet_tau20to25,
                            sf_vsjet_tau25to30, sf_vsjet_tau30to35,
                            sf_vsjet_tau35to40, sf_vsjet_tau40toInf,
                            correctionset](const float &pt) {
        double sf = 1.;
        Logger::get("TauIDvsJet_lt_SF_embedding")
            ->debug("ID {} - pt {}, wp {} "
                    "sf_vsjet_tau20to25 {}, sf_vsjet_tau25to30 {}, "
                    "sf_vsjet_tau30to35{}, sf_vsjet_tau35to40 {}, "
                    "sf_vsjet_tau40toInf {},",
                    correctionset, pt, wp, sf_vsjet_tau20to25,
                    sf_vsjet_tau25to30, sf_vsjet_tau30to35, sf_vsjet_tau35to40,
                    sf_vsjet_tau40toInf);
        if (pt >= 20.0 && pt < 25.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau20to25, wp});
        } else if (pt >= 25.0 && pt < 30.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau25to30, wp});
        } else if (pt >= 30.0 && pt < 35.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau30to35, wp});
        } else if (pt >= 35.0 && pt < 40.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau35to40, wp});
        } else if (pt >= 40.0 && pt < 10000.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau40toInf, wp});
        } else {
            sf = 1.;
        }
        Logger::get("TauIDvsJet_lt_SF_embedding")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {pt});
    return df1;
}

/**
 * @brief Function used to evaluate vsJets tau id scale factors in the tt
channel with the correctionlib for tauembedded samples

 * @param df The input dataframe
 * @param correctionManager The CorrectionManager object
 * @param decaymode decay mode of the tau
 * @param wp working point of the ID cut
 * @param sf_vsjet_tauDM0 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM1 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM10 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM11 id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param correctionset name of the correction set containing the tau id scale
factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_vsJet_tt_embedding(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correctionManager,
    const std::string &pt,
    const std::string &decaymode, const std::string &wp,
    const std::string &sf_vsjet_tauDM0, const std::string &sf_vsjet_tauDM1,
    const std::string &sf_vsjet_tauDM10, const std::string &sf_vsjet_tauDM11,
    const std::string &id_output, const std::string &sf_file,
    const std::string &correctionset) {

    Logger::get("TauIDvsJet_tt_SF_embedding")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("TauIDvsJet_tt_SF_embedding")
        ->debug("ID - Name {}", correctionset);
    auto evaluator = correctionManager.loadCorrection(sf_file, correctionset);
    auto idSF_calculator = [evaluator, wp, sf_vsjet_tauDM0, sf_vsjet_tauDM1,
                            sf_vsjet_tauDM10, sf_vsjet_tauDM11,
                            correctionset](const float &pt, const int &decaymode) {
        double sf = 1.;
        Logger::get("TauIDvsJet_tt_SF_embedding")
            ->debug("ID {} - decaymode {}, wp {}, pt {} "
                    "sf_vsjet_tauDM0 {}, sf_vsjet_tauDM1 {}, "
                    "sf_vsjet_tauDM10{}, sf_vsjet_tauDM11 {}, ",
                    correctionset, decaymode, wp, pt, sf_vsjet_tauDM0,
                    sf_vsjet_tauDM1, sf_vsjet_tauDM10, sf_vsjet_tauDM11);
        if (decaymode == 0) {
            sf = evaluator->evaluate({pt, decaymode, sf_vsjet_tauDM0, wp});
        } else if (decaymode == 1) {
            sf = evaluator->evaluate({pt, decaymode, sf_vsjet_tauDM1, wp});
        } else if (decaymode == 10) {
            sf = evaluator->evaluate({pt, decaymode, sf_vsjet_tauDM10, wp});
        } else if (decaymode == 11) {
            sf = evaluator->evaluate({pt, decaymode, sf_vsjet_tauDM11, wp});
        } else {
            sf = 1.;
        }
        Logger::get("TauIDvsJet_tt_SF_embedding")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {pt, decaymode});
    return df1;
}
} // namespace tau

namespace embedding {
/**
 * @brief Function used to readout the embedding selection trigger scalefactors
 *
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param pt_1 the pt of the leading generator particle in the event. This
 * corresponds to the leading muon selected by the embedding selection
 * @param eta_1 the eta of the leading generator particle in the event. This
 * corresponds to the leading muon selected by the embedding selection
 * @param pt_2 the pt of the subleading generator particle in the event. This
 * corresponds to the subleading muon selected by the embedding selection
 * @param eta_2 the eta of the subleading generator particle in the event. This
 * corresponds to the subleading muon selected by the embedding selection
 * @param output name of the output column
 * @param sf_file path to the correctionlib file containing the scale factor
 * @param idAlgorithm name of the scale factor in the correctionlib file
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode
selection_trigger(ROOT::RDF::RNode df,
                  correctionManager::CorrectionManager &correctionManager,
                  const std::string &pt_1, const std::string &eta_1,
                  const std::string &pt_2, const std::string &eta_2,
                  const std::string &output, const std::string &sf_file,
                  const std::string &idAlgorithm) {

    Logger::get("EmbeddingSelectionTriggerSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator = correctionManager.loadCorrection(sf_file, idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator](const float &pt_1, const float &eta_1, const float &pt_2,
                    const float &eta_2) {
            Logger::get("EmbeddingSelectionTriggerSF")
                ->debug(" pt_1 {}, eta_1 {}, pt_2 {}, eta_2 {}", pt_1, eta_1,
                        pt_2, eta_2);
            double sf = 1.;
            sf = evaluator->evaluate(
                {pt_1, std::abs(eta_1), pt_2, std::abs(eta_2)});
            Logger::get("EmbeddingSelectionTriggerSF")->debug("sf {}", sf);
            return sf;
        },
        {pt_1, eta_1, pt_2, eta_2});
    return df1;
}
/**
 * @brief Function used to readout the embedding selection trigger scalefactors
 * WARNING: This function is deprecated, please use the new function with
 * CorrectionManager
 * @param df the input dataframe
 * @param pt_1 the pt of the leading generator particle in the event. This
 * corresponds to the leading muon selected by the embedding selection
 * @param eta_1 the eta of the leading generator particle in the event. This
 * corresponds to the leading muon selected by the embedding selection
 * @param pt_2 the pt of the subleading generator particle in the event. This
 * corresponds to the subleading muon selected by the embedding selection
 * @param eta_2 the eta of the subleading generator particle in the event. This
 * corresponds to the subleading muon selected by the embedding selection
 * @param output name of the output column
 * @param sf_file path to the correctionlib file containing the scale factor
 * @param idAlgorithm name of the scale factor in the correctionlib file
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode
selection_trigger(ROOT::RDF::RNode df, const std::string &pt_1,
                  const std::string &eta_1, const std::string &pt_2,
                  const std::string &eta_2, const std::string &output,
                  const std::string &sf_file, const std::string &idAlgorithm) {
    Logger::get("EmbeddingSelectionTriggerSF")
        ->warn("Not using CorrectionManager is deprecated, this function will "
               "be removed in a future release, switch to the new function "
               "using CorrectionManager");
    Logger::get("EmbeddingSelectionTriggerSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator](const float &pt_1, const float &eta_1, const float &pt_2,
                    const float &eta_2) {
            Logger::get("EmbeddingSelectionTriggerSF")
                ->debug(" pt_1 {}, eta_1 {}, pt_2 {}, eta_2 {}", pt_1, eta_1,
                        pt_2, eta_2);
            double sf = 1.;
            sf = evaluator->evaluate(
                {pt_1, std::abs(eta_1), pt_2, std::abs(eta_2)});
            Logger::get("EmbeddingSelectionTriggerSF")->debug("sf {}", sf);
            return sf;
        },
        {pt_1, eta_1, pt_2, eta_2});
    return df1;
}
/**
 * @brief Function used to readout the embedding selection trigger scalefactors.
 *
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param pt the pt of the generator particle in the event. This corresponds to
 * one of the muons selected by the embedding selection
 * @param eta the eta of the generator particle in the event. This corresponds
 * to one of the muons selected by the embedding selection
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode
selection_id(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correctionManager,
             const std::string &pt, const std::string &eta,
             const std::string &output, const std::string &sf_file,
             const std::string &idAlgorithm) {
    Logger::get("EmbeddingSelectionIDSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator = correctionManager.loadCorrection(sf_file, idAlgorithm);
    auto df1 =
        df.Define(output,
                  [evaluator](const float &pt, const float &eta) {
                      Logger::get("EmbeddingSelectionIDSF")
                          ->debug(" pt {}, eta {},", pt, eta);
                      double sf = 1.;
                      sf = evaluator->evaluate({pt, std::abs(eta)});
                      Logger::get("EmbeddingSelectionIDSF")->debug("sf {}", sf);
                      return sf;
                  },
                  {pt, eta});
    return df1;
}
/**
 * @brief Function used to readout the embedding selection trigger scalefactors.
 * WARNING: This function is deprecated, please use the new function with
 * CorrectionManager
 * @param df the input dataframe
 * @param pt the pt of the generator particle in the event. This corresponds to
 * one of the muons selected by the embedding selection
 * @param eta the eta of the generator particle in the event. This corresponds
 * to one of the muons selected by the embedding selection
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode selection_id(ROOT::RDF::RNode df, const std::string &pt,
                              const std::string &eta, const std::string &output,
                              const std::string &sf_file,
                              const std::string &idAlgorithm) {
    Logger::get("EmbeddingSelectionIDSF")
        ->warn("Not using CorrectionManager is deprecated, this function will "
               "be removed in a future release, switch to the new function "
               "using CorrectionManager");

    Logger::get("EmbeddingSelectionIDSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 =
        df.Define(output,
                  [evaluator](const float &pt, const float &eta) {
                      Logger::get("EmbeddingSelectionIDSF")
                          ->debug(" pt {}, eta {},", pt, eta);
                      double sf = 1.;
                      sf = evaluator->evaluate({pt, std::abs(eta)});
                      Logger::get("EmbeddingSelectionIDSF")->debug("sf {}", sf);
                      return sf;
                  },
                  {pt, eta});
    return df1;
}
/**
 * @brief Function used to readout SFs from the muon scale factor measurements
 *
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param pt the pt of the muon
 * @param eta the eta of the muon
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param correctiontype the type of the correction. Use `emb` for embedding and
 * `mc` for monte carlo
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @param extrapolation_factor The extrapolation factor to be used for the scale
 * factor, defaults to 1.
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode
muon_sf(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correctionManager,
        const std::string &pt, const std::string &eta,
        const std::string &output, const std::string &sf_file,
        const std::string correctiontype, const std::string &idAlgorithm,
        const float &extrapolation_factor = 1.0) {

    Logger::get("EmbeddingMuonSF")->debug("Correction - Name {}", idAlgorithm);
    auto evaluator = correctionManager.loadCorrection(sf_file, idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator, correctiontype, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("EmbeddingMuonSF")
                ->debug(" pt {}, eta {}, correctiontype {}, extrapolation "
                        "factor {}",
                        pt, eta, correctiontype, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, std::abs(eta), correctiontype});
            Logger::get("EmbeddingMuonSF")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function used to readout SFs from the muon scale factor measurements
 * WARNING: This function is deprecated, please use the new function with
 * CorrectionManager
 * @param df the input dataframe
 * @param pt the pt of the muon
 * @param eta the eta of the muon
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param correctiontype the type of the correction. Use `emb` for embedding and
 * `mc` for monte carlo
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @param extrapolation_factor The extrapolation factor to be used for the scale
 * factor, defaults to 1.
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode muon_sf(ROOT::RDF::RNode df, const std::string &pt,
                         const std::string &eta, const std::string &output,
                         const std::string &sf_file,
                         const std::string correctiontype,
                         const std::string &idAlgorithm,
                         const float &extrapolation_factor = 1.0) {
    Logger::get("EmbeddingMuonSF")
        ->warn("Not using CorrectionManager is deprecated, this function will "
               "be removed in a future release, switch to the new function "
               "using CorrectionManager");

    Logger::get("EmbeddingMuonSF")->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator, correctiontype, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("EmbeddingMuonSF")
                ->debug(" pt {}, eta {}, correctiontype {}, extrapolation "
                        "factor {}",
                        pt, eta, correctiontype, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, std::abs(eta), correctiontype});
            Logger::get("EmbeddingMuonSF")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function used to readout SFs from the electron scale factor
 * measurements
 *
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param pt the pt of the electron
 * @param eta the eta of the electron
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param correctiontype the type of the correction. Use `emb` for embedding and
 * `mc` for monte carlo
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @param extrapolation_factor The extrapolation factor to be used for the scale
 * factor, defaults to 1.
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode
electron_sf(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correctionManager,
            const std::string &pt, const std::string &eta,
            const std::string &output, const std::string &sf_file,
            const std::string correctiontype, const std::string &idAlgorithm,
            const float &extrapolation_factor = 1.0) {

    Logger::get("EmbeddingElectronSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator = correctionManager.loadCorrection(sf_file, idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator, correctiontype, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("EmbeddingElectronSF")
                ->debug(" pt {}, eta {}, correctiontype {}, extrapolation "
                        "factor {}",
                        pt, eta, correctiontype, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, eta, correctiontype});
            Logger::get("EmbeddingElectronSF")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function used to readout SFs from the electron scale factor
 * measurements
 * WARNING: This function is deprecated, please use the new function with
 * CorrectionManager
 * @param df the input dataframe
 * @param pt the pt of the electron
 * @param eta the eta of the electron
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param correctiontype the type of the correction. Use `emb` for embedding and
 * `mc` for monte carlo
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @param extrapolation_factor The extrapolation factor to be used for the scale
 * factor, defaults to 1.
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode electron_sf(ROOT::RDF::RNode df, const std::string &pt,
                             const std::string &eta, const std::string &output,
                             const std::string &sf_file,
                             const std::string correctiontype,
                             const std::string &idAlgorithm,
                             const float &extrapolation_factor = 1.0) {
    Logger::get("EmbeddingElectronSF")
        ->warn("Not using CorrectionManager is deprecated, this function will "
               "be removed in a future release, switch to the new function "
               "using CorrectionManager");
    Logger::get("EmbeddingElectronSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator, correctiontype, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("EmbeddingElectronSF")
                ->debug(" pt {}, eta {}, correctiontype {}, extrapolation "
                        "factor {}",
                        pt, eta, correctiontype, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, eta, correctiontype});
            Logger::get("EmbeddingElectronSF")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function to evaluate the di-tau trigger or etau/mutau cross trigger
 * scale factor for embedded events from a xpog file
 *
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param pt the name of the column containing the tau pt variable
 * @param decaymode the name of the column containing the tau decay mode
 * variable
 * @param output name of the scale factor column
 * @param wp the name of the the tau id working point VVVLoose-VVTight
 * @param sf_file path to the file with the tau trigger scale factors
 * @param type the type of the tau trigger, available are "ditau", "etau",
 * "mutau", "ditauvbf"
 * @param corrtype name of the tau trigger correction type, available are
 * "eff_data", "eff_mc", "sf"
 * @param syst name of the systematic variation, options are "nom", "up", "down"
 * @return ROOT::RDF::RNode a new dataframe containing the new sf column
 */

ROOT::RDF::RNode
ditau_trigger_sf(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correctionManager,
                 const std::string &pt, const std::string &decaymode,
                 const std::string &output, const std::string &wp,
                 const std::string &sf_file, const std::string &type,
                 const std::string &corrtype, const std::string &syst) {

    Logger::get("ditau_trigger")
        ->debug("Setting up function for di-tau trigger sf");
    Logger::get("ditau_trigger")
        ->debug("correction type {}, file {}", corrtype, sf_file);
    // tauTriggerSF is the only correction set in the file for now, might change
    // with official sf release -> change into additional input parameter
    auto evaluator = correctionManager.loadCorrection(sf_file, "tauTriggerSF");
    Logger::get("ditau_trigger")->debug("WP {} - trigger type {}", wp, type);
    auto trigger_sf_calculator = [evaluator, wp, type, corrtype,
                                  syst](const float &pt, const int &decaymode) {
        float sf = 1.;
        Logger::get("ditau_trigger")
            ->debug("decaymode {}, pt {}", decaymode, pt);
        if (pt > 0) {
            if (decaymode == 0 || decaymode == 1 || decaymode == 10 ||
                decaymode == 11) {
                sf = evaluator->evaluate(
                    {pt, decaymode, type, wp, corrtype, syst});
            } else {
                sf = evaluator->evaluate({pt, -1, type, wp, corrtype, syst});
            }
        }
        Logger::get("ditau_trigger")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(output, trigger_sf_calculator, {pt, decaymode});
    return df1;
}
/**
 * @brief Function to evaluate the di-tau trigger or etau/mutau cross trigger
 * scale factor for embedded events from a xpog file
 * WARNING: This function is deprecated, please use the new function with
 * CorrectionManager
 * @param df the input dataframe
 * @param pt the name of the column containing the tau pt variable
 * @param decaymode the name of the column containing the tau decay mode
 * variable
 * @param output name of the scale factor column
 * @param wp the name of the the tau id working point VVVLoose-VVTight
 * @param sf_file path to the file with the tau trigger scale factors
 * @param type the type of the tau trigger, available are "ditau", "etau",
 * "mutau", "ditauvbf"
 * @param corrtype name of the tau trigger correction type, available are
 * "eff_data", "eff_mc", "sf"
 * @param syst name of the systematic variation, options are "nom", "up", "down"
 * @return ROOT::RDF::RNode a new dataframe containing the new sf column
 */

ROOT::RDF::RNode
ditau_trigger_sf(ROOT::RDF::RNode df, const std::string &pt,
                 const std::string &decaymode, const std::string &output,
                 const std::string &wp, const std::string &sf_file,
                 const std::string &type, const std::string &corrtype,
                 const std::string &syst) {
    Logger::get("ditau_trigger")
        ->warn("Not using CorrectionManager is deprecated, this function will "
               "be removed in a future release, switch to the new function "
               "using CorrectionManager");
    Logger::get("ditau_trigger")
        ->debug("Setting up function for di-tau trigger sf");
    Logger::get("ditau_trigger")
        ->debug("correction type {}, file {}", corrtype, sf_file);
    // tauTriggerSF is the only correction set in the file for now, might change
    // with official sf release -> change into additional input parameter
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at("tauTriggerSF");
    Logger::get("ditau_trigger")->debug("WP {} - trigger type {}", wp, type);
    auto trigger_sf_calculator = [evaluator, wp, type, corrtype,
                                  syst](const float &pt, const int &decaymode) {
        float sf = 1.;
        Logger::get("ditau_trigger")
            ->debug("decaymode {}, pt {}", decaymode, pt);
        if (pt > 0) {
            if (decaymode == 0 || decaymode == 1 || decaymode == 10 ||
                decaymode == 11) {
                sf = evaluator->evaluate(
                    {pt, decaymode, type, wp, corrtype, syst});
            } else {
                sf = evaluator->evaluate({pt, -1, type, wp, corrtype, syst});
            }
        }
        Logger::get("ditau_trigger")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(output, trigger_sf_calculator, {pt, decaymode});
    return df1;
}
} // namespace embedding
} // namespace scalefactor
#endif /* GUARD_SCALEFACTORS_H */