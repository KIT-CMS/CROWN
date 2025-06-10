#ifndef GUARD_MUONS_H
#define GUARD_MUONS_H

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/RoccoR.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "RooFunctor.h"
#include "RooWorkspace.h"
#include "correction.h"

namespace physicsobject {
namespace muon {

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
 * @brief This function calculates muon ID scale factors (SFs) for a single
 * muon dependening on its pseudorapidity (\f$\eta\f$) and transverse momentum
 * (\f$p_T\f$). The scale factors are loaded from a correctionlib file using a
 * specified scale factor name and variation.
 *
 * Recommendations by MuonPOG:
 * - [Run2](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_6_1)
 * - [Run3](https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_6_2)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * muon scale factor file
 * @param outputname name of the output column containing the ID scale factor
 * @param pt name of the column containing the transverse momentum of a muon
 * @param eta name of the column containing the pseudorapidity of a muon
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon scale factor for the ID correction,
 * e.g. "NUM_MediumID_DEN_TrackerMuons"
 * @param variation name the scale factor variation, "nominal" for the nominal
 * scale factor and "systup"/"systdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Id(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correction_manager,
                    const std::string &outputname, const std::string &pt,
                    const std::string &eta, const std::string &sf_file,
                    const std::string &sf_name, const std::string &variation) {
    Logger::get("physicsobject::muon::scalefactor::Id")
        ->debug("Setting up functions for muon id sf");
    Logger::get("physicsobject::muon::scalefactor::Id")
        ->debug("ID - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, variation](const float &pt, const float &eta) {
            Logger::get("physicsobject::muon::scalefactor::Id")
                ->debug("ID - pt {}, eta {}", pt, eta);
            double sf = 1.;
            // check to prevent muons with default values due to tau energy
            // correction shifts below good tau pt selection
            if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                sf = evaluator->evaluate({std::abs(eta), pt, variation});
            }
            return sf;
        },
        {pt, eta});
    return df1;
}

/**
 * @brief This function calculates muon iso scale factors (SFs) for a single
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
 * @param outputname name of the output column containing the iso scale factor
 * @param pt name of the column containing the transverse momentum of a muon
 * @param eta name of the column containing the pseudorapidity of a muon
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon scale factor for the iso correction,
 * e.g. "NUM_TightRelIso_DEN_MediumID"
 * @param variation name the scale factor variation, "nominal" for the nominal
 * scale factor and "systup"/"systdown" for the up/down variation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode Iso(ROOT::RDF::RNode df,
                     correctionManager::CorrectionManager &correction_manager,
                     const std::string &outputname, const std::string &pt,
                     const std::string &eta, const std::string &sf_file,
                     const std::string &sf_name, const std::string &variation) {
    Logger::get("physicsobject::muon::scalefactor::Iso")
        ->debug("Setting up functions for muon iso sf");
    Logger::get("physicsobject::muon::scalefactor::Iso")
        ->debug("ISO - Name {}", sf_name);
    auto evaluator = correction_manager.loadCorrection(sf_file, sf_name);
    auto df1 = df.Define(
        outputname,
        [evaluator, variation](const float &pt, const float &eta) {
            Logger::get("physicsobject::muon::scalefactor::Iso")
                ->debug("ISO - pt {}, eta {}", pt, eta);
            double sf = 1.;
            // check to prevent muons with default values due to tau energy
            // correction shifts below good tau pt selection
            if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                sf = evaluator->evaluate({std::abs(eta), pt, variation});
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
 */
ROOT::RDF::RNode
Trigger(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname, const std::string &pt,
        const std::string &eta, const std::string &sf_file,
        const std::string &sf_name, const std::string &variation) {
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
 * @brief Function used to evaluate id scale factors for muons from
 * a RooWorkspace.
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param id_output name of the id scale factor column
 * @param workspace_name path to the Rooworkspace
 * @param id_functor_name name of the function from the workspace
 * @param id_arguments arguments of the function
 *
 * @return a new dataframe containing the new column
 *
 * @warning This function is deprecated and should not be used!
 */
ROOT::RDF::RNode Id_rooworkspace(ROOT::RDF::RNode df, const std::string &pt,
                                 const std::string &eta,
                                 const std::string &id_output,
                                 const std::string &workspace_name,
                                 const std::string &id_functor_name,
                                 const std::string &id_arguments) {

    Logger::get("physicsobject::muon::scalefactor::Id_rooworkspace")
        ->debug("Setting up functions for muon sf");
    Logger::get("physicsobject::muon::scalefactor::Id_rooworkspace")
        ->debug("ID - Function {} // argset {}", id_functor_name, id_arguments);

    const std::shared_ptr<RooFunctorThreadsafe> id_function =
        loadFunctor(workspace_name, id_functor_name, id_arguments);
    auto df1 =
        utility::EvaluateWorkspaceFunction(df, id_output, id_function, pt, eta);
    return df1;
}

/**
 * @brief Function used to evaluate iso scale factors for muons from
 * a RooWorkspace.
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param iso muon iso
 * @param iso_output name of the iso scale factor column
 * @param workspace_name path to the Rooworkspace
 * @param iso_functor_name name of the function from the workspace
 * @param iso_arguments arguments of the function
 *
 * @return a new dataframe containing the new column
 *
 * @warning This function is deprecated and should not be used!
 */
ROOT::RDF::RNode Iso_rooworkspace(ROOT::RDF::RNode df, const std::string &pt,
                                  const std::string &eta,
                                  const std::string &iso,
                                  const std::string &iso_output,
                                  const std::string &workspace_name,
                                  const std::string &iso_functor_name,
                                  const std::string &iso_arguments) {

    Logger::get("physicsobject::muon::scalefactor::Iso_rooworkspace")
        ->debug("Setting up functions for muon sf");
    Logger::get("physicsobject::muon::scalefactor::Iso_rooworkspace")
        ->debug("Iso - Function {} // argset {}", iso_functor_name,
                iso_arguments);

    const std::shared_ptr<RooFunctorThreadsafe> iso_function =
        loadFunctor(workspace_name, iso_functor_name, iso_arguments);
    auto df1 = utility::EvaluateWorkspaceFunction(df, iso_output, iso_function,
                                                  pt, eta, iso);
    return df1;
}
} // end namespace scalefactor
} // end namespace muon
} // end namespace physicsobject
#endif /* GUARD_MUONS_H */