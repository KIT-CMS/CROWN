#ifndef GUARD_FATJETS_H
#define GUARD_FATJETS_H

#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <algorithm>

namespace physicsobject {
namespace fatjet {
namespace quantity {

/**
 * @brief Applies jet identification criteria based on JSON-defined jet ID corrections for AK8 jets.
 *
 * This function loads jet ID definitions from correctionlib JSON files for the specified
 * jet collection and evaluates two sets of criteria:
 *  - Tight ID
 *  - Tight Lepton Veto ID
 *
 * It uses these evaluations to assign a jet ID code to each jet in the input dataframe:
 *  - 6 : passes Tight and Tight Lepton Veto IDs
 *  - 2 : passes Tight ID but fails Tight Lepton Veto ID
 *  - 0 : fails Tight ID
 * (Ref. https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV#Recommendations_for_the_13_6_AN1)
 *
 * The jet ID is returned as a vector of int, compatible with NanoAOD v9 conventions.
 *
 * @param df Input ROOT RDataFrame containing jet variables
 * @param correction_manager correction manager responsible for loading the
 * correction scale uncertainty patch file
 * @param outputname Name of the new column to hold the computed jet ID flags
 * @param jet_eta Name of the branch for jet pseudorapidity
 * @param jet_chHEF Name of the branch for charged hadron energy fraction
 * @param jet_neHEF Name of the branch for neutral hadron energy fraction
 * @param jet_chEmEF Name of the branch for charged electromagnetic energy fraction
 * @param jet_neEmEF Name of the branch for neutral electromagnetic energy fraction
 * @param jet_chMult Name of the branch for charged multiplicity
 * @param jet_neMult Name of the branch for neutral multiplicity
 * @param jet_id_file Path to the jet ID JSON file containing correction definitions
 * @param jet_name Prefix of the jet collection used to select the appropriate corrections
 *
 * @return a RDataFrame with the new jet ID column appended
 */

ROOT::RDF::RNode 
ID(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &jet_eta,
        const std::string &jet_chHEF,
        const std::string &jet_neHEF,
        const std::string &jet_chEmEF,
        const std::string &jet_neEmEF,
        const std::string &jet_chMult,
        const std::string &jet_neMult,
        const std::string &jet_id_file,
        const std::string &jet_name) {

    // Load the jet ID correction file with Tight criteria
    auto tightID = 
        correction_manager.loadCorrection(jet_id_file, jet_name + "_Tight");

    // Load the jet ID correction file with TightLeptonVeto criteria
    auto tightLepVetoID = 
        correction_manager.loadCorrection(jet_id_file, jet_name + "_TightLeptonVeto"); 

    auto compute_jet_id = [tightID, tightLepVetoID](const ROOT::RVec<float> &eta,
                                                    const ROOT::RVec<float> &chHEF,
                                                    const ROOT::RVec<float> &neHEF,
                                                    const ROOT::RVec<float> &chEmEF,
                                                    const ROOT::RVec<float> &neEmEF,
                                                    const ROOT::RVec<UChar_t> &chMult,
                                                    const ROOT::RVec<UChar_t> &neMult) {

        size_t nJets = eta.size();
        ROOT::RVec<int> jetId(nJets); 
        for (size_t i = 0; i < nJets; ++i) {
            UChar_t mult = chMult.at(i) + neMult.at(i);
            bool passTight = false, passTightLepVeto = false;

            passTight = (tightID->evaluate(
                {eta.at(i), chHEF.at(i), neHEF.at(i), chEmEF.at(i),
                neEmEF.at(i), chMult.at(i), neMult.at(i), mult}
            ) > 0.5);

            passTightLepVeto = (tightLepVetoID->evaluate(
                {eta.at(i), chHEF.at(i), neHEF.at(i), chEmEF.at(i),
                neEmEF.at(i), chMult.at(i), neMult.at(i), mult}
            ) > 0.5);

            if (passTight && passTightLepVeto) jetId[i] = 6;
            else if (passTight && !passTightLepVeto) jetId[i] = 2;
            else jetId[i] = 0;
        }
        return jetId;
    };

    auto df1 = df.Define(outputname, compute_jet_id,
                     {jet_eta, jet_chHEF, jet_neHEF,
                      jet_chEmEF, jet_neEmEF,
                      jet_chMult, jet_neMult});
    return df1;
}

/**
 * @brief This function calculates a discriminator score from two ParticleNet 
 * tagger outputs (a signal X and the QCD background). The signal can e.g. be 
 * \f$X\rightarrow bb\f$ or \f$X\rightarrow cc\f$ etc. 
 * 
 * The score is computed using the formula: 
 * \f[
 * Score = \frac{P(X)}{P(X) + P(QCD)}
 * \f]
 *
 * @note This function is mainly needed when working with `nanoAODv9`, in newer 
 * versions this ratio is already included as a branch.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the XvsQCD score
 * @param pNet_X_decay name of the column containing the ParticleNet score for 
 * the signal process (\f$X\rightarrow ...\f$)
 * @param pNet_QCD name of the column containing the ParticleNet score for the 
 * QCD background
 * @param fatjet_collection name of the column containing a collection (vector) 
 * of good fatjet indices
 * @param position position of the fatjet in the collection (vector) that should 
 * be used for the score calculation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
ParticleNet_XvsQCD(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pNet_X_decay, const std::string &pNet_QCD,
                     const std::string &fatjet_collection, const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<float> &X_tagger,
                                const ROOT::RVec<float> &QCD_tagger,
                                const ROOT::RVec<int> &fatjets) {
                        float X_decay = default_float;
                        float QCD = default_float;
                        float X_vs_QCD = default_float;
                        if (position >= 0) {
                            const int index = fatjets.at(position);
                            if (index >= 0) {
                                X_decay = X_tagger.at(index);
                                QCD = QCD_tagger.at(index);
                                X_vs_QCD = X_decay / (X_decay + QCD);
                            }
                        }
                        return X_vs_QCD;
                     },
                     {pNet_X_decay, pNet_QCD, fatjet_collection});
}

/**
 * @brief This function calculates the ratio of two N-subjettiness variables 
 * of one fatjet, typically used for discriminating jet substructure:
 * \f[
 * \tau_{ratio} = \frac{\tau_N}{\tau_{N-1}}
 * \f]
 * For example, \f$\tau_{21} = \tau_2 / \tau_1\f$ for 2-prong vs. 1-prong decay 
 * identification. The resulting ratio is always between 0 and 1. If it is close 
 * to 0, the jet is more likely to be a N-prong jet, while a value close to 1 
 * indicates a (N-1)-prong jet.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the N-subjettiness ratio
 * @param tau_N name of the column containing the N-subjettiness variable for 
 * the numerator with number `N`
 * @param tau_Nm1 name of the column containing the N-subjettiness variable for 
 * the denominator with number `N-1`
 * @param fatjet_collection name of the column containing a collection (vector) 
 * of good fatjet indices
 * @param position position of the fatjet in the collection (vector) that should 
 * be used for the ratio calculation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
NsubjettinessRatio(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tau_N, const std::string &tau_Nm1,
                    const std::string &fatjet_collection, const int &position) {
    return df.Define(outputname,
                    [position](const ROOT::RVec<float> &tau_N,
                               const ROOT::RVec<float> &tau_Nm1,
                               const ROOT::RVec<int> &fatjets) {
                        float nsubjet_N = default_float;
                        float nsubjet_Nm1 = default_float;
                        float ratio = default_float;
                        if (position >= 0) { 
                            const int index = fatjets.at(position);
                            if (index >= 0) {
                                nsubjet_N = tau_N.at(index);
                                nsubjet_Nm1 = tau_Nm1.at(index);
                                ratio = nsubjet_N / nsubjet_Nm1;
                            }
                        }
                        return ratio;
                    },
                    {tau_N, tau_Nm1, fatjet_collection});
}
} // end namespace quantity
} // end namespace fatjet
} // end namespace physicsobject
#endif /* GUARD_FATJETS_H */
