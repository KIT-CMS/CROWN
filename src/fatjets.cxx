#ifndef GUARDFATJETS_H
#define GUARDFATJETS_H

#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include <typeinfo>

namespace fatjet {
/// Function to find a fatjet which matches to the leading b-jet from a b-jet
/// pair. The match is done with a deltaR criterium.
///
/// \param[in] df the input dataframe
/// \param[out] output_name the name of the selected fatjet (index)
/// \param[in] good_fatjet_collection name of the collection with the indices of
/// good fatjets \param[in] fatjet_pt name of the fatjet pts \param[in]
/// fatjet_eta name of the fatjet etas \param[in] fatjet_phi name of the fatjet
/// phis \param[in] fatjet_mass name of the fatjet masses \param[in] bpair_p4_1
/// four vector of the leading b-jet from the selected b-jet pair \param[in]
/// deltaRmax maximum required distance in dR between b-jet and a fatjet
/// candidate
///
/// \return a dataframe containing the new mask
// auto FindFatjetMatchingBjet(ROOT::RDF::RNode df, const std::string &output_name,
//                             const std::string &good_fatjet_collection,
//                             const std::string &fatjet_pt,
//                             const std::string &fatjet_eta,
//                             const std::string &fatjet_phi,
//                             const std::string &fatjet_mass,
//                             const std::string &bpair_p4_1,
//                             const float &deltaRmax) {
//     Logger::get("fatjet::FatjetMatchingToBjet")->debug("Setting up algorithm");
//     auto df1 = df.Define(
//         output_name,
//         [deltaRmax](const ROOT::RVec<int> &good_fatjet_collection,
//                     const ROOT::RVec<float> &fatjet_pt,
//                     const ROOT::RVec<float> &fatjet_eta,
//                     const ROOT::RVec<float> &fatjet_phi,
//                     const ROOT::RVec<float> &fatjet_mass,
//                     const ROOT::Math::PtEtaPhiMVector &bpair_p4_1) {
//             ROOT::RVec<int> selected_fatjet = {-1};
//             if ((good_fatjet_collection.size() > 0) && (bpair_p4_1.pt() > 0)) {
//                 Logger::get("fatjet::FatjetMatchingToBjet")
//                     ->debug("Running algorithm on at least one good fatjet");
//                 for (auto &index : good_fatjet_collection) {
//                     ROOT::Math::PtEtaPhiMVector fatjet_candidate =
//                         ROOT::Math::PtEtaPhiMVector(
//                             fatjet_pt.at(index), fatjet_eta.at(index),
//                             fatjet_phi.at(index), fatjet_mass.at(index));
//                     Logger::get("fatjet::FatjetMatchingToBjet")
//                         ->debug("{} fatjet candidate vector: {}", index,
//                                 fatjet_candidate);
//                     if ((ROOT::Math::VectorUtil::DeltaR(
//                              bpair_p4_1, fatjet_candidate) < deltaRmax)) {
//                         selected_fatjet = {static_cast<int>(index)};
//                         Logger::get("fatjet::FatjetMatchingToBjet")
//                             ->debug("Final fatjet {}", selected_fatjet[0]);
//                         break;
//                     }
//                 }
//             }
//             return selected_fatjet;
//         },
//         {good_fatjet_collection, fatjet_pt, fatjet_eta, fatjet_phi, fatjet_mass,
//          bpair_p4_1});
//     return df1;
// }
/// Function to find a fatjet with the highest particleNet X(tm) vs QCD score 
///
/// \param[in] df the input dataframe
/// \param[out] output_name the name of the selected fatjet (index)
/// \param[in] good_fatjet_collection name of the collection with the indices of
/// good fatjets \param[in] fatjet_pNet_Xtm name of the variable with the Xtm particleNet scores 
/// \param[in] fatjet_pNet_QCD name of the variable with the QCD particleNet scores 
///
/// \return a dataframe containing the new mask
auto FindXtmFatjet(ROOT::RDF::RNode df, const std::string &output_name,
                            const std::string &good_fatjet_collection,
                            const std::string &fatjet_pNet_XtmVsQCD) {
    Logger::get("fatjet::FindXtmFatjet")->debug("Setting up algorithm");
    auto df1 = df.Define(
        output_name,
        [](const ROOT::RVec<int> &good_fatjet_collection,
           const ROOT::RVec<float> &XmtVsQCD_tagger) {
            ROOT::RVec<int> selected_fatjet = {-1};
            float highest_pNet_value = default_float;
            if ((good_fatjet_collection.size() > 0)) {
                Logger::get("fatjet::FindXtmFatjet")
                    ->debug("Running algorithm on at least one good fatjet");
                // float Xtm = default_float;
                // float QCD = default_float;
                float Xtm_vs_QCD = default_float;
                for (auto &index : good_fatjet_collection) {
                    // Xtm = Xtm_tagger.at(index);
                    // QCD = QCD_tagger.at(index);
                    Xtm_vs_QCD = XmtVsQCD_tagger.at(index);
                    if (Xtm_vs_QCD > highest_pNet_value) {
                        highest_pNet_value = Xtm_vs_QCD;
                        selected_fatjet = {static_cast<int>(index)};
                    }
                }
                Logger::get("fatjet::FindXtmFatjet")
                            ->debug("Final fatjet {}", selected_fatjet[0]);
            }
            return selected_fatjet;
        },
        {good_fatjet_collection, fatjet_pNet_XtmVsQCD});
    return df1;
}
// This function flags events, where a suitable particle pair is found.
/// A pair is considered suitable, if a PairSelectionAlgo (like
/// ditau_pairselection::mutau::PairSelectionAlgo) returns indices, that are
/// not -1. Events, where any of the particle indices is -1 are vetoed
/// by this filter.
///
/// \param df The input dataframe
/// \param flagname The name of the generated flag column
/// \param pairname The name of the column, containing the indices of the
/// particles in the particle quantity vectors.
/// \returns a dataframe with the
/// new flag
ROOT::RDF::RNode flagGoodFatjets(ROOT::RDF::RNode df, const std::string &flagname,
                               const std::string &fatjetname) {
    using namespace ROOT::VecOps;
    return df.Define(
        flagname,
        // [](const ROOT::RVec<int> &fatjets) { return fatjets.size(); },
        [](const ROOT::RVec<int> &fatjets) {bool flagging  = (fatjets.size()<1) ? 0 :1; return flagging; },

        {fatjetname});
}
} // end namespace fatjet

namespace quantities {
namespace fatjet {
/// Function to writeout the value of the softdrop mass for a fatjet.
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] m_softdrop name of the column that contains the softdrop mass of
/// the fatjets
/// \param[in] fatjetcollection name of the vector that contains fatjet indices
/// of the fatjets belonging to the collection, its length constitutes the
/// output quantity \param position The position in the fatjet collection
/// vector, which is used to store the index of the particle in the particle
/// quantity vectors.
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode msoftdrop(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &m_softdrop,
                           const std::string &fatjetcollection,
                           const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<float> &softdrop_masses,
                                const ROOT::RVec<int> &fatjetcollection) {
                         float mass = default_float;
                         try {
                             const int index = fatjetcollection.at(position);
                             mass = softdrop_masses.at(index);
                         } catch (const std::out_of_range &e) {
                         }
                         return mass;
                     },
                     {m_softdrop, fatjetcollection});
}
/// Function to writeout the value of the particleNet Xtm vs QCD tagger for a
/// fatjet.
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] pNet_Xtm name of the column that contains the particleNet raw
/// value for X->tm of the fatjets \param[in] pNet_QCD name of the column that
/// contains the particleNet raw value for QCD of the fatjets \param[in]
/// fatjetcollection name of the vector that contains fatjet indices of the
/// fatjets belonging to the collection, its length constitutes the output
/// quantity \param position The position in the fatjet collection vector, which
/// is used to store the index of the particle in the particle quantity vectors.
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode
particleNet_XtmVsQCD(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pNet_XtmVsQCD, 
                     const std::string &fatjetcollection, const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<float> &XtmVsQCD_tagger,
                                const ROOT::RVec<int> &fatjetcollection) {
                         float Xtm_vs_QCD = default_float;
                         try {
                             const int index = fatjetcollection.at(position);
                             Xtm_vs_QCD = XtmVsQCD_tagger.at(index);
                         } catch (const std::out_of_range &e) {
                         }
                         return Xtm_vs_QCD;
                     },
                     {pNet_XtmVsQCD, fatjetcollection});
}
/// Function to writeout the N- over N-1-prong n-subjettiness ratio for a
/// fatjet. This variable is should discriminate between fatjets with N subjets
/// and fatjets with N-1 subjets. Lower values mean the fatjet is more N-prong
/// like, higher values more N-1-prong like.
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] tauN name of the column that contains the n-subjettiness for
/// N-prong of the fatjets \param[in] tauNm1 name of the column that contains
/// the n-subjettiness for N-1-prong of the fatjets \param[in] fatjetcollection
/// name of the vector that contains fatjet indices of the fatjets belonging to
/// the collection, its length constitutes the output quantity \param position
/// The position in the fatjet collection vector, which is used to store the
/// index of the particle in the particle quantity vectors.
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode
nsubjettiness_ratio(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tauN, const std::string &tauNm1,
                    const std::string &fatjetcollection, const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<float> &nsubjettiness_N,
                                const ROOT::RVec<float> &nsubjettiness_Nm1,
                                const ROOT::RVec<int> &fatjetcollection) {
                         float nsubjet_N = default_float;
                         float nsubjet_Nm1 = default_float;
                         float ratio = default_float;
                         try {
                             const int index = fatjetcollection.at(position);
                             nsubjet_N = nsubjettiness_N.at(index);
                             nsubjet_Nm1 = nsubjettiness_Nm1.at(index);
                             ratio = nsubjet_N / nsubjet_Nm1;
                         } catch (const std::out_of_range &e) {
                         }
                         return ratio;
                     },
                     {tauN, tauNm1, fatjetcollection});
}

ROOT::RDF::RNode
subjet_1(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &subjetIDx1,
                    const std::string &fatjetcollection, const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<short> &subjetIDx1,
                                const ROOT::RVec<int> &fatjetcollection) {
                        //  short subjet1 = default_short;
                        ROOT::RVec<short> subjet1;
                         try {
                             const int index = fatjetcollection.at(position);
                             subjet1 = {subjetIDx1.at(index)};
                         } catch (const std::out_of_range &e) {
                         }
                         return subjet1;
                     },
                     {subjetIDx1, fatjetcollection});
}

// ROOT::RDF::RNode trigger_mu_in_fatjet(ROOT::RDF::RNode df,
//                                  const std::string &outputname,
//                                  const std::string &fatjet_p4,
//                                  const std::string &muon_pt,
//                                  const std::string &muon_eta,
//                                  const std::string &muon_phi,
//                                  const std::string &muon_mass) {

//     auto match_muon_with_fatjet = [](const ROOT::RVec<float> &fatjet_p4,
//                                     const ROOT::RVec<float> &muon_pts,
//                                     const ROOT::RVec<float> &muon_etas,
//                                     const ROOT::RVec<float> &muon_phis,
//                                     const ROOT::RVec<float> &muon_masses) {

//          for (unsigned int i = 0; i < muon_pts.size(); i++) {

//             ROOT::Math::PtEtaPhiMVector muon_p4(muon_pts.at(i), muon_etas.at(i),
//                     muon_phis.at(i), muon_masses.at(i));

//             if (ROOT::Math::VectorUtil::DeltaR(fatjet_p4, muon_p4) < 0.8) {
//                 Logger::get("fatjet::trigger_mu_in_fatjet")
//                     ->debug("Muon found in fatjet with pt {} and eta {}",
//                             muon_pts.at(i), muon_etas.at(i));
//                 return 1;
//             }
//             else {
//                 Logger::get("fatjet::trigger_mu_in_fatjet")
//                 ->debug("Muon NOT found in fatjet with pt {} and eta {}",
//                     muon_pts.at(i), muon_etas.at(i));
//                 return 0;
//             }
//          }                                                    
                                 
// };

//     auto df1 =
//         df.Define(outputname, match_muon_with_fatjet,
//                   {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass});
//     return df1;
//     }
} // end namespace fatjet
} // end namespace quantities
#endif /* GUARDFATJETS_H */