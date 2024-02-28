#ifndef GUARD_GENPARTICLES_H
#define GUARD_GENPARTICLES_H

#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "bitset"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
typedef std::bitset<20> IntBits;

enum class GenMatchingCode : int {
    NONE = -1,
    IS_ELE_PROMPT = 1,
    IS_MUON_PROMPT = 2,
    IS_ELE_FROM_TAU = 3,
    IS_MUON_FROM_TAU = 4,
    IS_TAU_HAD_DECAY = 5,
    IS_FAKE = 6,
    IS_ELE_FROM_W = 7,
    IS_MUON_FROM_W = 8
};

namespace genmatching {

namespace tau {
/**
 * @brief implementation of the genmatching used in tau analyses, based
 * on
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * implementation. The genmatches are represented as integer flags:
    Decaytype         | Value
    ------------------|-------
    NONE              | -1,
    IS_ELE_PROMPT     | 1,
    IS_MUON_PROMPT    | 2,
    IS_ELE_FROM_TAU   | 3,
    IS_MUON_FROM_TAU  | 4,
    IS_TAU_HAD_DECAY  | 5,
    IS_FAKE           | 6

    The genmatch is caluclated individually for each taupair candidate.

 *
 * @param df The dataframe to be used for the genmatching.
 * @param outputname The name of the output column.
 * @param hadronicGenTaus The name of the column containing the hadronicGenTaus
 from genmatching::tau::hadronicGenTaus
 * @param genparticles_pdgid The name of the column containing the pdgids of the
 genparticles
 * @param genparticles_statusFlag The name of the column containing the
 statusFlags of the genparticles
 * @param genparticles_pt The name of the column containing the pT of the
 genparticles
 * @param genparticles_eta The name of the column containing the eta of the
 genparticles
 * @param genparticles_phi The name of the column containing the phi of the
 genparticles
 * @param genparticles_mass The name of the column containing the mass of the
 genparticles
 * @param lepton_p4 The name of the column containing the p4 of the lepton to be
 genmatched
 * @return a dataframe with the genmatching as a column named outputname
 */
ROOT::RDF::RNode genmatching(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &hadronicGenTaus,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &lepton_p4) {
    auto match_lepton = [](const std::vector<int> &hadronicGenTaus,
                           const ROOT::RVec<int> &pdgids,
                           const ROOT::RVec<unsigned short> &status_flags,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &lepton_p4) {
        // find closest lepton fulfilling the requirements
        float min_delta_r = 9999;
        int closest_genparticle_index = 0;
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            int pdgid = std::abs(pdgids.at(i));
            // check
            // 1. if there is a gen electron or muon close to the lepton
            // 2. that the genparticle pt is larger than 8 GeV
            // 3. the genparticle is isPrompt (statusbit 0) or
            // isDirectPromptTauDecayProduct (statusbit 5)
            bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
                              IntBits(status_flags.at(i)).test(5));
            if ((pdgid == 11 || pdgid == 13) && pts.at(i) > 8 && statusbit) {
                ROOT::Math::PtEtaPhiMVector gen_p4(pts.at(i), etas.at(i),
                                                   phis.at(i), masses.at(i));
                float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(gen_p4, lepton_p4);
                if (delta_r < min_delta_r) {
                    closest_genparticle_index = i;
                    min_delta_r = delta_r;
                }
            }
        }
        Logger::get("genmatching::tau::genmatching")
            ->debug("closest genlepton {} // DeltaR {}",
                    closest_genparticle_index, min_delta_r);
        // now loop trough the gentaus and check, if they are closer to the
        // lepton than the closest lepton genparticle
        for (auto hadronicGenTau : hadronicGenTaus) {
            // check if the hadronicGenTau is closer to the lepton than the
            // closest lepton genparticle
            ROOT::Math::PtEtaPhiMVector hadronicGenTau_p4(
                pts.at(hadronicGenTau), etas.at(hadronicGenTau),
                phis.at(hadronicGenTau), masses.at(hadronicGenTau));
            float gentau_delta_r =
                ROOT::Math::VectorUtil::DeltaR(hadronicGenTau_p4, lepton_p4);
            // the decay is considered a hadronic decay (statusbit 5) if
            // 1. the hadronicGenTau pt is larger than 15 GeV
            // 2. the delta_r is smaller than 0.2
            // 3. the delta_r is smaller than the closest lepton genparticle
            // delta_r
            if (hadronicGenTau_p4.Pt() > 15 && gentau_delta_r < 0.2 &&
                gentau_delta_r < min_delta_r) {
                // statusbit 5 is hadronic tau decay
                Logger::get("genmatching::tau::genmatching")
                    ->debug(
                        "found hadronicGenTau closer than closest lepton: {}",
                        gentau_delta_r);
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_TAU_HAD_DECAY");
                return (int)GenMatchingCode::IS_TAU_HAD_DECAY;
            }
        }
        // if it is not a hadronic decay, check if the lepton is close
        // enough (deltaR < 0.2)
        int closest_pdgid = std::abs(pdgids.at(closest_genparticle_index));
        if (min_delta_r < 0.2) {
            bool prompt =
                IntBits(status_flags.at(closest_genparticle_index)).test(0);
            bool from_tau =
                IntBits(status_flags.at(closest_genparticle_index)).test(5);
            if (closest_pdgid == 11 && prompt) {
                // statusbit 1 is prompt electron
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_PROMPT");
                return (int)GenMatchingCode::IS_ELE_PROMPT;
            }
            if (closest_pdgid == 13 && prompt) {
                // statusbit 2 is prompt muon
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_PROMPT");
                return (int)GenMatchingCode::IS_MUON_PROMPT;
            }
            if (closest_pdgid == 11 && from_tau) {
                // statusbit 3 is electron from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_FROM_TAU");
                return (int)GenMatchingCode::IS_ELE_FROM_TAU;
            }
            if (closest_pdgid == 13 && from_tau) {
                // statusbit 4 is muon from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_FROM_TAU");
                return (int)GenMatchingCode::IS_MUON_FROM_TAU;
            }
        }
        // if no genlepton was found within the deltaR < 0.2, return fake
        // (statusbit 6)
        Logger::get("genmatching::tau::genmatching")->debug("IS_FAKE");
        return (int)GenMatchingCode::IS_FAKE;
    };

    auto df1 =
        df.Define(outputname, match_lepton,
                  {hadronicGenTaus, genparticles_pdgid, genparticles_statusFlag,
                   genparticles_pt, genparticles_eta, genparticles_phi,
                   genparticles_mass, lepton_p4});
    return df1;
}


// ROOT::RDF::RNode fatjet_lepton_genmatching(ROOT::RDF::RNode df, const std::string &outputname,
//                              const std::string &hadronicGenTaus,
//                              const std::string &genparticles_pdgid,
//                              const std::string &genparticles_statusFlag,
//                              const std::string &genparticles_pt,
//                              const std::string &genparticles_eta,
//                              const std::string &genparticles_phi,
//                              const std::string &genparticles_mass,
//                              const std::string &lepton_p4) {
//     auto match_lepton = [](const std::vector<int> &hadronicGenTaus,
//                            const ROOT::RVec<int> &pdgids,
//                            const ROOT::RVec<unsigned short> &status_flags,
//                            const ROOT::RVec<float> &pts,
//                            const ROOT::RVec<float> &etas,
//                            const ROOT::RVec<float> &phis,
//                            const ROOT::RVec<float> &masses,
//                            const ROOT::Math::PtEtaPhiMVector &lepton_p4) {
//         // find closest lepton fulfilling the requirements
//         float min_delta_r = 0.8;
//         int closest_genparticle_index = 0;

//         int r_value ;
//         for (unsigned int i = 0; i < pdgids.size(); i++) {
//             int pdgid = std::abs(pdgids.at(i));
            
//             bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
//                               IntBits(status_flags.at(i)).test(5));

//             bool prompt =
//                 IntBits(status_flags.at(closest_genparticle_index)).test(0);
//             bool from_tau =
//                 IntBits(status_flags.at(closest_genparticle_index)).test(5);
//             if ( pts.at(i) > 8 && statusbit) {
//                 ROOT::Math::PtEtaPhiMVector gen_p4(pts.at(i), etas.at(i),
//                                                    phis.at(i), masses.at(i));
//                 float delta_r =
//                     ROOT::Math::VectorUtil::DeltaR(gen_p4, lepton_p4);
//                 if (delta_r < 0.8) {
//                     closest_genparticle_index = i;
//                     min_delta_r = delta_r;
//                     int closest_pdgid = std::abs(pdgids.at(closest_genparticle_index));
//                     if (closest_pdgid == 11 && prompt) {
//                         // statusbit 1 is prompt electron
//                         Logger::get("genmatching::tau::fatjet_genmatching")
//                             ->info("IS_ELE_PROMPT");
//                         // return 1;
//                         r_value = 1;
//                     }
//                     if (closest_pdgid == 13 && prompt) {
//                         // statusbit 2 is prompt muon
//                         Logger::get("genmatching::tau::fatjet_genmatching")
//                             ->info("IS_MUON_PROMPT");
//                         // return 2;
//                         r_value = 2;
//                     }
//                     if (closest_pdgid == 11 && from_tau) {
//                         // statusbit 3 is electron from tau
//                         Logger::get("genmatching::tau::fatjet_genmatching")
//                             ->info("IS_ELE_FROM_TAU");
//                         // return 3;
//                         r_value = 3;
//                     }
//                     if (closest_pdgid == 13 && from_tau) {
//                         // statusbit 4 is muon from tau
//                         Logger::get("genmatching::tau::fatjet_genmatching")
//                             ->info("IS_MUON_FROM_TAU");
//                         // return 4;
//                         r_value = 4;
//                     }
                    
//                 }else{
//                     Logger::get("genmatching::tau::fatjet_genmatching")
//                     ->info("Genmatching Muon from Tau is NOT found at DeltaR {}",delta_r);
//                     // return 6;
//                     r_value = 6;
//                 }
//             }
//             return r_value;
//         }
//            Logger::get("genmatching::tau::fatjet_genmatching")
//                     ->debug("End of code r_value {}", r_value);
                    
//     };

//     auto df1 =
//         df.Define(outputname, match_lepton,
//                   {hadronicGenTaus, genparticles_pdgid, genparticles_statusFlag,
//                    genparticles_pt, genparticles_eta, genparticles_phi,
//                    genparticles_mass, lepton_p4});
//     return df1;
// }

ROOT::RDF::RNode fatjet_lepton_genmatching(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &lepton_p4) {
    auto match_lepton = [](const ROOT::RVec<int> &pdgids,
                           const ROOT::RVec<unsigned short> &status_flags,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &lepton_p4) {

        float min_delta_r = 0.8;
        int closest_genparticle_index = 0;
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            int pdgid = std::abs(pdgids.at(i));
            // check
            // 1. if there is a gen electron or muon close to the lepton
            // 2. that the genparticle pt is larger than 8 GeV
            // 3. the genparticle is isPrompt (statusbit 0) or
            // isDirectPromptTauDecayProduct (statusbit 5)
            bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
                              IntBits(status_flags.at(i)).test(5));
            if ((pdgid == 11 || pdgid == 13) && pts.at(i) > 8 && statusbit) {
                ROOT::Math::PtEtaPhiMVector gen_p4(pts.at(i), etas.at(i),
                                                   phis.at(i), masses.at(i));
                float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(gen_p4, lepton_p4);
                if (delta_r < min_delta_r) {
                    closest_genparticle_index = i;
                    min_delta_r = delta_r;
                }
            }
        }
        Logger::get("genmatching::tau::genmatching")
            ->info("closest genlepton {} // DeltaR {}",
                    closest_genparticle_index, min_delta_r);

        int closest_pdgid = std::abs(pdgids.at(closest_genparticle_index));

        if (min_delta_r < 0.8) {

            bool prompt =
                IntBits(status_flags.at(closest_genparticle_index)).test(0);
            bool from_tau =
                IntBits(status_flags.at(closest_genparticle_index)).test(5);

            if (closest_pdgid == 11 && prompt) {
                // statusbit 1 is prompt electron
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_PROMPT");
                return 1;
            }
            if (closest_pdgid == 13 && prompt) {
                // statusbit 2 is prompt muon
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_PROMPT");
                return 2;
            }
            if (closest_pdgid == 11 && from_tau) {
                // statusbit 3 is electron from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_FROM_TAU");
                return 3;
            }
            if (closest_pdgid == 13 && from_tau) {
                // statusbit 4 is muon from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_FROM_TAU");
                return 4;
            }
            Logger::get("genmatching::tau::genmatching")
                    ->debug("Lepton genmatching is not found in dR {}",min_delta_r);
            return 6;

        }
    };

        auto df1 =
        df.Define(outputname, match_lepton,
                  { genparticles_pdgid, genparticles_statusFlag,
                   genparticles_pt, genparticles_eta, genparticles_phi,
                   genparticles_mass, lepton_p4});
    return df1;

    }

ROOT::RDF::RNode fatjet_tau_had_genmatching(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &hadronicGenTaus,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &lepton_p4) {
    auto match_lepton = [](const std::vector<int> &hadronicGenTaus,
                           const ROOT::RVec<int> &pdgids,
                           const ROOT::RVec<unsigned short> &status_flags,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &lepton_p4) {
        // find closest lepton fulfilling the requirements
        float min_delta_r = 9999;
        int closest_genparticle_index = 0;

        Logger::get("genmatching::tau::fatjet_genmatching")
                        ->info("\n &&&& start genmatchine  {}");

            int gen_value = -1;


            if (hadronicGenTaus.size() != 0) {

                
                Logger::get("genmatching::tau::fatjet_genmatching")
                        ->info("\n Gentaus with non zero size  {}", hadronicGenTaus.size());
                return gen_value;
                for (auto hadronicGenTau : hadronicGenTaus) {

                    Logger::get("genmatching::tau::fatjet_genmatching")
                        ->info("Gentaus with non zero size index {}", hadronicGenTau);
                    if (hadronicGenTau != -1) {

         


                        if (pts.at(hadronicGenTau) == 0 || etas.at(hadronicGenTau) == 0 || phis.at(hadronicGenTau) == 0 || masses.at(hadronicGenTau) == 0){
                            Logger::get("genmatching::tau::fatjet_genmatching")
                            ->info("Gentaus with some 0 quantities pt {}, eta {}, phi {}, mass {}", pts.at(hadronicGenTau), etas.at(hadronicGenTau), phis.at(hadronicGenTau), masses.at(hadronicGenTau));
                            
                            Logger::get("genmatching::tau::fatjet_genmatching") ->info("Genmatching Tau is NOT found");

                            // return gen_value;
                        }
                        else{

                         Logger::get("genmatching::tau::fatjet_genmatching")
                            ->info("Gentaus with some NON 0 quantities pt {}, eta {}, phi {}, mass {}", pts.at(hadronicGenTau), etas.at(hadronicGenTau), phis.at(hadronicGenTau), masses.at(hadronicGenTau));

                        ROOT::Math::PtEtaPhiMVector hadronicGenTau_p4(
                        pts.at(hadronicGenTau), etas.at(hadronicGenTau),
                        phis.at(hadronicGenTau), masses.at(hadronicGenTau));
                        float gentau_delta_r =
                            ROOT::Math::VectorUtil::DeltaR(hadronicGenTau_p4, lepton_p4);

                        if ( gentau_delta_r < 0.8 ) {


                            gen_value = 5;
                            // return gen_value;
                        }
                    //   gen_value = 5;
                    //     return gen_value;
                    }
                    }
                }

                }

            Logger::get("genmatching::tau::fatjet_genmatching") ->info("Genmatching value is {}", gen_value);

            return gen_value;

    };

    auto df1 =
        df.Define(outputname, match_lepton,
                  {hadronicGenTaus, genparticles_pdgid, genparticles_statusFlag,
                   genparticles_pt, genparticles_eta, genparticles_phi,
                   genparticles_mass, lepton_p4});

    return df1;
}

/**
 * @brief implementation of the genmatching used in tau analyses, based
 * on
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * implementation. The genmatches are represented as integer flags:
    Decaytype         | Value
    ------------------|-------
    NONE              | -1,
    IS_ELE_PROMPT     | 1,
    IS_MUON_PROMPT    | 2,
    IS_ELE_FROM_TAU   | 3,
    IS_MUON_FROM_TAU  | 4,
    IS_TAU_HAD_DECAY  | 5,
    IS_FAKE           | 6
    IS_ELE_FROM_W     | 7,
    IS_MUON_FROM_W    | 8

    The genmatch is caluclated individually for each taupair candidate.

 *
 * @param df The dataframe to be used for the genmatching.
 * @param outputname The name of the output column.
 * @param hadronicGenTaus The name of the column containing the hadronicGenTaus
 from genmatching::tau::hadronicGenTaus
 * @param genparticles_pdgid The name of the column containing the pdgids of the
 genparticles
 * @param genparticles_statusFlag The name of the column containing the
 statusFlags of the genparticles
 * @param genparticles_pt The name of the column containing the pT of the
 genparticles
 * @param genparticles_eta The name of the column containing the eta of the
 genparticles
 * @param genparticles_phi The name of the column containing the phi of the
 genparticles
 * @param genparticles_mass The name of the column containing the mass of the
 genparticles
 * @param genparticle_motheridx The name of the column containing the mother
 indices of the genparticles
 * @param genparticles_status The name of the column containing the status of
the genparticles
 * @param lepton_p4 The name of the column containing the p4 of the lepton to be
 genmatched
 * @return a dataframe with the genmatching as a column named outputname
 */
ROOT::RDF::RNode genmatching_wh(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &hadronicGenTaus, const std::string &genparticles_pdgid,
    const std::string &genparticles_statusFlag,
    const std::string &genparticles_pt, const std::string &genparticles_eta,
    const std::string &genparticles_phi, const std::string &genparticles_mass,
    const std::string &genparticle_motheridx,
    const std::string &genparticles_status, const std::string &lepton_p4) {
    auto match_lepton = [](const std::vector<int> &hadronicGenTaus,
                           const ROOT::RVec<int> &pdgids,
                           const ROOT::RVec<int> &mother_idx,
                           const ROOT::RVec<int> &status_flags,
                           const ROOT::RVec<int> &status,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &lepton_p4) {
        // find closest lepton fulfilling the requirements
        float min_delta_r = 9999;
        int closest_genparticle_index = 0;
        int closest_genparticle_mother_pdgid = 0;
        int closest_genparticle_mother_statusFlag = 0;
        int closest_genparticle_mother_status = 0;
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            int pdgid = std::abs(pdgids.at(i));
            // check
            // 1. if there is a gen electron or muon close to the lepton
            // 2. that the genparticle pt is larger than 8 GeV
            // 3. the genparticle is isPrompt (statusbit 0) or
            // isDirectPromptTauDecayProduct (statusbit 5)
            bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
                              IntBits(status_flags.at(i)).test(5));
            if ((pdgid == 11 || pdgid == 13) && pts.at(i) > 8 && statusbit) {
                ROOT::Math::PtEtaPhiMVector gen_p4(pts.at(i), etas.at(i),
                                                   phis.at(i), masses.at(i));
                float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(gen_p4, lepton_p4);
                if (delta_r < min_delta_r) {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("pdgids {}, status {}, status_flags {}, "
                                "mother_idx {}",
                                pdgids, status, status_flags, mother_idx);
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("mother_idx {}, pdgids {}, status {}, "
                                "status_flags {}",
                                mother_idx.at(i), pdgids.at(i),
                                status_flags.at(i), status.at(i));
                    if (mother_idx.at(i) == -1) {
                        closest_genparticle_index = i;
                        closest_genparticle_mother_pdgid = pdgids.at(i);
                        closest_genparticle_mother_status = status.at(i);
                        closest_genparticle_mother_statusFlag =
                            status_flags.at(i);
                        min_delta_r = delta_r;
                    } else {
                        closest_genparticle_index = i;
                        closest_genparticle_mother_pdgid =
                            pdgids.at(mother_idx.at(i));
                        closest_genparticle_mother_status =
                            status.at(mother_idx.at(i));
                        closest_genparticle_mother_statusFlag =
                            status_flags.at(mother_idx.at(i));
                        min_delta_r = delta_r;
                    }
                }
            }
        }
        Logger::get("genmatching::tau::genmatching")
            ->debug("closest genlepton {} // DeltaR {}",
                    closest_genparticle_index, min_delta_r);
        // now loop through the gentaus and check, if they are closer to the
        // lepton than the closest lepton genparticle
        for (auto hadronicGenTau : hadronicGenTaus) {
            // check if the hadronicGenTau is closer to the lepton than the
            // closest lepton genparticle
            ROOT::Math::PtEtaPhiMVector hadronicGenTau_p4(
                pts.at(hadronicGenTau), etas.at(hadronicGenTau),
                phis.at(hadronicGenTau), masses.at(hadronicGenTau));
            float gentau_delta_r =
                ROOT::Math::VectorUtil::DeltaR(hadronicGenTau_p4, lepton_p4);
            // the decay is considered a hadronic decay (statusbit 5) if
            // 1. the hadronicGenTau pt is larger than 15 GeV
            // 2. the delta_r is smaller than 0.2
            // 3. the delta_r is smaller than the closest lepton genparticle
            // delta_r
            if (hadronicGenTau_p4.Pt() > 15 && gentau_delta_r < 0.2 &&
                gentau_delta_r < min_delta_r) {
                // statusbit 5 is hadronic tau decay
                Logger::get("genmatching::tau::genmatching")
                    ->debug(
                        "found hadronicGenTau closer than closest lepton: {}",
                        gentau_delta_r);
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_TAU_HAD_DECAY");
                return (int)GenMatchingCode::IS_TAU_HAD_DECAY;
            }
        }
        // if it is not a hadronic decay, check if the lepton is close
        // enough (deltaR < 0.2)
        int closest_pdgid = std::abs(pdgids.at(closest_genparticle_index));
        if (min_delta_r < 0.2) {
            bool prompt =
                IntBits(status_flags.at(closest_genparticle_index)).test(0);
            bool from_tau =
                IntBits(status_flags.at(closest_genparticle_index)).test(5);
            if (closest_pdgid == 11 && prompt) {
                // statusbit 1 is prompt electron
                if (abs(closest_genparticle_mother_pdgid) == 24) {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_ELE_FROM_W");
                    return (int)GenMatchingCode::IS_ELE_FROM_W;
                } else {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_ELE_PROMPT");
                    return (int)GenMatchingCode::IS_ELE_PROMPT;
                }
            }
            if (closest_pdgid == 13 && prompt) {
                // statusbit 2 is prompt muon
                if (abs(closest_genparticle_mother_pdgid) == 24) {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_MUON_FROM_W");
                    return (int)GenMatchingCode::IS_MUON_FROM_W;
                } else {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_MUON_PROMPT");
                    return (int)GenMatchingCode::IS_MUON_PROMPT;
                }
            }
            if (closest_pdgid == 11 && from_tau) {
                // statusbit 3 is electron from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_FROM_TAU");
                return (int)GenMatchingCode::IS_ELE_FROM_TAU;
            }
            if (closest_pdgid == 13 && from_tau) {
                // statusbit 4 is muon from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_FROM_TAU");
                return (int)GenMatchingCode::IS_MUON_FROM_TAU;
            }
        }
        // if no genlepton was found within the deltaR < 0.2, return fake
        // (statusbit 6)
        Logger::get("genmatching::tau::genmatching")->debug("IS_FAKE");
        return (int)GenMatchingCode::IS_FAKE;
    };

    auto df1 = df.Define(
        outputname, match_lepton,
        {hadronicGenTaus, genparticles_pdgid, genparticle_motheridx,
         genparticles_statusFlag, genparticles_status, genparticles_pt,
         genparticles_eta, genparticles_phi, genparticles_mass, lepton_p4});
    return df1;
}
/**
 * @brief function to find all hadronicGenTaus needed for the genmatching, based
 * on
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * implementation. Loop trough all genparticles and check, if a genparticle,
 that is prompt, without any leptonic daughters can be found. If this
 genparticle has a neutrino as daughter, the hadronicGenTau is added to the list
 of hadronicGenTaus.
 *
 * @param df The dataframe to be extended with the hadronicGenTaus
 * @param outputname The name of the column to be added to the dataframe
 * @param genparticles_pdgid The name of the column containing the pdgids of the
 genparticles
 * @param genparticles_statusFlag The name of the column containing the
 statusflags of the genparticles
 * @param genparticles_motherid The name of the column containing the motherids
 of the genparticles
 * @return a new dataframe with the hadronicGenTaus added to the dataframe

 */

ROOT::RDF::RNode hadronicGenTaus(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &genparticles_pdgid,
                                 const std::string &genparticles_statusFlag,
                                 const std::string &genparticles_motherid) {

    auto gentaus = [](const ROOT::RVec<int> &pdgids,
                      const ROOT::RVec<unsigned short> &status_flags,
                      const ROOT::RVec<short> &mother_index) {
        // set default values for the output
        std::vector<int> hadronicGenTaus;
        if (pdgids.size() == 0) {
            hadronicGenTaus.push_back(-1);
            return hadronicGenTaus;
        }
        // debug printout
        for (int i = 0; i < pdgids.size(); i++) {
            Logger::get("genmatching::tau::genpair")
                ->debug("genparticle index {}, pdgid: {}, status_flags {}, "
                        "mother_index {}",
                        i, pdgids.at(i), status_flags.at(i),
                        mother_index.at(i));
        }
        // now loop though the genparticles and find the taus
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            // check if the particle is a tau
            if (std::abs(pdgids.at(i)) == 15) {
                // check if the particle is a stable one
                bool prompt = IntBits(status_flags.at(i)).test(0);
                if (prompt) {
                    Logger::get("genmatching::tau::genpair")
                        ->debug("Found prompt tau: {}", i);
                    // find all daughters of the tau by checking which particles
                    // have the tauindex i as mother
                    std::vector<int> daughters;
                    for (unsigned int j = 0; j < mother_index.size(); j++) {
                        if (mother_index.at(j) == i) {
                            daughters.push_back(j);
                            Logger::get("genmatching::tau::genpair")
                                ->debug("daughters of {} : {}", i, j);
                        }
                    }
                    // check if the tau has at least one daughter
                    if (daughters.size() > 0) {
                        // check if the tau has at least one daughter that is a
                        // tau
                        bool hasTauDaughter = false;
                        bool hasLeptonDaughter = false;
                        for (unsigned int j = 0; j < daughters.size(); j++) {
                            int daughter_pdgid =
                                std::abs(pdgids.at(daughters.at(j)));
                            Logger::get("genmatching::tau::genpair")
                                ->debug("daughter {} : pdgid {}", j,
                                        daughter_pdgid);
                            if (daughter_pdgid == 15) {
                                hasTauDaughter = true;
                            }
                            if (daughter_pdgid == 11 || daughter_pdgid == 13) {
                                hasLeptonDaughter = true;
                            }
                        }
                        // in this case, it is not the correct tau, continue
                        if (hasTauDaughter || hasLeptonDaughter) {
                            continue;
                        }
                        // now run again through the daughters and check if
                        // there is a neutrino, if this is the case, we have the
                        // correct gentau
                        for (unsigned int j = 0; j < daughters.size(); j++) {
                            int daughter_pdgid =
                                std::abs(pdgids.at(daughters.at(j)));
                            if (daughter_pdgid == 12 || daughter_pdgid == 14 ||
                                daughter_pdgid == 16) {
                                Logger::get("genmatching::tau::genpair")
                                    ->debug("gentau found: {}", i);
                                hadronicGenTaus.push_back(i);
                            }
                        }
                    }
                }
            }
        }
        Logger::get("genmatching::tau::genpair")
            ->debug("found {} hadronic hadronicGenTaus",
                    hadronicGenTaus.size());
        for (int i = 0; i < hadronicGenTaus.size(); i++) {
            Logger::get("genmatching::tau::genpair")
                ->debug("hadronicGenTaus {} : {}", i, hadronicGenTaus.at(i));
        }
        return hadronicGenTaus;
    };

    auto df1 = df.Define(
        outputname, gentaus,
        {genparticles_pdgid, genparticles_statusFlag, genparticles_motherid});
    return df1;
}

ROOT::RDF::RNode trigger_mu_in_fatjet(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &fatjet_p4,
                                 const std::string &muon_pt,
                                 const std::string &muon_eta,
                                 const std::string &muon_phi,
                                 const std::string &muon_mass) {

    auto match_muon_with_fatjet = [](const ROOT::Math::PtEtaPhiMVector &fatjet_p4,
                                    const ROOT::RVec<float> &muon_pts,
                                    const ROOT::RVec<float> &muon_etas,
                                    const ROOT::RVec<float> &muon_phis,
                                    const ROOT::RVec<float> &muon_masses) {

         for (unsigned int i = 0; i < muon_pts.size(); i++) {

            ROOT::Math::PtEtaPhiMVector muon_p4(muon_pts.at(i), muon_etas.at(i),
                    muon_phis.at(i), muon_masses.at(i));
            float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(muon_p4, fatjet_p4);

            if (delta_r < 0.8 && muon_pts.at(i) > 50 ){
                Logger::get("fatjet::trigger_mu_in_fatjet")
                    ->debug("Muon found in fatjet with pt {} and eta {}",
                            muon_pts.at(i), muon_etas.at(i));
                return 1;
            }
            else {
                Logger::get("fatjet::trigger_mu_in_fatjet")
                ->debug("Muon NOT found in fatjet with pt {} and eta {}",
                    muon_pts.at(i), muon_etas.at(i));
                return 0;
            }
         }                                                    
                                 
};

    auto df1 =
        df.Define(outputname, match_muon_with_fatjet,
                  {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass});
    return df1;
    }

} // end namespace tau

} // end namespace genmatching

#endif /* GUARD_GENPARTICLES_H */