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
    IS_FAKE = 6
};

namespace genmatching {

namespace tau {
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
                           const ROOT::RVec<int> &status_flags,
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
            // 1. if there is a gen electron or muos close to the lepton
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
            ->debug("closted genelepton {} // DeltaR {}",
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
/**
 * @brief function to find all hadronicGenTaus needed for the genmatching, based
 * on
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * implementation
 *
 * @param df
 * @param outputname
 * @param genparticles_pdgid
 * @param genparticles_statusFlag
 * @param genparticles_motherid
 * @return ROOT::RDF::RNode
 */

ROOT::RDF::RNode hadronicGenTaus(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &genparticles_pdgid,
                                 const std::string &genparticles_statusFlag,
                                 const std::string &genparticles_motherid) {

    auto gentaus = [](const ROOT::RVec<int> &pdgids,
                      const ROOT::RVec<int> &status_flags,
                      const ROOT::RVec<int> &mother_index) {
        // set default values for the output
        std::vector<int> hadronicGenTaus;
        if (pdgids.size() == 0) {
            hadronicGenTaus.push_back(-1);
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

} // end namespace tau

} // end namespace genmatching

#endif /* GUARD_GENPARTICLES_H */