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

enum class MatchingGenTauCode : int {
    NONE = -1,
    IS_ELE_PROMPT = 1,
    IS_MUON_PROMPT = 2,
    IS_ELE_FROM_TAU = 3,
    IS_MUON_FROM_TAU = 4,
    IS_TAU_HAD_DECAY = 5,
    IS_FAKE = 6,
    IS_ELE_PROMPT_FROM_W = 7,
    IS_MUON_PROMPT_FROM_W = 8
};

namespace genparticles {
namespace tau {

/**
 * @brief This function finds all hadronic generator-level taus needed for the 
 * matching to the reconstructed hadronic taus, based on the implementation from
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 *
 * The procedure is go trough all genparticles and check, if a genparticle,
 * that is prompt, without any leptonic daughters can be found. If this
 * genparticle has a neutrino as daughter, this genparticle is indentified as
 * a hadronic generator-level tau (GenTau). 
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the hadronic tau 
 * indices vector
 * @param genparticles_pdg_id name of the column containing the PDG IDs of the
 * genparticles
 * @param genparticles_status_flags name of the column containing the status 
 * flags of the genparticles, e.g. isPrompt, isHardProcess, isLastCopy, ...
 * @param genparticles_mother_index name of the column containing the mother 
 * particle indices of the genparticles
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode HadronicGenTaus(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &genparticles_pdg_id,
                                 const std::string &genparticles_status_flags,
                                 const std::string &genparticles_mother_index) {

    auto gentaus = [](const ROOT::RVec<int> &pdg_ids,
                      const ROOT::RVec<int> &status_flags,
                      const ROOT::RVec<int> &mother_indices) {
        // set default values for the output
        std::vector<int> hadGenTaus;
        if (pdg_ids.size() == 0) {
            hadGenTaus.push_back(-1);
            return hadGenTaus;
        }
        // debug printout
        for (int i = 0; i < pdg_ids.size(); i++) {
            Logger::get("genparticles::tau::HadronicGenTaus")
                ->debug("genparticle index {}, pdgid: {}, status_flags {}, "
                        "mother_index {}",
                        i, pdg_ids.at(i), status_flags.at(i),
                        mother_indices.at(i));
        }
        // now loop though the genparticles and find the taus
        for (unsigned int i = 0; i < pdg_ids.size(); i++) {
            // check if the particle is a tau
            if (std::abs(pdg_ids.at(i)) == 15) {
                // check if the particle is a stable one
                bool prompt = IntBits(status_flags.at(i)).test(0);
                if (prompt) {
                    Logger::get("genparticles::tau::HadronicGenTaus")
                        ->debug("Found prompt tau: {}", i);
                    // find all daughters of the tau by checking which particles
                    // have the tauindex i as mother
                    std::vector<int> daughters;
                    for (unsigned int j = 0; j < mother_indices.size(); j++) {
                        if (mother_indices.at(j) == i) {
                            daughters.push_back(j);
                            Logger::get("genparticles::tau::HadronicGenTaus")
                                ->debug("daughters of {} : {}", i, j);
                        }
                    }
                    // check if the tau has at least one daughter
                    if (daughters.size() > 0) {
                        // check if the tau has at least one daughter that is a
                        // tau -> indicates a leptonic decay
                        bool hasTauDaughter = false;
                        bool hasLeptonDaughter = false;
                        for (unsigned int j = 0; j < daughters.size(); j++) {
                            int daughter_pdgid =
                                std::abs(pdg_ids.at(daughters.at(j)));
                            Logger::get("genparticles::tau::HadronicGenTaus")
                                ->debug("daughter {} : pdgid {}", j,
                                        daughter_pdgid);
                            if (daughter_pdgid == 15) {
                                hasTauDaughter = true;
                            }
                            if (daughter_pdgid == 11 || daughter_pdgid == 13) {
                                hasLeptonDaughter = true;
                            }
                        }
                        // in this case, the tau decayed leptonically or is not 
                        // in its final state, therefore, it is not the correct 
                        // tau, continue
                        if (hasTauDaughter || hasLeptonDaughter) {
                            continue;
                        }
                        // now run again through the daughters and check if
                        // there is a neutrino, if this is the case, we have the
                        // correct gentau
                        for (unsigned int j = 0; j < daughters.size(); j++) {
                            int daughter_pdgid =
                                std::abs(pdg_ids.at(daughters.at(j)));
                            if (daughter_pdgid == 12 || daughter_pdgid == 14 ||
                                daughter_pdgid == 16) {
                                Logger::get("genparticles::tau::HadronicGenTaus")
                                    ->debug("gentau found: {}", i);
                                hadGenTaus.push_back(i);
                            }
                        }
                    }
                }
            }
        }
        Logger::get("genparticles::tau::HadronicGenTaus")
            ->debug("found {} hadronic GenTaus",
                    hadGenTaus.size());
        for (int i = 0; i < hadGenTaus.size(); i++) {
            Logger::get("genparticles::tau::HadronicGenTaus")
                ->debug("hadronicGenTaus {} : {}", i, hadGenTaus.at(i));
        }
        return hadGenTaus;
    };
    auto df1 = df.Define(
        outputname, gentaus,
        {genparticles_pdg_id, genparticles_status_flags, genparticles_mother_index});
    return df1;
}

/**
 * @brief This function determines the true origin of a reconstructed hadronic tau 
 * by matching it to generator-level particles. The implementation is based on 
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * 
 * The matching is represented by integer flags:
 *   Decaytype             | Value
 *   ----------------------|-------
 *   IS_ELE_PROMPT         | 1
 *   IS_MUON_PROMPT        | 2
 *   IS_ELE_FROM_TAU       | 3
 *   IS_MUON_FROM_TAU      | 4
 *   IS_TAU_HAD_DECAY      | 5
 *   IS_FAKE (not matched) | 6
 *
 * The matching logic is as follows:
 * 1. For each reconstructed tau, first, the closest "prompt" or "from tau decay"
 *    generator-level electron or muon with \f$p_T\f$ > 8 GeV is found. The distance 
 *    (\f$\Delta R\f$) to this lepton is saved.
 * 2. Next, an iteration is done through pre-identified generator-level hadronic taus. 
 *    If a gen. tau with \f$p_T\f$ > 15 GeV is found within a cone of \f$\Delta R\f$ < 0.2 
 *    of the reco. tau, and it is closer than the closest electron/muon found in step 1, 
 *    the match is classified as `IS_TAU_HAD_DECAY` (value: 5).
 * 3. If no such hadronic tau is found, it re-evaluates the closest electron/muon from 
 *    step 1. If this lepton is within \f$\Delta R\f$ < 0.2 of the reco. tau, the match is 
 *    classified based on the lepton's identity and origin:
 *    - Prompt electron (other): `IS_ELE_PROMPT` (value: 1)
 *    - Prompt muon (other): `IS_MUON_PROMPT` (value: 2)
 *    - Electron from a tau decay: `IS_ELE_FROM_TAU` (value: 3)
 *    - Muon from a tau decay: `IS_MUON_FROM_TAU` (value: 4)
 * 4. If nothing of the above is matched, the reco. tau is classified as a `IS_FAKE` 
 *    (value: 6).
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the gen. matching value
 * @param hadronic_gen_taus name of the column containing the hadronic gen. tau indices
 * found with `genparticles::tau::HadronicGenTaus`
 * @param genparticles_pdg_id name of the column containing the PDG IDs of the
 * genparticles
 * @param genparticles_status_flags name of the column containing the status 
 * flags of the genparticles, e.g. isPrompt, isHardProcess, 
 * isLastCopy, ...
 * @param genparticles_pt name of the column containing the \f$p_T\f$ of the
 * genparticles
 * @param genparticles_eta name of the column containing the \f$\eta\f$ of the
 * genparticles
 * @param genparticles_phi name of the column containing the \f$\phi\f$ of the
 * genparticles
 * @param genparticles_mass name of the column containing the mass of the
 * genparticles
 * @param reco_had_tau name of the column containing the Lorentz vector of the 
 * reconstructed hadronic tau lepton
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode GenMatching(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &hadronic_gen_taus, const std::string &genparticles_pdg_id,
    const std::string &genparticles_status_flags, const std::string &genparticles_pt, 
    const std::string &genparticles_eta, const std::string &genparticles_phi, 
    const std::string &genparticles_mass, const std::string &reco_had_tau) {
    auto match_tau = [](const std::vector<int> &had_gen_taus,
                           const ROOT::RVec<int> &pdg_ids,
                           const ROOT::RVec<int> &status_flags,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &reco_had_tau) {
        // find closest lepton fulfilling the requirements
        float min_delta_r = 9999;
        int closest_genparticle_index = 0;

        Logger::get("genparticles::tau::GenMatching")
            ->debug("pdg_ids {}, status_flags {}", pdg_ids, status_flags);

        for (unsigned int i = 0; i < pdg_ids.size(); i++) {
            int pdgid = std::abs(pdg_ids.at(i));
            // check
            // 1. if there is a gen electron or muon close to the lepton
            // 2. that the genparticle pt is larger than 8 GeV
            // 3. the genparticle is isPrompt (statusbit 0) or
            // isDirectPromptTauDecayProduct (statusbit 5)
            bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
                              IntBits(status_flags.at(i)).test(5));
            if ((pdgid == 11 || pdgid == 13) && pts.at(i) > 8 && statusbit) {
                ROOT::Math::PtEtaPhiMVector probe_gen_tau(pts.at(i), etas.at(i),
                                                          phis.at(i), masses.at(i));
                float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(probe_gen_tau, reco_had_tau);
                if (delta_r < min_delta_r) {
                    Logger::get("genparticles::tau::GenMatching")
                        ->debug("pdg_ids {}, status_flags {}",
                                pdg_ids.at(i), status_flags.at(i));
                    closest_genparticle_index = i;
                    min_delta_r = delta_r;
                }
            }
        }
        Logger::get("genparticles::tau::GenMatching")
            ->debug("closest genlepton {} // DeltaR {}",
                    closest_genparticle_index, min_delta_r);
        // now loop through the gentaus and check, if they are closer to the
        // lepton than the closest lepton genparticle
        for (auto gen_tau : had_gen_taus) {
            // check if the gen_tau is closer to the lepton than the
            // closest lepton genparticle
            ROOT::Math::PtEtaPhiMVector probe_gen_tau(
                pts.at(gen_tau), etas.at(gen_tau),
                phis.at(gen_tau), masses.at(gen_tau));
            float gen_tau_delta_r =
                ROOT::Math::VectorUtil::DeltaR(probe_gen_tau, reco_had_tau);
            // the decay is considered a hadronic decay (statusbit 5) if
            // 1. the hadronic gen. tau pt is larger than 15 GeV
            // 2. the delta_r is smaller than 0.2
            // 3. the delta_r is smaller than the closest lepton genparticle
            // delta_r
            if (probe_gen_tau.Pt() > 15 && gen_tau_delta_r < 0.2 &&
                gen_tau_delta_r < min_delta_r) {
                // statusbit 5 is hadronic tau decay
                Logger::get("genparticles::tau::GenMatching")
                    ->debug(
                        "found hadronic gen. tau closer than closest lepton: {}",
                        gen_tau_delta_r);
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_TAU_HAD_DECAY");
                return (int)MatchingGenTauCode::IS_TAU_HAD_DECAY;
            }
        }
        // if it is not a hadronic decay, check if the lepton is close
        // enough (deltaR < 0.2)
        int closest_pdgid = std::abs(pdg_ids.at(closest_genparticle_index));
        if (min_delta_r < 0.2) {
            bool prompt =
                IntBits(status_flags.at(closest_genparticle_index)).test(0);
            bool from_tau =
                IntBits(status_flags.at(closest_genparticle_index)).test(5);
            if (closest_pdgid == 11 && prompt) {
                // statusbit 1 is prompt electron
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_ELE_PROMPT");
                return (int)MatchingGenTauCode::IS_ELE_PROMPT;
            }
            if (closest_pdgid == 13 && prompt) {
                // statusbit 2 is prompt muon
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_MUON_PROMPT");
                return (int)MatchingGenTauCode::IS_MUON_PROMPT;
            }
            if (closest_pdgid == 11 && from_tau) {
                // statusbit 3 is electron from tau
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_ELE_FROM_TAU");
                return (int)MatchingGenTauCode::IS_ELE_FROM_TAU;
            }
            if (closest_pdgid == 13 && from_tau) {
                // statusbit 4 is muon from tau
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_MUON_FROM_TAU");
                return (int)MatchingGenTauCode::IS_MUON_FROM_TAU;
            }
        }
        // if no genlepton was found within the deltaR < 0.2, return fake
        // (statusbit 6)
        Logger::get("genparticles::tau::GenMatching")->debug("IS_FAKE");
        return (int)MatchingGenTauCode::IS_FAKE;
    };
    auto df1 = df.Define(
        outputname, match_tau,
        {hadronic_gen_taus, genparticles_pdg_id, genparticles_status_flags, 
         genparticles_pt, genparticles_eta, genparticles_phi, genparticles_mass, reco_had_tau});
    return df1;
}

/**
 * @brief This function determines the true origin of a reconstructed hadronic tau 
 * by matching it to generator-level particles. The implementation is based on 
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * 
 * @note This function additionally matches if the prompt electron/muon decayed from a W boson.
 *
 * The matching is represented by integer flags:
 *   Decaytype             | Value
 *   ----------------------|-------
 *   IS_ELE_PROMPT         | 1
 *   IS_MUON_PROMPT        | 2
 *   IS_ELE_FROM_TAU       | 3
 *   IS_MUON_FROM_TAU      | 4
 *   IS_TAU_HAD_DECAY      | 5
 *   IS_FAKE (not matched) | 6
 *   IS_ELE_PROMPT_FROM_W  | 7
 *   IS_MUON_PROMPT_FROM_W | 8
 *
 * The matching logic is as follows:
 * 1. For each reconstructed tau, first, the closest "prompt" or "from tau decay"
 *    generator-level electron or muon with \f$p_T\f$ > 8 GeV is found. The distance 
 *    (\f$\Delta R\f$) to this lepton is saved.
 * 2. Next, an iteration is done through pre-identified generator-level hadronic taus. 
 *    If a gen. tau with \f$p_T\f$ > 15 GeV is found within a cone of \f$\Delta R\f$ < 0.2 
 *    of the reco. tau, and it is closer than the closest electron/muon found in step 1, 
 *    the match is classified as `IS_TAU_HAD_DECAY` (value: 5).
 * 3. If no such hadronic tau is found, it re-evaluates the closest electron/muon from 
 *    step 1. If this lepton is within \f$\Delta R\f$ < 0.2 of the reco. tau, the match is 
 *    classified based on the lepton's identity and origin:
 *    - Prompt electron (from W boson): `IS_ELE_PROMPT_FROM_W` (value: 7)
 *    - Prompt electron (other): `IS_ELE_PROMPT` (value: 1)
 *    - Prompt muon (from W boson): `IS_MUON_PROMPT_FROM_W` (value: 8)
 *    - Prompt muon (other): `IS_MUON_PROMPT` (value: 2)
 *    - Electron from a tau decay: `IS_ELE_FROM_TAU` (value: 3)
 *    - Muon from a tau decay: `IS_MUON_FROM_TAU` (value: 4)
 * 4. If nothing of the above is matched, the reco. tau is classified as a `IS_FAKE` 
 *    (value: 6).
 *
 * @param df input dataframe
 * @param outputname name of the output column containing the gen. matching value
 * @param hadronic_gen_taus name of the column containing the hadronic gen. tau indices
 * found with `genparticles::tau::HadronicGenTaus`
 * @param genparticles_pdg_id name of the column containing the PDG IDs of the
 * genparticles
 * @param genparticles_status_flags name of the column containing the status 
 * flags of the genparticles, e.g. isPrompt, isHardProcess, 
 * isLastCopy, ...
 * @param genparticles_mother_index name of the column containing the mother 
 * particle indices of the genparticles
 * @param genparticles_pt name of the column containing the \f$p_T\f$ of the
 * genparticles
 * @param genparticles_eta name of the column containing the \f$\eta\f$ of the
 * genparticles
 * @param genparticles_phi name of the column containing the \f$\phi\f$ of the
 * genparticles
 * @param genparticles_mass name of the column containing the mass of the
 * genparticles
 * @param reco_had_tau name of the column containing the Lorentz vector of the 
 * reconstructed hadronic tau lepton
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode GenMatching(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &hadronic_gen_taus, const std::string &genparticles_pdg_id,
    const std::string &genparticles_status_flags, const std::string &genparticles_mother_index,
    const std::string &genparticles_pt, const std::string &genparticles_eta, 
    const std::string &genparticles_phi, const std::string &genparticles_mass, 
    const std::string &reco_had_tau) {
    auto match_tau = [](const std::vector<int> &had_gen_taus,
                           const ROOT::RVec<int> &pdg_ids,
                           const ROOT::RVec<int> &status_flags,
                           const ROOT::RVec<int> &mother_idx,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &reco_had_tau) {
        // find closest lepton fulfilling the requirements
        float min_delta_r = 9999;
        int closest_genparticle_index = 0;
        int closest_genparticle_mother_pdgid = 0;
        int closest_genparticle_mother_statusFlag = 0;

        Logger::get("genparticles::tau::GenMatching")
            ->debug("pdg_ids {}, status_flags {}, mother_idx {}",
                    pdg_ids, status_flags, mother_idx);

        for (unsigned int i = 0; i < pdg_ids.size(); i++) {
            int pdgid = std::abs(pdg_ids.at(i));
            // check
            // 1. if there is a gen electron or muon close to the lepton
            // 2. that the genparticle pt is larger than 8 GeV
            // 3. the genparticle is isPrompt (statusbit 0) or
            // isDirectPromptTauDecayProduct (statusbit 5)
            bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
                              IntBits(status_flags.at(i)).test(5));
            if ((pdgid == 11 || pdgid == 13) && pts.at(i) > 8 && statusbit) {
                ROOT::Math::PtEtaPhiMVector probe_gen_tau(pts.at(i), etas.at(i),
                                                          phis.at(i), masses.at(i));
                float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(probe_gen_tau, reco_had_tau);
                if (delta_r < min_delta_r) {
                    Logger::get("genparticles::tau::GenMatching")
                        ->debug("mother_idx {}, pdg_ids {}, status_flags {}",
                                mother_idx.at(i), pdg_ids.at(i), status_flags.at(i));
                    if (mother_idx.at(i) == -1) {
                        // mother index of -1 means that there is no mother particle
                        // usually this happens for the particles of the initial pp
                        // collision 
                        closest_genparticle_index = i;
                        closest_genparticle_mother_pdgid = pdg_ids.at(i);
                        closest_genparticle_mother_statusFlag =
                            status_flags.at(i);
                        min_delta_r = delta_r;
                    } else {
                        closest_genparticle_index = i;
                        closest_genparticle_mother_pdgid =
                            pdg_ids.at(mother_idx.at(i));
                        closest_genparticle_mother_statusFlag =
                            status_flags.at(mother_idx.at(i));
                        min_delta_r = delta_r;
                    }
                }
            }
        }
        Logger::get("genparticles::tau::GenMatching")
            ->debug("closest genlepton {} // DeltaR {}",
                    closest_genparticle_index, min_delta_r);
        // now loop through the gentaus and check, if they are closer to the
        // lepton than the closest lepton genparticle
        for (auto gen_tau : had_gen_taus) {
            // check if the gen_tau is closer to the lepton than the
            // closest lepton genparticle
            ROOT::Math::PtEtaPhiMVector probe_gen_tau(
                pts.at(gen_tau), etas.at(gen_tau),
                phis.at(gen_tau), masses.at(gen_tau));
            float gen_tau_delta_r =
                ROOT::Math::VectorUtil::DeltaR(probe_gen_tau, reco_had_tau);
            // the decay is considered a hadronic decay (statusbit 5) if
            // 1. the hadronic gen. tau pt is larger than 15 GeV
            // 2. the delta_r is smaller than 0.2
            // 3. the delta_r is smaller than the closest lepton genparticle
            // delta_r
            if (probe_gen_tau.Pt() > 15 && gen_tau_delta_r < 0.2 &&
                gen_tau_delta_r < min_delta_r) {
                // statusbit 5 is hadronic tau decay
                Logger::get("genparticles::tau::GenMatching")
                    ->debug(
                        "found hadronic gen. tau closer than closest lepton: {}",
                        gen_tau_delta_r);
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_TAU_HAD_DECAY");
                return (int)MatchingGenTauCode::IS_TAU_HAD_DECAY;
            }
        }
        // if it is not a hadronic decay, check if the lepton is close
        // enough (deltaR < 0.2)
        int closest_pdgid = std::abs(pdg_ids.at(closest_genparticle_index));
        if (min_delta_r < 0.2) {
            bool prompt =
                IntBits(status_flags.at(closest_genparticle_index)).test(0);
            bool from_tau =
                IntBits(status_flags.at(closest_genparticle_index)).test(5);
            if (closest_pdgid == 11 && prompt) {
                // statusbit 7 is prompt electron from W boson
                // statusbit 1 is prompt electron
                if (abs(closest_genparticle_mother_pdgid) == 24) {
                    Logger::get("genparticles::tau::GenMatching")
                        ->debug("IS_ELE_PROMPT_FROM_W");
                    return (int)MatchingGenTauCode::IS_ELE_PROMPT_FROM_W;
                } else {
                    Logger::get("genparticles::tau::GenMatching")
                        ->debug("IS_ELE_PROMPT");
                    return (int)MatchingGenTauCode::IS_ELE_PROMPT;
                }
            }
            if (closest_pdgid == 13 && prompt) {
                // statusbit 8 is prompt muon from W boson
                // statusbit 2 is prompt muon
                if (abs(closest_genparticle_mother_pdgid) == 24) {
                    Logger::get("genparticles::tau::GenMatching")
                        ->debug("IS_MUON_PROMPT_FROM_W");
                    return (int)MatchingGenTauCode::IS_MUON_PROMPT_FROM_W;
                } else {
                    Logger::get("genparticles::tau::GenMatching")
                        ->debug("IS_MUON_PROMPT");
                    return (int)MatchingGenTauCode::IS_MUON_PROMPT;
                }
            }
            if (closest_pdgid == 11 && from_tau) {
                // statusbit 3 is electron from tau
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_ELE_FROM_TAU");
                return (int)MatchingGenTauCode::IS_ELE_FROM_TAU;
            }
            if (closest_pdgid == 13 && from_tau) {
                // statusbit 4 is muon from tau
                Logger::get("genparticles::tau::GenMatching")
                    ->debug("IS_MUON_FROM_TAU");
                return (int)MatchingGenTauCode::IS_MUON_FROM_TAU;
            }
        }
        // if no genlepton was found within the deltaR < 0.2, return fake
        // (statusbit 6)
        Logger::get("genparticles::tau::GenMatching")->debug("IS_FAKE");
        return (int)MatchingGenTauCode::IS_FAKE;
    };
    auto df1 = df.Define(
        outputname, match_tau,
        {hadronic_gen_taus, genparticles_pdg_id, genparticles_status_flags,
         genparticles_mother_index, genparticles_pt, genparticles_eta, 
         genparticles_phi, genparticles_mass, reco_had_tau});
    return df1;
}
} // end namespace tau
} // end namespace genparticles
#endif /* GUARD_GENPARTICLES_H */
