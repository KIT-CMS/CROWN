#ifndef GUARD_TRIGGERS_H
#define GUARD_TRIGGERS_H

#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "bitset"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include <regex>

typedef std::bitset<30> IntBits;

namespace trigger {

/**
 * @brief Function used to try and match a object with a trigger object. An
object is successfully matched, if they overlap within the given deltaR cone
(matchDeltaR), the pt and eta thresholds are met, and the triggerobject id and
the triggerobject bit are set to the correct value.

For the trigger objects, two additional quantities are stored, the id of the
particle and some further information on the trigger path encoded in a bitmap.
The triggerobject ID is encoded as follows:
    @code
    1 = Jet
    2 = MET
    3 = HT
    4 = MHT
    6 = FatJet
    11 = Electron (PixelMatched e/gamma)
    22 = Photon (PixelMatch-vetoed e/gamma)
    13 = Muon
    15 = Tau
    @endcode

The additional trigger information is encoded as follows. A more detailed
description can be found in the corresponding NanoAOD producer
[link](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L17)


Electrons                           | Value | Bit (value used in the config)
------------------------------------|-------|-------
ignore matching                     |  -    | -1
no match                            |  0    | -
CaloIdL_TrackIdL_IsoVL              |  1    | 0
1e (WPTight)                        |  2    | 1
1e (WPLoose)                        |  4    | 2
OverlapFilter PFTau                 |  8    | 3
2e                                  |  16   | 4
1e-1mu                              |  32   | 5
1e-1tau                             |  64   | 6
3e                                  |  128  | 7
2e-1mu                              |  256  | 8
1e-2mu                              |  512  | 9
1e (32_L1DoubleEG_AND_L1SingleEGOr) |  1024 | 10

Muons               | Value | Bit (value used in the config)
--------------------|-------|-------
ignore matching     |  -    | -1
no match            |  0    | -
TrkIsoVVL           |  1    | 0
Iso                 |  2    | 1
OverlapFilter PFTau |  4    | 2
1mu                 |  8    | 3
2mu                 |  16   | 4
1mu-1e              |  32   | 5
1mu-1tau            |  64   | 6
3mu                 |  128  | 7
2mu-1e              |  256  | 8
1mu-2e              |  512  | 9

Taus                 | Value | Bit (value used in the config)
---------------------|-------|-------
ignore matching      |  -    | -1
no match             |  0    | -
LooseChargedIso      |  1    | 0
MediumChargedIso     |  2    | 1
TightChargedIso      |  4    | 2
TightID OOSC photons |  8    | 3
HPS                  |  16   | 4
single-tau + tau+MET |  32   | 5
di-tau               |  64   | 6
e-tau                |  128  | 7
mu-tau               |  256  | 8
VBF+di-tau for Tau   |  512  | 9

jets                                    | Value | Bit (value used in the config)
----------------------------------------|-------|-------
no match                                |  0    | -
VBF cross-cleaned from loose iso PFTau  |  1    | 0
 * @param particle the `ROOT::Math::PtEtaPhiMVector` vector of the object to
match
 * @param triggerobject_pts `ROOT::RVec<float>` of trigger object pts
 * @param triggerobject_etas `ROOT::RVec<float>` of trigger object etas
 * @param triggerobject_phis `ROOT::RVec<float>` of trigger object phis
 * @param triggerobject_bits `ROOT::RVec<float>` of trigger object bits.
Depending on the trigger object id, this bitmap has a different meaning as
listed in th table above.
 * @param triggerobject_ids `ROOT::RVec<float>` of trigger object ids. The ID
encoding is explained above.
 * @param matchDeltaR The maximum deltaR value used for the match
 * @param pt_cut pt cut value on the trigger object pt
 * @param eta_cut eta cut value on the trigger object eta
 * @param trigger_particle_id_cut the triggerobject id required
 * @param triggerbit_cut the triggerobject filter bit position required
 * @return true, if all criteria are met, false otherwise
 */

bool matchParticle(const ROOT::Math::PtEtaPhiMVector &particle,
                   ROOT::RVec<float> &triggerobject_pts,
                   ROOT::RVec<float> &triggerobject_etas,
                   ROOT::RVec<float> &triggerobject_phis,
                   ROOT::RVec<int> &triggerobject_bits,
                   ROOT::RVec<int> &triggerobject_ids, const float &matchDeltaR,
                   const float &pt_cut, const float &eta_cut,
                   const int &trigger_particle_id_cut,
                   const int &triggerbit_cut) {
    Logger::get("CheckTriggerMatch")->warn("Checking Triggerobjects");
    Logger::get("CheckTriggerMatch")
        ->warn("Total number of triggerobjects: {}", triggerobject_pts.size());
    for (std::size_t idx = 0; idx < triggerobject_pts.size(); ++idx) {
        Logger::get("CheckTriggerMatch")->warn("Triggerobject Nr. {}", idx);
        Logger::get("CheckTriggerMatch")
            ->warn("bit Value: {}", IntBits(triggerobject_bits[idx]));
        Logger::get("CheckTriggerMatch")
            ->warn("bit Value: {}", triggerobject_bits[idx]);
        auto triggerobject = ROOT::Math::RhoEtaPhiVectorF(
            0, triggerobject_etas[idx], triggerobject_phis[idx]);
        // We check the deltaR match as well as that the pt and eta of the
        // triggerobject are above the given thresholds
        bool deltaR = ROOT::Math::VectorUtil::DeltaR(triggerobject, particle) <
                      matchDeltaR;
        // if we don't want to do any matching here, the triggerbut_cut value is
        // -1
        Logger::get("CheckTriggerMatch")
            ->warn("bit Value: {}", triggerobject_bits[idx]);
        bool bit = (triggerbit_cut == -1) ||
                   (IntBits(triggerobject_bits[idx]).test(triggerbit_cut));
        bool id = triggerobject_ids[idx] == trigger_particle_id_cut;
        bool pt = particle.pt() > pt_cut;
        bool eta = abs(particle.eta()) < eta_cut;
        Logger::get("CheckTriggerMatch")
            ->warn("-------------------------------------------------------");
        Logger::get("CheckTriggerMatch")->warn("deltaR Check: {}", deltaR);
        Logger::get("CheckTriggerMatch")
            ->warn("deltaR Value: {}",
                    ROOT::Math::VectorUtil::DeltaR(triggerobject, particle));
        Logger::get("CheckTriggerMatch")->warn("id Check: {}", id);
        Logger::get("CheckTriggerMatch")
            ->warn("id Value: {}", triggerobject_ids[idx]);
        Logger::get("CheckTriggerMatch")->warn("bit Check: {}", bit);
        Logger::get("CheckTriggerMatch")
            ->warn("bit Value: {}", IntBits(triggerobject_bits[idx]));
        Logger::get("CheckTriggerMatch")->warn("pt Check: {}", pt);
        Logger::get("CheckTriggerMatch")
            ->warn("pt Value: {}", triggerobject_pts[idx]);
        Logger::get("CheckTriggerMatch")->warn("eta Check: {}", eta);
        Logger::get("CheckTriggerMatch")
            ->warn("eta Value: {}", triggerobject_etas[idx]);
        Logger::get("CheckTriggerMatch")
            ->warn("-------------------------------------------------------");
        if (deltaR && bit && id && pt && eta) {
            // remove the matching object from the object vectors so it cant be
            // matched by the next particle as well (if there is one)
            triggerobject_ids.erase(triggerobject_ids.begin() + idx);
            triggerobject_bits.erase(triggerobject_bits.begin() + idx);
            triggerobject_pts.erase(triggerobject_pts.begin() + idx);
            triggerobject_etas.erase(triggerobject_etas.begin() + idx);
            triggerobject_phis.erase(triggerobject_phis.begin() + idx);
            return true;
        }
    }
    return false;
};
/**
 * @brief Function to generate a trigger flag based on an hlt path and trigger
 * object matching for the given object. This relies on the
 * trigger::matchParticle function which does the matching test.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle_p4 `ROOT::Math::PtEtaPhiMVector` of the object to be checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param hltpath name of the hlt path to be checked, this can be a valid regex.
 * If more than one matching HLT path is found, the function will throw an
 * exception, if no matching HLT path is found, the function will return a
 * dataframe with a flag of false for all entries
 * @param pt_cut minimal pt value for the triggerobject
 * @param eta_cut maximal pt value for the triggerobject
 * @param trigger_particle_id_cut trigger id value the triggerobject has to
 * match (details can be found in the documentation of trigger::matchParticle)
 * @param triggerbit_cut trigger bit value the triggerobject has to match
 * (details can be found in the documentation of trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger flag column
 */

ROOT::RDF::RNode GenerateSingleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &hltpath, const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold) {

    auto triggermatch =
        [DeltaR_threshold, pt_cut, eta_cut, trigger_particle_id_cut,
         triggerbit_cut](bool hltpath,
                         const ROOT::Math::PtEtaPhiMVector &particle_p4,
                         ROOT::RVec<int> triggerobject_bits,
                         ROOT::RVec<int> triggerobject_ids,
                         ROOT::RVec<float> triggerobject_pts,
                         ROOT::RVec<float> triggerobject_etas,
                         ROOT::RVec<float> triggerobject_phis) {
            Logger::get("GenerateSingleTriggerFlag")->warn("Checking Trigger");
            bool result = false;
            bool match_result = false;
            if (hltpath) {
                Logger::get("CheckTriggerMatch")
                    ->warn("Checking Triggerobject match with particles ....");
                match_result = matchParticle(
                    particle_p4, triggerobject_pts, triggerobject_etas,
                    triggerobject_phis, triggerobject_bits, triggerobject_ids,
                    DeltaR_threshold, pt_cut, eta_cut, trigger_particle_id_cut,
                    triggerbit_cut);
            }
            result = hltpath & match_result;
            Logger::get("GenerateSingleTriggerFlag")
                ->warn("---> HLT Match: {}", hltpath);
            Logger::get("GenerateSingleTriggerFlag")
                ->warn("---> Total Match: {}", match_result);
            Logger::get("GenerateSingleTriggerFlag")
                ->warn("--->>>> result: {}", result);
            return result;
        };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hltpath_regex = std::regex(hltpath);
    // loop over all available trigger names and check if the hltpath is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hltpath_regex)) {
            Logger::get("GenerateSingleTriggerFlag")
                ->warn("Found matching trigger: {}", trigger);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("GenerateSingleTriggerFlag")
            ->info("No matching trigger for {} found, returning false for "
                   "trigger flag {}",
                   hltpath, triggerflag_name);
        auto df1 = df.Define(triggerflag_name, []() { return false; });
        return df1;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("GenerateSingleTriggerFlag")
            ->warn("More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        Logger::get("GenerateSingleTriggerFlag")
            ->warn("Found matching trigger: {}", matched_trigger_names[0]);
        auto df1 =
            df.Define(triggerflag_name, triggermatch,
                      {matched_trigger_names[0], particle_p4,
                       triggerobject_bits, triggerobject_id, triggerobject_pt,
                       triggerobject_eta, triggerobject_phi});
        return df1;
    }
}

/**
 * @brief Function to generate a trigger flag based on an hlt path and trigger
 * object matching for the given object. This relies on the
 * trigger::matchParticle function which does the matching test. The
 * implementation is similar to the trigger::GenerateSingleTriggerFlag function,
 * but here two objects have to be matched for a passing flag.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle1_p4 `ROOT::Math::PtEtaPhiMVector` of the first object to be
 * checked
 * @param particle2_p4 `ROOT::Math::PtEtaPhiMVector` of the second object to be
 * checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param hltpath name of the hlt path to be checked, this can be a valid regex.
 * If more than one matching HLT path is found, the function will throw an
 * exception, if no matching HLT path is found, the function will return a
 * dataframe with a flag of false for all entries
 * @param p1_pt_cut minimal pt value for the triggerobject matching the first
 * object
 * @param p2_pt_cut minimal pt value for the triggerobject matching the second
 * object
 * @param p1_eta_cut maximal pt value for the triggerobject matching the first
 * object
 * @param p2_eta_cut maximal pt value for the triggerobject matching the second
 * object
 * @param p1_trigger_particle_id_cut trigger id value the triggerobject matching
 * the first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_trigger_particle_id_cut trigger id value the triggerobject matching
 * the second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p1_triggerbit_cut trigger bit value the triggerobject matching the
 * first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_triggerbit_cut trigger bit value the triggerobject matching the
 * second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger flag column
 *
 */

ROOT::RDF::RNode GenerateDoubleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &hltpath,
    const float &p1_pt_cut, const float &p2_pt_cut, const float &p1_eta_cut,
    const float &p2_eta_cut, const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold) {

    auto triggermatch = [DeltaR_threshold, p1_pt_cut, p2_pt_cut, p1_eta_cut,
                         p2_eta_cut, p1_trigger_particle_id_cut,
                         p2_trigger_particle_id_cut, p1_triggerbit_cut,
                         p2_triggerbit_cut](
                            bool hltpath,
                            const ROOT::Math::PtEtaPhiMVector &particle1_p4,
                            const ROOT::Math::PtEtaPhiMVector &particle2_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        Logger::get("GenerateDoubleTriggerFlag")->warn("Checking Trigger");
        bool result = false;
        bool match_result_p1 = false;
        bool match_result_p2 = false;
        if (hltpath) {
            Logger::get("GenerateDoubleTriggerFlag")
                ->warn("Checking Triggerobject match with particles ....");
            Logger::get("GenerateDoubleTriggerFlag")->warn("First particle");
            match_result_p1 = matchParticle(
                particle1_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p1_pt_cut, p1_eta_cut,
                p1_trigger_particle_id_cut, p1_triggerbit_cut);
            Logger::get("GenerateDoubleTriggerFlag")->warn("Second particle");
            match_result_p2 = matchParticle(
                particle2_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p2_pt_cut, p2_eta_cut,
                p2_trigger_particle_id_cut, p2_triggerbit_cut);
        }
        result = hltpath & match_result_p1 & match_result_p2;
        Logger::get("GenerateDoubleTriggerFlag")
            ->warn("---> HLT Match: {}", hltpath);
        Logger::get("GenerateDoubleTriggerFlag")
            ->warn("---> Total Match P1: {}", match_result_p1);
        Logger::get("GenerateDoubleTriggerFlag")
            ->warn("---> Total Match P2: {}", match_result_p2);
        Logger::get("GenerateDoubleTriggerFlag")
            ->warn("--->>>> result: {}", result);
        return result;
    };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hltpath_regex = std::regex(hltpath);
    // loop over all available trigger names and check if the hltpath is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hltpath_regex)) {
            Logger::get("GenerateDoubleTriggerFlag")
                ->warn("Found matching trigger: {}", trigger);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("GenerateDoubleTriggerFlag")
            ->info("No matching trigger for {} found, returning false for "
                   "trigger flag {}",
                   hltpath, triggerflag_name);
        auto df1 = df.Define(triggerflag_name, []() { return false; });
        return df1;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("GenerateDoubleTriggerFlag")
            ->warn("More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        Logger::get("GenerateDoubleTriggerFlag")
            ->warn("Found matching trigger: {}", matched_trigger_names[0]);
        auto df1 =
            df.Define(triggerflag_name, triggermatch,
                      {matched_trigger_names[0], particle1_p4, particle2_p4,
                       triggerobject_bits, triggerobject_id, triggerobject_pt,
                       triggerobject_eta, triggerobject_phi});
        return df1;
    }
}

/**
 * @brief Function to generate a trigger flag based on a trigger
 * object matching for two given objects. This relies on the
 * trigger::matchParticle function which does the matching test. The
 * implementation is similar to the trigger::GenerateDoubleTriggerFlag function,
 * but here, no hlt path is required, only the matching of trigger objects is
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle1_p4 `ROOT::Math::PtEtaPhiMVector` of the first object to be
 * checked
 * @param particle2_p4 `ROOT::Math::PtEtaPhiMVector` of the second object to be
 * checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param p1_pt_cut minimal pt value for the triggerobject matching the first
 * object
 * @param p2_pt_cut minimal pt value for the triggerobject matching the second
 * object
 * @param p1_eta_cut maximal pt value for the triggerobject matching the first
 * object
 * @param p2_eta_cut maximal pt value for the triggerobject matching the second
 * object
 * @param p1_trigger_particle_id_cut trigger id value the triggerobject matching
 * the first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_trigger_particle_id_cut trigger id value the triggerobject matching
 * the second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p1_triggerbit_cut trigger bit value the triggerobject matching the
 * first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_triggerbit_cut trigger bit value the triggerobject matching the
 * second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger flag column
 *
 */


ROOT::RDF::RNode MatchDoubleTriggerObject(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi,
    const float &p1_pt_cut, const float &p2_pt_cut, const float &p1_eta_cut,
    const float &p2_eta_cut, const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold) {

    auto triggermatch = [DeltaR_threshold, p1_pt_cut, p2_pt_cut, p1_eta_cut,
                         p2_eta_cut, p1_trigger_particle_id_cut,
                         p2_trigger_particle_id_cut, p1_triggerbit_cut,
                         p2_triggerbit_cut](
                            const ROOT::Math::PtEtaPhiMVector &particle1_p4,
                            const ROOT::Math::PtEtaPhiMVector &particle2_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        bool match_result_p1 = false;
        bool match_result_p2 = false;
        Logger::get("MatchDoubleTriggerObject")
            ->warn("Checking Triggerobject match with particles ....");
        Logger::get("MatchDoubleTriggerObject")->warn("First particle");
        match_result_p1 = matchParticle(
            particle1_p4, triggerobject_pts, triggerobject_etas,
            triggerobject_phis, triggerobject_bits, triggerobject_ids,
            DeltaR_threshold, p1_pt_cut, p1_eta_cut,
            p1_trigger_particle_id_cut, p1_triggerbit_cut);
        Logger::get("MatchDoubleTriggerObject")->warn("Second particle");
        match_result_p2 = matchParticle(
            particle2_p4, triggerobject_pts, triggerobject_etas,
            triggerobject_phis, triggerobject_bits, triggerobject_ids,
            DeltaR_threshold, p2_pt_cut, p2_eta_cut,
            p2_trigger_particle_id_cut, p2_triggerbit_cut);
        bool result = match_result_p1 & match_result_p2;
        Logger::get("MatchDoubleTriggerObject")
            ->warn("---> Total Match P1: {}", match_result_p1);
        Logger::get("MatchDoubleTriggerObject")
            ->warn("---> Total Match P2: {}", match_result_p2);
        Logger::get("MatchDoubleTriggerObject")
            ->warn("--->>>> result: {}", result);
        return result;
    };
    auto df1 =
        df.Define(triggerflag_name, triggermatch,
                  {particle1_p4, particle2_p4, triggerobject_bits, triggerobject_id,
                   triggerobject_pt, triggerobject_eta, triggerobject_phi});
    return df1;
}

/**
 * @brief Function to generate a trigger flag based on a trigger
 * object match only. This relies on the
 * trigger::matchParticle function which does the matching test.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle_p4 `ROOT::Math::PtEtaPhiMVector` of the object to be checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param pt_cut minimal pt value for the triggerobject
 * @param eta_cut maximal pt value for the triggerobject
 * @param trigger_particle_id_cut trigger id value the triggerobject has to
 * match (details can be found in the documentation of trigger::matchParticle)
 * @param triggerbit_cut trigger bit value the triggerobject has to match
 * (details can be found in the documentation of trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger object match flag column
 */

ROOT::RDF::RNode MatchSingleTriggerObject(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold) {

    auto triggermatch = [DeltaR_threshold, pt_cut, eta_cut,
                         trigger_particle_id_cut, triggerbit_cut](
                            const ROOT::Math::PtEtaPhiMVector &particle_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        Logger::get("MatchSingleTriggerObject")->warn("Checking Trigger");
        Logger::get("MatchSingleTriggerObject")
            ->warn("Checking Triggerobject match with particles ....");
        bool match_result =
            matchParticle(particle_p4, triggerobject_pts, triggerobject_etas,
                          triggerobject_phis, triggerobject_bits,
                          triggerobject_ids, DeltaR_threshold, pt_cut, eta_cut,
                          trigger_particle_id_cut, triggerbit_cut);
        Logger::get("MatchSingleTriggerObject")
            ->warn("--->>>> match_result: {}", match_result);
        return match_result;
    };
    auto df1 =
        df.Define(triggerflag_name, triggermatch,
                  {particle_p4, triggerobject_bits, triggerobject_id,
                   triggerobject_pt, triggerobject_eta, triggerobject_phi});
    return df1;
}


} // end namespace trigger
#endif /* GUARD_TRIGGERS_H */