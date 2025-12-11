#ifndef GUARD_TRIGGERS_H
#define GUARD_TRIGGERS_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "bitset"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <regex>

typedef std::bitset<30> IntBits;

namespace trigger {

/**
 * @brief This function tries to match an object with a trigger object. An
 * object is successfully matched, if they overlap within a given \f$\Delta R\f$
 * cone and the \f$p_T\f$ / \f$|\eta|\f$ thresholds are met. Further, the trigger
 * object ID and the trigger object filter bit are set to the correct value.
 *
 * For the trigger objects, two additional quantities are stored, the ID of the
 * particle and some further information on the trigger path encoded in a bitmap.
 * The trigger object ID is encoded as follows:
 * @code
 *   1 = Jet
 *   2 = MET
 *   3 = HT
 *   4 = MHT
 *   6 = FatJet
 *   11 = Electron (PixelMatched e/gamma)
 *   22 = Photon (PixelMatch-vetoed e/gamma)
 *   13 = Muon
 *   15 = Tau
 * @endcode
 *
 * A more detailed description about the encoded additional trigger information
 * can be found in the corresponding NanoAOD producer
 * [triggerObjects_cff.py](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py)
 * It is recommended to look in the CMSSW branch that was used to produce the
 * NanoAOD version you are working with and not the master branch.
 * As an example, for Run2 (NanoAODv9) for electrons:
 *
 * Electrons                           | Value | Bit (value used in the config)
 * ------------------------------------|-------|-------
 * ignore matching                     |  -    | -1
 * no match                            |  0    | -
 * CaloIdL_TrackIdL_IsoVL              |  1    | 0
 * 1e (WPTight)                        |  2    | 1
 * 1e (WPLoose)                        |  4    | 2
 * OverlapFilter PFTau                 |  8    | 3
 * 2e                                  |  16   | 4
 * 1e-1mu                              |  32   | 5
 * 1e-1tau                             |  64   | 6
 * 3e                                  |  128  | 7
 * 2e-1mu                              |  256  | 8
 * 1e-2mu                              |  512  | 9
 * 1e (32_L1DoubleEG_AND_L1SingleEGOr) |  1024 | 10
 *
 * @param particle Lorentz vector of the object/particle that should be
 * matched to a trigger object
 * @param triggerobject_pts vector of trigger object \f$p_T\f$ values
 * @param triggerobject_etas vector of trigger object \f$\eta\f$ values
 * @param triggerobject_phis vector of trigger object \f$\phi\f$ values
 * @param triggerobject_ids vector of trigger object IDs
 * @param triggerobject_filterbits vector of trigger object filter bits.
 * Depending on the trigger object ID (e.g. electron, muon, ...), this
 * bitmap has a different meaning.
 * @param pt_threshold \f$p_T\f$ threshold value the tested object has to
 * exceed in order to be considered a match
 * @param eta_threshold \f$\eta\f$ threshold value the tested object has to
 * be below in order to be considered a match
 * @param trigger_particle_id_value trigger object ID value that should be
 * tested
 * @param trigger_bit_value trigger object filter bit position that should be
 * tested. If no bit matching is desired, set this value to -1.
 * @param deltaR_threshold maximum \f$\Delta R\f$ value between the trigger
 * object and the tested object in order to be considered a match
 *
 * @return boolean value indicating whether a match was found
 */
bool matchParticle(const ROOT::Math::PtEtaPhiMVector &particle,
                   ROOT::RVec<float> &triggerobject_pts,
                   ROOT::RVec<float> &triggerobject_etas,
                   ROOT::RVec<float> &triggerobject_phis,
                   ROOT::RVec<int> &triggerobject_ids,
                   ROOT::RVec<int> &triggerobject_filterbits,
                   const float &pt_threshold, const float &eta_threshold,
                   const int &trigger_particle_id_value, const int &trigger_bit_value,
                   const float &deltaR_threshold) {
    Logger::get("trigger::matchParticle")->debug("Checking Triggerobjects");
    Logger::get("trigger::matchParticle")
        ->debug("Total number of triggerobjects: {}", triggerobject_pts.size());
    for (std::size_t idx = 0; idx < triggerobject_pts.size(); ++idx) {
        Logger::get("trigger::matchParticle")->debug("Triggerobject Nr. {}", idx);
        auto triggerobject = ROOT::Math::RhoEtaPhiVectorF(
            0, triggerobject_etas[idx], triggerobject_phis[idx]);
        // We check the deltaR match as well as that the pt and eta of the
        // triggerobject are above the given thresholds
        bool deltaR = ROOT::Math::VectorUtil::DeltaR(triggerobject, particle) <
                      deltaR_threshold;
        bool pt = particle.pt() > pt_threshold;
        bool eta = abs(particle.eta()) < eta_threshold;
        // if you don't want to do bit matching here, the trigger_bit_value
        // has to be set to -1
        bool bit = (trigger_bit_value == -1) ||
                   (IntBits(triggerobject_filterbits[idx]).test(trigger_bit_value));
        bool id = triggerobject_ids[idx] == trigger_particle_id_value;

        Logger::get("trigger::matchParticle")
            ->debug("-------------------------------------------------------");

        Logger::get("trigger::matchParticle")
            ->debug("deltaR_threshold: {}, Check: {}", deltaR_threshold, deltaR);
        Logger::get("trigger::matchParticle")
            ->debug("deltaR value: {}",
                ROOT::Math::VectorUtil::DeltaR(triggerobject, particle));

        Logger::get("trigger::matchParticle")
            ->debug("trigger_particle_id_value: {}, Check: {}",
                trigger_particle_id_value, id);
        Logger::get("trigger::matchParticle")
            ->debug("id value: {}", triggerobject_ids[idx]);

        Logger::get("trigger::matchParticle")
            ->debug("trigger_bit_value: {}, Check: {}", trigger_bit_value, bit);
        Logger::get("trigger::matchParticle")
            ->debug("bit value: {}", IntBits(triggerobject_filterbits[idx]));

        Logger::get("trigger::matchParticle")
            ->debug("pt_threshold: {}, Check: {}", pt_threshold, pt);
        Logger::get("trigger::matchParticle")
            ->debug("pt value (trg): {}, pt value (reco): {}",
                triggerobject_pts[idx], particle.pt());

        Logger::get("trigger::matchParticle")
            ->debug("eta_threshold: {}, Check: {}", eta_threshold, eta);
        Logger::get("trigger::matchParticle")
            ->debug("eta (trg) value: {}, eta (reco) value: {}",
                triggerobject_etas[idx], abs(particle.eta()));

        Logger::get("trigger::matchParticle")
            ->debug("-------------------------------------------------------");
        if (deltaR && bit && id && pt && eta) {
            // remove the matching object from the object vectors so it cannot be
            // matched by the next particle as well (if there is one)
            triggerobject_pts.erase(triggerobject_pts.begin() + idx);
            triggerobject_etas.erase(triggerobject_etas.begin() + idx);
            triggerobject_phis.erase(triggerobject_phis.begin() + idx);
            triggerobject_ids.erase(triggerobject_ids.begin() + idx);
            triggerobject_filterbits.erase(triggerobject_filterbits.begin() + idx);
            return true;
        }
    }
    return false;
};
/**
 * @brief This function generates a trigger flag based on an HLT path and
 * trigger object matching for a selected object. This relies on the
 * `trigger::matchParticle` function which does the object to trigger
 * object matching test.
 *
 * @note This function is defined for single object triggers only. For double
 * object triggers be referred to the `trigger::DoubleObjectFlag` functions.
 *
 * @note If more than one matching HLT path is found, the function will
 * throw an exception, if no matching HLT path is found, the function will
 * return a dataframe with a flag with false for all entries.
 *
 * @param df input dataframe
 * @param outputname name of the output flag
 * @param particle name of the column containing the Lorentz vector of the
 * object/particle that should be matched to a trigger object
 * @param triggerobject_pt name of the column containing the trigger object
 * \f$p_T\f$ values
 * @param triggerobject_eta name of the column containing the trigger object
 * \f$\eta\f$ values
 * @param triggerobject_phi name of the column containing the trigger object
 * \f$\phi\f$ values
 * @param triggerobject_id name of the column containing the trigger object
 * IDs
 * @param triggerobject_filterbit name of the column containing the trigger
 * object filter bits. Depending on the trigger object ID (e.g. electron,
 * muon, ...), this bitmap has a different meaning.
 * @param hlt_path name of the column containing the HLT path to be checked,
 * this can be a valid regex.
 * @param pt_threshold \f$p_T\f$ threshold value the tested object has to
 * exceed in order to be considered a match
 * @param eta_threshold \f$\eta\f$ threshold value the tested object has to
 * be below in order to be considered a match
 * @param trigger_particle_id_value trigger object ID value that should be
 * tested
 * @param trigger_bit_value trigger object filter bit position that should be
 * tested. If no bit matching is desired, set this value to -1.
 * @param deltaR_threshold maximum \f$\Delta R\f$ value between the trigger
 * object and the tested object in order to be considered a match
 *
 * @return a new dataframe containing the trigger flag column
 */
ROOT::RDF::RNode SingleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &triggerobject_id, const std::string &triggerobject_filterbit,
    const std::string &hlt_path, const float &pt_threshold, const float &eta_threshold,
    const int &trigger_particle_id_value, const int &trigger_bit_value,
    const float &deltaR_threshold) {
    // In nanoAODv12 the type of trigger object ID was changed to UShort_t
    // For v9 compatibility a type casting is applied
    auto [df1, triggerobject_id_column] = utility::Cast<ROOT::RVec<UShort_t>, ROOT::RVec<Int_t>>(
            df, triggerobject_id+"_v12", "ROOT::VecOps::RVec<UShort_t>", triggerobject_id);
    // In nanoAODv15 the type of trigger object ID was changed to ULong64_t
    // For v9 and v12 compatibility a type casting is applied
    auto [df2, triggerobject_filterbit_column] = utility::Cast<ROOT::RVec<ULong64_t>, ROOT::RVec<Int_t>>(
            df1, triggerobject_filterbit+"_v15", "ROOT::VecOps::RVec<ULong64_t>", triggerobject_filterbit);

    auto trigger_matching = [pt_threshold, eta_threshold, trigger_particle_id_value,
                             trigger_bit_value, deltaR_threshold](
                                bool hlt_path_match,
                                const ROOT::Math::PtEtaPhiMVector &particle,
                                ROOT::RVec<float> triggerobject_pts,
                                ROOT::RVec<float> triggerobject_etas,
                                ROOT::RVec<float> triggerobject_phis,
                                ROOT::RVec<UShort_t> triggerobject_ids_v12,
                                ROOT::RVec<ULong64_t> triggerobject_filterbits_v15) {
        auto triggerobject_ids = static_cast<ROOT::RVec<int>>(triggerobject_ids_v12);
        auto triggerobject_filterbits = static_cast<ROOT::RVec<int>>(triggerobject_filterbits_v15);

        bool match_result = false;
        if (hlt_path_match) {
            Logger::get("trigger::SingleObjectFlag")
                ->debug("Checking triggerobject match with particle ....");
            match_result = matchParticle(
                particle, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_ids, triggerobject_filterbits,
                pt_threshold, eta_threshold, trigger_particle_id_value,
                trigger_bit_value, deltaR_threshold);
        }
        bool result = hlt_path_match & match_result;
        Logger::get("trigger::SingleObjectFlag")
            ->debug("---> HLT Matching: {}", hlt_path_match);
        Logger::get("trigger::SingleObjectFlag")
            ->debug("---> Object Matching: {}", match_result);
        Logger::get("trigger::SingleObjectFlag")
            ->debug("--->>>> result: {}", result);
        return match_result;
    };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hlt_path_regex = std::regex(hlt_path);
    // loop over all available trigger names and check if the HLT path is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hlt_path_regex)) {
            Logger::get("trigger::SingleObjectFlag")
                ->debug("Found matching trigger: {}", trigger);
            Logger::get("trigger::SingleObjectFlag")
                ->debug("For HLT path: {}", hlt_path);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("trigger::SingleObjectFlag")
            ->debug("No matching trigger for {} found, returning false for "
                   "trigger flag {}", hlt_path, outputname);
        auto df3 = df2.Define(outputname, []() { return false; });
        return df3;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("trigger::SingleObjectFlag")
            ->debug("More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        auto df3 =
            df2.Define(outputname, trigger_matching,
                      {matched_trigger_names[0], particle, triggerobject_pt,
                       triggerobject_eta, triggerobject_phi,
                       triggerobject_id_column, triggerobject_filterbit_column});
        return df3;
    }
}

/**
 * @brief This function generates a trigger flag based on trigger object
 * matching for a selected object. This relies on the `trigger::matchParticle`
 * function which does the object to trigger object matching test.
 *
 * @note This function is defined for single object triggers only. For double
 * object triggers be referred to the `trigger::DoubleObjectFlag` functions.
 *
 * @param df input dataframe
 * @param outputname name of the output flag
 * @param particle name of the column containing the Lorentz vector of the
 * object/particle that should be matched to a trigger object
 * @param triggerobject_pt name of the column containing the trigger object
 * \f$p_T\f$ values
 * @param triggerobject_eta name of the column containing the trigger object
 * \f$\eta\f$ values
 * @param triggerobject_phi name of the column containing the trigger object
 * \f$\phi\f$ values
 * @param triggerobject_id name of the column containing the trigger object
 * IDs
 * @param triggerobject_filterbit name of the column containing the trigger
 * object filter bits. Depending on the trigger object ID (e.g. electron,
 * muon, ...), this bitmap has a different meaning.
 * @param pt_threshold \f$p_T\f$ threshold value the tested object has to
 * exceed in order to be considered a match
 * @param eta_threshold \f$\eta\f$ threshold value the tested object has to
 * be below in order to be considered a match
 * @param trigger_particle_id_value trigger object ID value that should be
 * tested
 * @param trigger_bit_value trigger object filter bit position that should be
 * tested. If no bit matching is desired, set this value to -1.
 * @param deltaR_threshold maximum \f$\Delta R\f$ value between the trigger
 * object and the tested object in order to be considered a match
 *
 * @return a new dataframe containing the trigger flag column
 */
ROOT::RDF::RNode SingleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &triggerobject_id, const std::string &triggerobject_filterbit,
    const float &pt_threshold, const float &eta_threshold,
    const int &trigger_particle_id_value, const int &trigger_bit_value,
    const float &deltaR_threshold) {
    // In nanoAODv12 the type of trigger object ID was changed to UShort_t
    // For v9 compatibility a type casting is applied
    auto [df1, triggerobject_id_column] = utility::Cast<ROOT::RVec<UShort_t>, ROOT::RVec<Int_t>>(
            df, triggerobject_id+"_v12", "ROOT::VecOps::RVec<UShort_t>", triggerobject_id);
    // In nanoAODv15 the type of trigger object ID was changed to ULong64_t
    // For v9 and v12 compatibility a type casting is applied
    auto [df2, triggerobject_filterbit_column] = utility::Cast<ROOT::RVec<ULong64_t>, ROOT::RVec<Int_t>>(
            df1, triggerobject_filterbit+"_v15", "ROOT::VecOps::RVec<ULong64_t>", triggerobject_filterbit);

    auto trigger_matching = [pt_threshold, eta_threshold, trigger_particle_id_value,
                             trigger_bit_value, deltaR_threshold](
                                const ROOT::Math::PtEtaPhiMVector &particle,
                                ROOT::RVec<float> triggerobject_pts,
                                ROOT::RVec<float> triggerobject_etas,
                                ROOT::RVec<float> triggerobject_phis,
                                ROOT::RVec<UShort_t> triggerobject_ids_v12,
                                ROOT::RVec<ULong64_t> triggerobject_filterbits_v15) {
        auto triggerobject_ids = static_cast<ROOT::RVec<int>>(triggerobject_ids_v12);
        auto triggerobject_filterbits = static_cast<ROOT::RVec<int>>(triggerobject_filterbits_v15);
        Logger::get("trigger::SingleObjectFlag")
            ->debug("Checking triggerobject match with particle ....");
        bool match_result = matchParticle(
            particle, triggerobject_pts, triggerobject_etas,
            triggerobject_phis, triggerobject_ids, triggerobject_filterbits,
            pt_threshold, eta_threshold, trigger_particle_id_value,
            trigger_bit_value, deltaR_threshold);
        Logger::get("trigger::SingleObjectFlag")
            ->debug("---> Matching: {}", match_result);
        return match_result;
    };

    auto df3 = df2.Define(outputname, trigger_matching,
                   {particle, triggerobject_pt, triggerobject_eta, triggerobject_phi,
                    triggerobject_id_column, triggerobject_filterbit_column});
    return df3;
}

/**
 * @brief This function generates a trigger flag based on an HLT path and
 * trigger object matching for the selected objects. This relies on the
 * `trigger::matchParticle` function which does the object to trigger
 * object matching test.
 *
 * @note This function is defined for double object triggers only. For single
 * object triggers be referred to the `trigger::SingleObjectFlag` functions.
 *
 * @note If more than one matching HLT path is found, the function will
 * throw an exception, if no matching HLT path is found, the function will
 * return a dataframe with a flag with false for all entries.
 *
 * @param df input dataframe
 * @param outputname name of the output flag
 * @param particle_1 name of the column containing the Lorentz vector of the
 * first object/particle that should be matched to a trigger object
 * @param particle_2 name of the column containing the Lorentz vector of the
 * second object/particle that should be matched to a trigger object
 * @param triggerobject_pt name of the column containing the trigger object
 * \f$p_T\f$ values
 * @param triggerobject_eta name of the column containing the trigger object
 * \f$\eta\f$ values
 * @param triggerobject_phi name of the column containing the trigger object
 * \f$\phi\f$ values
 * @param triggerobject_id name of the column containing the trigger object
 * IDs
 * @param triggerobject_filterbit name of the column containing the trigger
 * object filter bits. Depending on the trigger object ID (e.g. electron,
 * muon, ...), this bitmap has a different meaning.
 * @param hlt_path name of the column containing the HLT path to be checked,
 * this can be a valid regex.
 * @param pt_threshold_1 \f$p_T\f$ threshold value the first tested object
 * has to exceed in order to be considered a match
 * @param pt_threshold_2 \f$p_T\f$ threshold value the second tested object
 * has to exceed in order to be considered a match
 * @param eta_threshold_1 \f$\eta\f$ threshold value the first tested object
 * has to be below in order to be considered a match
 * @param eta_threshold_2 \f$\eta\f$ threshold value the second tested object
 * has to be below in order to be considered a match
 * @param trigger_particle_id_value_1 trigger object ID value that should be
 * tested for the first object
 * @param trigger_particle_id_value_2 trigger object ID value that should be
 * tested for the second object
 * @param trigger_bit_value_1 trigger object filter bit position that should be
 * tested for the first object. If no bit matching is desired, set this value
 * to -1.
 * @param trigger_bit_value_2 trigger object filter bit position that should be
 * tested for the second object. If no bit matching is desired, set this value
 * to -1.
 * @param deltaR_threshold maximum \f$\Delta R\f$ value between the trigger
 * object and the tested object in order to be considered a match
 *
 * @return a new dataframe containing the trigger flag column
 */
ROOT::RDF::RNode DoubleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle_1, const std::string &particle_2,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit, const std::string &hlt_path,
    const float &pt_threshold_1, const float &pt_threshold_2,
    const float &eta_threshold_1, const float &eta_threshold_2,
    const int &trigger_particle_id_value_1, const int &trigger_particle_id_value_2,
    const int &trigger_bit_value_1, const int &trigger_bit_value_2,
    const float &deltaR_threshold) {
    // In nanoAODv12 the type of trigger object ID was changed to UShort_t
    // For v9 compatibility a type casting is applied
    auto [df1, triggerobject_id_column] = utility::Cast<ROOT::RVec<UShort_t>, ROOT::RVec<Int_t>>(
            df, triggerobject_id+"_v12", "ROOT::VecOps::RVec<UShort_t>", triggerobject_id);
    // In nanoAODv15 the type of trigger object ID was changed to ULong64_t
    // For v9 and v12 compatibility a type casting is applied
    auto [df2, triggerobject_filterbit_column] = utility::Cast<ROOT::RVec<ULong64_t>, ROOT::RVec<Int_t>>(
            df1, triggerobject_filterbit+"_v15", "ROOT::VecOps::RVec<ULong64_t>", triggerobject_filterbit);

    auto trigger_matching = [pt_threshold_1, pt_threshold_2, eta_threshold_1,
                             eta_threshold_2, trigger_particle_id_value_1,
                             trigger_particle_id_value_2, trigger_bit_value_1,
                             trigger_bit_value_2, deltaR_threshold](
                                bool hlt_path_match,
                                const ROOT::Math::PtEtaPhiMVector &particle_1,
                                const ROOT::Math::PtEtaPhiMVector &particle_2,
                                ROOT::RVec<float> triggerobject_pts,
                                ROOT::RVec<float> triggerobject_etas,
                                ROOT::RVec<float> triggerobject_phis,
                                ROOT::RVec<UShort_t> triggerobject_ids_v12,
                                ROOT::RVec<ULong64_t> triggerobject_filterbits_v15) {
        auto triggerobject_ids = static_cast<ROOT::RVec<int>>(triggerobject_ids_v12);
        auto triggerobject_filterbits = static_cast<ROOT::RVec<int>>(triggerobject_filterbits_v15);

        bool match_result_p1 = false;
        bool match_result_p2 = false;
        if (hlt_path_match) {
            Logger::get("trigger::DoubleObjectFlag")
                ->debug("Checking triggerobject match with particle ....");
            Logger::get("trigger::DoubleObjectFlag")->debug("First particle");
            match_result_p1 = matchParticle(
                particle_1, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_ids, triggerobject_filterbits,
                pt_threshold_1, eta_threshold_1, trigger_particle_id_value_1,
                trigger_bit_value_1, deltaR_threshold);
            Logger::get("trigger::DoubleObjectFlag")->debug("Second particle");
            match_result_p2 = matchParticle(
                particle_2, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_ids, triggerobject_filterbits,
                pt_threshold_2, eta_threshold_2, trigger_particle_id_value_2,
                trigger_bit_value_2, deltaR_threshold);
        }

        bool result = hlt_path_match & match_result_p1 & match_result_p2;
        Logger::get("trigger::SingleObjectFlag")
            ->debug("---> HLT Matching: {}", hlt_path_match);
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("---> First Object Matching: {}", match_result_p1);
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("---> Second Object Matching: {}", match_result_p2);
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("--->>>> result: {}", result);
        return result;
    };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hlt_path_regex = std::regex(hlt_path);
    // loop over all available trigger names and check if the HLT path is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hlt_path_regex)) {
            Logger::get("trigger::DoubleObjectFlag")
                ->debug("Found matching trigger: {}", trigger);
            Logger::get("trigger::DoubleObjectFlag")
                ->debug("For HLT path: {}", hlt_path);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("No matching trigger for {} found, returning false for "
                   "trigger flag {}",
                   hlt_path, outputname);
        auto df3 = df2.Define(outputname, []() { return false; });
        return df3;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        auto df3 =
            df2.Define(outputname, trigger_matching,
                      {matched_trigger_names[0], particle_1, particle_2,
                       triggerobject_pt, triggerobject_eta, triggerobject_phi,
                       triggerobject_id_column, triggerobject_filterbit_column});
        return df3;
    }
}

/**
 * @brief This function generates a trigger flag based on the trigger object
 * matching for the selected objects. This relies on the `trigger::matchParticle`
 * function which does the object to trigger object matching test.
 *
 * @note This function is defined for double object triggers only. For single
 * object triggers be referred to the `trigger::SingleObjectFlag` functions.
 *
 * @param df input dataframe
 * @param outputname name of the output flag
 * @param particle_1 name of the column containing the Lorentz vector of the
 * first object/particle that should be matched to a trigger object
 * @param particle_2 name of the column containing the Lorentz vector of the
 * second object/particle that should be matched to a trigger object
 * @param triggerobject_pt name of the column containing the trigger object
 * \f$p_T\f$ values
 * @param triggerobject_eta name of the column containing the trigger object
 * \f$\eta\f$ values
 * @param triggerobject_phi name of the column containing the trigger object
 * \f$\phi\f$ values
 * @param triggerobject_id name of the column containing the trigger object
 * IDs
 * @param triggerobject_filterbit name of the column containing the trigger
 * object filter bits. Depending on the trigger object ID (e.g. electron,
 * muon, ...), this bitmap has a different meaning.
 * @param pt_threshold_1 \f$p_T\f$ threshold value the first tested object
 * has to exceed in order to be considered a match
 * @param pt_threshold_2 \f$p_T\f$ threshold value the second tested object
 * has to exceed in order to be considered a match
 * @param eta_threshold_1 \f$\eta\f$ threshold value the first tested object
 * has to be below in order to be considered a match
 * @param eta_threshold_2 \f$\eta\f$ threshold value the second tested object
 * has to be below in order to be considered a match
 * @param trigger_particle_id_value_1 trigger object ID value that should be
 * tested for the first object
 * @param trigger_particle_id_value_2 trigger object ID value that should be
 * tested for the second object
 * @param trigger_bit_value_1 trigger object filter bit position that should be
 * tested for the first object. If no bit matching is desired, set this value
 * to -1.
 * @param trigger_bit_value_2 trigger object filter bit position that should be
 * tested for the second object. If no bit matching is desired, set this value
 * to -1.
 * @param deltaR_threshold maximum \f$\Delta R\f$ value between the trigger
 * object and the tested object in order to be considered a match
 *
 * @return a new dataframe containing the trigger flag column
 */
ROOT::RDF::RNode DoubleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle_1, const std::string &particle_2,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit,
    const float &pt_threshold_1, const float &pt_threshold_2,
    const float &eta_threshold_1, const float &eta_threshold_2,
    const int &trigger_particle_id_value_1, const int &trigger_particle_id_value_2,
    const int &trigger_bit_value_1, const int &trigger_bit_value_2,
    const float &deltaR_threshold) {
    // In nanoAODv12 the type of trigger object ID was changed to UShort_t
    // For v9 compatibility a type casting is applied
    auto [df1, triggerobject_id_column] = utility::Cast<ROOT::RVec<UShort_t>, ROOT::RVec<Int_t>>(
            df, triggerobject_id+"_v12", "ROOT::VecOps::RVec<UShort_t>", triggerobject_id);
    // In nanoAODv15 the type of trigger object ID was changed to ULong64_t
    // For v9 and v12 compatibility a type casting is applied
    auto [df2, triggerobject_filterbit_column] = utility::Cast<ROOT::RVec<ULong64_t>, ROOT::RVec<Int_t>>(
            df1, triggerobject_filterbit+"_v15", "ROOT::VecOps::RVec<ULong64_t>", triggerobject_filterbit);

    auto trigger_matching = [pt_threshold_1, pt_threshold_2, eta_threshold_1,
                             eta_threshold_2, trigger_particle_id_value_1,
                             trigger_particle_id_value_2, trigger_bit_value_1,
                             trigger_bit_value_2, deltaR_threshold](
                                const ROOT::Math::PtEtaPhiMVector &particle_1,
                                const ROOT::Math::PtEtaPhiMVector &particle_2,
                                ROOT::RVec<float> triggerobject_pts,
                                ROOT::RVec<float> triggerobject_etas,
                                ROOT::RVec<float> triggerobject_phis,
                                ROOT::RVec<UShort_t> triggerobject_ids_v12,
                                ROOT::RVec<ULong64_t> triggerobject_filterbits_v15) {
        auto triggerobject_ids = static_cast<ROOT::RVec<int>>(triggerobject_ids_v12);
        auto triggerobject_filterbits = static_cast<ROOT::RVec<int>>(triggerobject_filterbits_v15);
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("Checking triggerobject match with particle ....");
        Logger::get("trigger::DoubleObjectFlag")->debug("First particle");
        bool match_result_p1 = matchParticle(
            particle_1, triggerobject_pts, triggerobject_etas,
            triggerobject_phis, triggerobject_ids, triggerobject_filterbits,
            pt_threshold_1, eta_threshold_1, trigger_particle_id_value_1,
            trigger_bit_value_1, deltaR_threshold);
        Logger::get("trigger::DoubleObjectFlag")->debug("Second particle");
        bool match_result_p2 = matchParticle(
            particle_2, triggerobject_pts, triggerobject_etas,
            triggerobject_phis, triggerobject_ids, triggerobject_filterbits,
            pt_threshold_2, eta_threshold_2, trigger_particle_id_value_2,
            trigger_bit_value_2, deltaR_threshold);
        bool result = match_result_p1 & match_result_p2;
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("---> Matching p1: {}", match_result_p1);
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("---> Matching p2: {}", match_result_p2);
        Logger::get("trigger::DoubleObjectFlag")
            ->debug("--->>>> result: {}", result);
        return result;
    };

    auto df3 =
        df2.Define(outputname, trigger_matching,
                    {particle_1, particle_2, triggerobject_pt, triggerobject_eta,
                    triggerobject_phi, triggerobject_id_column, triggerobject_filterbit_column});
    return df3;
}

/**
 * @brief This function generates a new column containing the prescale value
 * for a trigger given run and lumiblock. It is read from an external JSON file
 * containing prescale values for a specific trigger.
 *
 * @note The JSON files can be produced with the `brilcalc` tool from LUM POG
 * which needs to be run on lxplus. More details can be found here in the
 * [documentation](https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html#brilcalctrg)
 * of `brilcalc`.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading
 * the prescale trigger values
 * @param outputname name of the output column containing the prescale value
 * @param hlt_path name of the column containing the HLT path
 * @param run name of the column containing the run number
 * @param lumiblock name of the column containing the event lumiblock
 * @param prescale_file relative path to the JSON containing the prescale values
 *
 * @return a new dataframe containing the prescale value column
 */
ROOT::RDF::RNode GetPrescaleValues(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &hlt_path,
    const std::string &run, const std::string &lumiblock,
    const std::string &prescale_file) {

    Logger::get("trigger::GetPrescaleValues")
        ->debug("reading json from {}", prescale_file);
    const nlohmann::json prescale_json =
        *correction_manager.loadjson(prescale_file);

    auto get_prescale = [prescale_json](const Bool_t hlt, const UInt_t run,
                                        const UInt_t lumiblock) {
        int prescale = -1;
        Logger::get("trigger::GetPrescaleValues")
            ->debug("run: {}, lumiblock: {}", run, lumiblock);

        if (hlt == false) {
            prescale = -2;
            Logger::get("trigger::GetPrescaleValues")
                ->debug("no HLT hit, prescale value: {}", prescale);
            return prescale;
        }

        const std::string s_run = std::to_string(run);
        if (prescale_json.find(s_run) != prescale_json.end()) {
            Logger::get("trigger::GetPrescaleValues")
                ->debug("found run in JSON ...");

            unsigned highest_lumi = 1;
            for (auto &[i_key, i_value] : prescale_json[s_run].items()) {
                Logger::get("trigger::GetPrescaleValues")
                    ->debug("... checking lumi {}, prescale {} ...",
                            std::stoi(i_key), int(i_value));
                if (lumiblock >= std::stoi(i_key)) {
                    if (std::stoi(i_key) >= highest_lumi) {
                        highest_lumi = std::stoi(i_key);
                        prescale = i_value;
                        Logger::get("trigger::GetPrescaleValues")
                            ->debug("... assigning prescale value: {}",
                                    prescale);
                    }
                }
            }
        } else {
            prescale = -3;
            Logger::get("trigger::GetPrescaleValues")
                ->debug(
                    "could not find run and lumi in JSON, prescale value: {}",
                    prescale);
        }
        return prescale;
    };
    auto df1 = df.Define(outputname, get_prescale, {hlt_path, run, lumiblock});
    return df1;
}
} // end namespace trigger
#endif /* GUARD_TRIGGERS_H */
