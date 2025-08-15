#ifndef GUARD_EVENT_H
#define GUARD_EVENT_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include <nlohmann/json.hpp>
#include <openssl/sha.h>

namespace event {
namespace filter {

/**
 * @brief This function applies a filter to the input dataframe using a Golden
 * JSON file, which contains a mapping of valid run-luminosity pairs. The
 * dataframe is filtered by checking if the run and luminosity values for each
 * row match the entries in the Golden JSON. Rows with invalid run-luminosity
 * pairs are removed.
 *
 * The Golden JSON files are taken from the CMS recommendations.
 *
 * Run2: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3 (not
 * added yet)
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading
 * the Golden JSON
 * @param filtername name of the filter to be applied (used in the dataframe
 * report)
 * @param run name of the run column
 * @param luminosity name of the luminosity column
 * @param json_path path to the Golden JSON file
 *
 * @return a filtered dataframe
 */
ROOT::RDF::RNode
GoldenJSON(ROOT::RDF::RNode df,
           correctionManager::CorrectionManager &correction_manager,
           const std::string &filtername, const std::string &run,
           const std::string &luminosity, const std::string &json_path) {
    nlohmann::json golden_json = *correction_manager.loadjson(json_path);
    auto jsonFilterlambda = [golden_json](UInt_t run, UInt_t luminosity) {
        bool matched = false;
        // check if the run exists
        if (golden_json.find(std::to_string(run)) != golden_json.end()) {
            // now loop over all luminosity blocks and check if the event is
            // valid
            for (auto &luminosity_range : golden_json[std::to_string(run)]) {
                if (luminosity >= luminosity_range[0] &&
                    luminosity <= luminosity_range[1]) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                Logger::get("event::filter::GoldenJSON")
                    ->debug("Run {} / luminosity {} not in json file", run,
                            luminosity);
            }
        }
        return matched;
    };
    return df.Filter(jsonFilterlambda, {run, luminosity}, filtername);
}
} // end namespace filter

namespace quantity {

/**
 * @brief This function defines a new column in the dataframe with seeds for a
 * random number generator for each event.
 * 
 * The seed value for each event is calculated by concatenating event index
 * variables and a seed value to `{seed}_{lumi}_{run}_{event}`. From that, a
 * SHA256 hash is calculated. The first four bytes of the hash are then used
 * to create a 32-bit unsigned integer, which serves as the event seed.
 * 
 * @param df input dataframe
 * @param outputname name of the new column containing the generated event seeds
 * @param lumi name of the column containing the luminosity block number
 * @param run name of the column containing the run number
 * @param event name of the column containing the event number
 * @param master_seed master seed value to be added to the hash used for event seed
 * generation
 *
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode
GenerateSeed(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &lumi,
    const std::string &run,
    const std::string &event,
    const UInt_t &master_seed = 42
) {
   
    auto generate_seed = [master_seed] (
        const unsigned int &lumi,
        const unsigned int &run,
        const unsigned long long &event
    ) {
        // string for setting the seed value
        const std::string seed_string = std::to_string(master_seed) + "_" + std::to_string(lumi) + "_" + std::to_string(run) + "_" + std::to_string(event);

        // create a SHA256 has from the seed string
        unsigned char hash[SHA256_DIGEST_LENGTH];
        SHA256(reinterpret_cast<const unsigned char*>(seed_string.c_str()), seed_string.size(), hash);

        // use the first for bits of the hash to create a 32-bit unsigned integer as seed
        unsigned int event_seed = 0;
        for (int i = 0; i < 4; ++i) {
            event_seed = (event_seed << 8) | hash[i];
        }

        return event_seed;
    };

    return df.Define(
        outputname,
        generate_seed,
        {lumi, run, event}
    );
}

} // end namespace quantity


namespace lhe {

/**
 * @brief Extract the flavor of a leptonic Drell-Yan process from the LHE information.
 * This function checks if exactly two final-state LHE leptons (electron, muon, tau) of
 * the same flavor are present in the event. If so, the PDG ID of the lepton that the
 * Z boson or the boson decays into is returned. If no such pair is found, the function
 * returns -1.
 * 
 * Possible values of the returned PDG ID are:
 * - 11 for a decay into an electron pair,
 * - 13 for a decay into a muon pair,
 * - 15 for a decay into a tau pair.
 * 
 * @param df input dataframe
 * @param outputname name of the new column containing the PDG ID of the decay particles
 * @param lhe_pdgid name of the column containing the PDG ID values of the LHE particles
 * @param lhe_status name of the column containing the status of the LHE particles
 * 
 * @return a new dataframe with the decay flavor column added.
 */
ROOT::RDF::RNode DrellYanDecayFlavor(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &lhe_pdgid,
    const std::string &lhe_status
) {
    auto flavor_flag = [] (
        const ROOT::RVec<int> &lhe_pdgid,
        const ROOT::RVec<int> &lhe_status
    ) {
        int decay_flavor_pdgid = -1; 
        const std::vector<int> decay_flavors = {11, 13, 15};
        for (const auto &decay_flavor : decay_flavors) {
            if (ROOT::VecOps::Sum((lhe_status == 1) && (abs(lhe_pdgid) == decay_flavor)) == 2) {
                decay_flavor_pdgid = decay_flavor;
                break;
            }
        }
        return decay_flavor_pdgid;
    };

    return df.Define(
        outputname,
        flavor_flag,
        {lhe_pdgid, lhe_status}
    );
}

} // end namespace lhe
} // end namespace event

#endif /* GUARD_EVENT_H */