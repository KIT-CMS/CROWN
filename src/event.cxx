#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include <cmath>
#include <exception>
#include <nlohmann/json.hpp>
#include <openssl/sha.h>

namespace event {
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
 * @param master_seed master seed value to be added to the hash used for event
 * seed generation
 *
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode GenerateSeed(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &lumi, const std::string &run,
                              const std::string &event,
                              const UInt_t &master_seed = 42) {

    auto generate_seed = [master_seed](const unsigned int &lumi,
                                       const unsigned int &run,
                                       const unsigned long long &event) {
        // string for setting the seed value
        const std::string seed_string =
            std::to_string(master_seed) + "_" + std::to_string(lumi) + "_" +
            std::to_string(run) + "_" + std::to_string(event);

        // create a SHA256 has from the seed string
        unsigned char hash[SHA256_DIGEST_LENGTH];
        SHA256(reinterpret_cast<const unsigned char *>(seed_string.c_str()),
               seed_string.size(), hash);

        // use the first for bits of the hash to create a 32-bit unsigned
        // integer as seed
        unsigned int event_seed = 0;
        for (int i = 0; i < 4; ++i) {
            event_seed = (event_seed << 8) | hash[i];
        }

        return event_seed;
    };

    return df.Define(outputname, generate_seed, {lumi, run, event});
}

/**
 * @brief Weight corresponds to the `event::filter::HashFilter` filter.
 *
 * This filter is unphysical and has been introduced to create independent
 * event samples for different eras. Thus, when the sample is filtered, the
 * lower number of selected MC events needs to be compensated with a weight,
 * which scales the sample back to the correct normalization. The weight is
 * obtained by inverting the size of the range `[u_min, u_max)` (with both
 * numbers between 0 and 1), from which the events have been drawn:
 *
 * ```cpp
 * 1 / (u_max - u_min)
 * ```
 *
 * This function can be used to split up a MC sample into several subsamples
 * for different data-taking eras.
 *
 * @param df The input data frame
 * @param output Name of the output weight
 * @param u_min Minimal value of the considered interval for the event passing
 * the filter; must be a value between 0 and 1
 * @param u_max Maximal value of the considered interval for the event passing
 * the filter; must be a value between 0 and 1
 *
 * @return the data frame with the new column
 *
 * @throws std::runtime_error if `u_min` or `u_max` are outside the `[0, 1]`
 * interval, or if `u_min >= u_max`.
 */
ROOT::RDF::RNode SampleSplitWeight(ROOT::RDF::RNode df,
                                   const std::string &output,
                                   const double &u_min, const double &u_max) {
    // Define the logger name
    const std::string logger_name = "event::quantity::SampleSplitWeight";

    // Catch misconfiguration of u_min and u_max
    if (u_min < 0 || u_min > 1 || u_max < 0 || u_max > 1) {
        const std::string msg =
            "u_min and u_max must be numbers between 0 and 1";
        Logger::get(logger_name)->critical(msg);
        throw std::runtime_error(msg);
    }
    if (u_min >= u_max) {
        const std::string msg = "u_min must be strictly smaller than u_max";
        Logger::get(logger_name)->critical(msg);
        throw std::runtime_error(msg);
    }

    return df.Define(output,
                     [u_min, u_max]() { return (float)1 / (u_max - u_min); });
}

} // end namespace quantity

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

/**
 * @brief Keep only events whose hash falls into a given fraction of the hash
 * range, e.g., to split up a sample into mutually exclusive ones.
 *
 * The event hash can be generated using the `event::quantity::GenerateSeed`
 * function. The output is a `uint32_t`. The range `[u_min, u_max)` defines
 * which events are selected, with `u_min` and `u_max` being numbers between 0
 * and 1. These boundaries are scaled to the `uint32_t` value range. The event
 * passes the filter if `event_hash` is within the scaled interval:
 *
 * ```
 * u_min * 2^32 <= event_hash < u_max * 2^32
 * ```
 *
 * This function can be used to split up a MC sample into several subsamples
 * for different data-taking eras.
 *
 * @param df The input data frame
 * @param filtername Name of the filter
 * @param event_hash Name of the column containing the event hash values
 * @param u_min Minimal value of the considered interval for the event passing
 * the filter; will be mapped to uint32 range; must be a value between 0 and 1
 * @param u_max Maximal value of the considered interval for the event passing
 * the filter; will be mapped to uint32 range; must be a value between 0 and 1
 *
 * @return the filtered data frame
 *
 * @throws std::runtime_error if `u_min` or `u_max` are outside the `[0, 1]`
 * interval, or if `u_min >= u_max`.
 */
ROOT::RDF::RNode HashFilter(ROOT::RDF::RNode df, const std::string &filtername,
                            const std::string &event_hash, const double &u_min,
                            const double &u_max) {
    // Define the logger name
    const std::string logger_name = "event::filter::HashFilter";

    // Catch misconfiguration of u_min and u_max
    if (u_min < 0 || u_min > 1 || u_max < 0 || u_max > 1) {
        const std::string msg =
            "u_min and u_max must be numbers between 0 and 1";
        Logger::get(logger_name)->critical(msg);
        throw std::runtime_error(msg);
    }
    if (u_min >= u_max) {
        const std::string msg = "u_min must be strictly smaller than u_max";
        Logger::get(logger_name)->critical(msg);
        throw std::runtime_error(msg);
    }

    // Define lower and upper bounds on the hash. As it is a uint32, the
    // largest possible value is 2**32 - 1
    auto hash_min = uint64_t(u_min * std::pow(2, 32));
    auto hash_max = uint64_t(u_max * std::pow(2, 32));

    auto filter = [logger_name, hash_min,
                   hash_max](const unsigned int &event_hash) {
        // Keep the event if the hash is inside the [hash_min, hash_max)
        // interval
        Logger::get(logger_name)
            ->debug("Got hash {}, check if it is in [{}, {})", event_hash,
                    hash_min, hash_max);
        return (hash_min <= event_hash) && (event_hash < hash_max);
    };

    return df.Filter(filter, {event_hash}, filtername);
}

} // end namespace filter
} // end namespace event
