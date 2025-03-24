#ifndef GUARDEVENT_H
#define GUARDEVENT_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include <nlohmann/json.hpp>

namespace event {

/**
 * @brief Function to filter events based on their `run` and `luminosity` block
 * values from the golden json. The json files are taken from the CMS
 * recommendations.
 *
 * Run2: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
 *
 * Run3: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3 (not
 * added yet)
 *
 * @param df input dataframe
 * @param correctionManager correction manager to be used to read the golden
 * json
 * @param filtername name of the filter, used in the dataframe report
 * @param run column containing the run value
 * @param luminosity column containing the luminosity block value
 * @param json_path path to the golden json file containing all valid
 * runs and luminosity blocks
 * @return a filtered dataframe
 */
ROOT::RDF::RNode
GoldenJSONFilter(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correctionManager,
                 const std::string &filtername, const std::string &run,
                 const std::string &luminosity, const std::string &json_path) {
    nlohmann::json golden_json = *correctionManager.loadjson(json_path);
    auto jsonFilterlambda = [golden_json](UInt_t run, UInt_t luminosity) {
        bool matched = false;
        // check if the run exists
        if (golden_json.find(std::to_string(run)) != golden_json.end()) {
            // now loop over all luminosity blocks and check if the event is
            // valid
            for (auto &luminosityrange : golden_json[std::to_string(run)]) {
                if (luminosity >= luminosityrange[0] &&
                    luminosity <= luminosityrange[1]) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                Logger::get("GoldenJSONFilter")
                    ->debug("Run {} / luminosity {} not in json file", run,
                            luminosity);
            }
        }
        return matched;
    };
    return df.Filter(jsonFilterlambda, {run, luminosity}, filtername);
}
} // namespace event

#endif /* GUARDEVENT_H */