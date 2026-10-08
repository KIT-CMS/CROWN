#ifndef GUARD_EVENT_H
#define GUARD_EVENT_H

#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include <nlohmann/json.hpp>
#include <openssl/sha.h>
#include <array>
#include <memory>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace event {
namespace quantity {

namespace {

std::string SampleNormalizationParseNickFromPath(const std::string &sample_id) {
    auto tree_sep = sample_id.rfind('/');
    std::string filename =
        (tree_sep == std::string::npos) ? sample_id : sample_id.substr(0, tree_sep);
    std::vector<std::string> parts;
    size_t pos = 0;
    size_t next;
    while ((next = filename.find('/', pos)) != std::string::npos) {
        parts.push_back(filename.substr(pos, next - pos));
        pos = next + 1;
    }
    parts.push_back(filename.substr(pos));
    if (parts.size() < 3) {
        throw std::runtime_error(
            "event::quantity::SampleNormalization: path '" + filename +
            "' has too few segments to contain a nick (expected "
            ".../era/nick/scope/nick_N.root)");
    }
    return parts[parts.size() - 3];
}

double SampleNormalizationLookupField(const nlohmann::json &norm_table,
                                      const std::string &nick,
                                      const std::string &field) {
    if (!norm_table.contains(nick)) {
        Logger::get("event::quantity::SampleNormalization")
            ->error("nick '{}' not found in the normalization table, "
                    "sample_database is missing an entry for this ",
                    nick);
        throw std::runtime_error(
            "event::quantity::SampleNormalization: unknown nick " + nick);
    }
    return norm_table.at(nick).at(field).get<double>();
}

} // namespace

/**
 * @brief This function defines three per-file constant columns --
 * `xsec_output` (cross section), `ngen_weight_output` (1 / number of
 * generated events) and `genweight_output` (effective normalization factor
 * accounting for negative generator weights) -- by looking up this file's
 * sample nick in a `nick -> {xsec, nevents, generator_weight}` JSON table.
 *
 * The sample nick is parsed from the input file's path at runtime via
 * `ROOT::RDF::RSampleInfo` (path convention: `.../{era}/{nick}/{scope}/
 * {nick}_{N}.root`), so one compiled executable correctly normalizes every
 * nick contained in its input sample. An unknown nick throws rather than
 * silently defaulting.
 *
 * @param df input dataframe
 * @param correctionManager correction manager responsible for loading the
 * normalization JSON table
 * @param xsec_output name of the new column containing the cross section
 * @param ngen_weight_output name of the new column containing 1/nevents
 * @param genweight_output name of the new column containing the effective
 * generator-weight normalization factor
 * @param norm_table_path path to the `nick -> {xsec, nevents,
 * generator_weight}` JSON lookup table
 *
 * @return a dataframe with the three new columns
 */
ROOT::RDF::RNode
SampleNormalization(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &xsec_output,
                    const std::string &ngen_weight_output,
                    const std::string &genweight_output,
                    const std::string &norm_table_path) {
    nlohmann::json norm_table = *correctionManager.loadjson(norm_table_path);

    // crossSectionPerEventWeight -- the raw xsec (pb)
    auto df1 = df.DefinePerSample(
        xsec_output,
        [norm_table](unsigned int /*slot*/, const ROOT::RDF::RSampleInfo &id) {
            return SampleNormalizationLookupField(
                norm_table, SampleNormalizationParseNickFromPath(id.AsString()),
                "xsec");
        });
    // numberGeneratedEventsWeight -- 1/nevents
    auto df2 = df1.DefinePerSample(
        ngen_weight_output,
        [norm_table](unsigned int /*slot*/, const ROOT::RDF::RSampleInfo &id) {
            return 1.0 / SampleNormalizationLookupField(
                             norm_table,
                             SampleNormalizationParseNickFromPath(id.AsString()),
                             "nevents");
        });
    // negative_events_fraction -- effective normalization factor
    auto df3 = df2.DefinePerSample(
        genweight_output,
        [norm_table](unsigned int /*slot*/, const ROOT::RDF::RSampleInfo &id) {
            return SampleNormalizationLookupField(
                norm_table, SampleNormalizationParseNickFromPath(id.AsString()),
                "generator_weight");
        });
    return df3;
}

ROOT::RDF::RNode
STXSNormalization(ROOT::RDF::RNode df,
                  correctionManager::CorrectionManager &correctionManager,
                  const std::string &family_output,
                  const std::string &norm_output,
                  const std::string &scale_up_output,
                  const std::string &scale_down_output,
                  const std::string &norm_table_path,
                  const std::string &htxs_bin_column,
                  const std::string &lhe_scale_up_column,
                  const std::string &lhe_scale_down_column) {
    struct STXSTable {
        int family;
        std::unordered_map<int, std::array<double, 3>> bins;
    };
    nlohmann::json norm_table = *correctionManager.loadjson(norm_table_path);
    auto tables = std::make_shared<std::vector<STXSTable>>();
    auto index = std::make_shared<std::unordered_map<std::string, int>>();
    for (const auto &[nick, entry] : norm_table.items()) {
        if (!entry.contains("stxs")) {
            continue;
        }
        STXSTable table{entry.at("stxs").at("family").get<int>(), {}};
        for (const auto &[bin, factors] : entry.at("stxs").at("bins").items()) {
            table.bins[std::stoi(bin)] = {factors.at(0).get<double>(),
                                          factors.at(1).get<double>(),
                                          factors.at(2).get<double>()};
        }
        (*index)[nick] = tables->size();
        tables->push_back(table);
    }
    const std::string index_column = family_output + "_nick_index";
    auto df1 = df.DefinePerSample(
        index_column,
        [index](unsigned int /*slot*/, const ROOT::RDF::RSampleInfo &id) {
            auto it = index->find(
                SampleNormalizationParseNickFromPath(id.AsString()));
            return it == index->end() ? -1 : it->second;
        });
    auto lookup = [tables](int idx, int bin) -> const std::array<double, 3> * {
        if (idx < 0) {
            return nullptr;
        }
        auto it = tables->at(idx).bins.find(bin);
        return it == tables->at(idx).bins.end() ? nullptr : &it->second;
    };
    auto df2 = df1.Define(
        family_output,
        [tables](int idx) { return idx < 0 ? 0 : tables->at(idx).family; },
        {index_column});
    auto df3 = df2.Define(
        norm_output,
        [lookup](int idx, int bin) {
            auto factors = lookup(idx, bin);
            return factors ? (*factors)[0] : 1.0;
        },
        {index_column, htxs_bin_column});
    auto df4 = df3.Define(
        scale_up_output,
        [lookup](int idx, int bin, float scale) {
            auto factors = lookup(idx, bin);
            return factors ? scale / (*factors)[1] : scale;
        },
        {index_column, htxs_bin_column, lhe_scale_up_column});
    return df4.Define(
        scale_down_output,
        [lookup](int idx, int bin, float scale) {
            auto factors = lookup(idx, bin);
            return factors ? scale / (*factors)[2] : scale;
        },
        {index_column, htxs_bin_column, lhe_scale_down_column});
}

ROOT::RDF::RNode
STXSLheScale(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correctionManager,
             const std::string &output, const std::string &table_path,
             const std::string &variation, const std::string &family_column,
             const std::string &htxs_bin_column,
             const std::string &scale_up_column,
             const std::string &scale_down_column) {
    if (variation == "nominal") {
        return df.Define(output, []() { return 1.0; });
    }
    const auto split = variation.rfind('_');
    const std::string nuisance = variation.substr(0, split);
    const bool use_up = variation.substr(split + 1) == "up";
    nlohmann::json table = *correctionManager.loadjson(table_path);
    auto members = std::make_shared<std::set<std::pair<int, int>>>();
    for (const auto &pair : table.at(nuisance)) {
        members->insert({pair.at(0).get<int>(), pair.at(1).get<int>()});
    }
    return df.Define(
        output,
        [members, use_up](int family, int bin, double up, double down) {
            return members->count({family, bin}) ? (use_up ? up : down) : 1.0;
        },
        {family_column, htxs_bin_column, scale_up_column, scale_down_column});
}

/**
 * @brief This function creates a new column with `sign(genWeight) /
 * negative_fraction`. This is the standard normalization for MC generators
 * that produce negative event weights (e.g. amc@NLO, POWHEG): using only the
 * sign of the generator weight, divided by the sample's effective
 * normalization factor (`1 - 2 * (fraction of negative-weight events)`),
 * normalizes by the *effective* number of events instead of the raw event
 * count, so negative-weight events correctly dilute the yield rather than
 * being ignored or double-penalized.
 *
 * @param df input dataframe
 * @param outputname name of the new column
 * @param genweight_quantity name of the column containing the generator
 * weight (`Float_t`, as stored in NanoAOD)
 * @param negative_fraction_quantity name of the column containing the
 * sample's effective normalization factor (`1 - 2 * negative-weight
 * fraction`)
 *
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode
NormalizedGenWeightSign(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &genweight_quantity,
                        const std::string &negative_fraction_quantity) {
    return df.Define(
        outputname,
        [](const float &genWeight, const double &negative_fraction) {
            double sign = (genWeight < 0) ? -1.0 : 1.0;
            return sign / negative_fraction;
        },
        {genweight_quantity, negative_fraction_quantity});
}

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
} // end namespace filter
} // end namespace event

#endif /* GUARD_EVENT_H */