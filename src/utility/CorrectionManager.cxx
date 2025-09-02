#include "../../include/utility/CorrectionManager.hxx"
#include "../../include/utility/Logger.hxx"
#include <fstream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
/// namespace used for the CorrectionManager
namespace correctionManager {
/**
 * @brief Load a correctionlib correction from a file and return it. This
 * function works for most corrections found in correctionlib files (the default
 * type for corrections is correction::Correction) . If the requested correction
 * is already loaded, it will return a pointer to the already loaded correction.
 * If not, the correction will be loaded and stored in the CorrectionManager.
 *
 * @param filePath The path to the correctionlib file
 * @param corrName The name of the correction to load
 * @return const correction::Correction*
 */
const correction::Correction *
CorrectionManager::loadCorrection(const std::string &filePath,
                                  const std::string &corrName) {
    auto filePath_it = correction_map.find(filePath);
    if (filePath_it == correction_map.end()) {
        Logger::get("CorrectionManager")
            ->debug("CorrectionFile {} not loaded yet, adding {} to the "
                    "CorrectionManager...", 
                    filePath, corrName);
        auto result = correction_map.emplace(
            filePath,
            std::unordered_map<
                std::string, std::shared_ptr<const correction::Correction>>());
        filePath_it = result.first;
    }

    auto &corr_name = filePath_it->second;
    auto corrName_it = corr_name.find(corrName);
    if (corrName_it == corr_name.end()) {
        auto correctionset = correction::CorrectionSet::from_file(filePath);
        corrName_it =
            corr_name.emplace(corrName, correctionset->at(corrName)).first;
        Logger::get("CorrectionManager")
            ->debug("Created correction for model - name: {} - {}", filePath,
                    corrName);
        n_corrections++;
    } else {
        Logger::get("CorrectionManager")
            ->debug("Correction already exists for file - name: {} - "
                    "{}, using preloaded correction",
                    filePath, corrName);
    }
    return corrName_it->second.get();
};
/**
 * @brief Load a correctionlib correction from a file and return it. This
 * function works correction::CompoundCorrection objects. If the requested
 * CompoundCorrection is already loaded, it will return a pointer to the already
 * loaded CompoundCorrection. If not, the CompoundCorrection will be loaded and
 * stored in the CorrectionManager.
 *
 * @param filePath The path to the correctionlib file
 * @param corrName The name of the correction to load
 * @return const correction::Correction*
 */
const correction::CompoundCorrection *
CorrectionManager::loadCompoundCorrection(const std::string &filePath,
                                          const std::string &corrName) {
    auto filePath_it = correctionCompound_map.find(filePath);
    if (filePath_it == correctionCompound_map.end()) {
        Logger::get("CorrectionManager")
            ->debug("Compound CorrectionFile {} not loaded yet, adding it to the "
                    "CorrectionManager...",
                    filePath);
        auto result = correctionCompound_map.emplace(
            filePath,
            std::unordered_map<
                std::string,
                std::shared_ptr<const correction::CompoundCorrection>>());
        n_corrections++;
        filePath_it = result.first;
    }

    auto &corr_name = filePath_it->second;
    auto corrName_it = corr_name.find(corrName);
    if (corrName_it == corr_name.end()) {
        auto correctionset = correction::CorrectionSet::from_file(filePath);
        corrName_it =
            corr_name.emplace(corrName, correctionset->compound().at(corrName))
                .first;
        Logger::get("CorrectionManager")
            ->debug("Created CompoundCorrection for model - name: {} - {}",
                    filePath, corrName);
        n_corrections++;
    } else {
        Logger::get("CorrectionManager")
            ->debug("CompoundCorrection already exists for file - name: {} - "
                    "{}, using preloaded CompoundCorrection",
                    filePath, corrName);
    }
    return corrName_it->second.get();
};
/**
 * @brief Load a json file from a file and return it. This
 * function works works for all json files. If the requested
 * json file is already loaded, it will return a pointer to the already
 * loaded json. If not, the json will be loaded and
 * stored in the CorrectionManager.
 *
 * @param filePath The path to the json file
 * @return const nlohmann::json
 */
const nlohmann::json *CorrectionManager::loadjson(const std::string &filePath) {
    auto json_it = json_map.find(filePath);
    if (json_it == json_map.end()) {
        Logger::get("CorrectionManager")
            ->debug("Json file {} not loaded yet, adding it to the "
                    "CorrectionManager...",
                    filePath);
        std::ifstream json_file(filePath);
        auto result = json_map.emplace(
            filePath,
            std::make_shared<nlohmann::json>(nlohmann::json::parse(json_file)));
        n_corrections++;
        json_it = result.first;
    }
    return json_it->second.get();
}

/**
 * @brief Report the number of corrections managed by the CorrectionManager
 */
void CorrectionManager::report() {
    Logger::get("CorrectionManager")
        ->info("CorrectionManager manages {} different corrections for this "
               "production",
               n_corrections);
};

} // namespace correctionManager