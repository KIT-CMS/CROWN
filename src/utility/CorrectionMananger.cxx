#include <memory>
#include <string>
#include "../../include/utility/CorrectionManager.hxx"
#include "../../include/utility/Logger.hxx"

const correction::Correction *
CorrectionManager::loadCorrection(const std::string &filePath,
                                  const std::string &corrName) {
    auto filePath_it = correction_map.find(filePath);
    if (filePath_it == correction_map.end()) {
        Logger::get("CorrectionManager")
            ->debug("CorrectionFile {} not loaded yet, adding it...", filePath);
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
const correction::CompoundCorrection *
CorrectionManager::loadCompoundCorrection(const std::string &filePath,
                                          const std::string &corrName) {
    auto filePath_it = correctionCompound_map.find(filePath);
    if (filePath_it == correctionCompound_map.end()) {
        Logger::get("CorrectionManager")
            ->debug("CorrectionFile {} not loaded yet, adding it...", filePath);
        auto result = correctionCompound_map.emplace(
            filePath,
            std::unordered_map<
                std::string,
                std::shared_ptr<const correction::CompoundCorrection>>());
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
            ->debug("Correction already exists for file - name: {} - "
                   "{}, using preloaded CompoundCorrection",
                   filePath, corrName);
    }
    return corrName_it->second.get();
};
void CorrectionManager::report() {
    Logger::get("CorrectionManager")
        ->info("CorrectionManager manages {} different corrections for this production", n_corrections);
};