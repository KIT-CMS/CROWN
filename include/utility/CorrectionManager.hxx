#ifndef GUARD_CORRECTION_MANAGER
#define GUARD_CORRECTION_MANAGER

#include "correction.h"
#include <memory>
#include <unordered_map>

namespace correctionManager {
class CorrectionManager {
  public:
    CorrectionManager(){};
    ~CorrectionManager() {
        correction_map.clear();
        correctionCompound_map.clear();
    }
    const correction::Correction *loadCorrection(const std::string &filePath,
                                                 const std::string &corrName);
    const correction::CompoundCorrection *
    loadCompoundCorrection(const std::string &filePath,
                           const std::string &corrName);

    void report();

  private:
    std::unordered_map<
        std::string,
        std::unordered_map<std::string,
                           std::shared_ptr<const correction::Correction>>>
        correction_map;
    std::unordered_map<
        std::string,
        std::unordered_map<
            std::string, std::shared_ptr<const correction::CompoundCorrection>>>
        correctionCompound_map;
    int n_corrections = 0;
};
} // namespace correctionManager
#endif /* GUARD_CORRECTION_MANAGER */