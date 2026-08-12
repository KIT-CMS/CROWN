#ifndef GUARD_TAUVARIATIONS_H
#define GUARD_TAUVARIATIONS_H

#include <regex>
#include <format>
#include <vector>
#include <functional>
#include "correction.h"


namespace physicsobject {

namespace tau {

namespace scalefactor {

const int DEFAULT_DECAY_MODE = -10;
const float DEFAULT_PT_MIN = -10.0;
const float DEFAULT_PT_MAX = -10.0;

// Class to handle custom variations of tau ID vs jet scale factors
class TauIDVsJetVariation {
public:
    // Constructor
    TauIDVsJetVariation(const std::string&);

    // Wrap correction::Correction::evaluate function
    std::function<double (const std::vector<correction::Variable::Type>&)> wrap_evaluate(const correction::Correction*) const;

private:
    // Name of the variation passed to the constructor
    std::string custom_variation_;

    // Name of the variation in the correction file
    std::string variation_;

    // Flags for custom variations and for imposed selections on DM and pt
    bool has_dm_selection_;
    bool has_pt_selection_;

    // Decay mode and pt selection values for the variation
    int decay_mode_ = DEFAULT_DECAY_MODE;
    float pt_min_ = DEFAULT_PT_MIN;
    float pt_max_ = DEFAULT_PT_MAX;

    // Get index of a variable in the list of inputs of a correction::Correction
    size_t get_variable_index(const correction::Correction*, const std::string&) const;
};

} // end namespace scalefactor

} // end namespace tau

} // end namespace physicsobject


#endif // GUARD_TAUVARIATIONS_H
