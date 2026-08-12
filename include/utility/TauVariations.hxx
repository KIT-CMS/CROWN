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

// Class to handle custom variations of tau ID vs jet scale factors
class TauIDVsJetVariation {
public:
    // Constructor
    TauIDVsJetVariation(const std::string&);

    // Wrap correction::Correction::evaluate function
    std::function<double (const std::vector<Variable::Type>&)> wrap_evaluate(const correction::Correction*) const;

private:
    // Name of the variation in the correction file
    std::string variation_;

    // Decay mode and pt selection values for the variation
    int decay_mode_ = -10;
    float pt_min_ = -10.0;
    float pt_max_ = -10.0;

    // Flags for custom variations and for imposed selections on DM and pt
    bool has_dm_selection_ = false;
    bool has_pt_selection_ = false;

    // Get index of a variable in the list of inputs of a correction::Correction
    size_t TauIDVsJetVariation::get_variable_index(const correction::Correction&, const std::string&) const;
};

} // end namespace scalefactor

} // end namespace tau

} // end namespace physicsobject


#endif // GUARD_TAUVARIATIONS_H
