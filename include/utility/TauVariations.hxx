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

// Map generator-level match types to their corresponding integer values
enum struct GenType : std::vector<int> {
    NONE = {},         // no generator-level match information
    GEN_ELE = {1, 3},  // genuine electron
    GEN_MU = {2, 4},   // genuine muon
    GEN_TAU = {5}      // genuine tau
};

// Map eta regions of electrons to ranges of absolute eta values
enum struct EtaRegion : std::vector<float> {
    NONE = {},            // no eta region information
    BARREL = {0.0, 1.5},  // barrel electron
    ENDCAP = {1.5, 2.5}   // endcap electron
};

// Encapsulate logic for generator-level match restrictions
class GenMatchRestriction {
public:
    GenMatchRestriction();
    GenMatchRestriction(const GenType &gen_type);
    bool is_selected(const int &gen_match);

private:
    bool restrict_;
    std::vector<int> gen_matches_;
}

// Encapsulate logic for decay mode restrictions
class DecayModeRestriction {
public:
    DecayModeRestriction();
    DecayModeRestriction(const std::vector<int> &decay_modes);
    bool is_selected(const int &decay_mode);

private:
    bool restrict_;
    std::vector<int> decay_modes_;
}

// Encapsulate logic for pt restrictions
class PtRestriction {
public:
    PtRestriction();
    PtRestriction(const float &, const float &);
    PtRestriction(const float &);

private:
    bool restrict_;
    std::pair<float, float> pt_range_;
}

// Encapsulate logic for eta restrictions
class EtaRestriction {
public:
    EtaRestriction();
    EtaRestriction(const std::pair<float, float> &);
    EtaRestriction(const float &, const float &);
    EtaRestriction(const EtaRegion &);
    bool is_selected(const float &);

private:
    bool restrict_;
    std::pair<float, float> abs_eta_range_;
}

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
