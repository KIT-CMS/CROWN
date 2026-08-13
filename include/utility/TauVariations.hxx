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

// Map genuine object types to indices of generator-level match values
enum struct GenType {
    GEN_ELE,  // genuine electron
    GEN_MU,   // genuine muon
    GEN_TAU   // genuine tau 
};
const std::unordered_map<GenType, std::vector<int>> GEN_TYPES = {
    {GenType::GEN_ELE, {1, 3}},
    {GenType::GEN_MU, {2, 4}},
    {GenType::GEN_TAU, {5}}
};

// Map eta regions of electrons to ranges of absolute eta values
enum struct EtaRegion {
    BARREL,  // barrel electron
    ENDCAP   // endcap electron
};
const std::unordered_map<EtaRegion, std::pair<float, float>> ETA_REGIONS = {
    {EtaRegion::BARREL, {0.0, 1.5}},
    {EtaRegion::ENDCAP, {1.5, 2.5}}
};

// Encapsulate logic for generator-level match restrictions
class GenMatchRestriction {
public:
    GenMatchRestriction();
    GenMatchRestriction(const std::vector<int>&);
    bool is_active() const;
    bool is_selected(const int &) const;
    std::string repr() const;

private:
    bool restrict_;
    std::vector<int> gen_matches_;
};

// Encapsulate logic for decay mode restrictions
class DecayModeRestriction {
public:
    DecayModeRestriction();
    DecayModeRestriction(const std::vector<int>&);
    DecayModeRestriction(const int&);
    bool is_active() const;
    bool is_selected(const int &) const;
    std::string repr() const;

private:
    bool restrict_;
    std::vector<int> decay_modes_;
};

// Encapsulate logic for pt restrictions
class PtRestriction {
public:
    PtRestriction();
    PtRestriction(const float &, const float &);
    PtRestriction(const float &);
    bool is_active() const;
    bool is_selected(const float &) const;
    std::string repr() const;

private:
    bool restrict_;
    std::pair<float, float> pt_range_;
};

// Encapsulate logic for eta restrictions
class EtaRestriction {
public:
    EtaRestriction();
    EtaRestriction(const std::pair<float, float> &);
    EtaRestriction(const float &, const float &);
    bool is_active() const;
    bool is_selected(const float &) const;
    std::string repr() const;

private:
    bool restrict_;
    std::pair<float, float> abs_eta_range_;
};

// Class to handle custom variations of tau ID vs jet scale factors
class TauIDVsJetVariation {
public:
    // Constructor
    TauIDVsJetVariation(const std::string&);

    // Wrap correction::Correction::evaluate function
    std::function<double (const std::vector<correction::Variable::Type>&)> wrap_evaluate(const correction::Correction*) const;

private:
    // Name of the variation
    std::string variation_;

    // Name of the variation accessed in the correction file
    std::string cfile_variation_;

    // Flag that indicates whether the variation is a custom variation
    bool is_custom_variation_;

    // Restrictions for generator-level match, decay mode, pt, and eta
    GenMatchRestriction gen_match_restriction_;
    DecayModeRestriction decay_mode_restriction_;
    PtRestriction pt_restriction_;
    EtaRestriction eta_restriction_;

    // Match variation to custom variation pattern and extract parameters
    std::pair<bool, std::unordered_map<std::string, std::string>> match_custom_variation(const std::string&) const;

    // Get index of a variable in the list of inputs of a correction::Correction
    size_t get_variable_index(const correction::Correction*, const std::string&, const size_t&) const;

    // Throw an exception for variable with out-of-range index
    void throw_variable_out_of_range(const std::string&, const size_t&) const;
};

} // end namespace scalefactor

} // end namespace tau

} // end namespace physicsobject


#endif // GUARD_TAUVARIATIONS_H
