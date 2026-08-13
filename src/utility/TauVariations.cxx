#include <regex>
#include <format>
#include <vector>
#include <functional>
#include "correction.h"

#include "../../include/utility/Logger.hxx"
#include "../../include/utility/TauVariations.hxx"


namespace physicsobject {

namespace tau {

namespace scalefactor {

// -----------------------------------------------------------------------------
// physicsobject::tau::scalefactor::GenMatchRestriction
// -----------------------------------------------------------------------------

GenMatchRestriction::GenMatchRestriction() : restrict_(false), gen_matches_({}) {}

GenMatchRestriction::GenMatchRestriction(const std::vector<int> &gen_matches) : restrict_(true), gen_matches_(gen_matches) {
    // Validate that only allowed generator-level match values are passed
    for (const auto &gen_match : gen_matches_) {
        if (
            gen_match != 1 && gen_match != 2 && gen_match != 3
            && gen_match != 4 && gen_match != 5
        ) {
            auto msg = std::format(
                "Invalid generator-level match value: {}. Allowed values are "
                "1, 2, 3, 4, and 5.", gen_match
            );
            throw std::invalid_argument(msg);
        }
    }
}

bool GenMatchRestriction::is_active() const {
    return restrict_;
}

bool GenMatchRestriction::is_selected(const int &gen_match) const {
    if (!restrict_) {
        return true;
    }
    return (
        std::find(gen_matches_.begin(), gen_matches_.end(), gen_match) != gen_matches_.end()
    );
}

std::string GenMatchRestriction::repr() const {
    std::string joined_gen_matches = "";
    for (size_t i = 0; i < gen_matches_.size(); ++i) {
        if (i > 0) {
            joined_gen_matches += ", ";
        }
        joined_gen_matches += std::to_string(gen_matches_[i]);
    }
    return std::format(
        "GenMatchRestriction(is_active={}, gen_matches=[{}])",
        restrict_, joined_gen_matches
    );
}

// -----------------------------------------------------------------------------
// physicsobject::tau::scalefactor::DecayModeRestriction
// -----------------------------------------------------------------------------

DecayModeRestriction::DecayModeRestriction() : restrict_(false), decay_modes_(std::vector<int>()) {}

DecayModeRestriction::DecayModeRestriction(const std::vector<int> &decay_modes) : restrict_(true), decay_modes_(decay_modes) {
    // Validate that only allowed decay mode values are passed
    for (const auto &decay_mode : decay_modes_) {
        if (
            decay_mode != 0 && decay_mode != 1 && decay_mode != 10
            && decay_mode != 11
        ) {
            auto msg = std::format(
                "Invalid decay mode value: {}. Allowed values are 0, 1, 10, and "
                "11.", decay_mode
            );
            throw std::invalid_argument(msg);
        }
    }
}

DecayModeRestriction::DecayModeRestriction(const int &decay_mode)
    : restrict_(true), decay_modes_({decay_mode}) {
    // Validate that only allowed decay mode values are passed
    if (
        decay_mode != 0 && decay_mode != 1 && decay_mode != 10
        && decay_mode != 11
    ) {
        auto msg = std::format(
            "Invalid decay mode value: {}. Allowed values are 0, 1, 10, and "
            "11.", decay_mode
        );
        throw std::invalid_argument(msg);
    }
}

bool DecayModeRestriction::is_active() const {
    return restrict_;
}

bool DecayModeRestriction::is_selected(const int &decay_mode) const {
    if (!restrict_) {
        return true;
    }
    return (
        std::find(decay_modes_.begin(), decay_modes_.end(), decay_mode)
        != decay_modes_.end()
    );
}

std::string DecayModeRestriction::repr() const {
    std::string joined_decay_modes = "";
    for (size_t i = 0; i < decay_modes_.size(); ++i) {
        if (i > 0) {
            joined_decay_modes += ", ";
        }
        joined_decay_modes += std::to_string(decay_modes_[i]);
    }
    return std::format(
        "DecayModeRestriction(is_active={}, decay_modes=[{}])",
        restrict_, joined_decay_modes
    );
}

// -----------------------------------------------------------------------------
// physicsobject::tau::scalefactor::PtRestriction
// -----------------------------------------------------------------------------

PtRestriction::PtRestriction() : restrict_(false), pt_range_(std::make_pair(-10.f, -10.f)) {}

PtRestriction::PtRestriction(const float &pt_min, const float &pt_max) : restrict_(true), pt_range_(std::make_pair(pt_min, pt_max)) {
    if (pt_range_.first >= pt_range_.second) {
        auto msg = std::format(
            "Invalid pt range: [{}, {}). The lower bound must be smaller than "
            "the upper bound.", pt_range_.first, pt_range_.second
        );
        throw std::invalid_argument(msg);
    }
}

PtRestriction::PtRestriction(const float &pt_min)
    : restrict_(true), pt_range_(std::make_pair(pt_min, std::numeric_limits<float>::infinity())) {
    // No validation needed since infinity is always greater than any finite pt_min
}

bool PtRestriction::is_active() const {
    return restrict_;
}

bool PtRestriction::is_selected(const float &pt) const {
    if (!restrict_) {
        return true;
    }
    return pt >= pt_range_.first && pt < pt_range_.second;
}

std::string PtRestriction::repr() const {
    return std::format(
        "PtRestriction(is_active={}, pt_range=[{}, {}))",
        restrict_, pt_range_.first, pt_range_.second
    );
}

// -----------------------------------------------------------------------------
// physicsobject::tau::scalefactor::EtaRestriction
// -----------------------------------------------------------------------------

EtaRestriction::EtaRestriction() : restrict_(false), abs_eta_range_(std::make_pair(-10.f, -10.f)) {}

EtaRestriction::EtaRestriction(const std::pair<float, float> &abs_eta_range) : restrict_(true), abs_eta_range_(abs_eta_range) {
    if (abs_eta_range_.first >= abs_eta_range_.second) {
        auto msg = std::format(
            "Invalid eta range: [{}, {}). The lower bound must be smaller than "
            "the upper bound.", abs_eta_range_.first, abs_eta_range_.second
        );
        throw std::invalid_argument(msg);
    }
}

EtaRestriction::EtaRestriction(const float &abs_eta_min, const float &abs_eta_max)
    : restrict_(true), abs_eta_range_(std::make_pair(abs_eta_min, abs_eta_max)) {
    if (abs_eta_range_.first >= abs_eta_range_.second) {
        auto msg = std::format(
            "Invalid eta range: [{}, {}). The lower bound must be smaller than "
            "the upper bound.", abs_eta_range_.first, abs_eta_range_.second
        );
        throw std::invalid_argument(msg);
    }
}

bool EtaRestriction::is_active() const {
    return restrict_;
}

bool EtaRestriction::is_selected(const float &eta) const {
    if (!restrict_) {
        return true;
    }
    return abs(eta) >= abs_eta_range_.first && abs(eta) < abs_eta_range_.second;
}

std::string EtaRestriction::repr() const {
    return std::format(
        "EtaRestriction(is_active={}, abs_eta_range=[{}, {}))",
        restrict_, abs_eta_range_.first, abs_eta_range_.second
    );
}

// ----------------------------------------------------------------------------
// physicsobject::tau::scalefactor::TauIDVsJetVariation
// ----------------------------------------------------------------------------

// --- public ------------------------------------------------------------------

/**
 * @brief Construct a new `TauIDVsJetVariation`.
 *
 * The object is constructed by parsing the `variation` string. If the string
 * matches a pattern for a custom tau ID variation, selections are parsed from
 * the variation and stored in the object. The variation passed to the
 * correction is just the direction of the variation, either "up" or "down",
 * while a selection on `pt` and `eta` is imposed. If the selection is not
 * passed, the nominal value (obtained with the `"nom"` variation) is used.
 
 * If the string does not match the pattern, the variation is stored as-is. 
 * In this case, the `"variation"` must correspond to an available value in the
 * correction file used to evaluate the scale factor. The value is then used for
 * the evaluation of the scale factor. No selections are imposed.
 *
 * The `variation` string is matched with the regular expression
 * `"(up|down)_custom(_dm(0|1|10|11))?(_pt(\\d+)to(\\d+))?"`.
 *
 * - the first group captures the direction of the variation, either "up" or
 *   "down"
 *
 * - the second group captures the optional decay mode selection. Allowed values
 *   for the decay mode are 0, 1, 10, and 11. If this group is not matched, no
 *   decay mode selection takes place.
 *
 * - the third group captures the optional \f$p_{\text{T}}\f$ selection. The
 *   lower and upper value of the considered \f$p_{\text{T}}\f$ range are
 *   captured from the values before and after `"to"`. These numbers must
 *   represent unsigned integers. If the group is not matched, no
 *   \f$p_{\text{T}}\f$ selection takes place.
 *
 * Example:
 *
 * The variation string `"up_custom_dm10_pt20to40"` translates into the
 * following pseudocode for the evaluation of the scale factor:
 *
 * ```cpp
 * double pt, eta;
 * int decay_mode, gen_match;
 * std::string wp, vsele_wp, sf_dependence;
 * double sf;
 * if (pt >= 20 && pt < 40 && decay_mode == 10) {
 *     sf = correction->evaluate(pt, decay_mode, gen_match, wp, vsele_wp, "up", sf_dependence);
 * else {
 *     sf = correction->evaluate(pt, decay_mode, gen_match, wp, vsele_wp, "nom", sf_dependence);
 * }
 * ```
 *
 * For non-matching variation string `"down_syst_alleras", the scale factor is
 * directly evaluated with the correction's evaluate function:
 *
 * ```cpp
 * double pt, eta;
 * int decay_mode, gen_match;
 * std::string wp, vsele_wp, sf_dependence;
 * double sf;
 * sf = correction->evaluate(pt, decay_mode, gen_match, wp, vsele_wp, "down_syst_alleras", sf_dependence);
 * ```
 *
 * @param variation Name of the tau ID vs jets scale factor variation
 */
TauIDVsJetVariation::TauIDVsJetVariation(const std::string &variation) :
    gen_match_restriction_(GenMatchRestriction()),
    decay_mode_restriction_(DecayModeRestriction()),
    pt_restriction_(PtRestriction()),
    eta_restriction_(EtaRestriction()) {

    // Set the variation name
    variation_ = variation;
    Logger::get("TauIDVsJetVariation")->debug(
        "Handling tau ID vs jet variation {}", variation_
    );

    // Match the variation to the custom variation pattern
    auto r = match_custom_variation(variation);
    bool is_custom = r.first;
    std::unordered_map<std::string, std::string> match_results = r.second;

    // Set attribute flag that indicates whether the variation is a custom
    // variation
    is_custom_variation_ = is_custom;

    if (is_custom) {
        // The variation in the correction file is just the total shift in the
        // direction declared in the custom variation string
        cfile_variation_ = match_results["direction"];

        // Set generator-level match restriction if matched in the custom
        // variation string
        if (match_results.find("gen_match") != match_results.end()) {
            std::string gen_object_str = match_results["gen_match"];
            auto gen_matches = std::vector<int>();
            if (gen_object_str == "genEle") {
                gen_matches = GEN_TYPES.at(GenType::GEN_ELE);
            } else if (gen_object_str == "genMu") {
                gen_matches = GEN_TYPES.at(GenType::GEN_MU);
            } else if (gen_object_str == "genTau") {
                gen_matches = GEN_TYPES.at(GenType::GEN_TAU);
            }
            gen_match_restriction_ = GenMatchRestriction(gen_matches);
        }

        // Set decay mode restriction if matched in the custom variation string
        if (match_results.find("decay_mode") != match_results.end()) {
            decay_mode_restriction_ = DecayModeRestriction(
                std::stoi(match_results["decay_mode"])
            );
        }

        // Set pt restriction if matched in the custom variation string
        if (
            match_results.find("pt_min") != match_results.end()
            && match_results.find("pt_max") != match_results.end()
        ) {
            float pt_min = std::stof(match_results["pt_min"]);
            if (match_results["pt_max"] == "Inf") {
                pt_restriction_ = PtRestriction(pt_min);
            } else {
                float pt_max = std::stof(match_results["pt_max"]);
                pt_restriction_ = PtRestriction(pt_min, pt_max);
            }
        }

        // Set the eta restriction
        if (match_results.find("eta_region") != match_results.end()) {
            auto eta_region_str = match_results["eta_region"];
            auto abs_eta_range = std::pair<float, float>();
            if (eta_region_str == "barrel") {
                abs_eta_range = ETA_REGIONS.at(EtaRegion::BARREL);
            } else if (eta_region_str == "endcap") {
                abs_eta_range = ETA_REGIONS.at(EtaRegion::ENDCAP);
            } else if (eta_region_str == "wheel1") {
                abs_eta_range = ETA_REGIONS.at(EtaRegion::WHEEL_1);
            } else if (eta_region_str == "wheel2") {
                abs_eta_range = ETA_REGIONS.at(EtaRegion::WHEEL_2);
            } else if (eta_region_str == "wheel3") {
                abs_eta_range = ETA_REGIONS.at(EtaRegion::WHEEL_3);
            } else if (eta_region_str == "wheel4") {
                abs_eta_range = ETA_REGIONS.at(EtaRegion::WHEEL_4);
            } else if (eta_region_str == "wheel5") {
                abs_eta_range = ETA_REGIONS.at(EtaRegion::WHEEL_5);
            }
            eta_restriction_ = EtaRestriction(abs_eta_range);
        }

        // Raise an exception if all optional groups are missing
        if (
            !gen_match_restriction_.is_active()
            && !decay_mode_restriction_.is_active()
            && !pt_restriction_.is_active()
            && !eta_restriction_.is_active()
        ) {
            auto msg = std::format(
                "Invalid custom variation string: {}. The string must contain "
                "at least one of the optional groups for generator-level "
                "match, decay mode, pt, or eta selection.", variation
            );
            throw std::invalid_argument(msg);
        }

    } else {
        // The variation name is assumed to be accessible from the correction
        // file. Therefore, both variation names are set to the same value and
        // no restrictions are defined.
        cfile_variation_ = variation;
    }

    // Final debug output to show the stored values of the variation and
    // selections
    Logger::get("TauIDVsJetVariation")->debug(
        "Set up variation handler with the following values:"
    );
    Logger::get("TauIDVsJetVariation")->debug(
        "  correction file variation: {}", cfile_variation_);

    Logger::get("TauIDVsJetVariation")->debug(
        "  gen match restriction:     {}", gen_match_restriction_.repr());
    Logger::get("TauIDVsJetVariation")->debug(
        "  decay mode restriction:    {}", decay_mode_restriction_.repr());
    Logger::get("TauIDVsJetVariation")->debug(
        "  pt restriction:            {}", pt_restriction_.repr());
    Logger::get("TauIDVsJetVariation")->debug(
        "  eta restriction:           {}", eta_restriction_.repr());
}

/**
 * @brief Wrap the `correction::Correction::evaluate` function to allow for
 * custom tau ID vs jet scale factor variations.
 *
 * The returned wrapper function has the same signature as the
 * `correction::Correction::evaluate` method.
 *
 * @param evaluator Pointer to the `correction::Correction` object used to
 * evaluate
 *
 * @return Function that takes a list of inputs and returns the scale factor.
 */
std::function<double (const std::vector<correction::Variable::Type>&)>
TauIDVsJetVariation::wrap_evaluate(const correction::Correction *evaluator)
const {
    // Get indices of the variables in the list of inputs of the correction::Correction
    size_t gen_index = get_variable_index(evaluator, "genmatch", -1);
    size_t decay_mode_index = get_variable_index(evaluator, "dm", -1);
    size_t pt_index = get_variable_index(evaluator, "pt", -1);
    size_t eta_index = get_variable_index(evaluator, "eta", -1);
    size_t syst_index = get_variable_index(evaluator, "syst", -1);

    // Define the evaluate wrapper function
    auto wrapper = [
        is_custom_variation = this->is_custom_variation_,
        cfile_variation = this->cfile_variation_,
        gen_match_restriction = this->gen_match_restriction_,
        decay_mode_restriction = this->decay_mode_restriction_,
        pt_restriction = this->pt_restriction_,
        eta_restriction = this->eta_restriction_,
        evaluator,
        gen_index,
        decay_mode_index,
        pt_index,
        eta_index,
        syst_index
    ] (const std::vector<correction::Variable::Type> &values) {
        // If no custom selections are imposed, just evaluate the correction
        // factor using the provided values
        if (!is_custom_variation) {
            Logger::get("TauIDVsJetVariation")->debug("Default evaluation of correction");
            return evaluator->evaluate(values);
        }

        // For custom selections, the input values for the evaluate function
        // need to be manipulated manually.
        Logger::get("TauIDVsJetVariation")->debug(
            "Custom evaluation of correction for custom variation");

        // Set default values for selection inputs
        int gen_match = -10;
        int decay_mode = -10;
        double pt = -10.;
        double eta = -10.;

        // Set the values if they are needed for imposed restrictions
        if (gen_match_restriction.is_active()) {
            if (gen_index == static_cast<size_t>(-1)) {
                auto msg = std::format(
                    "Variable genmatch not found in the list of inputs of the correction. "
                    "Please check the variable name and ensure it is present in the "
                    "correction inputs."
                );
                throw std::out_of_range(msg);
            }
            gen_match = std::get<int>(values[gen_index]);
        }
        if (decay_mode_restriction.is_active()) {
            if (decay_mode_index == static_cast<size_t>(-1)) {
                auto msg = std::format(
                    "Variable dm not found in the list of inputs of the correction. "
                    "Please check the variable name and ensure it is present in the "
                    "correction inputs."
                );
                throw std::out_of_range(msg);
            }
            decay_mode = std::get<int>(values[decay_mode_index]);
        }
        if (pt_restriction.is_active()) {
            if (pt_index == static_cast<size_t>(-1)) {
                auto msg = std::format(
                    "Variable pt not found in the list of inputs of the correction. "
                    "Please check the variable name and ensure it is present in the "
                    "correction inputs."
                );
                throw std::out_of_range(msg);
            }
            pt = std::get<double>(values[pt_index]);
        }
        if (eta_restriction.is_active()) {
            if (eta_index == static_cast<size_t>(-1)) {
                auto msg = std::format(
                    "Variable eta not found in the list of inputs of the correction. "
                    "Please check the variable name and ensure it is present in the "
                    "correction inputs."
                );
                throw std::out_of_range(msg);
            }
            eta = std::get<double>(values[eta_index]);
        }

        // Print debug output for the values of the selection inputs
        Logger::get("TauIDVsJetVariation")->debug(
            "Checking selections for");
        Logger::get("TauIDVsJetVariation")->debug(
            "  gen_match        {}", gen_match);
        Logger::get("TauIDVsJetVariation")->debug(
            "  decay_mode       {}", decay_mode);
        Logger::get("TauIDVsJetVariation")->debug(
            "  pt               {}", pt);
        Logger::get("TauIDVsJetVariation")->debug(
            "  eta              {}", eta);

        // Check whether the event passes selections if restrictions are imposed
        auto is_selected = (
            gen_match_restriction.is_selected(gen_match)
            && decay_mode_restriction.is_selected(decay_mode)
            && pt_restriction.is_selected(pt)
            && eta_restriction.is_selected(eta)
        );

        Logger::get("TauIDVsJetVariation")->debug(
            "Selection results in selection status {}", is_selected);

        // If the event is marked as selected, evaluate the correction
        // factor with the correct variation direction. If the selection is not
        // passed, evaluate with the nominal variation
        std::vector<correction::Variable::Type> values_copy = values;
        if (is_selected) {
            values_copy[syst_index] = cfile_variation;
        } else {
            values_copy[syst_index] = "nom";
        }
        Logger::get("TauIDVsJetVariation")->debug(
            "Evaluating correction with variation {}",
            std::get<std::string>(values_copy[syst_index])
        );

        return evaluator->evaluate(values_copy);
    };

    return wrapper;
}

// --- private -----------------------------------------------------------------

std::pair<bool, std::unordered_map<std::string, std::string>> TauIDVsJetVariation::match_custom_variation(
    const std::string &custom_variation
) const {
    // Define regular expression that catches custom variation definitions
    auto custom_pattern = std::regex(
        "(up|down)_custom(_(genEle|genMu|genTau))?(_dm(0|1|10|11))?(_pt(\\d+)to(\\d+|Inf))?(_(barrel|endcap|wheel[1-5]))?",
        std::regex_constants::ECMAScript
    );
    Logger::get("TauIDVsJetVariation")->debug(
        "Parsing tau ID vs jet variation: {}", custom_variation
    );

    bool matched = false;
    std::unordered_map<std::string, std::string> results;
    std::smatch matches;
    if (std::regex_match(custom_variation, matches, custom_pattern)) {
        // Set flag that the regex pattern has been matched
        matched = true;

        // Capture matched groups if they are not empty, add them to the results
        // map
        results.insert({"direction", matches[1].str()});
        if (!matches[2].str().empty()) {
            results.insert({"gen_match", matches[3].str()});
        }
        if (!matches[4].str().empty()) {
            results.insert({"decay_mode", matches[5].str()});
        }
        if (!matches[6].str().empty()) {
            results.insert({"pt_min", matches[7].str()});
            results.insert({"pt_max", matches[8].str()});
        }
        if (!matches[9].str().empty()) {
            results.insert({"eta_region", matches[10].str()});
        }
    }

    return std::make_pair(matched, results);
}

size_t TauIDVsJetVariation::get_variable_index(
    const correction::Correction *evaluator, const std::string &name, const size_t &default_index
) const {
    // Go through the list of the evaluator's inputs and find the index of the
    // variable with the given name
    for (size_t i = 0; i < evaluator->inputs().size(); ++i) {
        if (evaluator->inputs()[i].name() == name) {
            return i;
        }
    }

    // If the variable is not found, return the default index
    return default_index;
}

void TauIDVsJetVariation::throw_variable_out_of_range(const std::string &name, const size_t &index) const {
    if (index == -1) {
        auto msg = std::format(
            "Variable {} not found in the list of inputs of the correction. "
            "Please check the variable name and ensure it is present in the "
            "correction inputs.", name
        );
        throw std::out_of_range(msg);
    }
}

} // end namespace scalefactor

} // end namespace tau

} // end namespace physicsobject
