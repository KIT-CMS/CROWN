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
 * The variation string `"up_custom_dm10_pt20to40"` would translate into the the
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
TauIDVsJetVariation::TauIDVsJetVariation(const std::string &variation) {
    // Store the input variation's name
    custom_variation_ = variation;

    // Define regular expression that catches custom variation definitions
    auto custom_pattern = std::regex(
        "(up|down)_custom(_dm(0|1|10|11))?(_pt(\\d+)to(\\d+))?",
        std::regex_constants::ECMAScript
    );
    Logger::get("TauIDVsJetVariation")->debug(
        "Parsing tau ID vs jet variation: {}", variation
    );

    std::smatch matches;
    if (std::regex_match(variation, matches, custom_pattern)) {
        // Set variation in the correction file and DM and pt selection values
        // if regex pattern has been matched

        // Capture matched groups
        std::string direction = matches[1].str();
        std::string dm_str = matches[3].str();
        std::string pt_min_str = matches[5].str();
        std::string pt_max_str = matches[6].str();

        // If corresponding groups have been matched, set the selection flags to
        // true
        has_dm_selection_ = matches[2].matched;
        has_pt_selection_ = matches[4].matched;

        // Raise an exception if both optional groups are missing
        if (!has_dm_selection_ && !has_pt_selection_) {
            auto msg = std::format(
                "Invalid custom variation string: {}. The string must contain "
                "at least one of the optional groups for decay mode or pt "
                "selection.", variation
            );
            throw std::invalid_argument(msg);
        }

        // Set private member variables based on matched groups or default
        // values
        variation_ = direction;
        if (has_dm_selection_) {
            decay_mode_ = std::stoi(dm_str);
        }
        if (has_pt_selection_) {
            pt_min_ = std::stof(pt_min_str);
            pt_max_ = std::stof(pt_max_str);
        }
    } else {
        // If the custom_pattern regular expression is not matched, set the
        // variation to the input string and do not restrict decay mode or pt
        // range (i.e., set default values). In this case, the variation will
        // be uses as-is on the correction file.

        // Set main attributes
        variation_ = variation;

        // Do not impose selections on decay mode or pt range
        has_dm_selection_ = false;
        has_pt_selection_ = false;
    }

    // Final debug output to show the stored values of the variation and
    // selections
    Logger::get("TauIDVsJetVariation")->debug(
        "Set up variation with the following values:"
    );
    Logger::get("TauIDVsJetVariation")->debug(
        "  correction file variation: {}", variation_);
    Logger::get("TauIDVsJetVariation")->debug(
        "  decay mode selection:      {}", has_dm_selection_);
    if (has_dm_selection_) {
        Logger::get("TauIDVsJetVariation")->debug(
        "  decay mode selected:       {}", decay_mode_);
    }
    Logger::get("TauIDVsJetVariation")->debug(
        "  pt selection:              {}", has_pt_selection_);
    if (has_pt_selection_) {
        Logger::get("TauIDVsJetVariation")->debug(
        "  pt range selected:         [{}, {})", pt_min_, pt_max_);
    }
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
    // Get indices of pt and decay mode in the evaluate function inputs
    size_t pt_index = get_variable_index(evaluator, "pt");
    size_t dm_index = get_variable_index(evaluator, "dm");
    size_t variation_index = get_variable_index(evaluator, "syst");

    // Capture selection flags and values for the wrapper function
    auto custom_variation = custom_variation_;
    auto has_dm_selection = has_dm_selection_;
    auto has_pt_selection = has_pt_selection_;
    auto sel_decay_mode = decay_mode_;
    auto sel_pt_min = pt_min_;
    auto sel_pt_max = pt_max_;

    // Capture correction evaluation to evaluate the correction file scale
    // factors 
    const std::string correction_variation = variation_;

    // Define the evaluate wrapper function
    auto wrapper = [
        evaluator,
        pt_index,
        dm_index,
        variation_index,
        custom_variation,
        has_dm_selection,
        has_pt_selection,
        sel_decay_mode,
        sel_pt_min,
        sel_pt_max,
        correction_variation
    ] (const std::vector<correction::Variable::Type> &values) {
        // If no custom selections are imposed, just evaluate the correction
        // factor using the provided values
        if (!has_dm_selection && !has_pt_selection) {
            Logger::get("TauIDVsJetVariation")->debug("Default evaluation of correction");
            return evaluator->evaluate(values);
        }

        // For custom selections, the input values for the evaluate function
        // need to be manipulated manually.
        Logger::get("TauIDVsJetVariation")->debug(
            "Custom evaluation of correction for custom variation {}",
            custom_variation);

        // Get pt and decay mode values
        auto pt = std::get<double>(values[pt_index]);
        auto decay_mode = std::get<int>(values[dm_index]);

        Logger::get("TauIDVsJetVariation")->debug(
            "Checking selections for pt {}, decay mode {}",
            pt, decay_mode
        );

        // Check whether the event passes the decay mode and pt selections
        auto dm_selected = decay_mode == sel_decay_mode;
        auto pt_selected = pt >= sel_pt_min && pt < sel_pt_max;
        auto selected = (
            (has_dm_selection && dm_selected) || !has_dm_selection
        ) && (
            (has_pt_selection && pt_selected) || !has_pt_selection
        );
        Logger::get("TauIDVsJetVariation")->debug(
            "Selection results for custom variation");
        Logger::get("TauIDVsJetVariation")->debug(
            "  decay mode selection {}", dm_selected
        );
        Logger::get("TauIDVsJetVariation")->debug(
            "  pt selection         {}", pt_selected
        );
        Logger::get("TauIDVsJetVariation")->debug(
            "  overall selection    {}", selected
        );

        // If the event is marked as selected, evaluate the correction
        // factor with the correct variation direction. If the selection is not
        // passed, evaluate with the nominal variation
        std::vector<correction::Variable::Type> values_copy = values;
        if (selected) {
            values_copy[variation_index] = correction_variation;
        } else {
            values_copy[variation_index] = "nom";
        }
        Logger::get("TauIDVsJetVariation")->debug(
            "Evaluating correction with variation {}",
            std::get<std::string>(values_copy[variation_index])
        );

        return evaluator->evaluate(values_copy);
    };

    return wrapper;
}

// --- private -----------------------------------------------------------------

size_t TauIDVsJetVariation::get_variable_index(
    const correction::Correction *evaluator, const std::string &name
) const {
    // Go through the list of the evaluator's inputs and find the index of the
    // variable with the given name
    for (size_t i = 0; i < evaluator->inputs().size(); ++i) {
        if (evaluator->inputs()[i].name() == name) {
            return i;
        }
    }

    // Raise an exception if the variable has not been found
    auto msg = std::format(
        "Variable name {} not found in evaluator inputs", name
    );
    throw std::out_of_range(msg);
}

} // end namespace scalefactor

} // end namespace tau

} // end namespace physicsobject
