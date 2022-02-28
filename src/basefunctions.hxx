#ifndef GUARDBASEFUCTIONS_H
#define GUARDBASEFUCTIONS_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/Logger.hxx"
#include "utility/RooFunctorThreadsafe.hxx"
#include "utility/utility.hxx"
#include <nlohmann/json.hpp>

enum Channel { MT = 0, ET = 1, TT = 2, EM = 3 };

/// Namespace used for common basefunctions. Theses functions return a lambda
/// function to be used in a data frame define

namespace basefunctions {

/**
 * @brief Function to filter events based on their run and luminosity block
 * values
 *
 * @param df the dataframe to filter
 * @param json_path the path to the golden json file containing all valid runs
 * ans luminosity blocks
 * @param run th column containing the run value
 * @param luminosity the column containing the luminosity block value
 * @param filtername the name of the filter
 * @return a filtered dataframe
 */
auto JSONFilter(auto &df, const std::string &json_path, const std::string &run,
                const std::string &luminosity, const std::string &filtername) {
    std::ifstream i(json_path);
    nlohmann::json golden_json;
    i >> golden_json;
    auto jsonFilterlambda = [golden_json](UInt_t run, UInt_t luminosity) {
        bool matched = false;
        // check if the run exists
        if (golden_json.find(std::to_string(run)) != golden_json.end()) {
            // now loop over all luminosity blocks and check if the event is
            // valid
            for (auto &luminosityrange : golden_json[std::to_string(run)]) {
                if (luminosity >= luminosityrange[0] &&
                    luminosity <= luminosityrange[1]) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                Logger::get("JSONFilter")
                    ->debug("Run {} / luminosity {} not in json file", run,
                            luminosity);
            }
        }
        return matched;
    };
    return df.Filter(jsonFilterlambda, {run, luminosity}, filtername);
}

/// Function to add an input quantity under a different name
///
/// \param df the dataframe to add the quantity to
/// \param inputname name of the existing column
/// \param outputname name of the new column
///
/// \returns a dataframe with the new column

template <typename T>
auto rename(auto &df, const std::string &inputname,
            const std::string &outputname) {
    return df.Define(outputname, [](const T &q) { return q; }, {inputname});
}
/// Function to add a new quantity with a defined value
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column
/// \param value the value to be added
///
/// \returns a dataframe with the new column

template <typename T>
auto DefineQuantity(auto &df, const std::string &outputname, T const &value) {
    return df.Define(outputname, [value]() { return value; }, {});
}
/// This function filters events, where neither of the input flags is true.
/// This is used to filter events which do not pass an underlying requirement in
/// any systematic variation.
///
/// \param df The input dataframe
/// \param filtername The name of the filter, used in the Dataframe report
/// \param flags Parameter pack of column names that contain the considered
/// flags of type bool
///
/// \returns a filtered dataframe
template <class... Flags>
auto FilterFlagsAny(auto &df, const std::string &filtername,
                    const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Filter(
        ROOT::RDF::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return Any(flags); }),
        FlagList, filtername);
}

/// This function defines a flag being true if either of the input flags is
/// true.
///
/// \param df The input dataframe
/// \param outputflag The name of the new column
/// \param flags Parameter pack of column names that contain the considered
/// flags of type bool
///
/// \returns a dataframe containing the new column
template <class... Flags>
auto CombineFlagsAny(auto &df, const std::string &outputflag,
                     const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Define(
        outputflag,
        ROOT::RDF::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return Any(flags); }),
        FlagList);
}

/// Function to apply a selection on an integer quantity.
/// Returns true if the value is smaller than the given cut value
///
/// \param df The input dataframe
/// \param quantity The quantity the selection refers to
/// \param selection The values of the quantity that are accepted
/// \param filtername The name of the filter, used in the Dataframe report
///
/// \returns a filtered dataframe
template <typename T>
auto FilterIntSelection(auto &df, const std::string &quantity,
                        const std::vector<T> &selection,
                        const std::string &filtername) {
    return df.Filter(
        [selection](const T probe) {
            return std::find(selection.begin(), selection.end(), probe) !=
                   selection.end();
        },
        {quantity}, filtername);
}

/// Function to apply a maximal filter requirement to a quantity.
/// Returns true if the value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
auto FilterMax(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = values < cut;
        return mask;
    };
}

/// Function to apply a maximal filter requirement to a quantity.
/// Returns true if the absolute value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
auto FilterAbsMax(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = abs(values) < cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to a quantity.
/// Returns true if the value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
auto FilterMin(const float &cut) {
    // As in ROOT, for min we use >=
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = values >= cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to an integer quantity.
/// Returns true if the value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
auto FilterMinInt(const int &cut) {
    // As in ROOT, for min we use >=
    return [cut](const ROOT::RVec<int> &values) {
        ROOT::RVec<int> mask = values >= cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to a quantity.
/// Returns true if the absolute value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
auto FilterAbsMin(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = abs(values) >= cut;
        return mask;
    };
}

/// Function to combine two RVec Masks by multiplying the two RVec elementwise
///
/// \param mask_1 The first mask
/// \param mask_2 The second mask
///
/// \returns a lambda function which returns the multiplication of the two masks
auto MultiplyTwoMasks() {
    return [](const ROOT::RVec<Int_t> &mask_1,
              const ROOT::RVec<Int_t> &mask_2) { return mask_1 * mask_2; };
}

/// Function to filter the Tau ID in NanoAOD. The disciminator output is stored
/// as a bitmask. In order to filter based on a given WP, a bitwise comparison
/// between the given ID and the filtered index is performed. Return true if the
/// working point is passed by a tau.
///
/// \param index The bitmask index to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
auto FilterID(const int &index) {
    return [index](const ROOT::RVec<UChar_t> &IDs) {
        ROOT::RVec<int> mask;
        for (auto const ID : IDs) {
            if (index > 0)
                mask.push_back(std::min(1, int(ID & 1 << index - 1)));
            else
                mask.push_back(int(1));
        }
        return mask;
    };
}
/// Function to filter the Jet ID in NanoAOD. Ths values are stored bitwise, so
/// the integer value has to be decoded into binary, and then the value of the
/// index bit has to be compared.
///
/// \param index The bitmask index to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
auto FilterJetID(const int &index) {
    return [index](const ROOT::RVec<Int_t> &IDs) {
        ROOT::RVec<int> mask;
        for (auto const ID : IDs) {
            mask.push_back(std::min(1, (ID >> index) & 1));
        }
        return mask;
    };
}
/// Function to filter the Jet pileup ID in NanoAOD. This ID is applied on jets
/// below a given pt threshold. The jet pileup ID has 4 possible values:
/// 0==fail, 4==pass loose, 6==pass loose and medium, 7==pass loose, medium and
/// tight
///
/// \param PUindex The index to be used for comparison
/// \param PUptcut The pt threshold to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
auto FilterJetPUID(const int &PUindex, const float &PUptcut) {
    return [PUindex, PUptcut](const ROOT::RVec<Int_t> &PUIDs,
                              const ROOT::RVec<float> &jet_pts) {
        ROOT::RVec<int> tmp_mask1 = PUIDs >= PUindex;
        ROOT::RVec<int> tmp_mask2 = jet_pts >= PUptcut;
        ROOT::RVec<int> mask = (tmp_mask1 + tmp_mask2) > 0;
        return mask;
    };
}
/// Function to evaluate a `RooWorkspace` function and put the output into a new
/// dataframe column
///
/// \param[in] df The dataframe, where the new column should be added
/// \param[in] outputname name of the new column
/// \param[in] function A `RooFunctor` pointer, which has to be loaded from a
/// Roo Workspace \param[in] inputs a paramter pack containing all column names
/// needed to be able to evaluate the workspace function
///
/// \returns a dataframe with the newly defined output column
template <class... Inputs>
auto evaluateWorkspaceFunction(
    auto &df, const std::string &outputname,
    const std::shared_ptr<RooFunctorThreadsafe> &function,
    const Inputs &...inputs) {
    Logger::get("evaluateWorkspaceFunction")
        ->debug("Starting evaluation for {}", outputname);
    auto getValue = [function](const ROOT::RVec<float> &values) {
        Logger::get("evaluateWorkspaceFunction")
            ->debug("Type: {} ", typeid(function).name());
        std::vector<double> argvalues(values.begin(), values.end());
        auto result = function->eval(argvalues.data());
        Logger::get("evaluateWorkspaceFunction")->debug("result {}", result);
        return result;
    };
    std::vector<std::string> InputList;
    utility::appendParameterPackToVector(InputList, inputs...);
    const auto nInputs = sizeof...(Inputs);
    Logger::get("evaluateWorkspaceFunction")->debug("nInputs: {} ", nInputs);
    auto df1 = df.Define(
        outputname, utility::PassAsVec<nInputs, float>(getValue), InputList);
    // change back to ROOT::RDF as soon as fix is available
    return df1;
}
/// Helper function to recursively define columns for each entry of a vector
/// quantity
///
/// \param df Input dataframe
/// \param name name of the vector quantity
/// \param names vector of names for the new columns
/// \param idx index of the current recursion loop, should not be set outside
/// this function
///
/// \returns a lambda function to be used in RDF Define
template <typename T>
auto UnrollVectorQuantity(auto &df, const std::string &name,
                          const std::vector<std::string> &names,
                          const size_t &idx = 0) {
    if (idx >= names.size()) {
        return df;
    }
    auto df1 = df.Define(
        names.at(idx),
        [idx](const std::vector<T> &quantities) { return quantities.at(idx); },
        {name});
    return UnrollVectorQuantity<T>(df1, name, names, idx + 1);
}

} // namespace basefunctions

#endif /* GUARDBASEFUNCTIONS_H */
