#ifndef GUARDBASEFUCTIONS_H
#define GUARDBASEFUCTIONS_H

#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <nlohmann/json.hpp>

enum Channel { MT = 0, ET = 1, TT = 2, EM = 3 };

/// Namespace used for common basefunctions. Theses functions return a lambda
/// function to be used in a data frame define

namespace basefunctions {

// /**
//  * @brief Function to filter events based on their run and luminosity block
//  * values
//  *
//  * @param df the dataframe to filter
//  * @param json_path the path to the golden json file containing all valid
//  runs
//  * ans luminosity blocks
//  * @param run th column containing the run value
//  * @param luminosity the column containing the luminosity block value
//  * @param filtername the name of the filter
//  * @return a filtered dataframe
//  */
// auto JSONFilter(auto &df, const std::string &json_path, const std::string
// &run,
//                 const std::string &luminosity, const std::string &filtername)
//                 {
//     std::ifstream i(json_path);
//     nlohmann::json golden_json;
//     i >> golden_json;
//     auto jsonFilterlambda = [golden_json](UInt_t run, UInt_t luminosity) {
//         bool matched = false;
//         // check if the run exists
//         if (golden_json.find(std::to_string(run)) != golden_json.end()) {
//             // now loop over all luminosity blocks and check if the event is
//             // valid
//             for (auto &luminosityrange : golden_json[std::to_string(run)]) {
//                 if (luminosity >= luminosityrange[0] &&
//                     luminosity <= luminosityrange[1]) {
//                     matched = true;
//                     break;
//                 }
//             }
//             if (!matched) {
//                 Logger::get("JSONFilter")
//                     ->debug("Run {} / luminosity {} not in json file", run,
//                             luminosity);
//             }
//         }
//         return matched;
//     };
//     return df.Filter(jsonFilterlambda, {run, luminosity}, filtername);
// }

/// Function to add an input quantity under a different name
///
/// \param df the dataframe to add the quantity to
/// \param inputname name of the existing column
/// \param outputname name of the new column
///
/// \returns a dataframe with the new column

template <typename T>
ROOT::RDF::RNode rename(auto &df, const std::string &inputname,
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
ROOT::RDF::RNode DefineQuantity(auto &df, const std::string &outputname,
                                T const &value) {
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
ROOT::RDF::RNode FilterFlagsAny(auto &df, const std::string &filtername,
                                const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Filter(
        utility::PassAsVec<nFlags, bool>(
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
ROOT::RDF::RNode CombineFlagsAny(auto &df, const std::string &outputflag,
                                 const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Define(
        outputflag,
        utility::PassAsVec<nFlags, bool>(
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
ROOT::RDF::RNode FilterIntSelection(auto &df, const std::string &quantity,
                                    const std::vector<T> &selection,
                                    const std::string &filtername) {
    return df.Filter(
        [selection](const T probe) {
            return std::find(selection.begin(), selection.end(), probe) !=
                   selection.end();
        },
        {quantity}, filtername);
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
ROOT::RDF::RNode
evaluateWorkspaceFunction(auto &df, const std::string &outputname,
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
ROOT::RDF::RNode UnrollVectorQuantity(auto &df, const std::string &name,
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
