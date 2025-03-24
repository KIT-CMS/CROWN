#ifndef GUARDBASEFUCTIONS_H
#define GUARDBASEFUCTIONS_H

#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

namespace basefunctions {

/**
 * @brief Function to add an input column as a new column under a different
 * name.
 *
 * @param df input dataframe
 * @param outputname name of the new column
 * @param inputname name of the existing column
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Rename(ROOT::RDF::RNode df,
                               const std::string &outputname,
                               const std::string &inputname) {
    return df.Define(outputname, [](const T &q) { return q; }, {inputname});
}

/**
 * @brief Function to add a new column with a defined `value` of type `T`. The type
 * needs to be specified when calling this function with
 * `basefunctions::DefineQuantity<T>(...)`.
 *
 * @param df input dataframe
 * @param outputname name of the new column
 * @param value the value to be added
 *
 * @returns a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode DefineQuantity(ROOT::RDF::RNode df,
                                       const std::string &outputname,
                                       T const &value) {
    return df.Define(outputname, [value]() { return value; }, {});
}

/**
 * @brief This function creates a new column `outputname` with the negatives of
 * the values in the column `inputname`.
 *
 * Note that this function is implemented as a template, so specify the type `T`
 * of the objects in the input column when calling this function with
 * `basefunctions::Negative<T>(...)`. The type can also be a vector with
 * `ROOT::VecOps::RVec<T>`.
 *
 * @param df input dataframe
 * @param outputname name of the output column
 * @param inputname column name of the input column
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Negative(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &inputname) {
    auto negative = [](const T &input) { return -input; };
    return df.Define(outputname, negative, {inputname});
}

/**
 * @brief Function to evaluate a `RooWorkspace` function and put the output into
 * a new dataframe column.
 *
 * @param df input dataframe
 * @param outputname name of the new column
 * @param function a `RooFunctor` pointer, which has to be loaded from a Roo
 * Workspace
 * @param inputs a parameter pack containing all column names needed to be able
 * to evaluate the workspace function
 *
 * @return a dataframe with the new column
 */
template <class... Inputs>
inline ROOT::RDF::RNode
EvaluateWorkspaceFunction(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::shared_ptr<RooFunctorThreadsafe> &function,
                          const Inputs &...inputs) {
    Logger::get("EvaluateWorkspaceFunction")
        ->debug("Starting evaluation for {}", outputname);
    auto getValue = [function](const ROOT::RVec<float> &values) {
        Logger::get("EvaluateWorkspaceFunction")
            ->debug("Type: {} ", typeid(function).name());
        std::vector<double> argvalues(values.begin(), values.end());
        auto result = function->eval(argvalues.data());
        Logger::get("EvaluateWorkspaceFunction")->debug("result {}", result);
        return result;
    };
    std::vector<std::string> InputList;
    utility::appendParameterPackToVector(InputList, inputs...);
    const auto nInputs = sizeof...(Inputs);
    Logger::get("EvaluateWorkspaceFunction")->debug("nInputs: {} ", nInputs);
    auto df1 = df.Define(
        outputname, utility::PassAsVec<nInputs, float>(getValue), InputList);
    return df1;
}
} // namespace basefunctions

#endif /* GUARDBASEFUNCTIONS_H */
