#ifndef GUARD_BASEFUCTIONS_H
#define GUARD_BASEFUCTIONS_H

#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

namespace basefunctions {

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

#endif /* GUARD_BASEFUNCTIONS_H */
