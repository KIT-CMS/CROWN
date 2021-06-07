#ifndef GUARDBASEFUCTIONS_H
#define GUARDBASEFUCTIONS_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "RooFunctor.h"
#include "utility/Logger.hxx"
#include "utility/utility.hxx"

enum Channel { MT = 0, ET = 1, TT = 2, EM = 3 };

/// Namespace used for common basefunctions. Theses functions return a lambda
/// function to be used in a data frame define

namespace basefunctions {

/// Function to apply a maximal filter requirement to a quantity.
/// Returns true if the value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
auto FilterMax(float cut) {
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
auto FilterAbsMax(float cut) {
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
auto FilterMin(float cut) {
    // As in ROOT, for min we use >=
    return [cut](const ROOT::RVec<float> &values) {
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
auto FilterAbsMin(float cut) {
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
auto FilterID(int index) {
    return [index](const ROOT::RVec<UChar_t> &IDs) {
        ROOT::RVec<int> mask;
        for (auto const ID : IDs) {
            mask.push_back(std::min(1, int(ID & 1 << index - 1)));
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
auto FilterJetID(int index) {
    return [index](const ROOT::RVec<Int_t> &IDs) {
        ROOT::RVec<int> mask;
        for (auto const ID : IDs) {
            mask.push_back(std::min(1, (ID >> index) & 1));
        }
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
auto evaluateWorkspaceFunction(auto df, const std::string &outputname,
                               const std::shared_ptr<RooFunctor> function,
                               const Inputs &... inputs) {
    Logger::get("evaluateWorkspaceFunction")
        ->debug("Starting evaluation for {}", outputname);
    auto getValue = [function](const ROOT::RVec<float> &values) {
        Logger::get("evaluateWorkspaceFunction")
            ->debug("Type: {} // nPar {} // nObs {}", typeid(function).name(),
                    function->nPar(), function->nObs());
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
} // namespace basefunctions

#endif /* GUARDBASEFUNCTIONS_H */
