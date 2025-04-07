#ifndef GUARD_UTILITY_H
#define GUARD_UTILITY_H

#include <cmath>
#include <string>
#include <utility> // make_index_sequence
#include <vector>

#include "../../include/utility/Logger.hxx"
#include "../../include/utility/RooFunctorThreadsafe.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

/// Namespace used for common utility functions.
namespace utility {

/**
 * @brief Function to check if two double values are approximately equal
 *
 * @param value1 first double value to compare
 * @param value2 second double value to compare
 * @param maxDelta maximum difference between the two values
 * @return true or false
 */
inline bool ApproxEqual(double value1, double value2, double maxDelta = 1e-5) {
    if (value1 == value2) {
        return true;
    } else {
        double delta = std::abs(value1 - value2);
        if ((value1 + value2) != 0) {
            delta *= (2.0 / std::abs(value1 + value2));
        }
        return (delta < maxDelta);
    }
}

/**
 * @brief This function extracts and returns the last element from a tuple 
 * containing a variable number of arguments. The function uses `std::get` 
 * to access the last element of the tuple, allowing you to retrieve the 
 * last argument without knowing its type or index beforehand.
 *
 * @param args input tuple containing the arguments
 * 
 * @return last element of the tuple
 */
template <typename... Args>
constexpr auto extractLastArgument(const std::tuple<Args...>& args) {
    constexpr size_t N = sizeof...(Args);
    // Extract last argument
    return std::get<N - 1>(args);  
}

/**
 * @brief Function to append a parameter pack to a vector
 *
 * @param v the vector to append to
 * @param parameter string to append
 */
inline void appendParameterPackToVector(std::vector<std::string> &v,
                                        const std::string &parameter) {
    v.push_back(parameter);
}

/**
 * @brief Function to append a parameter pack to a vector
 *
 * @tparam ParameterPack
 * @param v the vector to append to
 * @param parameter the string to append
 * @param pack the parameter pack to append
 */
template <class... ParameterPack>
inline void appendParameterPackToVector(std::vector<std::string> &v,
                                        const std::string &parameter,
                                        const ParameterPack &...pack) {
    v.push_back(parameter);
    appendParameterPackToVector(v, pack...);
}

/// \cond
template <typename I, typename T, typename F> class PassAsVecHelper;

template <std::size_t... N, typename T, typename F>
class PassAsVecHelper<std::index_sequence<N...>, T, F> {
    template <std::size_t Idx> using AlwaysT = T;
    typename std::decay<F>::type fFunc;

  public:
    PassAsVecHelper(F &&f) : fFunc(std::forward<F>(f)) {}
    auto operator()(AlwaysT<N>... args) -> decltype(fFunc({args...})) {
        return fFunc({args...});
    }
};

template <std::size_t N, typename T, typename F>
auto PassAsVec(F &&f) -> PassAsVecHelper<std::make_index_sequence<N>, T, F> {
    return utility::PassAsVecHelper<std::make_index_sequence<N>, T, F>(
        std::forward<F>(f));
}
/// \endcond

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
    appendParameterPackToVector(InputList, inputs...);
    const auto nInputs = sizeof...(Inputs);
    Logger::get("EvaluateWorkspaceFunction")->debug("nInputs: {} ", nInputs);
    auto df1 = df.Define(
        outputname, PassAsVec<nInputs, float>(getValue), InputList);
    return df1;
}
} // end namespace utility
#endif /* GUARD_UTILITY_H */