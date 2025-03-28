#ifndef GUARDUTILITY_H
#define GUARDUTILITY_H

#include <cmath>
#include <string>
#include <utility> // make_index_sequence
#include <vector>

#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
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
} // end namespace utility
#endif /* GUARDUTILITY_H */