#ifndef GUARDUTILITY_H
#define GUARDUTILITY_H

#include <cmath>
#include <string>
#include <utility> // make_index_sequence
#include <vector>
/// Namespace used for common utility functions.
namespace utility {
inline bool ApproxEqual(auto value1, auto value2, double maxDelta = 1e-5) {
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
inline void appendParameterPackToVector(std::vector<std::string> &v,
                                        const std::string &parameter) {
    v.push_back(parameter);
}
template <class... ParameterPack>
inline void appendParameterPackToVector(std::vector<std::string> &v,
                                        const std::string &parameter,
                                        const ParameterPack &...pack) {
    v.push_back(parameter);
    appendParameterPackToVector(v, pack...);
}

/// !!!! Remove once we can switch to Root 6.25, where fix is included
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
} // end namespace utility
#endif /* GUARDUTILITY_H */