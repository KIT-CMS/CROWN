/// Namespace used for common utility functions.

#include <cstddef>
#include <utility>

namespace utility {
bool ApproxEqual(auto value1, auto value2, double maxDelta = 1e-5) {
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

template <std::size_t ...Idxs>
auto helper(std::index_sequence<Idxs...>) {
    return [] (decltype(Idxs, 0)...masks) { return (masks *...); };
}

template <std::size_t N>
auto MaskCombiner() {
    return helper(std::make_index_sequence<N>{});
};


} // end namespace utility
