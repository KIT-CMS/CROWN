/// Namespace used for common utility functions.
#include <cmath>
#include <vector>
#include <string>
#include <utility> // make_index_sequence

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
void appendParameterPackToVector(std::vector<std::string> &v,
                                 const std::string &parameter) {
    v.push_back(parameter);
}
template <class... ParameterPack>
void appendParameterPackToVector(std::vector<std::string> &v,
                                 const std::string &parameter,
                                 const ParameterPack &...pack) {
    v.push_back(parameter);
    appendParameterPackToVector(v, pack...);
}
}