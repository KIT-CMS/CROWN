#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

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
/// \param index The bitmask index to be used for comparion
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

} // namespace basefunctions