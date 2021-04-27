#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

enum Channel { MT = 0, ET = 1, TT = 2, EM = 3 };

namespace basefunctions {

auto FilterMax(float cut) {
  return [cut](const ROOT::RVec<float> &values) {
    ROOT::RVec<int> mask = values < cut;
    return mask;
  };
}

auto FilterMin(float cut) {
  // As in ROOT, for min we use >=
  return [cut](const ROOT::RVec<float> &values) {
    ROOT::RVec<int> mask = values >= cut;
    return mask;
  };
}

auto FilterAbsMax(float cut) {
  return [cut](const ROOT::RVec<float> &values) {
    ROOT::RVec<int> mask = abs(values) < cut;
    return mask;
  };
}

auto FilterAbsMin(float cut) {
  return [cut](const ROOT::RVec<float> &values) {
    ROOT::RVec<int> mask = abs(values) >= cut;
    return mask;
  };
}

auto MultiplyTwoMasks() {
  return [](const ROOT::RVec<Int_t> &mask_1, const ROOT::RVec<Int_t> &mask_2) {
    return mask_1 * mask_2;
  };
}

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