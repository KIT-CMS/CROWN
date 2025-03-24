#ifndef GUARD_VECS_H
#define GUARD_VECS_H

#include "ROOT/RVec.hxx"

namespace vectoroperations {
float calculateMT(ROOT::Math::PtEtaPhiMVector &particle,
                  ROOT::Math::PtEtaPhiMVector &met);

/// Function to combine two RVec Masks by multiplying the two RVec elementwise
///
/// \param mask_1 The first mask
/// \param mask_2 The second mask
///
/// \returns a lambda function which returns the multiplication of the two masks
inline auto MultiplyTwoMasks() {
    return [](const ROOT::RVec<Int_t> &mask_1,
              const ROOT::RVec<Int_t> &mask_2) { return mask_1 * mask_2; };
}
} // end namespace vectoroperations
#endif /* GUARD_VECS_H */