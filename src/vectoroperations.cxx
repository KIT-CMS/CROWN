#ifndef GUARD_VECS_H
#define GUARD_VECS_H

#include <Math/Vector4D.h>
namespace vectoroperations {
/**
 * @brief calculate the transverse mass. The transverse mass is defined as:
 \f[
    m_{\text{T}} = \sqrt{2 \, p_{T} \, \text{E}_{T} \,
 (1-\cos\Delta\phi)}
 \f]
 where \f$p_{T}\f$ is the transverse momentum of particle \f$p\f$,
 \f$\text{E}_{T}\f$ is the missing transverse energy, and
 \f$\cos\Delta\phi\f$ is the cosine of the angle between the particle and
 missing transverse energy.
 *
 * @param particle lorentz vector of the particle
 * @param met lorentz vector of the missing transverse energy
 * @return the transverse mass of the particle
 */
float calculateMT(ROOT::Math::PtEtaPhiMVector &particle,
                  ROOT::Math::PtEtaPhiMVector &met) {
    return (float)sqrt(2 * particle.Pt() * met.Pt() *
                       (1. - cos(particle.Phi() - met.Phi())));
}
} // end namespace vectoroperations
#endif /* GUARD_VECS_H */