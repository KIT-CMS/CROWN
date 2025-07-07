#ifndef GUARD_GENPARTICLES_H
#define GUARD_GENPARTICLES_H

#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "bitset"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>

namespace genparticles {
namespace tau {

ROOT::RDF::RNode HadronicGenTaus(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &genparticles_pdg_id,
                                 const std::string &genparticles_status_flags,
                                const std::string &genparticles_mother_index);
ROOT::RDF::RNode GenMatching(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &hadronic_gen_taus,
                             const std::string &genparticles_pdg_id,
                             const std::string &genparticles_status_flags,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &reco_had_tau);
ROOT::RDF::RNode GenMatching(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &hadronic_gen_taus, const std::string &genparticles_pdg_id,
    const std::string &genparticles_status_flags, const std::string &genparticle_motheridx,
    const std::string &genparticles_pt, const std::string &genparticles_eta, 
    const std::string &genparticles_phi, const std::string &genparticles_mass, 
    const std::string &reco_had_tau);
} // end namespace tau
} // end namespace genparticles
#endif /* GUARD_GENPARTICLES_H */
