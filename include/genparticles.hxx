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

namespace genleptons {
ROOT::RDF::RNode GenLepton(ROOT::RDF::RNode df,
                           const std::string &genparticles_pt,
                           const std::string &genparticles_eta,
                           const std::string &genparticles_phi,
                           const std::string &genparticles_mass,
                           const std::string &genparticles_pdgid,
                           const std::string &genparticles_status,
                           const std::string &genparticles_statusFlag);
ROOT::RDF::RNode GenLeptonPreFSR(ROOT::RDF::RNode df,
                           const std::string &genparticles_pt,
                           const std::string &genparticles_eta,
                           const std::string &genparticles_phi,
                           const std::string &genparticles_mass,
                           const std::string &genparticles_pdgid,
                           const std::string &genparticles_status,
                           const std::string &genparticles_statusFlag);
ROOT::RDF::RNode GenDressedLepton(ROOT::RDF::RNode df,
                           const std::string &genparticles_pt,
                           const std::string &genparticles_eta,
                           const std::string &genparticles_phi,
                           const std::string &genparticles_mass,
                           const std::string &genparticles_pdgid,
                           const std::string &genparticles_hasTauAnc,
                           const std::string &genlep_p4_1,
                           const std::string &genlep_p4_2);
}

namespace genflag {
ROOT::RDF::RNode DYGenFlag(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const int &pdgId);
ROOT::RDF::RNode WGenFlag(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const int &pdgId);
}

namespace genmatching {
namespace tau {
ROOT::RDF::RNode genmatching(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genTaus,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &lepton_p4);
                             
ROOT::RDF::RNode genmatching_wh(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genTaus,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &genparticle_motheridx,
                             const std::string &genparticles_status,
                             const std::string &lepton_p4);

ROOT::RDF::RNode hadronicGenTaus(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &genparticles_pdgid,
                                 const std::string &genparticles_statusFlag,
                                 const std::string &genparticles_motherid);
} // end namespace tau

} // namespace genmatching
#endif /* GUARD_GENPARTICLES_H */