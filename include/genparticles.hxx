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

ROOT::RDF::RNode fatjet_lepton_genmatching(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &lepton_p4);


ROOT::RDF::RNode fatjet_tau_had_genmatching(ROOT::RDF::RNode df, const std::string &outputname,
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

ROOT::RDF::RNode trigger_mu_in_fatjet(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &fatjet_p4,
                                 const std::string &muon_pt,
                                 const std::string &muon_eta,
                                 const std::string &muon_phi,
                                 const std::string &muon_mass);
} // end namespace tau

} // namespace genmatching
#endif /* GUARD_GENPARTICLES_H */