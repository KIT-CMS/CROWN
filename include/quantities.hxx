#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "basefunctions.hxx"
#include "defaults.hxx"
#include "utility/Logger.hxx"
#include "vectoroperations.hxx"
#include <Math/Vector4D.h>
namespace quantities {
ROOT::RDF::RNode pt(auto &df, const std::string &outputname,
                    const std::string &inputvector);
ROOT::RDF::RNode eta(auto &df, const std::string &outputname,
                     const std::string &inputvector);
ROOT::RDF::RNode phi(auto &df, const std::string &outputname,
                     const std::string &inputvector);
ROOT::RDF::RNode mass(auto &df, const std::string &outputname,
                      const std::string &inputvector);
ROOT::RDF::RNode dxy(auto &df, const std::string &outputname,
                     const int &position, const std::string &pairname,
                     const std::string &dxycolumn);
ROOT::RDF::RNode dz(auto &df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &dzcolumn);
ROOT::RDF::RNode charge(auto &df, const std::string &outputname,
                        const int &position, const std::string &pairname,
                        const std::string &chargecolumn);
ROOT::RDF::RNode m_vis(auto &df, const std::string &outputname,
                       const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pt_vis(auto &df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pzetamissvis(auto &df, const std::string &outputname,
                              const std::string &p_1_p4,
                              const std::string &p_2_p4,
                              const std::string &met);
ROOT::RDF::RNode mTdileptonMET(auto &df, const std::string &outputname,
                               const std::string &p_1_p4,
                               const std::string &p_2_p4,
                               const std::string &met);
ROOT::RDF::RNode mT(auto &df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met);
ROOT::RDF::RNode pt_tt(auto &df, const std::string &outputname,
                       const std::string &p_1_p4, const std::string &p_2_p4,
                       const std::string &met);
ROOT::RDF::RNode pt_ttjj(auto &df, const std::string &outputname,
                         const std::string &p_1_p4, const std::string &p_2_p4,
                         const std::string &jet_1_p4,
                         const std::string &jet_2_p4, const std::string &met);
ROOT::RDF::RNode mt_tot(auto &df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4,
                        const std::string &met);
ROOT::RDF::RNode isolation(auto &df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &isolationcolumn);
ROOT::RDF::RNode pdgid(auto &df, const std::string &outputname,
                       const int &position, const std::string &pairname,
                       const std::string &pdgidcolumn);
ROOT::RDF::RNode NumberOfGoodLeptons(auto &df, const std::string &outputname,
                                     const std::string &goodleptons);
namespace tau {
ROOT::RDF::RNode decaymode(auto &df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &decaymodecolumn);
ROOT::RDF::RNode genmatch(auto &df, const std::string &outputname,
                          const int &position, const std::string &pairname,
                          const std::string &genmatchcolumn);
ROOT::RDF::RNode matching_jet_pt(auto &df, const std::string &outputname,
                                 const int &position,
                                 const std::string &pairname,
                                 const std::string &taujet_index,
                                 const std::string &jetpt_column);
ROOT::RDF::RNode matching_genjet_pt(auto &df, const std::string &outputname,
                                    const int &position,
                                    const std::string &pairname,
                                    const std::string &taujet_index,
                                    const std::string &genjet_index,
                                    const std::string &genjetpt_column);
ROOT::RDF::RNode TauIDFlag(auto &df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &nameID, const int &idxID);
} // end namespace tau
} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */