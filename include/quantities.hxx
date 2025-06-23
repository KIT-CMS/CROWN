#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "defaults.hxx"
#include "utility/Logger.hxx"
#include "vectoroperations.hxx"
#include <Math/Vector4D.h>
namespace quantities {
ROOT::RDF::RNode scalarPtSum(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &pt_1, const std::string &pt_2,
                             const std::string &pt_3);
ROOT::RDF::RNode m_vis(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pt_W(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pt_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode deltaPhi(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &p_1_p4, const std::string &p_2_p4);
ROOT::RDF::RNode deltaPhi_WH(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &p_1_p4,
                             const std::string &p_2_p4,
                             const std::string &p_3_p4);
ROOT::RDF::RNode
p4_fastmtt(ROOT::RDF::RNode df, const std::string &outputname,
           const std::string &pt_1, const std::string &pt_2,
           const std::string &eta_1, const std::string &eta_2,
           const std::string &phi_1, const std::string &phi_2,
           const std::string &mass_1, const std::string &mass_2,
           const std::string &met_pt, const std::string &met_phi,
           const std::string &met_cov_xx, const std::string &met_cov_xy,
           const std::string &met_cov_yy, const std::string &decay_mode_1,
           const std::string &decay_mode_2, const std::string &finalstate);
ROOT::RDF::RNode pt_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors);
ROOT::RDF::RNode pzetamissvis(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &p_1_p4,
                              const std::string &p_2_p4,
                              const std::string &met);
ROOT::RDF::RNode mTdileptonMET(ROOT::RDF::RNode df,
                               const std::string &outputname,
                               const std::string &p_1_p4,
                               const std::string &p_2_p4,
                               const std::string &met);
ROOT::RDF::RNode deltaR(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4);
ROOT::RDF::RNode mT(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met);
ROOT::RDF::RNode pt_tt(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &p_1_p4, const std::string &p_2_p4,
                       const std::string &met);
ROOT::RDF::RNode pt_ttjj(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::string &p_1_p4, const std::string &p_2_p4,
                         const std::string &jet_1_p4,
                         const std::string &jet_2_p4, const std::string &met);
ROOT::RDF::RNode pt_dijet(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &jet_1_p4,
                          const std::string &jet_2_p4);
ROOT::RDF::RNode jet_hemisphere(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &jet_1_p4,
                                const std::string &jet_2_p4);
ROOT::RDF::RNode mt_tot(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4,
                        const std::string &met);
ROOT::RDF::RNode mtt_coll_approx(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4,
                        const std::string &met);
namespace tau {
ROOT::RDF::RNode matching_jet_pt(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const int &position,
                                 const std::string &pairname,
                                 const std::string &taujet_index,
                                 const std::string &jetpt_column);
ROOT::RDF::RNode matching_genjet_pt(
    ROOT::RDF::RNode df, const std::string &outputname, const int &position,
    const std::string &pairname, const std::string &taujet_index,
    const std::string &genjet_index, const std::string &genjetpt_column);
ROOT::RDF::RNode TauIDFlag(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &nameID, const int &idxID);
} // end namespace tau
} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */
