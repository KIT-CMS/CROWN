#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "basefunctions.hxx"
#include "defaults.hxx"
#include "utility/Logger.hxx"
#include "vectoroperations.hxx"
#include <Math/Vector4D.h>
namespace quantities {
ROOT::RDF::RNode pt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &inputvector);
ROOT::RDF::RNode eta(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector);
ROOT::RDF::RNode phi(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector);
ROOT::RDF::RNode mass(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &inputvector);
ROOT::RDF::RNode dxy(ROOT::RDF::RNode df, const std::string &outputname,
                     const int &position, const std::string &pairname,
                     const std::string &dxycolumn);
ROOT::RDF::RNode dz(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &dzcolumn);
ROOT::RDF::RNode charge(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &position, const std::string &pairname,
                        const std::string &chargecolumn);
ROOT::RDF::RNode m_vis(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::vector<std::string> &inputvectors);
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
ROOT::RDF::RNode mt_tot(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4,
                        const std::string &met);
ROOT::RDF::RNode isolation(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &isolationcolumn);
ROOT::RDF::RNode pdgid(ROOT::RDF::RNode df, const std::string &outputname,
                       const int &position, const std::string &pairname,
                       const std::string &pdgidcolumn);
ROOT::RDF::RNode NumberOfGoodLeptons(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const std::string &goodleptons);
namespace tau {
ROOT::RDF::RNode decaymode(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &decaymodecolumn);
ROOT::RDF::RNode genmatch(ROOT::RDF::RNode df, const std::string &outputname,
                          const int &position, const std::string &pairname,
                          const std::string &genmatchcolumn);
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
namespace muon {
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &idcolumn);
ROOT::RDF::RNode is_global(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &globalflagcolumn);

ROOT::RDF::RNode rndm(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &objCollection,
    int seed
);

ROOT::RDF::RNode applyRoccoRData(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &filename,
    const int &position,
    const std::string &objCollection,
    const std::string &chargColumn,
    const std::string &ptColumn,
    const std::string &etaColumn,
    const std::string &phiColumn,
    int error_set,
    int error_member
);

ROOT::RDF::RNode applyRoccoRMC(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &filename,
    const int &position,
    const std::string &objCollection,
    const std::string &chargColumn,
    const std::string &ptColumn,
    const std::string &etaColumn,
    const std::string &phiColumn,
    const std::string &genPtColumn,
    const std::string &nTrackerLayersColumn,
    const std::string &rndmColumn,
    int error_set,
    int error_member
);

} // namespace muon
} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */