#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "defaults.hxx"
#include "utility/Logger.hxx"
#include <Math/Vector4D.h>

namespace quantities {
ROOT::RDF::RNode DeltaPhi(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &vector_1, const std::string &vector_2);
ROOT::RDF::RNode DeltaEta(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &vector_1, const std::string &vector_2);
ROOT::RDF::RNode DeltaR(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2);
ROOT::RDF::RNode PairHemisphere(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &vector_1,
                                const std::string &vector_2);
ROOT::RDF::RNode PzetaMissVis(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &vector_1,
                              const std::string &vector_2,
                              const std::string &vector_3);
ROOT::RDF::RNode TransverseMass(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2);
ROOT::RDF::RNode TransverseMass(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2,
                        const std::string &vector_3);
ROOT::RDF::RNode CollinearApproxMtt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2,
                        const std::string &vector_3);
ROOT::RDF::RNode
FastMtt(ROOT::RDF::RNode df, const std::string &outputname,
           const std::string &pt_1, const std::string &pt_2,
           const std::string &eta_1, const std::string &eta_2,
           const std::string &phi_1, const std::string &phi_2,
           const std::string &mass_1, const std::string &mass_2,
           const std::string &met_pt, const std::string &met_phi,
           const std::string &met_cov_xx, const std::string &met_cov_xy,
           const std::string &met_cov_yy, const std::string &decay_mode_1,
           const std::string &decay_mode_2, const std::string &finalstate);
ROOT::RDF::RNode JetMatching(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &jet_quantity,
                            const std::string &object_jet_index,
                            const std::string &object_index_vector,
                            const int &position);
ROOT::RDF::RNode GenJetMatching(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &genjet_quantity,
                            const std::string &jet_genjet_index,
                            const std::string &object_jet_index,
                            const std::string &object_index_vector,
                            const int &position);
ROOT::RDF::RNode deltaPhi_WH(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &vector_1,
                             const std::string &vector_2,
                             const std::string &vector_3);
ROOT::RDF::RNode pt_W(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::vector<std::string> &vectors);
} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */
