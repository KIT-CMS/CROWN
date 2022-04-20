#ifndef GUARD_REWEIGHTING_H
#define GUARD_REWEIGHTING_H

namespace reweighting {
ROOT::RDF::RNode puweights(auto &df, const std::string &weightname,
                           const std::string &truePUMean,
                           const std::string &filename,
                           const std::string &histogramname);
ROOT::RDF::RNode topptreweighting(auto &df, const std::string &weightname,
                                  const std::string &gen_pdgids,
                                  const std::string &gen_status,
                                  const std::string &gen_pt);
ROOT::RDF::RNode zPtMassReweighting(auto &df, const std::string &weightname,
                                    const std::string &gen_boson,
                                    const std::string &workspace_file,
                                    const std::string &functor_name,
                                    const std::string &argset);
} // namespace reweighting
#endif /* GUARD_REWEIGHTING_H */