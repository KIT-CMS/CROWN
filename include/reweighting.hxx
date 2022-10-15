#ifndef GUARD_REWEIGHTING_H
#define GUARD_REWEIGHTING_H

namespace reweighting {
ROOT::RDF::RNode puweights(ROOT::RDF::RNode df, const std::string &weightname,
                           const std::string &truePUMean,
                           const std::string &filename,
                           const std::string &histogramname);
ROOT::RDF::RNode puweights(ROOT::RDF::RNode df, const std::string &weightname,
                           const std::string &truePU,
                           const std::string &filename,
                           const std::string &eraname,
                           const std::string &variation);
ROOT::RDF::RNode topptreweighting(ROOT::RDF::RNode df,
                                  const std::string &weightname,
                                  const std::string &gen_pdgids,
                                  const std::string &gen_status,
                                  const std::string &gen_pt);
ROOT::RDF::RNode zPtMassReweighting(ROOT::RDF::RNode df,
                                    const std::string &weightname,
                                    const std::string &gen_boson,
                                    const std::string &workspace_file,
                                    const std::string &functor_name,
                                    const std::string &argset);
} // namespace reweighting
#endif /* GUARD_REWEIGHTING_H */