#ifndef GUARDLVECS_H
#define GUARDLVECS_H

namespace lorentzvectors {
ROOT::RDF::RNode buildparticle(ROOT::RDF::RNode df,
                               const std::vector<std::string> &quantities,
                               const std::string &outputname,
                               const int &position);
ROOT::RDF::RNode build(ROOT::RDF::RNode df,
                       const std::vector<std::string> &obj_quantities,
                       const int pairindex, const std::string &obj_p4_name);
ROOT::RDF::RNode buildMet(ROOT::RDF::RNode df, const std::string &met_pt,
                          const std::string &met_phi,
                          const std::string &outputname);

/// namespace used for mutau lorentzvectors
namespace mutau {
ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &pairname,
                       const std::vector<std::string> &muon_quantities,
                       const std::vector<std::string> &tau_quantities,
                       const std::string &muon_p4_name,
                       const std::string &tau_p4_name);
} // namespace mutau
} // namespace lorentzvectors
#endif /* GUARDLVECS_H */