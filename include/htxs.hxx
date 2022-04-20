#ifndef GUARDHTXS_H
#define GUARDHTXS_H

namespace htxs {
ROOT::RDF::RNode ggHNLLOWeights(auto &df, const std::string &weight_name,
                                const std::string &rootfilename,
                                const std::string &generator,
                                const std::string &htxs_pth,
                                const std::string &htxs_njets);
ROOT::RDF::RNode
ggH_WG1_uncertainties(auto &df, const std::vector<std::string> &weight_names,
                      const std::string &htxs_flag, const std::string &htxs_pth,
                      const std::string &htxs_njets);
ROOT::RDF::RNode
qqH_WG1_uncertainties(auto &df, const std::vector<std::string> &weight_names,
                      const std::string &htxs_flag, const size_t &idx = 0);
} // namespace htxs
#endif /* GUARDHTXS_H */