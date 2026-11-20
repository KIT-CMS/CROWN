#pragma once

#include <ROOT/RDataFrame.hxx>
#include <utility/CorrectionManager.hxx>

namespace stxs {
ROOT::RDF::RNode
ggHNNLOWeights(ROOT::RDF::RNode df, const std::string &weight_name,
               const std::string &rootfilename, const std::string &generator,
               const std::string &htxs_pth, const std::string &htxs_njets);
ROOT::RDF::RNode
ggH_WG1_uncertainties(ROOT::RDF::RNode df,
                      const std::vector<std::string> &weight_names,
                      const std::string &htxs_flag, const std::string &htxs_pth,
                      const std::string &htxs_njets);
ROOT::RDF::RNode
qqH_WG1_uncertainties(ROOT::RDF::RNode df,
                      const std::vector<std::string> &weight_names,
                      const std::string &htxs_flag, const size_t &idx = 0);
ROOT::RDF::RNode
Normalization(ROOT::RDF::RNode df,
              correctionManager::CorrectionManager &correction_manager,
              const std::string &family_output, const std::string &norm_output,
              const std::string &stxs_table_path, const std::string &nick,
              const std::string &htxs_bin);
ROOT::RDF::RNode
LHEScale(ROOT::RDF::RNode df,
         correctionManager::CorrectionManager &correction_manager,
         const std::string &output, const std::string &stxs_table_path,
         const std::string &variation, const std::string &nick,
         const std::string &htxs_bin, const std::string &lhe_scale_weights);
} // namespace stxs
