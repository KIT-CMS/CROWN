#ifndef GUARD_EMBEDDING_H
#define GUARD_EMBEDDING_H

namespace embedding {
namespace electron {

ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap);
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correction_manager,
             const std::string &corrected_pt, const std::string &pt,
             const std::string &eta, const std::string &es_file,
             const std::string &correction_name, 
             const std::string &variation_barrel,
             const std::string &variation_endcap);
} // namespace electron

namespace tau {

ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &decay_mode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11);
} // namespace tau
} // namespace embedding
#endif /* GUARD_EMBEDDING_H */