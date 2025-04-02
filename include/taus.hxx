#ifndef GUARD_TAUS_H
#define GUARD_TAUS_H

namespace physicsobject {
namespace tau {

ROOT::RDF::RNode
PtCorrectionMC_eleFake(ROOT::RDF::RNode df,
                       correctionManager::CorrectionManager &correction_manager,
                       const std::string &corrected_pt, 
                       const std::string &pt,
                       const std::string &eta, 
                       const std::string &decay_mode,
                       const std::string &gen_match, 
                       const std::string &es_file,
                       const std::string &correction_name, 
                       const std::string &id_algorithm,
                       const std::string &variation_dm0_barrel, 
                       const std::string &variation_dm1_barrel,
                       const std::string &variation_dm0_endcap, 
                       const std::string &variation_dm1_endcap);
ROOT::RDF::RNode
PtCorrectionMC_muFake(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correction_manager,
                      const std::string &corrected_pt, 
                      const std::string &pt,
                      const std::string &eta, 
                      const std::string &decay_mode,
                      const std::string &gen_match, 
                      const std::string &es_file,
                      const std::string &correction_name,
                      const std::string &id_algorithm, 
                      const std::string &variation);
ROOT::RDF::RNode
PtCorrectionMC_genuineTau(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correction_manager,
                    const std::string &corrected_pt, 
                    const std::string &pt,
                    const std::string &eta, 
                    const std::string &decay_mode,
                    const std::string &gen_match, 
                    const std::string &es_file,
                    const std::string &correction_name,
                    const std::string &id_algorithm, 
                    const std::string &variation_dm0,
                    const std::string &variation_dm1, 
                    const std::string &variation_dm10,
                    const std::string &variation_dm11);
} // namespace tau
} // namespace physicsobject
#endif /* GUARD_TAUS_H */