#ifndef GUARD_ELECTRONS_H
#define GUARD_ELECTRONS_H

namespace physicsobject {
namespace electron {

ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &pt,
               const std::string &eta, const std::string &seed_gain,
               const std::string &es_resolution_up,
               const std::string &es_resolution_down,
               const std::string &es_file,
               const std::string &es_name,
               const std::string &era,
               const std::string &variation);  // for Run 2
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &pt,
               const std::string &eta,
               const std::string &delta_eta_sc,
               const std::string &r9,
               const std::string &event_seed,
               const std::string &sf_file,
               const std::string &sf_name,
               const std::string &variation);  // for Run 3
ROOT::RDF::RNode
PtCorrectionData(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correction_manager,
               const std::string &outputname, const std::string &pt,
               const std::string &eta, const std::string &delta_eta_sc,
               const std::string &seed_gain,
               const std::string &r9, const std::string &run,
               const std::string &sf_file,
               const std::string &sf_name);  // for Run 3
ROOT::RDF::RNode VetoECALGap(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &eta,
                             const std::string &delta_eta_sc,
                             const float &end_ecal_barrel,
                             const float &start_ecal_endcap);
ROOT::RDF::RNode
CutInteractionPoint(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &eta, const std::string &delta_eta_sc,
                    const std::string &dxy, const std::string &dz,
                    const float &ecal_barrel_endcap_boundary,
                    const float &max_dxy_barrel, const float &max_dz_barrel,
                    const float &max_dxy_endcap, const float &max_dz_endcap);

namespace scalefactor {
ROOT::RDF::RNode Id(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &outputname, const std::string &pt,
                    const std::string &eta, const std::string &era,
                    const std::string &wp, const std::string &sf_file,
                    const std::string &sf_name, const std::string &variation);
ROOT::RDF::RNode Trigger(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correction_manager,
                    const std::string &outputname, const std::string &pt,
                    const std::string &eta, const std::string &era,
                    const std::string &path_id_name, const std::string &sf_file,
                    const std::string &sf_name, const std::string &variation);
ROOT::RDF::RNode Trigger(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correction_manager,
                    const std::string &outputname, const std::string &pt,
                    const std::string &eta, const std::string &trigger_flag,
                    const std::string &era, const std::string &path_id_name,
                    const std::string &sf_file, const std::string &sf_name,
                    const std::string &variation);
} // end namespace scalefactor
} // end namespace electron
} // end namespace physicsobject
#endif /* GUARD_ELECTRONS_H */
