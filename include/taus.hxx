#ifndef GUARD_TAUS_H
#define GUARD_TAUS_H

namespace physicsobject {
namespace tau {

ROOT::RDF::RNode
PtCorrectionMC_eleFake(ROOT::RDF::RNode df,
                       correctionManager::CorrectionManager &correction_manager,
                       const std::string &outputname, const std::string &pt,
                       const std::string &eta, const std::string &decay_mode,
                       const std::string &gen_match, const std::string &es_file,
                       const std::string &correction_name,
                       const std::string &id_algorithm,
                       const std::string &variation_dm0_barrel,
                       const std::string &variation_dm1_barrel,
                       const std::string &variation_dm0_endcap,
                       const std::string &variation_dm1_endcap);
ROOT::RDF::RNode
PtCorrectionMC_muFake(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correction_manager,
                      const std::string &outputname, const std::string &pt,
                      const std::string &eta, const std::string &decay_mode,
                      const std::string &gen_match, const std::string &es_file,
                      const std::string &correction_name,
                      const std::string &id_algorithm,
                      const std::string &variation);
ROOT::RDF::RNode PtCorrectionMC_genuineTau(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &pt,
    const std::string &eta, const std::string &decay_mode,
    const std::string &gen_match, const std::string &es_file,
    const std::string &correction_name, const std::string &id_algorithm,
    const std::string &variation_dm0, const std::string &variation_dm1,
    const std::string &variation_dm10, const std::string &variation_dm11);

namespace scalefactor {

ROOT::RDF::RNode
Id_vsJet_lt(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correction_manager,
            const std::string &outputname,
            const std::string &pt, const std::string &decay_mode,
            const std::string &gen_match, 
            const std::string &sf_file,
            const std::string &sf_name,
            const std::vector<int> &selected_dms,
            const std::string &wp, const std::string &vsele_wp,
            const std::string &sf_dependence,
            const std::string &sf_vsjet_tau30to35,
            const std::string &sf_vsjet_tau35to40,
            const std::string &sf_vsjet_tau40to500,
            const std::string &sf_vsjet_tau500to1000,
            const std::string &sf_vsjet_tau1000toinf);
ROOT::RDF::RNode Id_vsJet_tt(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, 
    const std::string &pt, const std::string &decay_mode,
    const std::string &gen_match, 
    const std::string &sf_file, const std::string &sf_name,
    const std::vector<int> &selected_dms,
    const std::string &wp, const std::string &vsele_wp,
    const std::string &sf_dependence,
    const std::string &sf_vsjet_tauDM0,
    const std::string &sf_vsjet_tauDM1, 
    const std::string &sf_vsjet_tauDM10,
    const std::string &sf_vsjet_tauDM11);
ROOT::RDF::RNode
Id_vsEle(ROOT::RDF::RNode df,
         correctionManager::CorrectionManager &correction_manager,
         const std::string &outputname,
         const std::string &eta,
         const std::string &gen_match, 
         const std::string &sf_file, const std::string &sf_name,
         const std::string &wp, 
         const std::string &sf_vsele_barrel,
         const std::string &sf_vsele_endcap);

ROOT::RDF::RNode
Id_vsMu(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &eta,
        const std::string &gen_match, 
        const std::string &sf_file,
        const std::string &sf_name,
        const std::string &wp,  
        const std::string &sf_vsmu_wheel1,
        const std::string &sf_vsmu_wheel2, 
        const std::string &sf_vsmu_wheel3,
        const std::string &sf_vsmu_wheel4, 
        const std::string &sf_vsmu_wheel5);
ROOT::RDF::RNode
Trigger(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &pt, const std::string &decay_mode, 
        const std::string &sf_file,
        const std::string &sf_name,
        const std::string &trigger_name, const std::string &wp,
        const std::string &corr_type, const std::string &variation);
} // end namespace scalefactor
} // end namespace tau
} // end namespace physicsobject
#endif /* GUARD_TAUS_H */