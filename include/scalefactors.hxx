#ifndef GUARD_SCALEFACTORS_H
#define GUARD_SCALEFACTORS_H

namespace scalefactor {
namespace muon {
ROOT::RDF::RNode id_rooworkspace(auto &df, const std::string &pt,
                                 const std::string &eta,
                                 const std::string &id_output,
                                 const std::string &workspace_name,
                                 const std::string &id_functor_name,
                                 const std::string &id_arguments);
ROOT::RDF::RNode iso_rooworkspace(auto &df, const std::string &pt,
                                  const std::string &eta,
                                  const std::string &iso,
                                  const std::string &iso_output,
                                  const std::string &workspace_name,
                                  const std::string &iso_functor_name,
                                  const std::string &iso_arguments);
ROOT::RDF::RNode id(auto &df, const std::string &pt, const std::string &eta,
                    const std::string &year_id, const std::string &variation,
                    const std::string &id_output, const std::string &sf_file,
                    const std::string &idAlgorithm);
ROOT::RDF::RNode iso(auto &df, const std::string &pt, const std::string &eta,
                     const std::string &year_id, const std::string &variation,
                     const std::string &iso_output, const std::string &sf_file,
                     const std::string &idAlgorithm);
} // namespace muon
namespace tau {

ROOT::RDF::RNode
id_vsJet_lt(auto &df, const std::string &pt, const std::string &decayMode,
            const std::string &genMatch, const std::vector<int> &selectedDMs,
            const std::string &wp, const std::string &sf_vsjet_tau30to35,
            const std::string &sf_vsjet_tau35to40,
            const std::string &sf_vsjet_tau40to500,
            const std::string &sf_vsjet_tau500to1000,
            const std::string &sf_vsjet_tau1000toinf,
            const std::string &sf_dependence, const std::string &id_output,
            const std::string &sf_file, const std::string &idAlgorithm);
ROOT::RDF::RNode
id_vsJet_tt(auto &df, const std::string &pt, const std::string &decayMode,
            const std::string &genMatch, const std::vector<int> &selectedDMs,
            const std::string &wp, const std::string &sf_vsjet_tauDM0,
            const std::string &sf_vsjet_tauDM1,
            const std::string &sf_vsjet_tauDM10,
            const std::string &sf_vsjet_tauDM11,
            const std::string &sf_dependence, const std::string &id_output,
            const std::string &sf_file, const std::string &idAlgorithm);
ROOT::RDF::RNode
id_vsEle(auto &df, const std::string &eta, const std::string &decayMode,
         const std::string &genMatch, const std::vector<int> &selectedDMs,
         const std::string &wp, const std::string &sf_vsele_barrel,
         const std::string &sf_vsele_endcap, const std::string &id_output,
         const std::string &sf_file, const std::string &idAlgorithm);
ROOT::RDF::RNode
id_vsMu(auto &df, const std::string &eta, const std::string &decayMode,
        const std::string &genMatch, const std::vector<int> &selectedDMs,
        const std::string &wp, const std::string &sf_vsmu_wheel1,
        const std::string &sf_vsmu_wheel2, const std::string &sf_vsmu_wheel3,
        const std::string &sf_vsmu_wheel4, const std::string &sf_vsmu_wheel5,
        const std::string &id_output, const std::string &sf_file,
        const std::string &idAlgorithm);
} // namespace tau

namespace electron {

ROOT::RDF::RNode id(auto &df, const std::string &pt, const std::string &eta,
                    const std::string &year_id, const std::string &wp,
                    const std::string &variation, const std::string &id_output,
                    const std::string &sf_file, const std::string &idAlgorithm);
} // namespace electron
} // namespace scalefactor
#endif /* GUARD_SCALEFACTORS_H */