#ifndef GUARD_SCALEFACTORS_H
#define GUARD_SCALEFACTORS_H

namespace scalefactor {
namespace tau {

ROOT::RDF::RNode
id_vsJet_lt_embedding(ROOT::RDF::RNode df,
                      correctionManager::CorrectionManager &correctionManager,
                      const std::string &pt, const std::string &wp,
                      const std::string &sf_vsjet_tau20to25,
                      const std::string &sf_vsjet_tau25to30,
                      const std::string &sf_vsjet_tau30to35,
                      const std::string &sf_vsjet_tau35to40,
                      const std::string &sf_vsjet_tau40toInf,
                      const std::string &id_output, const std::string &sf_file,
                      const std::string &correctionset);
ROOT::RDF::RNode id_vsJet_tt_embedding(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correctionManager,
    const std::string &pt,
    const std::string &decaymode, const std::string &wp,
    const std::string &sf_vsjet_tauDM0, const std::string &sf_vsjet_tauDM1,
    const std::string &sf_vsjet_tauDM10, const std::string &sf_vsjet_tauDM11,
    const std::string &id_output, const std::string &sf_file,
    const std::string &correctionset);
} // namespace tau

namespace jet {

ROOT::RDF::RNode
btagSF(ROOT::RDF::RNode df,
       correctionManager::CorrectionManager &correctionManager,
       const std::string &pt, const std::string &eta,
       const std::string &btag_discr, const std::string &flavor,
       const std::string &jet_mask, const std::string &bjet_mask,
       const std::string &jet_veto_mask, const std::string &variation,
       const std::string &sf_output, const std::string &sf_file,
       const std::string &corr_algorithm);
// Deprecated function without CorrectionManager
ROOT::RDF::RNode
btagSF(ROOT::RDF::RNode df, const std::string &pt, const std::string &eta,
       const std::string &btag_discr, const std::string &flavor,
       const std::string &jet_mask, const std::string &bjet_mask,
       const std::string &jet_veto_mask, const std::string &variation,
       const std::string &sf_output, const std::string &sf_file,
       const std::string &corr_algorithm);
} // namespace jet
namespace embedding {
ROOT::RDF::RNode
selection_trigger(ROOT::RDF::RNode df,
                  correctionManager::CorrectionManager &correctionManager,
                  const std::string &pt_1, const std::string &eta_1,
                  const std::string &pt_2, const std::string &eta_2,
                  const std::string &output, const std::string &sf_file,
                  const std::string &idAlgorithm);
// Deprecated function without CorrectionManager
ROOT::RDF::RNode
selection_trigger(ROOT::RDF::RNode df, const std::string &pt_1,
                  const std::string &eta_1, const std::string &pt_2,
                  const std::string &eta_2, const std::string &output,
                  const std::string &sf_file, const std::string &idAlgorithm);
ROOT::RDF::RNode
selection_id(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correctionManager,
             const std::string &pt, const std::string &eta,
             const std::string &output, const std::string &sf_file,
             const std::string &idAlgorithm);
// Deprecated function without CorrectionManager
ROOT::RDF::RNode selection_id(ROOT::RDF::RNode df, const std::string &pt,
                              const std::string &eta, const std::string &output,
                              const std::string &sf_file,
                              const std::string &idAlgorithm);
ROOT::RDF::RNode
muon_sf(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correctionManager,
        const std::string &pt, const std::string &eta,
        const std::string &output, const std::string &sf_file,
        const std::string correctiontype, const std::string &idAlgorithm,
        const float &extrapolation_factor = 1.0);
// Deprecated function without CorrectionManager
ROOT::RDF::RNode muon_sf(ROOT::RDF::RNode df, const std::string &pt,
                         const std::string &eta, const std::string &output,
                         const std::string &sf_file,
                         const std::string correctiontype,
                         const std::string &idAlgorithm,
                         const float &extrapolation_factor = 1.0);
ROOT::RDF::RNode
electron_sf(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correctionManager,
            const std::string &pt, const std::string &eta,
            const std::string &output, const std::string &sf_file,
            const std::string correctiontype, const std::string &idAlgorithm,
            const float &extrapolation_factor = 1.0);
// Deprecated function without CorrectionManager
ROOT::RDF::RNode electron_sf(ROOT::RDF::RNode df, const std::string &pt,
                             const std::string &eta, const std::string &output,
                             const std::string &sf_file,
                             const std::string correctiontype,
                             const std::string &idAlgorithm,
                             const float &extrapolation_factor = 1.0);
ROOT::RDF::RNode
ditau_trigger_sf(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correctionManager,
                 const std::string &pt, const std::string &decaymode,
                 const std::string &output, const std::string &wp,
                 const std::string &sf_file, const std::string &type,
                 const std::string &corrtype, const std::string &syst);
// Deprecated function without CorrectionManager
ROOT::RDF::RNode
ditau_trigger_sf(ROOT::RDF::RNode df, const std::string &pt,
                 const std::string &decaymode, const std::string &output,
                 const std::string &wp, const std::string &sf_file,
                 const std::string &type, const std::string &corrtype,
                 const std::string &syst);
} // namespace embedding
} // namespace scalefactor
#endif /* GUARD_SCALEFACTORS_H */