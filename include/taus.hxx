#ifndef GUARD_TAUS_H
#define GUARD_TAUS_H

namespace physicsobject {
namespace tau {
ROOT::RDF::RNode CutDecayModes(ROOT::RDF::RNode df, const std::string &maskname,
                               const std::string &tau_dms,
                               const std::vector<int> &SelectedDecayModes);
ROOT::RDF::RNode CutTauID(ROOT::RDF::RNode df, const std::string &maskname,
                          const std::string &nameID, const int &idxID);
ROOT::RDF::RNode
PtCorrection_eleFake(ROOT::RDF::RNode df,
                     correctionManager::CorrectionManager &correctionManager,
                     const std::string &corrected_pt, const std::string &pt,
                     const std::string &eta, const std::string &decayMode,
                     const std::string &genMatch, const std::string &sf_file,
                     const std::string &jsonESname,
                     const std::string &idAlgorithm,
                     const std::string &sf_dm0_b, const std::string &sf_dm1_b,
                     const std::string &sf_dm0_e, const std::string &sf_dm1_e);
// deprecated version without CorrectionManager
ROOT::RDF::RNode
PtCorrection_eleFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const std::string &decayMode, const std::string &genMatch,
                     const std::string &sf_file, const std::string &jsonESname,
                     const std::string &idAlgorithm,
                     const std::string &sf_dm0_b, const std::string &sf_dm1_b,
                     const std::string &sf_dm0_e, const std::string &sf_dm1_e);
ROOT::RDF::RNode
PtCorrection_muFake(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &corrected_pt, const std::string &pt,
                    const std::string &eta, const std::string &decayMode,
                    const std::string &genMatch, const std::string &sf_file,
                    const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &sf_es);
// deprecated version without CorrectionManager
ROOT::RDF::RNode
PtCorrection_muFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &sf_es);
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &decayMode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11);
ROOT::RDF::RNode
PtCorrection_genTau(ROOT::RDF::RNode df,
                    correctionManager::CorrectionManager &correctionManager,
                    const std::string &corrected_pt, const std::string &pt,
                    const std::string &eta, const std::string &decayMode,
                    const std::string &genMatch, const std::string &sf_file,
                    const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &DM0,
                    const std::string &DM1, const std::string &DM10,
                    const std::string &DM11);
// deprecated version without CorrectionManager
ROOT::RDF::RNode
PtCorrection_genTau(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &DM0,
                    const std::string &DM1, const std::string &DM10,
                    const std::string &DM11);
} // namespace tau
} // namespace physicsobject
#endif /* GUARD_TAUS_H */