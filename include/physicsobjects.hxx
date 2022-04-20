#ifndef GUARD_PHYSICSOBJECTS_H
#define GUARD_PHYSICSOBJECTS_H

namespace physicsobject {
ROOT::RDF::RNode CutPt(auto &df, const std::string &quantity,
                       const std::string &maskname, const float &ptThreshold);
ROOT::RDF::RNode CutEta(auto &df, const std::string &quantity,
                        const std::string &maskname, const float &EtaThreshold);
ROOT::RDF::RNode CutDz(auto &df, const std::string &quantity,
                       const std::string &maskname, const float &Threshold);
ROOT::RDF::RNode CutDxy(auto &df, const std::string &quantity,
                        const std::string &maskname, const float &Threshold);
template <class... Masks>
ROOT::RDF::RNode CombineMasks(auto &df, const std::string &maskname,
                              const Masks &...masks);
ROOT::RDF::RNode VetoCandInMask(auto &df, const std::string &outputmaskname,
                                const std::string &inputmaskname,
                                const std::string &ditaupair, const int index);
ROOT::RDF::RNode FilterMasks(auto &df, const std::string &maskname);
ROOT::RDF::RNode LeptonVetoFlag(auto &df, const std::string &outputname,
                                const std::string &vetomap);
ROOT::RDF::RNode ObjectMassCorrectionWithPt(auto &df,
                                            const std::string &corrected_mass,
                                            const std::string &raw_mass,
                                            const std::string &raw_pt,
                                            const std::string &corrected_pt);
ROOT::RDF::RNode CheckForDiLeptonPairs(
    auto &df, const std::string &output_flag, const std::string &leptons_pt,
    const std::string &leptons_eta, const std::string &leptons_phi,
    const std::string &leptons_mass, const std::string &leptons_charge,
    const std::string &leptons_mask, const float dR_cut);
namespace muon {
ROOT::RDF::RNode CutID(auto &df, const std::string &maskname,
                       const std::string &nameID);
ROOT::RDF::RNode CutIsolation(auto &df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold);
ROOT::RDF::RNode CutDecayModes(auto &df, const std::string &maskname,
                               const std::string &tau_dms,
                               const std::vector<int> &SelectedDecayModes);
ROOT::RDF::RNode CutTauID(auto &df, const std::string &maskname,
                          const std::string &nameID, const int &idxID);
ROOT::RDF::RNode PtCorrection_byValue(auto &df, const std::string &corrected_pt,
                                      const std::string &pt,
                                      const std::string &decayMode,
                                      const float &sf_dm0, const float &sf_dm1,
                                      const float &sf_dm10,
                                      const float &sf_dm11);
ROOT::RDF::RNode
PtCorrection(auto &df, const std::string &corrected_pt, const std::string &pt,
             const std::string &eta, const std::string &decayMode,
             const std::string &genMatch, const std::string &sf_file,
             const std::string &jsonESname, const std::string &idAlgorithm,
             const std::string &DM0, const std::string &DM1,
             const std::string &DM10, const std::string &DM11,
             const std::vector<int> &SelectedDMs);
} // namespace muon

namespace electron {

ROOT::RDF::RNode CutID(auto &df, const std::string &maskname,
                       const std::string &nameID);
ROOT::RDF::RNode CutCBID(auto &df, const std::string &maskname,
                         const std::string &nameID, const int &IDvalue);
ROOT::RDF::RNode CutIsolation(auto &df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold);
} // end namespace electron
} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */