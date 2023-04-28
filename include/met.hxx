#ifndef GUARDMET_H
#define GUARDMET_H

namespace met {

ROOT::RDF::RNode calculateGenBosonVector(
    ROOT::RDF::RNode df, const std::string &genparticle_pt,
    const std::string &genparticle_eta, const std::string &genparticle_phi,
    const std::string &genparticle_mass, const std::string &genparticle_id,
    const std::string &genparticle_status,
    const std::string &genparticle_statusflag, const std::string outputname,
    bool is_data);
ROOT::RDF::RNode genBosonMass(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &inputvector);
ROOT::RDF::RNode genBosonPt(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &inputvector);
ROOT::RDF::RNode genBosonEta(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &inputvector);
ROOT::RDF::RNode genBosonPhi(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &inputvector);
ROOT::RDF::RNode genBosonRapidity(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &inputvector);
ROOT::RDF::RNode DefineRecoilsDilep(ROOT::RDF::RNode df,
                                    const std::string &paralell_outputname,
                                    const std::string &perpendicular_outputname,
                                    const std::string &lepton1,
                                    const std::string &lepton2,
                                    const std::string &met_p4);
ROOT::RDF::RNode
DefineRecoilsSinglelep(ROOT::RDF::RNode df,
                       const std::string &paralell_outputname,
                       const std::string &perpendicular_outputname,
                       const std::string &lepton1, const std::string &met_p4);
ROOT::RDF::RNode
propagateLeptonsToMet(ROOT::RDF::RNode df, const std::string &met,
                      const std::string &p4_1_uncorrected,
                      const std::string &p4_2_uncorrected,
                      const std::string &p4_1, const std::string &p4_2,
                      const std::string &outputname, bool apply_propagation);
ROOT::RDF::RNode
propagateLeptonsToMet(ROOT::RDF::RNode df, const std::string &met,
                      const std::string &p4_1_uncorrected,
                      const std::string &p4_2_uncorrected,
                      const std::string &p4_3_uncorrected,
                      const std::string &p4_1, const std::string &p4_2, const std::string &p4_3,
                      const std::string &outputname, bool apply_propagation);
ROOT::RDF::RNode propagateLeptonsToMet(ROOT::RDF::RNode df,
                                       const std::string &met,
                                       const std::string &p4_1_uncorrected,
                                       const std::string &p4_1,
                                       const std::string &outputname,
                                       bool apply_propagation);
ROOT::RDF::RNode propagateJetsToMet(
    ROOT::RDF::RNode df, const std::string &met,
    const std::string &jet_pt_corrected, const std::string &jet_eta_corrected,
    const std::string &jet_phi_corrected, const std::string &jet_mass_corrected,
    const std::string &jet_pt, const std::string &jet_eta,
    const std::string &jet_phi, const std::string &jet_mass,
    const std::string &outputname, bool apply_propagation, float min_jet_pt);
ROOT::RDF::RNode applyRecoilCorrections(
    ROOT::RDF::RNode df, const std::string &met, const std::string &genmet,
    const std::string &jet_pt, const std::string &outputname,
    const std::string &recoilfile, const std::string &systematicsfile,
    bool applyRecoilCorrections, bool resolution, bool response, bool shiftUp,
    bool shiftDown, bool isWjets);
} // end namespace met
#endif /* GUARDMET_H */
