#ifndef GUARDFINALSTATE_H
#define GUARDFINALSTATE_H

namespace finalstate{

ROOT::RDF::RNode n_taus(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &tau_pt,
                                const std::string &tau_eta,
                                const std::string &tau_dz,
                                const std::string &tau_vs_jet,
                                const std::string &tau_vs_ele,
                                const std::string &tau_vs_mu,
                                const std::string &tau_dm,
                                const float &TauPtThreshold,
                                const float &TauEtaThreshold,
                                const float &TauDZThreshold,
                                const float &TauVsJetThreshold,
                                const float &TauVsMuThreshold,
                                const float &TauVsEleThreshold,
                                const std::vector<int> &SelectedDecayModes);

ROOT::RDF::RNode n_muons(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &muon_pt,
                                 const std::string &muon_eta,
                                 const float &MuPtThreshold,
                                 const float &MuEtaThreshold);

ROOT::RDF::RNode n_eles(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &ele_pt,
                                 const std::string &ele_eta,
                                 const float &ElePtThreshold,
                                 const float &EleEtaThreshold);

namespace mutau{

ROOT::RDF::RNode single_mu_in_fatjet_mutau(ROOT::RDF::RNode df,
                                            const std::string &outputname,
                                            const std::string &fatjet_p4,
                                            const std::string &muon_pt,
                                            const std::string &muon_eta,
                                            const std::string &muon_phi,
                                            const std::string &muon_mass,
                                            const float &pt_threshold,
                                            const float &eta_threshold);

ROOT::RDF::RNode single_mu_in_fatjet_mutau_deltaR(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const float &pt_threshold,
                                             const float &eta_threshold);

ROOT::RDF::RNode single_mu_in_fatjet_mutau_mu_pT(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const float &pt_threshold,
                                             const float &eta_threshold);

ROOT::RDF::RNode single_mu_in_fatjet_mutau_mu_eta(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const float &pt_threshold,
                                             const float &eta_threshold);

ROOT::RDF::RNode single_mu_in_fatjet_mutau_mu_iso(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const std::string &muon_iso,
                                             const float &pt_threshold,
                                             const float &eta_threshold);                                             
}

}

#endif /* GUARDFINALSTATE_H */