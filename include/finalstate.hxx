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

}

#endif /* GUARDFINALSTATE_H */