#ifndef GUARD_FATJETS_H
#define GUARD_FATJETS_H

namespace physicsobject {
namespace fatjet {
namespace quantity {


ROOT::RDF::RNode
ID(ROOT::RDF::RNode df,
        correctionManager::CorrectionManager &correction_manager,
        const std::string &outputname,
        const std::string &jet_eta,
        const std::string &jet_chHEF,
        const std::string &jet_neHEF,
        const std::string &jet_chEmEF,
        const std::string &jet_neEmEF,
        const std::string &jet_chMult,
        const std::string &jet_neMult,
        const std::string &jet_id_file,
        const std::string &jet_name);
ROOT::RDF::RNode
ParticleNet_XvsQCD(ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &pNet_X_decay,
    const std::string &pNet_QCD,
    const std::string &fatjet_collection,
    const int &position);
ROOT::RDF::RNode
NsubjettinessRatio(ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &tau_N,
    const std::string &tau_Nm1,
    const std::string &fatjet_collection,
    const int &position);
} // end namespace quantity
} // end namespace fatjet
} // end namespace physicsobject
#endif /* GUARD_FATJETS_H */
