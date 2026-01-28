#ifndef GUARD_FATJETS_H
#define GUARD_FATJETS_H

namespace physicsobject {
namespace fatjet {
namespace quantity {

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