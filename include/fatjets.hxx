#ifndef GUARDFATJETS_H
#define GUARDFATJETS_H

namespace fatjet {
ROOT::RDF::RNode FindFatjetMatchingBjet(
    ROOT::RDF::RNode df, const std::string &output_name,
    const std::string &good_fatjet_collection, const std::string &fatjet_pt,
    const std::string &fatjet_eta, const std::string &fatjet_phi,
    const std::string &fatjet_mass, const std::string &bpair_p4_1,
    const float &deltaRmax);
ROOT::RDF::RNode FindXbbFatjet(
    ROOT::RDF::RNode df, const std::string &output_name,
    const std::string &good_fatjet_collection, const std::string &fatjet_pNet_Xbb, const std::string &fatjet_pNet_QCD);
} // namespace fatjet

namespace quantities {
namespace fatjet {
ROOT::RDF::RNode msoftdrop(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &m_softdrop,
                           const std::string &fatjetcollection,
                           const int &position);
ROOT::RDF::RNode
particleNet_XbbvsQCD(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pNet_Xbb, const std::string &pNet_QCD,
                     const std::string &fatjetcollection, const int &position);
ROOT::RDF::RNode
nsubjettiness_ratio(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tauN, const std::string &tauNm1,
                    const std::string &fatjetcollection, const int &position);
} // end namespace fatjet
} // end namespace quantities
#endif /* GUARDFATJETS_H */