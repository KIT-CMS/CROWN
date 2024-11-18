#ifndef GUARDFATJETS_H
#define GUARDFATJETS_H

namespace fatjet {
// ROOT::RDF::RNode FindFatjetMatchingBjet(
//     ROOT::RDF::RNode df, const std::string &output_name,
//     const std::string &good_fatjet_collection, const std::string &fatjet_pt,
//     const std::string &fatjet_eta, const std::string &fatjet_phi,
//     const std::string &fatjet_mass, const std::string &bpair_p4_1,
//     const float &deltaRmax);
ROOT::RDF::RNode FindXtmFatjet(
    ROOT::RDF::RNode df, const std::string &output_name,
    const std::string &good_fatjet_collection, const std::string &fatjet_pNet_XtmVsQCD);

ROOT::RDF::RNode FindXteFatjet(
    ROOT::RDF::RNode df, const std::string &output_name,
    const std::string &good_fatjet_collection, const std::string &fatjet_pNet_XteVsQCD);

ROOT::RDF::RNode FindXttFatjet(
    ROOT::RDF::RNode df, const std::string &output_name,
    const std::string &good_fatjet_collection, const std::string &fatjet_pNet_XttVsQCD);

ROOT::RDF::RNode flagGoodFatjets(ROOT::RDF::RNode df, const std::string &flagname,
                               const std::string &fatjetname);
} // namespace fatjet

namespace quantities {
namespace fatjet {
ROOT::RDF::RNode msoftdrop(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &m_softdrop,
                           const std::string &fatjetcollection,
                           const int &position);
ROOT::RDF::RNode
particleNet_XtmVsQCD(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pNet_XtmVsQCD,
                     const std::string &fatjetcollection, const int &position);
ROOT::RDF::RNode
particleNet_XteVsQCD(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pNet_XteVsQCD,
                     const std::string &fatjetcollection, const int &position);
ROOT::RDF::RNode
particleNet_XttVsQCD(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pNet_XttVsQCD,
                     const std::string &fatjetcollection, const int &position);
ROOT::RDF::RNode
nsubjettiness_ratio(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tauN, const std::string &tauNm1,
                    const std::string &fatjetcollection, const int &position);
ROOT::RDF::RNode
subjet_1(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &subjetIDx1,
                    const std::string &fatjetcollection, const int &position);

ROOT::RDF::RNode delta_R_gentau_fatjet(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &fatjet_p4,
                                 const std::string &tau_pt,
                                 const std::string &tau_eta,
                                 const std::string &tau_phi,
                                 const std::string &tau_mass);

ROOT::RDF::RNode m_inv_ditau(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &fatjet_pt,
                                 const std::string &met_fatjet_pt,
                                 const std::string &fatjet_mass); 
} // end namespace fatjet
} // end namespace quantities
#endif /* GUARDFATJETS_H */