#ifndef GUARD_LORENTZVECTOR_H
#define GUARD_LORENTZVECTOR_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/utility.hxx"
#include "defaults.hxx"
#include <Math/Vector4D.h>

namespace lorentzvector {

/**
 * @brief This function constructs a vectorial sum of an arbitrary number of 
 * Lorentz vectors. If one of the Lorentz vectors is not well defined (has 
 * default values), the function returns a default Lorentz vector.
 *
 * @tparam Lorentzvectors type of the input column
 * @param df input dataframe
 * @param outputname name of the output column containing the summed Lorentz vector
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <class... Lorentzvectors>
ROOT::RDF::RNode Combine(ROOT::RDF::RNode df, const std::string &outputname,
                            const Lorentzvectors &...LVs) {
    std::vector<std::string> LV_list;
    utility::appendParameterPackToVector(LV_list, LVs...);
    const auto n_LVs = sizeof...(Lorentzvectors);
    using namespace ROOT::VecOps;
    return df.Define(
        outputname,
        utility::PassAsVec<n_LVs, ROOT::Math::PtEtaPhiMVector>(
            [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &LVs) {
                for (int i = 0; i < LVs.size(); i++) {
                    auto vector = LVs.at(i);
                    if (vector.pt() < 0.) {
                        return default_lorentzvector;
                    }
                }
                return Sum(LVs, ROOT::Math::PtEtaPhiMVector());
            }),
        LV_list);
}
ROOT::RDF::RNode Build(ROOT::RDF::RNode df,
                       const std::string &outputname,
                       const std::string &pt,
                       const std::string &eta,
                       const std::string &phi,
                       const std::string &mass,
                       const std::string &index_vector,
                       const int position);
ROOT::RDF::RNode Build(ROOT::RDF::RNode df,
                       const std::string &outputname,
                       const std::string &pt,
                       const std::string &eta,
                       const std::string &phi,
                       const std::string &mass,
                       const int index);
ROOT::RDF::RNode Build(ROOT::RDF::RNode df,
                       const std::string &outputname,
                       const std::string &pt,
                       const std::string &eta,
                       const std::string &phi,
                       const std::string &mass);
ROOT::RDF::RNode BuildMET(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &met_pt, const std::string &met_phi);
ROOT::RDF::RNode BuildCollection(ROOT::RDF::RNode df,
                                   const std::string &outputname,
                                   const std::string &pt,
                                   const std::string &eta,
                                   const std::string &phi,
                                   const std::string &mass,
                                   const std::string &object_mask);                          
ROOT::RDF::RNode BuildCollection(ROOT::RDF::RNode df,
                                    const std::string &outputname,
                                    const std::string &pt,
                                    const std::string &eta,
                                    const std::string &phi,
                                    const std::string &mass);
ROOT::RDF::RNode Scale(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::string &vector,
                         const float &scalefactor);
ROOT::RDF::RNode GetPt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &vector);
ROOT::RDF::RNode GetEta(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &vector);
ROOT::RDF::RNode GetPhi(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &vector);
ROOT::RDF::RNode GetMass(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &vector);
ROOT::RDF::RNode GetEnergy(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &vector);

/// namespace used for mutau lorentzvectors
namespace mutau {
ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &pairname,
                       const std::vector<std::string> &muon_quantities,
                       const std::vector<std::string> &tau_quantities,
                       const std::string &muon_p4_name,
                       const std::string &tau_p4_name);
} // namespace mutau

} // end namespace lorentzvector
#endif /* GUARD_LORENTZVECTOR_H */
