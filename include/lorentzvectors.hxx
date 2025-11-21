#ifndef GUARD_LORENTZVECTORS_H
#define GUARD_LORENTZVECTORS_H

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
 * @tparam Lorentzvectors variadic template parameter pack representing the 
 * Lorentz vector columns
 * @param df input dataframe
 * @param outputname name of the output column containing the summed Lorentz vector
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode Sum(ROOT::RDF::RNode df, const std::string &outputname,
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
/**
 * @brief This function constructs a vectorial sum of an arbitrary number of 
 * Lorentz vectors (can also be only one) and returns its transverse momentum 
 * (\f$p_T\f$). If one of the Lorentz vectors is not well defined (has default 
 * values), the function returns a default value.
 *
 * @tparam Lorentzvectors variadic template parameter pack representing the 
 * Lorentz vector columns
 * @param df input dataframe
 * @param outputname name of the output column containing the \f$p_T\f$ of the 
 * summed Lorentz vectors
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode GetPt(ROOT::RDF::RNode df, const std::string &outputname,
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
                        return default_float;
                    }
                }
                ROOT::Math::PtEtaPhiMVector new_LV = Sum(LVs, ROOT::Math::PtEtaPhiMVector());
                return (float)new_LV.pt();
            }),
        LV_list);
}

/**
 * @brief This function constructs a vectorial sum of an arbitrary number of 
 * Lorentz vectors (can also be only one) and returns its pseudorapodity 
 * \f$\eta\f$. If one of the Lorentz vectors is not well defined (has default 
 * values), the function returns a default value.
 *
 * @tparam Lorentzvectors variadic template parameter pack representing the 
 * Lorentz vector columns
 * @param df input dataframe
 * @param outputname name of the output column containing the pseudorapidity 
 * \f$\eta\f$ of the summed Lorentz vectors
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode GetEta(ROOT::RDF::RNode df, const std::string &outputname,
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
                        return default_float;
                    }
                }
                ROOT::Math::PtEtaPhiMVector new_LV = Sum(LVs, ROOT::Math::PtEtaPhiMVector());
                return (float)new_LV.eta();
            }),
        LV_list);
}

/**
 * @brief This function constructs a vectorial sum of an arbitrary number of 
 * Lorentz vectors (can also be only one) and returns its azimuthal angle 
 * \f$\phi\f$. If one of the Lorentz vectors is not well defined (has default 
 * values), the function returns a default value.
 *
 * @tparam Lorentzvectors variadic template parameter pack representing the 
 * Lorentz vector columns
 * @param df input dataframe
 * @param outputname name of the output column containing the azimuthal angle 
 * \f$\phi\f$ of the summed Lorentz vectors
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode GetPhi(ROOT::RDF::RNode df, const std::string &outputname,
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
                        return default_float;
                    }
                }
                ROOT::Math::PtEtaPhiMVector new_LV = Sum(LVs, ROOT::Math::PtEtaPhiMVector());
                return (float)new_LV.phi();
            }),
        LV_list);
}

/**
 * @brief This function constructs a vectorial sum of an arbitrary number of 
 * Lorentz vectors (can also be only one) and returns its invariant mass. If 
 * one of the Lorentz vectors is not well defined (has default values), the 
 * function returns a default value.
 *
 * @tparam Lorentzvectors variadic template parameter pack representing the 
 * Lorentz vector columns
 * @param df input dataframe
 * @param outputname name of the output column containing the invariant mass 
 * of the summed Lorentz vectors
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode GetMass(ROOT::RDF::RNode df, const std::string &outputname,
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
                        return default_float;
                    }
                }
                ROOT::Math::PtEtaPhiMVector new_LV = Sum(LVs, ROOT::Math::PtEtaPhiMVector());
                return (float)new_LV.mass();
            }),
        LV_list);
}

/**
 * @brief This function constructs a vectorial sum of an arbitrary number of 
 * Lorentz vectors (can also be only one) and returns its energy. If one of 
 * the Lorentz vectors is not well defined (has default values), the function 
 * returns a default value.
 *
 * @tparam Lorentzvectors variadic template parameter pack representing the 
 * Lorentz vector columns
 * @param df input dataframe
 * @param outputname name of the output column containing the energy of the 
 * summed Lorentz vectors
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode GetEnergy(ROOT::RDF::RNode df, const std::string &outputname,
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
                        return default_float;
                    }
                }
                ROOT::Math::PtEtaPhiMVector new_LV = Sum(LVs, ROOT::Math::PtEtaPhiMVector());
                return (float)new_LV.energy();
            }),
        LV_list);
}

/**
 * @brief This function constructs a vectorial sum of an arbitrary number of 
 * Lorentz vectors (can also be only one) and returns its rapidity \f$y\f$. 
 * If one of the Lorentz vectors is not well defined (has default values), the
 * function returns a default value.
 *
 * @tparam Lorentzvectors variadic template parameter pack representing the 
 * Lorentz vector columns
 * @param df input dataframe
 * @param outputname name of the output column containing the rapidity \f$y\f$ 
 * of the summed Lorentz vectors
 * @param LVs Parameter pack of column names that contain the considered
 * Lorentz vectors, must be of type `ROOT::Math::PtEtaPhiMVector`
 *
 * @return a dataframe containing the new column
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode GetRapidity(ROOT::RDF::RNode df, const std::string &outputname,
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
                        return default_float;
                    }
                }
                ROOT::Math::PtEtaPhiMVector new_LV = Sum(LVs, ROOT::Math::PtEtaPhiMVector());
                return (float)new_LV.Rapidity();
            }),
        LV_list);
}
} // end namespace lorentzvector
#endif /* GUARD_LORENTZVECTORS_H */
