#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "../include/SVFit/FastMTT.hxx"
#include "../include/SVFit/MeasuredTauLepton.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>


namespace quantities {

/**
 * @brief This function calculates the spatial distance in the x-y-plane 
 * (\f$\Delta\phi\f$) between two Lorentz vectors.
 *
 * @note For the calculation the `ROOT::Math::VectorUtil::DeltaPhi()` 
 * function is used which already takes care of the periodicity of the 
 * azimuthal angle.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the \f$\Delta\phi\f$ 
 * value
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode DeltaPhi(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &vector_1,
                          const std::string &vector_2) {
    auto calculate_deltaPhi = [](ROOT::Math::PtEtaPhiMVector &p4_1,
                                 ROOT::Math::PtEtaPhiMVector &p4_2) {
        if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
            return default_float;
        return (float)ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_2);
    };
    return df.Define(outputname, calculate_deltaPhi, {vector_1, vector_2});
}

/**
 * @brief This function calculates the spatial distance in the 
 * \f$\eta\f$-\f$\phi\f$-plane (\f$\Delta R\f$) between two Lorentz vectors.
 * It is defined as
 * \f[ \Delta R = \sqrt{(\eta_1 - \eta_2)^2 + (\phi_1 - \phi_2)^2} \f]
 * where \f$\eta_1\f$ and \f$\phi_1\f$ are from the first Lorentz vector and \f$\eta_2\f$
 * and \f$\phi_2\f$ are from the second Lorentz vector.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the \f$\Delta R\f$ 
 * value
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode DeltaR(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2) {
    auto calculate_deltaR = [](ROOT::Math::PtEtaPhiMVector &p4_1,
                               ROOT::Math::PtEtaPhiMVector &p4_2) {
        if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
            return default_float;
        return (float)ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2);
    };
    return df.Define(outputname, calculate_deltaR, {vector_1, vector_2});
}

/**
 * @brief This function checks the hemisphere of a pair of particles. If both
 * particles are in the same hemisphere (both positive/negative \f$\eta\f$), 
 * the quantity is set to `1`, otherwise it is set to `0`.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the hemisphere value
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode PairHemisphere(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &vector_1,
                                const std::string &vector_2) {
    auto calculate_hemisphere = [](ROOT::Math::PtEtaPhiMVector &p4_1,
                                   ROOT::Math::PtEtaPhiMVector &p4_2) {
        if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
            return default_int;
        return (int)(p4_1.Eta() * p4_2.Eta() > 0);
    };
    return df.Define(outputname, calculate_hemisphere,
                     {vector_1, vector_2});
}

/**
 * @brief This function calculates the quantity `pZetaMissVis` from the two leptons
 * in the event and the MET vector. The variable is defined as:
 * \f[
 *    D_\zeta = p_\zeta^\text{miss} - 0.85 p_\zeta^\text{vis}
 *   \qquad;
 *   p_\zeta^\text{miss} = \vec{p}_\text{T}^\text{miss} \cdot \hat{\zeta}
 *   \qquad;
 *   p_\zeta^\text{vis} = (\vec{p}_\text{T}^{p_1} + \vec{p}_\text{T}^{p_2}) \cdot
 *   \hat{\zeta} 
 * \f] 
 * where \f$\vec{p}_\text{T}^{p_{1,2}}\f$ corresponds to the transverse momentum 
 * vector of the first (second) lepton and \f$\hat{\zeta}\f$ to the bisectional 
 * direction between the two leptons in the transverse plane.
 *
 * For more information check: D. Jang, “Search for MSSM Higgs decaying to tau pairs 
 * in pp collision at √s=1.96 TeV at CDF”. PhD thesis, Rutgers University, 2006. 
 * FERMILAB-THESIS-2006-11.
 *
 * @param df the input dataframe
 * @param outputname the name of the new column containing the PzetaMissVis value
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 * @param vector_3 name of the column containing the third Lorentz vector (MET vector)
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode PzetaMissVis(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &vector_1,
                              const std::string &vector_2,
                              const std::string &vector_3) {
    float alpha = 0.85;
    auto calculate_pzetamissvis = [alpha](ROOT::Math::PtEtaPhiMVector &p4_1,
                                          ROOT::Math::PtEtaPhiMVector &p4_2,
                                          ROOT::Math::PtEtaPhiMVector &p4_met) {
        auto met_3dvec = p4_met.Vect();
        met_3dvec.SetZ(0.0);
        // calculate zeta for the delepton system
        auto p1_norm = p4_1.Vect().Unit();
        auto p2_norm = p4_2.Vect().Unit();
        p1_norm.SetZ(0.0);
        p2_norm.SetZ(0.0);
        p1_norm = p1_norm.Unit();
        p2_norm = p2_norm.Unit();
        auto zeta = (p1_norm + p2_norm).Unit();

        auto dileptonsystem = p4_1.Vect() + p4_2.Vect();
        dileptonsystem.SetZ(0);
        auto pzetaVis = dileptonsystem.Dot(zeta);
        return met_3dvec.Dot(zeta) - (alpha * pzetaVis);
    };
    return df.Define(outputname, calculate_pzetamissvis, {vector_1, vector_2, vector_3});
}

/**
 * @brief This function calculates the transverse mass \f$m_T\f$ of a two particle 
 * system, where both particles are massless. The transverse mass is defined as:
 *
 * \f[
 *    m_{T} = \sqrt{2 \cdot p_{T,1} \cdot p_{T,2} \cdot
 * (1-\cos(\Delta\phi))}
 * \f]
 *
 * where \f$\Delta\phi\f$ is the azimuthal angle between the two particles.
 *
 * @note The transverse mass is usually used to estimate the mass of the W boson
 * based on a lepton (particle 1) and the missing transverse energy as the 
 * neutrino (particle 2). 
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the \f$m_T\f$ value
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode TransverseMass(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2) {
    auto calculate_MT = [](ROOT::Math::PtEtaPhiMVector &p4_1,
                           ROOT::Math::PtEtaPhiMVector &p4_2) {
        if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
            return default_float;
        return (float)sqrt(2 * p4_1.Pt() * p4_2.Pt() *
            (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_2))));
    };
    return df.Define(outputname, calculate_MT, {vector_1, vector_2});
}

/**
 * @brief This function calculates the transverse mass \f$m_T\f$ of a three 
 * particle system, where first the two first particles are summed up and 
 * used to calculate the transverse mass with the third particle. 
 * The transverse mass is defined as:
 *
 * \f[
 *    m_{T} = \sqrt{2 \cdot p_{T,1+2} \cdot p_{T,3} \cdot
 * (1-\cos(\Delta\phi))}
 * \f]
 *
 * where \f$\Delta\phi\f$ is the azimuthal angle between the third particle 
 * and the summed Lorentz vector of the first two particles.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the \f$m_T\f$ value
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 * @param vector_3 name of the column containing the third Lorentz vector
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode TransverseMass(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2, 
                        const std::string &vector_3) {
    auto calculate_MT = [](ROOT::Math::PtEtaPhiMVector &p4_1,
                           ROOT::Math::PtEtaPhiMVector &p4_2,
                           ROOT::Math::PtEtaPhiMVector &p4_3) {
        if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0 || p4_3.pt() < 0.0)
            return default_float;
        ROOT::Math::PtEtaPhiMVector sum_p4 = p4_1 + p4_2;
        return (float)sqrt(2 * sum_p4.Pt() * p4_3.Pt() *
            (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(sum_p4, p4_3))));
    };
    return df.Define(outputname, calculate_MT, {vector_1, vector_2, vector_3});
}

/**
 * @brief This function calculates the total transverse mass. This is usually
 * used to estimate the Higgs to \f$\tau\tau\f$ decay where multiple neutrinos
 * are involved and estimated via the MET vector. The total transverse mass is 
 * defined as:
 * \f[
 *   m_{T}^{tot} = \sqrt{m_{T}^2(p_{1},E_{T}^{miss}) +
 *    m_{T}^2(p_{2},E_{T}^{miss}) + m_{T}^2(p_{1},p_2) } 
 * \f]
 * where \f$ m_{T}^2 \f$ is the transverse mass, \f$ p_{1}\f$ and \f$ p_{2}\f$ 
 * are the lepton Lorentzvectors and \f$E_{T}^{miss}\f$ is the missing energy.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the total transverse mass
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 * @param vector_3 name of the column containing the third Lorentz vector (usually
 * the missing transverse energy vector)
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode TotalTransverseMass(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2,
                        const std::string &vector_3) {
    auto calculate_mt_tot = [](ROOT::Math::PtEtaPhiMVector &p4_1,
                               ROOT::Math::PtEtaPhiMVector &p4_2,
                               ROOT::Math::PtEtaPhiMVector &p4_met) {
        const float mt_1 = sqrt(2 * p4_1.Pt() * p4_met.Pt() *
            (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_met))));
        const float mt_2 = sqrt(2 * p4_2.Pt() * p4_met.Pt() *
            (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(p4_2, p4_met))));
        const float mt_mix = sqrt(2 * p4_1.Pt() * p4_2.Pt() *
            (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(p4_1, p4_2))));
        return (float)sqrt(mt_1 * mt_1 + mt_2 * mt_2 + mt_mix * mt_mix);
    };
    return df.Define(outputname, calculate_mt_tot, {vector_1, vector_2, vector_3});
}

/**
 * @brief This function calculates collinear mass approximation. It is defined through 
 * two equations, assuming, that neutrinos of \f$\tau\f$ decays fly into the same 
 * direction as visible decay products:
 * \f[
 *   p(\tau_{i}) = (1 + x_{\tau_{i}^{vis}}) \cdot p(\tau_{i}^{vis}), \qquad i = 1,2
 * \f]
 * where \f$ p(...)\f$ represents the Lorentz vectors of the \f$\tau\f$ leptons
 * \f$\tau_{1}\f$ and \f$\tau_{2}\f$. The fractions \f$x_{\tau_{i}^{vis}}\f$ are the 
 * additional amount of neutrino contributions, relative to the visible decay products. 
 * This means, the missing transverse energy vector \f$\vec{p}_{T}^{miss}\f$ can be 
 * computed as follows:
 * \f[
 *   \vec{p}_{T}^{miss} = x_{\tau_{1}^{vis}} \cdot \vec{p}_{T}(\tau_{1}^{vis}) 
 *    + x_{\tau_{2}^{vis}} \cdot \vec{p}_{T}(\tau_{2}^{vis})
 * \f]
 * This set of equations in turn allows to determine the values \f$x_{\tau_{i}^{vis}}\f$. 
 * Example for \f$i=1\f$:
 * \f[
 *   x_{\tau_{1}^{vis}} = \frac{p_{T}^{miss}}{p_{T}(\tau_{1}^{vis})} \cdot 
 *    \frac{\sin(\phi_{\tau_{2}^{vis}} - \phi_{miss})}{\sin(\phi_{\tau_{2}^{vis}} 
 *    - \phi_{\tau_{1}^{vis}})}
 * \f]
 * The collinear mass approximation is then computed from the sum of the full \f$\tau\f$ 
 * Lorentz vectors \f$p(\tau_{i})\f$.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the approximanted collinear mass
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 * @param vector_3 name of the column containing the third Lorentz vector (MET vector)
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode CollinearApproxMtt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &vector_1, const std::string &vector_2,
                        const std::string &vector_3) {
    auto calculate_mtt = [](ROOT::Math::PtEtaPhiMVector &p4_1,
                            ROOT::Math::PtEtaPhiMVector &p4_2,
                            ROOT::Math::PtEtaPhiMVector &p4_met) {
        // Require valid pt values
        if (p4_1.Pt() < 0. || p4_2.Pt() < 0. || p4_met.Pt() > 0.)
            return default_float;

        // Calculate the phi difference between the two visible particles
        const float delta_phi = p4_2.Phi() - p4_1.Phi();
        // Avoid division by zero by checking the sine of the phi difference
        if (std::fabs(sin(delta_phi)) < 1e-6)
            return default_float;

        const float x_1 = p4_met.Pt() / p4_1.Pt() * sin(p4_2.Phi() - p4_met.Phi()) / sin(delta_phi);
        const float x_2 = p4_met.Pt() / p4_2.Pt() * sin(p4_1.Phi() - p4_met.Phi()) / (-sin(delta_phi));
        ROOT::Math::PtEtaPhiMVector coll_lorentz = (1. + x_1) * p4_1 + (1. + x_2) * p4_2;
        return (float)coll_lorentz.mass();
    };
    return df.Define(outputname, calculate_mtt, {vector_1, vector_2, vector_3});
}

/**
 * @brief This function calculates the FastMTT Lorentz vector as an estimate 
 * for H\f$(\tau\tau)\f$ based on the information from both reconstructed 
 * leptons and the reconstructed MET. The implementation is based on
 * https://github.com/SVfit/ClassicSVfit/tree/fastMTT_19_02_2019
 *
 * @param df The dataframe to add the quantity to
 * @param outputname name of the new column containing the Lorentz vector
 * @param pt_1 the name of the column containing the pt of the first particle
 * @param pt_2 the name of the column containing the pt of the second particle
 * @param eta_1  the name of the column containing the eta of the first particle
 * @param eta_2 the name of the column containing the eta of the second particle
 * @param phi_1 the name of the column containing the phi of the first particle
 * @param phi_2 the name of the column containing the phi of the second particle
 * @param mass_1 the name of the column containing the mass of the first
 * particle
 * @param mass_2 the name of the column containing the mass of the second
 * particle
 * @param met_pt the name of the column containing the met pt
 * @param met_phi the name of the column containing the met phi
 * @param met_cov_xx the name of the column containing the met covariance xx
 * @param met_cov_xy the name of the column containing the met covariance xy
 * @param met_cov_yy the name of the column containing the met covariance yy
 * @param decay_mode_1 the name of the column containing the decay mode of the
 * first particle
 * @param decay_mode_2 the name of the column containing the decay mode of the
 * second particle
 * @param finalstate the final state of the ditaudecay. Supported are "mt",
 * "et", "tt", "em"
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode
FastMtt(ROOT::RDF::RNode df, const std::string &outputname,
           const std::string &pt_1, const std::string &pt_2,
           const std::string &eta_1, const std::string &eta_2,
           const std::string &phi_1, const std::string &phi_2,
           const std::string &mass_1, const std::string &mass_2,
           const std::string &met_pt, const std::string &met_phi,
           const std::string &met_cov_xx, const std::string &met_cov_xy,
           const std::string &met_cov_yy, const std::string &decay_mode_1,
           const std::string &decay_mode_2, const std::string &finalstate) {
    // initialize the FastMTT algorithm
    auto calculate_fast_mtt =
        [finalstate](const float &pt_1, const float &pt_2, const float &eta_1,
                     const float &eta_2, const float &phi_1, const float &phi_2,
                     const float &mass_1, const float &mass_2,
                     const float &met_pt, const float &met_phi,
                     const float &met_cov_xx, const float &met_cov_xy,
                     const float &met_cov_yy, const int &decay_mode_1,
                     const int &decay_mode_2) {
            std::vector<fastmtt::MeasuredTauLepton> measuredTauLeptons;
            TMatrixD covMET(2, 2);
            covMET[0][0] = met_cov_xx;
            covMET[1][0] = met_cov_xy;
            covMET[0][1] = met_cov_xy;
            covMET[1][1] = met_cov_yy;
            // build the met lorentz vector
            ROOT::Math::PtEtaPhiMVector met(met_pt, 0.0, met_phi, 0.0);
            // set the decay modes according to the final state
            auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
            auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
            int dm_1, dm_2;
            if (finalstate == "mt") {
                dm_1 = -1;
                dm_2 = decay_mode_2;
                auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToMuDecay;
                auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
            } else if (finalstate == "et") {
                dm_1 = -1;
                dm_2 = decay_mode_2;
                auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToElecDecay;
                auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToHadDecay;
            } else if (finalstate == "tt") {
                dm_1 = decay_mode_1;
                dm_2 = decay_mode_2;
            } else if (finalstate == "em") {
                dm_1 = -1;
                dm_2 = -1;
                auto decay_obj_1 = fastmtt::MeasuredTauLepton::kTauToElecDecay;
                auto decay_obj_2 = fastmtt::MeasuredTauLepton::kTauToMuDecay;
            } else {
                Logger::get("FastMTT")->error(
                    "Final state {} not supported by FastMTT", finalstate);
                return (ROOT::Math::PtEtaPhiMVector)LorentzVector();
            }
            measuredTauLeptons.push_back(fastmtt::MeasuredTauLepton(
                decay_obj_1, pt_1, eta_1, phi_1, mass_1, dm_1));
            measuredTauLeptons.push_back(fastmtt::MeasuredTauLepton(
                decay_obj_2, pt_2, eta_2, phi_2, mass_2, dm_2));
            FastMTT FastMTTAlgo;
            ROOT::Math::PtEtaPhiMVector result =
                FastMTTAlgo.run(measuredTauLeptons, met.X(), met.Y(), covMET);
            // ROOT::Math::PtEtaPhiMVector result(_result.Pt(), _result.Eta(),
            //                                    _result.Phi(), _result.M());
            Logger::get("FastMTT")->debug("FastMTT result: {}", result.M());
            return result;
        };
    return df.Define(outputname, calculate_fast_mtt,
                     {pt_1, pt_2, eta_1, eta_2, phi_1, phi_2, mass_1, mass_2,
                      met_pt, met_phi, met_cov_xx, met_cov_xy, met_cov_yy,
                      decay_mode_1, decay_mode_2});
}

/**
 * @brief This function calculates the deltaPhi between the lepton from a W boson
 * and the visible Higgs boson decay products.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the deltaR value
 * @param vector_1 name of the column containing the first Lorentz vector
 * @param vector_2 name of the column containing the second Lorentz vector
 * @param vector_3 name of the column containing the third Lorentz vector
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode deltaPhi_WH(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &vector_1,
                             const std::string &vector_2,
                             const std::string &vector_3) {
    auto calculate_deltaPhi = [](ROOT::Math::PtEtaPhiMVector &vector_1,
                                 ROOT::Math::PtEtaPhiMVector &vector_2,
                                 ROOT::Math::PtEtaPhiMVector &vector_3) {
        auto const dileptonsystem = vector_2 + vector_3;
        return ROOT::Math::VectorUtil::DeltaPhi(vector_1, dileptonsystem);
    };
    return df.Define(outputname, calculate_deltaPhi, {vector_1, vector_2, vector_3});
}

/**
 * @brief This function estimates the pt of the W boson from the visible lepton 
 * Lorentz vector, the MET Lorentz vector and the neutrino Lorentz vector 
 * component from the Higgs system.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the pt value
 * @param vectors vector with three names of the columns containing the
 * required Lorentz vectors
 *
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode pt_W(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::vector<std::string> &vectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2,
           const ROOT::Math::PtEtaPhiMVector &p4_3) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0 || p4_3.pt() < 0.0)
                return default_float;
            auto const w_p4 = p4_1 + p4_2 - p4_3;
            return (float)w_p4.Pt();
        },
        vectors);
}

/// namespace for tau specific quantities
namespace tau {

/// Function to writeout the pt of the reco jet associated with the given
/// tau.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association
/// between the tau and the reco jet \param jetpt_column name of the column
/// containing the recojet pt values
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode matching_jet_pt(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const int &position,
                                 const std::string &pairname,
                                 const std::string &taujet_index,
                                 const std::string &jetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<float> &jetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         return jetpt.at(jetindex, default_float);
                     },
                     {pairname, taujet_index, jetpt_column});
}
/// Function to writeout the pt of the gen jet associated with the reco jet,
/// which is associated with the given tau. \code
///  Tau --> recoJet --> GenJet
///   \endcode
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association
/// between the tau and the reco jet \param genjet_index name of the column
/// containing the association between the reco jet and the gen jet \param
/// genjetpt_column name of the column containing the genJet pt values
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode matching_genjet_pt(
    ROOT::RDF::RNode df, const std::string &outputname, const int &position,
    const std::string &pairname, const std::string &taujet_index,
    const std::string &genjet_index, const std::string &genjetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<int> &genjets,
                                const ROOT::RVec<float> &genjetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         const int genjetindex = genjets.at(jetindex, -1);
                         return genjetpt.at(genjetindex, default_float);
                     },
                     {pairname, taujet_index, genjet_index, genjetpt_column});
}
/// Function to writeout a flag if a tau passes a specific tau id cut. The
/// particle is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the flag
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param nameID name of the ID column in the NanoAOD
/// \param idxID bitvalue of the WP the has to be passed
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode TauIDFlag(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &nameID, const int &idxID) {
    return df.Define(
        outputname,
        [position, idxID](const ROOT::RVec<int> &pair,
                          const ROOT::RVec<UChar_t> &IDs) {
            Logger::get("tauIDFlag")
                ->debug(
                    "position tau in pair {}, pair {}, id bit {}, vsjet ids {}",
                    position, pair, idxID, IDs);
            const int index = pair.at(position);
            const int ID = IDs.at(index, default_int);
            if (ID != default_int)
                return std::min(1, int(ID & 1 << (idxID - 1)));
            else
                return int(ID);
        },
        {pairname, nameID});
}
} // end namespace tau
} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */
