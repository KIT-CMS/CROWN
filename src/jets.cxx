#ifndef GUARDJETS_H
#define GUARDJETS_H

#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TRandom3.h"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>

namespace jet {
/// Function to veto jets overlapping with tau candidates
///
/// \param[in] df the input dataframe
/// \param[out] output_col the name of the produced mask \param[in] jet_eta name
/// of the jet etas \param[in] jet_phi name of the jet phis \param[in] p4_1 four
/// vector of the first tau candidate \param[in] p4_2 four vector of the second
/// tau candidate \param[in] deltaRmin minimum required distance in dR between
/// jets and tau candidates
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
VetoOverlappingJets(ROOT::RDF::RNode df, const std::string &output_col,
                    const std::string &jet_eta, const std::string &jet_phi,
                    const std::string &p4_1, const std::string &p4_2,
                    const float &deltaRmin) {
    auto df1 = df.Define(
        output_col,
        [deltaRmin](const ROOT::RVec<float> &jet_eta,
                    const ROOT::RVec<float> &jet_phi,
                    const ROOT::Math::PtEtaPhiMVector &p4_1,
                    ROOT::Math::PtEtaPhiMVector &p4_2) {
            Logger::get("VetoOverlappingJets")->debug("Checking jets");
            ROOT::RVec<int> mask(jet_eta.size(), 1);
            for (std::size_t idx = 0; idx < mask.size(); ++idx) {
                ROOT::Math::RhoEtaPhiVectorF jet(0, jet_eta.at(idx),
                                                 jet_phi.at(idx));
                Logger::get("VetoOverlappingJets")
                    ->debug("Jet:  Eta: {} Phi: {} ", jet.Eta(), jet.Phi());
                Logger::get("VetoOverlappingJets")
                    ->debug("Letpon 1 {}:  Eta: {} Phi: {}, Pt{}", p4_1,
                            p4_1.Eta(), p4_1.Phi(), p4_1.Pt());
                Logger::get("VetoOverlappingJets")
                    ->debug("Lepton 2 {}:  Eta: {} Phi: {}, Pt{}", p4_2,
                            p4_2.Eta(), p4_2.Phi(), p4_2.Pt());
                auto deltaR_1 = ROOT::Math::VectorUtil::DeltaR(jet, p4_1);
                auto deltaR_2 = ROOT::Math::VectorUtil::DeltaR(jet, p4_2);
                Logger::get("VetoOverlappingJets")
                    ->debug("DeltaR 1 {}", deltaR_1);
                Logger::get("VetoOverlappingJets")
                    ->debug("DeltaR 2 {}", deltaR_2);
                mask[idx] = (deltaR_1 > deltaRmin && deltaR_2 > deltaRmin);
            }
            return mask;
        },
        {jet_eta, jet_phi, p4_1, p4_2});
    return df1;
}

/// Function to determine pt order of jets
///
/// \param[in] df the input dataframe
/// \param[out] output_col the name of the produced mask \param[in] jet_pt name
/// of the jet pts \param[in] jetmask name of the mask marking all valid jets to
/// be considered
///
/// \return a dataframe containing a list of jet indices sorted by pt
ROOT::RDF::RNode OrderJetsByPt(ROOT::RDF::RNode df,
                               const std::string &output_col,
                               const std::string &jet_pt,
                               const std::string &jetmask_name) {
    auto df1 = df.Define(
        output_col,
        [output_col, jetmask_name](const ROOT::RVec<int> &jetmask,
                                   const ROOT::RVec<float> &jet_pt) {
            Logger::get("OrderJetsByPt")
                ->debug("Ordering good jets from {} by pt, output stored in {}",
                        jetmask_name, output_col);
            Logger::get("OrderJetsByPt")->debug("Jetpt before {}", jet_pt);
            Logger::get("OrderJetsByPt")->debug("Mask {}", jetmask);
            auto good_jets_pt =
                ROOT::VecOps::Where(jetmask > 0, jet_pt, (float)0.);
            Logger::get("OrderJetsByPt")->debug("Jetpt after {}", good_jets_pt);
            // we have to convert the result into an RVec of ints since argsort
            // gives back an unsigned long vector
            auto temp = ROOT::VecOps::Intersect(
                ROOT::VecOps::Argsort(good_jets_pt,
                                      [](double x, double y) { return x > y; }),
                ROOT::VecOps::Nonzero(good_jets_pt));
            Logger::get("OrderJetsByPt")->debug("jet Indices {}", temp);
            ROOT::RVec<int> result(temp.size());
            std::transform(temp.begin(), temp.end(), result.begin(),
                           [](unsigned long int x) { return (int)x; });
            Logger::get("OrderJetsByPt")->debug("jet Indices int {}", result);
            return result;
        },
        {jetmask_name, jet_pt});
    return df1;
}
} // end namespace jet

namespace physicsobject {
namespace jet {

/// Function to cut jets based on the jet ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID, const int &idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterJetID(idxID), {nameID});
    return df1;
}
/// Function to cut jets based on the jet pileup ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed \param[in] jet_pt name of the
/// input jet pts \param[in] jet_pt_cut threshold for the input jet pts
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutPUID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const std::string &jet_pt,
                         const int &idxID, const float &jet_pt_cut) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterJetPUID(idxID, jet_pt_cut),
                  {nameID, jet_pt});
    return df1;
}

/// Function to shift and smear jet pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_jet_pt the name of the shifted and smeared jet pts
/// \param[in] jet_pt name of the input jet pts
/// \param[in] jet_eta name of the jet etas
/// \param[in] jet_phi name of the jet phis
/// \param[in] gen_jet_pt name of the gen jet pts
/// \param[in] gen_jet_eta name of the gen jet etas
/// \param[in] gen_jet_phi name of the gen jet phis
/// \param[in] rho name of the pileup density
/// \param[in] energy_shift_sources vector of JEC unc source names to be applied
/// in one group
/// \param[in] energy_shift_state parameter to control jet energy
/// scale shift: 0 - nominal; 1 - Up; -1 - Down
/// \param[in] energy_reso_shift parameter to control jet energy resolution
/// shift: 0 - nominal; 1 - Up; -1 - Down
///
/// \return a dataframe containing the modified jet pts
ROOT::RDF::RNode
JetPtCorrection(ROOT::RDF::RNode df, const std::string &corrected_jet_pt,
                const std::string &jet_pt, const std::string &jet_eta,
                const std::string &jet_phi, const std::string &gen_jet_pt,
                const std::string &gen_jet_eta, const std::string &gen_jet_phi,
                const std::string &rho,
                const std::vector<std::string> &energy_shift_sources,
                const int &energy_shift_state, const int &energy_reso_shift) {
    // configure readout tool here and capture it in lambda below
    // dummy lambdas for now:
    std::vector<std::function<float(float, float)>> JetEnergyShiftSources;
    for (const auto &source : energy_shift_sources) {
        JetEnergyShiftSources.push_back(
            [&](float x1, float x2) { return 0.05; });
    }
    auto JetEnergyResolution = [](const float pt, const float eta,
                                  const float rho) {
        return 1.0 / std::sqrt(pt);
    };
    auto JetEnergyResolutionSF = [](const float pt, const float eta,
                                    const int jer_shift) {
        return 1.02 + 0.01 * jer_shift;
    };
    // lambda run with dataframe
    auto JetEnergyCorrectionLambda =
        [JetEnergyShiftSources, JetEnergyResolution, JetEnergyResolutionSF,
         energy_shift_state, energy_reso_shift](
            const ROOT::RVec<float> &pt_values,
            const ROOT::RVec<float> &eta_values,
            const ROOT::RVec<float> &phi_values,
            const ROOT::RVec<float> &gen_pt_values,
            const ROOT::RVec<float> &gen_eta_values,
            const ROOT::RVec<float> &gen_phi_values, const float &rho_value) {
            ROOT::RVec<float> pt_values_corrected;
            for (int i = 0; i < pt_values.size(); i++) {
                float pt_scale_shift = 0.0;
                // jet energy scale should already be corrected
                // apply uncertainty shifts related to the energy scale
                if (energy_shift_state != 0.0) {
                    // Single source keeps sign
                    if (JetEnergyShiftSources.size() == 1) {
                        pt_scale_shift = energy_shift_state *
                                         JetEnergyShiftSources.at(0)(
                                             pt_values.at(i), eta_values.at(i));
                    } else {
                        for (const auto &source : JetEnergyShiftSources) {
                            pt_scale_shift += std::pow(
                                source(pt_values.at(i), eta_values.at(i)), 2.0);
                        }
                        pt_scale_shift =
                            energy_shift_state * std::sqrt(pt_scale_shift);
                    }
                }
                pt_values_corrected.push_back(pt_values.at(i) + pt_scale_shift);
                Logger::get("JetEnergyResolution")
                    ->debug("JE scale: Shifting jet pt from {} to {} ",
                            pt_values.at(i), pt_values_corrected.at(i));
                // apply jet energy smearing - hybrid method
                float reso = JetEnergyResolution(pt_values_corrected.at(i),
                                                 eta_values.at(i), rho_value);
                float resoSF =
                    JetEnergyResolutionSF(pt_values_corrected.at(i),
                                          eta_values.at(i), energy_reso_shift);
                // search for gen jet
                ROOT::Math::RhoEtaPhiVectorF jet(pt_values_corrected.at(i),
                                                 eta_values.at(i),
                                                 phi_values.at(i));
                float genjetpt = -1.0;
                Logger::get("JetEnergyResolution")
                    ->debug("Going to smear jet:  Eta: {} Phi: {} ", jet.Eta(),
                            jet.Phi());
                double min_dR = std::numeric_limits<double>::infinity();
                for (int j = 0; j < gen_pt_values.size(); j++) {
                    ROOT::Math::RhoEtaPhiVectorF genjet(gen_pt_values.at(j),
                                                        gen_eta_values.at(j),
                                                        gen_phi_values.at(j));
                    Logger::get("JetEnergyResolution")
                        ->debug("Checking gen Jet:  Eta: {} Phi:", genjet.Eta(),
                                genjet.Phi());
                    auto deltaR =
                        ROOT::Math::VectorUtil::DeltaR(jet, genjet) < 0.4;
                    if (deltaR > min_dR)
                        continue;
                    if (deltaR < 0.4 &&
                        std::abs(pt_values_corrected.at(i) -
                                 gen_pt_values.at(j)) >
                            3.0 * reso * pt_values_corrected.at(i)) {
                        min_dR = deltaR;
                        genjetpt = gen_pt_values.at(j);
                    }
                }
                if (genjetpt > 0.0) { // matched gen jet
                    Logger::get("JetEnergyResolution")
                        ->debug("Found gen jet for hybrid smearing method");
                    pt_values_corrected.at(i) +=
                        (resoSF - 1.0) * (pt_values_corrected.at(i) - genjetpt);
                } else {
                    Logger::get("JetEnergyResolution")
                        ->debug(
                            "No gen jet found. Applying stochastic smearing.");
                    TRandom3 randm = TRandom3(
                        static_cast<int>((eta_values.at(i) + 5) * 1000) * 1000 +
                        static_cast<int>((phi_values.at(i) + 4) * 1000) +
                        10000);
                    double shift =
                        randm.Gaus(0, reso) *
                        std::sqrt(std::max(resoSF * resoSF - 1, 0.0f));
                    pt_values_corrected.at(i) *= std::max(0.0, 1.0 + shift);
                }
                Logger::get("JetEnergyResolution")
                    ->debug("Shifting jet pt from {} to {} ",
                            pt_values.at(i) + pt_scale_shift,
                            pt_values_corrected.at(i));
                // if (pt_values_corrected.at(i)>15.0), this
                // correction should be propagated to MET
                // (requirement for type I corrections)
            }
            return pt_values_corrected;
        };
    auto df1 = df.Define(
        corrected_jet_pt, JetEnergyCorrectionLambda,
        {jet_pt, jet_eta, jet_phi, gen_jet_pt, gen_jet_eta, gen_jet_phi, rho});
    return df1;
}

/// Function to select jets passing a ID requirement, using
/// basefunctions::FilterMin
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the rawID column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe
/// \param[in] idThreshold minimal ID value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutRawID(ROOT::RDF::RNode df, const std::string &quantity,
                          const std::string &maskname,
                          const float &idThreshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMin(idThreshold), {quantity});
    return df1;
}
} // end namespace jet
} // end namespace physicsobject

namespace quantities {
namespace jet {
/// Function to determine number of jets in a jet collection
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] jetcollection name of the vector that contains jet indices of the
/// jets belonging to the collection, its length constitutes the output quantity
///
/// \return a dataframe containing the number of jets in the jet collection
ROOT::RDF::RNode NumberOfJets(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &jetcollection) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &jetcollection) {
                         Logger::get("NumberOfJets")->debug("Counting jets");
                         Logger::get("NumberOfJets")
                             ->debug("NJets {}", jetcollection.size());
                         return (int)jetcollection.size();
                     },
                     {jetcollection});
}
/// Function to writeout the value of the btagger for a jet. The tag value is
/// identified via the discriminant of e.g. the DeepJet tagger from nanoAOD
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] btagcolumn name of the column that contains btag values of
/// the jets
/// \param[in] jetcollection name of the vector that contains jet indices of the
/// jets belonging to the collection, its length constitutes the output quantity
/// \param position The position in the jet collection vector, which is used to
/// store the index of the particle in the particle quantity vectors.
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode btagValue(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &btagcolumn,
                           const std::string &jetcollection,
                           const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<float> &btagvalues,
                                const ROOT::RVec<int> &jetcollection) {
                         float btagValue = default_float;
                         try {
                             const int index = jetcollection.at(position);
                             btagValue = btagvalues.at(index);
                         } catch (const std::out_of_range &e) {
                         }
                         return btagValue;
                     },
                     {btagcolumn, jetcollection});
}
} // end namespace jet
} // end namespace quantities
#endif /* GUARDJETS_H */