#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "basefunctions.hxx"
#include "correction.h"
#include "utility/utility.hxx"
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
/// Namespace containing function to apply cuts on physics objects. The
/// cut results are typically stored within a mask, which is represented by
/// an `ROOT::RVec<int>`.
///    \code
///    In the mask
///    1 --> cut is passed by the object
///    0 --> cut is not passed by the object
///    \endcode
/// multiple cuts can be combined by multiplying masks using
/// physicsobject::CombineMasks.
namespace physicsobject {
/// Function to select objects above a pt threshold, using
/// basefunctions::FilterMin
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the pt column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] ptThreshold minimal pt value
///
/// \return a dataframe containing the new mask
auto CutPt(auto &df, const std::string &quantity, const std::string &maskname,
           const float &ptThreshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMin(ptThreshold), {quantity});
    return df1;
}
/// Function to select objects blow an eta threshold, using
/// basefunctions::FilterAbsMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the eta column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] EtaThreshold maximal eta value
///
/// \return a dataframe containing the new mask
auto CutEta(auto &df, const std::string &quantity, const std::string &maskname,
            const float &EtaThreshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterAbsMax(EtaThreshold),
                         {quantity});
    return df1;
}
/// Function to select objects below an Dz threshold, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the Dz column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] Threshold maximal Dz value
///
/// \return a dataframe containing the new mask
auto CutDz(auto &df, const std::string &quantity, const std::string &maskname,
           const float &Threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterAbsMax(Threshold), {quantity});
    return df1;
}
/// Function to select objects below an Dxy threshold, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the Dxy column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] Threshold maximal Dxy value
///
/// \return a dataframe containing the new mask
auto CutDxy(auto &df, const std::string &quantity, const std::string &maskname,
            const float &Threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterAbsMax(Threshold), {quantity});
    return df1;
}
/// Function to combine a list of masks into a single mask. This is done be
/// multiplying all input masks
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe
/// \param[in] masks a parameter pack containing an arbitrary number of
/// `std::vector<std::string>` objects. Each string is the name of a mask to be
/// combined
///
/// \return a dataframe containing the new mask
template <class... Masks>
auto CombineMasks(auto &df, const std::string &maskname,
                  const Masks &...masks) {
    auto multiplyMasks = [](const ROOT::RVec<ROOT::RVec<int>> &x) {
        ROOT::RVec<int> result(x[0].size(), 1);
        for (auto &xx : x) {
            result *= xx;
        }
        return result;
    };
    // std::vector<std::string> MaskList{{masks...}}; does weird things in case
    // of two arguments in masks
    std::vector<std::string> MaskList;
    utility::appendParameterPackToVector(MaskList, masks...);
    const auto nMasks = sizeof...(Masks);
    return df.Define(
        maskname, ROOT::RDF::PassAsVec<nMasks, ROOT::RVec<int>>(multiplyMasks),
        MaskList);
}

/// Function to take a mask and create a new one where a tau candidate is set to
/// false
///
/// \param[in] df the input dataframe
/// \param[out] outputmaskname the name of the new mask to be added as column to
/// the dataframe \param[in] inputmaskname the name of the input mask \param[in]
/// ditaupair name of the column of the ditaupair \param[in] index index of the
/// tau candidate to be ignored by mask
///
/// \return a dataframe containing the new mask
auto VetoCandInMask(auto &df, const std::string &outputmaskname,
                    const std::string &inputmaskname,
                    const std::string &ditaupair, const int index) {
    return df.Define(outputmaskname,
                     [index, inputmaskname](const ROOT::RVec<int> &mask,
                                            const ROOT::RVec<int> &pair) {
                         Logger::get("VetoCandInMask")
                             ->debug("Vetoing the selected candidate (index "
                                     "{}) from the mask {}",
                                     index, inputmaskname);
                         auto newmask = mask;
                         if (pair.at(index) >= 0)
                             newmask.at(pair.at(index)) = 0;
                         return newmask;
                     },
                     {inputmaskname, ditaupair});
}

/// Function to filter events based on a mask. If the mask contains at least
/// one object, the event is filtered
///
///   \param df the input dataframe
///   \param maskname the name of the column containing the vetomap
///    to be used
///
///   \return a new df with the events filtered
auto FilterMasks(auto &df, const std::string &maskname) {
    auto df1 = df.Filter(
        [](const ROOT::RVec<Int_t> &mask) {
            auto result = Any(mask);
            return result;
        },
        {maskname});
    return df1;
}

/// Function to generate a veto based on a mask. If the mask contains at least
/// one object, the veto flag is set, meaning that this event contains at least
/// one object matching the requirements of the veto map.
///
/// \code
///  In the veto column
///  1 --> the event contains at least one object matching the veto map
///  0 --> the event does not contain any object matching the veto map
/// \endcode
///
/// \param df the input dataframe
/// \param outputname name of the new quantity containing the veto flags
/// \param vetomap the name of the column containing the vetomap to be used
///
/// \return a new df containing the veto flag column
auto LeptonVetoFlag(auto &df, const std::string &outputname,
                    const std::string &vetomap) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &mask) {
                         return ROOT::VecOps::Nonzero(mask).size() != 0;
                     },
                     {vetomap});
}

/// Function to correct object mass in alignment with object pt correction
///
/// \param[in] df the input dataframe
/// \param[out] corrected_mass the name of the corrected masses to be determined
/// \param[in] raw_mass name of the input mass \param[in] raw_pt name of the
/// uncorrected object pts \param[in] corrected_pt name of the corrected object
/// pts
///
/// \return a dataframe containing the modified object masses
auto ObjectMassCorrectionWithPt(auto &df, const std::string &corrected_mass,
                                const std::string &raw_mass,
                                const std::string &raw_pt,
                                const std::string &corrected_pt) {
    auto mass_correction_lambda =
        [](const ROOT::RVec<float> &mass_values,
           const ROOT::RVec<float> &pt_values,
           const ROOT::RVec<float> &corrected_pt_values) {
            ROOT::RVec<float> corrected_mass_values(mass_values.size());
            for (int i = 0; i < mass_values.size(); i++) {
                corrected_mass_values[i] = mass_values.at(i) *
                                           corrected_pt_values.at(i) /
                                           pt_values.at(i);
            }
            return corrected_mass_values;
        };
    auto df1 = df.Define(corrected_mass, mass_correction_lambda,
                         {raw_mass, raw_pt, corrected_pt});
    return df1;
}

/// Function to check whether at least one lepton pair is present
///
/// \param[in] df the input dataframe
/// \param[out] output_flag the name of the bool column that is created
/// \param[in] leptons_pt name of the input pt column of the lepton collection
/// \param[in] leptons_eta name of the input eta column of the lepton collection
/// \param[in] leptons_phi name of the input phi column of the lepton collection
/// \param[in] leptons_mass name of the input mass column of the lepton
/// collection \param[in] leptons_charge name of the input charge column of the
/// lepton collection \param[in] leptons_mask name of the input mask column of
/// the lepton collection that marks lepton to be taken into account \param[in]
/// dR_cut minimum required angular distance between the leptons
///
/// \return a dataframe containing the new bool column
auto CheckForDiLeptonPairs(
    auto &df, const std::string &output_flag, const std::string &leptons_pt,
    const std::string &leptons_eta, const std::string &leptons_phi,
    const std::string &leptons_mass, const std::string &leptons_charge,
    const std::string &leptons_mask, const float dR_cut) {
    auto pair_finder_lambda = [dR_cut](const ROOT::RVec<float> &pt_values,
                                       const ROOT::RVec<float> &eta_values,
                                       const ROOT::RVec<float> &phi_values,
                                       const ROOT::RVec<float> &mass_values,
                                       const ROOT::RVec<int> &charge_values,
                                       const ROOT::RVec<int> &mask) {
        const auto valid_lepton_indices = ROOT::VecOps::Nonzero(mask);
        for (auto it1 = valid_lepton_indices.begin();
             it1 != valid_lepton_indices.end(); it1++) {
            for (auto it2 = it1 + 1; it2 != valid_lepton_indices.end(); it2++) {
                if (charge_values.at(*it1) != charge_values.at(*it2)) {
                    auto p4_1 = ROOT::Math::PtEtaPhiMVector(
                        pt_values.at(*it1), eta_values.at(*it1),
                        phi_values.at(*it1), mass_values.at(*it1));
                    auto p4_2 = ROOT::Math::PtEtaPhiMVector(
                        pt_values.at(*it2), eta_values.at(*it2),
                        phi_values.at(*it2), mass_values.at(*it2));
                    if (ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2) >= dR_cut)
                        return true;
                }
            }
        }
        return false;
    };
    auto df1 = df.Define(output_flag, pair_finder_lambda,
                         {leptons_pt, leptons_eta, leptons_phi, leptons_mass,
                          leptons_charge, leptons_mask});
    return df1;
}

/// Muon specific functions
namespace muon {
/// Function to cut on muons based on the muon ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
auto CutID(auto &df, const std::string &maskname, const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
/// Function to cut on muons based on the muon isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
auto CutIsolation(auto &df, const std::string &maskname,
                  const std::string &isolationName, const float &Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}

} // end namespace muon
/// Tau specific functions
namespace tau {
/// Function to cut on taus based on the tau decay mode
///
/// \param[in] df the input dataframe
/// \param[in] tau_dms name of the column with tau decay modes
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] SelectedDecayModes a `std::vector<int>` containing the
/// decay modes, that should pass the cut
///
/// \return a dataframe containing the new mask
auto CutDecayModes(auto &df, const std::string &maskname,
                   const std::string &tau_dms,
                   const std::vector<int> &SelectedDecayModes) {
    auto df1 = df.Define(
        maskname,
        [SelectedDecayModes](const ROOT::RVec<Int_t> &decaymodes) {
            ROOT::RVec<int> mask;
            for (auto n : decaymodes) {
                mask.push_back(int(std::find(SelectedDecayModes.begin(),
                                             SelectedDecayModes.end(),
                                             n) != SelectedDecayModes.end()));
            }
            return mask;
        },
        {tau_dms});
    return df1;
}
/// Function to cut taus based on the tau ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
auto CutTauID(auto &df, const std::string &maskname, const std::string &nameID,
              const int &idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterID(idxID), {nameID});
    return df1;
}
/// Function to correct tau pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt \param[in] decayMode
/// \param[in] sf_dm0 scale factor to be applied to taus with decay mode 0
/// \param[in] sf_dm1 scale factor to be applied to other 1 prong taus
/// \param[in] sf_dm10 scale factor to be applied to taus with decay mode 10
/// \param[in] sf_dm11 scale factor to be applied to other 3 prong taus
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
auto PtCorrection_byValue(auto &df, const std::string &corrected_pt,
                          const std::string &pt, const std::string &decayMode,
                          const float &sf_dm0, const float &sf_dm1,
                          const float &sf_dm10, const float &sf_dm11) {
    auto tau_pt_correction_lambda =
        [sf_dm0, sf_dm1, sf_dm10, sf_dm11](const ROOT::RVec<float> &pt_values,
                                           const ROOT::RVec<int> &decay_modes) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (decay_modes.at(i) == 0)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm0;
                else if (decay_modes.at(i) > 0 && decay_modes.at(i) < 5)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm1;
                else if (decay_modes.at(i) == 10)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm10;
                else if (decay_modes.at(i) > 10 && decay_modes.at(i) < 15)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm11;
                else
                    corrected_pt_values[i] = pt_values.at(i);
            }
            return corrected_pt_values;
        };
    auto df1 =
        df.Define(corrected_pt, tau_pt_correction_lambda, {pt, decayMode});
    return df1;
}
/// Function to correct tau pt with correctionlib
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] DM0 variation decay mode 0
/// \param[in] DM1 variation decay mode 1
/// \param[in] DM10 variation decay mode 10
/// \param[in] DM11 variation decay mode 11
/// DM values: "nom","up","down"
///
/// \return a dataframe containing the new mask
auto PtCorrection(auto &df, const std::string &corrected_pt,
                  const std::string &pt, const std::string &eta,
                  const std::string &decayMode, const std::string &genMatch,
                  const std::string &sf_file, const std::string &jsonESname,
                  const std::string &idAlgorithm, const std::string &DM0,
                  const std::string &DM1, const std::string &DM10,
                  const std::string &DM11,
                  const std::vector<int> &SelectedDMs) {
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda =
        [evaluator, idAlgorithm, DM0, DM1, DM10, DM11,
         SelectedDMs](const ROOT::RVec<float> &pt_values,
                      const ROOT::RVec<float> &eta_values,
                      const ROOT::RVec<int> &decay_modes,
                      const ROOT::RVec<UChar_t> &genmatch) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM0});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM1});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 10) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM10});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 11) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM11});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
                Logger::get("tauEnergyCorrection")
                    ->debug("tau pt before {}, tau pt after {}",
                            pt_values.at(i), corrected_pt_values.at(i));
            }
            return corrected_pt_values;
        };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
} // end namespace tau

namespace electron {
/// Function to cut electrons based on the electron MVA ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
auto CutID(auto &df, const std::string &maskname, const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
} /// Function to cut jets based on the cut based electron ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// IDvalue value of the WP the has to be passed
///
/// \return a dataframe containing the new mask
auto CutCBID(auto &df, const std::string &maskname, const std::string &nameID,
             const int &IDvalue) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMinInt(IDvalue), {nameID});
    return df1;
}
/// Function to cut electrons based on the electron isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
auto CutIsolation(auto &df, const std::string &maskname,
                  const std::string &isolationName, const float &Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}

} // end namespace electron

} // namespace physicsobject