#include "ROOT/RDataFrame.hxx"
#include "basefunctions.hxx"
#include "utility/utility.hxx"
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
/// Namespace containing function to apply filters on physics objects. The
/// filter results are typically stored within a mask, which is represented by
/// an `ROOT::RVec<int>`.
///    \code
///    In the mask
///    1 --> filter is passed by the object
///    0 --> filter is not passed by the object
///    \endcode
/// multiple filters can be combined by multiplying masks using
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
auto CutPt(auto df, const std::string quantity, const std::string maskname,
           const float ptThreshold) {
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
auto CutEta(auto df, const std::string quantity, const std::string maskname,
            const float EtaThreshold) {
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
auto CutDz(auto df, const std::string quantity, const std::string maskname,
           const float Threshold) {
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
auto CutDxy(auto df, const std::string quantity, const std::string maskname,
            const float Threshold) {
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
auto CombineMasks(auto df, const std::string &maskname,
                  const Masks &... masks) {
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

/// Function to filter events based on a mask. If the mask contains at least
/// one object, the event is filtered
///
///   \param df the input dataframe
///   \param maskname the name of the column containing the vetomap
///    to be used
///
///   \return a new df with the events filtered
auto FilterMasks(auto df, const std::string maskname) {
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
auto LeptonVetoFlag(auto df, const std::string &outputname,
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
auto ObjectMassCorrectionWithPt(auto df, const std::string corrected_mass,
                                const std::string raw_mass,
                                const std::string raw_pt,
                                const std::string corrected_pt) {
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

/// Muon specific functions
namespace muon {
/// Function to filter muons based on the muon ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
auto FilterID(auto df, const std::string maskname, const std::string nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
/// Function to filter muons based on the muon isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
auto FilterIsolation(auto df, const std::string maskname,
                     const std::string isolationName, const float Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}

} // end namespace muon
/// Tau specific functions
namespace tau {
/// Function to filter taus based on the tau decay mode
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] SelectedDecayModes a `std::vector<int>` containing the
/// decay modes, that should pass the filter
///
/// \return a dataframe containing the new mask
auto FilterDecayModes(auto df, const std::string maskname,
                      const std::vector<int> SelectedDecayModes) {
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
        {"Tau_decayMode"});
    return df1;
}
/// Function to filter taus based on the tau ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
auto FilterTauID(auto df, const std::string maskname, const std::string nameID,
                 const int idxID) {
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
auto PtCorrection(auto df, const std::string corrected_pt, const std::string pt,
                  const std::string decayMode, const float sf_dm0,
                  const float sf_dm1, const float sf_dm10,
                  const float sf_dm11) {
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

} // end namespace tau

namespace electron {
/// Function to filter electrons based on the electron ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
auto FilterID(auto df, const std::string maskname, const std::string nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
/// Function to filter electrons based on the electron isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
auto FilterIsolation(auto df, const std::string maskname,
                     const std::string isolationName, const float Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}

} // end namespace electron

} // namespace physicsobject