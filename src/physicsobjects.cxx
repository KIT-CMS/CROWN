#ifndef GUARD_PHYSICSOBJECTS_H
#define GUARD_PHYSICSOBJECTS_H

#include "../include/basefilters.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

/** 
 * This namespace contains functions to apply cuts on physics objects. The
 * cut results are typically stored within a mask of the same length as the 
 * physics objects vector and is represented by a `ROOT::RVec<int>`.
 *    \code
 *    In the mask
 *    1 --> object passes a cut
 *    0 --> object does not pass a cut
 *    \endcode
 * Multiple cuts can be combined by multiplying masks using
 * physicsobject::CombineMasks.
 */
namespace physicsobject {

/**
 * @brief This function defines a mask in the dataframe that selects objects based on 
 * eta-dependent upper and lower thresholds for a given object `quantity`. The selection 
 * criteria differ between the barrel and endcap regions of the detector.
 *
 * @param df input dataframe
 * @param outputname name of the output mask column
 * @param eta name of the eta column of a physics object in the dataframe
 * @param quantity name of the quantity column to apply the selection on
 * @param barrel_endcap_boundary absolute eta boundary separating the 
 * barrel and endcap regions
 * @param lower_threshold_barrel lower threshold for selection in the barrel region
 * @param upper_threshold_barrel upper threshold for selection in the barrel region
 * @param lower_threshold_endcap lower threshold for selection in the endcap region
 * @param upper_threshold_endcap upper threshold for selection in the endcap region
 *
 * @return a dataframe containing the new mask as a column
 */
ROOT::RDF::RNode CutQunatityBarrelEndcap(ROOT::RDF::RNode df, 
                                         const std::string &outputname,
                                         const std::string &eta, 
                                         const std::string &quantity,
                                         const float &barrel_endcap_boundary, 
                                         const float &lower_threshold_barrel,
                                         const float &upper_threshold_barrel, 
                                         const float &lower_threshold_endcap,
                                         const float &upper_threshold_endcap) {
    auto lambda = [barrel_endcap_boundary, lower_threshold_barrel, 
                   upper_threshold_barrel, lower_threshold_endcap, 
                   upper_threshold_endcap](
                     const ROOT::RVec<float> &etas,
                     const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask =
            (
                (
                    (abs(etas) < barrel_endcap_boundary) && 
                    (values >= lower_threshold_barrel) &&
                    (values < upper_threshold_barrel)
                ) ||
                (
                    (abs(etas) >= barrel_endcap_boundary) &&
                    (values >= lower_threshold_endcap) &&
                    (values < upper_threshold_endcap)
                )
            );
        return mask;
    };

    auto df1 = df.Define(outputname, lambda, {eta, quantity});
    return df1;
}

/**
 * @brief This function updates the mask by setting the value at the specified index 
 * to zero, effectively vetoing the object at that position.
 *
 * @param df input dataframe
 * @param outputname name of the output mask column after vetoing the object
 * @param object_mask name of the mask column that contains object selection flags
 * @param index index of an object that specifies which object to veto
 *
 * @return a dataframe containing the new mask as a column
 */
ROOT::RDF::RNode VetoSingleObject(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &object_mask,
                                  const int index) {
    return df.Define(outputname,
        [index, object_mask](const ROOT::RVec<int> &mask) {
            Logger::get("VetoSingleObject")
                ->debug("Vetoing the object at index {} from the mask {}",
                        index, object_mask);
            auto new_mask = mask;
            new_mask.at(index) = 0;
            return new_mask;
        },
        {object_mask});
}

/**
 * @brief This function updates the mask by setting the value at the specified index 
 * to zero, effectively vetoing the object at that position.
 *
 * @param df input dataframe
 * @param outputname name of the output mask column after vetoing the object
 * @param object_mask name of the mask column that contains object selection flags
 * @param index_vector name of the column containing indices of objects
 * @param position position in the index vector that specifies which object to veto
 *
 * @return a dataframe containing the new mask as a column
 */
ROOT::RDF::RNode VetoSingleObject(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &object_mask,
                                  const std::string &index_vector,
                                  const int position) {
    return df.Define(outputname,
        [position, object_mask](const ROOT::RVec<int> &mask,
                            const ROOT::RVec<int> &index_vector) {
            Logger::get("VetoSingleObject")
                ->debug("Vetoing the object at index {} from the mask {}",
                        index_vector.at(position), object_mask);
            auto new_mask = mask;
            if (index_vector.at(position) >= 0)
                new_mask.at(index_vector.at(position)) = 0;
            return new_mask;
        },
        {object_mask, index_vector});
}

/**
 * @brief This function counts the number of selected objects in the provided object 
 * mask and stores the number.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the object count
 * @param object_mask name of the mask column indicating selected objects
 *
 * @return a dataframe with a new column
 */
ROOT::RDF::RNode Number(ROOT::RDF::RNode df, 
                        const std::string &outputname,
                        const std::string &object_mask) {
    return df.Define(outputname,
        [](const ROOT::RVec<int> &mask) {
            int count = ROOT::VecOps::Nonzero(mask).size();
            return count;
        },
        {object_mask});
}
/**
 * @brief This function checks how many objects are selected based on the provided 
 * mask and returns `true` if the count matches the specified number.
 *
 * @param df input dataframe
 * @param outputname name of the output flag column
 * @param object_mask name of the mask column indicating selected objects
 * @param number expected number of selected objects
 *
 * @return a dataframe with a new flag column
 */
ROOT::RDF::RNode NumberFlag(ROOT::RDF::RNode df, 
                            const std::string &outputname,
                            const std::string &object_mask, 
                            const int &number) {
    return df.Define(outputname,
        [number](const ROOT::RVec<int> &mask) {
            return ROOT::VecOps::Nonzero(mask).size() == number;
        },
        {object_mask});
}

/**
 * @brief This function checks if there is at least one selected object in the provided 
 * object mask. If at least one object is selected, the output flag is set to `true`, 
 * otherwise, it is set to `false`. This flag can be used to veto events that include 
 * a specific object.
 *
 * @param[in] df input dataframe
 * @param[in] outputname name of the output column storing the veto flag
 * @param[in] object_mask name of the mask column indicating selected objects
 *
 * @return a dataframe with a new flag column
 */
ROOT::RDF::RNode Veto(ROOT::RDF::RNode df,
                      const std::string &outputname,
                      const std::string &object_mask) {
    return df.Define(outputname,
        [](const ROOT::RVec<int> &mask) {
            return ROOT::VecOps::Nonzero(mask).size() != 0;
        },
        {object_mask});
}

/**
 * @brief This function extracts the indices of objects that are marked as selected 
 * in the given object mask.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the object indices
 * @param object_mask name of the mask column indicating selected objects
 *
 * @return a dataframe with a new column
 */
ROOT::RDF::RNode GetIndices(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &object_mask) {
    return df.Define(outputname,
        [](const ROOT::RVec<int> &mask) {
            ROOT::VecOps::RVec<int> indices = ROOT::VecOps::Nonzero(mask);
            Logger::get("GetIndices")->debug("indices = {}", indices);
            return indices;
        },
        {object_mask});
}

/**
 * @brief This function checks if any object in a given collection (`object_mask`) 
 * overlaps with a specified target object defined by its four-momentum within a 
 * user-defined `delta_r_cut` distance. If an object is found within this threshold, 
 * the function returns `true`, indicating an overlap with the target object. The 
 * result can be used to veto events that have overlapping objects.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the overlap veto flag
 * @param target_p4 name of the column containing the target four-momentum as a 
 * `ROOT::Math::PtEtaPhiMVector`
 * @param object_pt name of the column containing object transverse momenta
 * @param object_eta name of the column containing object pseudorapidities
 * @param object_phi name of the column containing object azimuthal angles
 * @param object_mass name of the column containing object masses
 * @param object_mask name of the column containing the mask of selected objects to 
 * compare with
 * @param delta_r_cut minimal deltaR distance allowed between objects and target to 
 * count as not overlapping
 *
 * @return a dataframe with a new flag column
 */
ROOT::RDF::RNode OverlapVeto(ROOT::RDF::RNode df, 
                             const std::string &outputname, 
                             const std::string &target_p4,
                             const std::string &object_pt,
                             const std::string &object_eta, 
                             const std::string &object_phi,
                             const std::string &object_mass, 
                             const std::string &object_mask, 
                             const float delta_r_cut) {
    auto veto_overlapping_target =
        [delta_r_cut](const ROOT::Math::PtEtaPhiMVector &p4,
                     const ROOT::RVec<float> &pts,
                     const ROOT::RVec<float> &etas,
                     const ROOT::RVec<float> &phis,
                     const ROOT::RVec<float> &masses,
                     const ROOT::RVec<int> &mask) {

            ROOT::RVec<int> selected_indices =
                ROOT::VecOps::Nonzero(mask);
            const auto selected_pts =
                ROOT::VecOps::Take(pts, selected_indices);
            const auto selected_etas =
                ROOT::VecOps::Take(etas, selected_indices);
            const auto selected_phis =
                ROOT::VecOps::Take(phis, selected_indices);
            const auto selected_masses =
                ROOT::VecOps::Take(masses, selected_indices);
            auto selected_p4s =
                ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
                    selected_pts, selected_etas, selected_phis, selected_masses);
            
            for (const auto &p4_test : selected_p4s) {
                if (ROOT::Math::VectorUtil::DeltaR(p4_test, p4) < delta_r_cut) {
                    return true;
                }
            }
            return false;
        };
    auto df1 = df.Define(outputname, veto_overlapping_target,
        {target_p4, object_pt, object_eta, object_phi, object_mass, object_mask});
    return df1;
}

/**
 * @brief This function modifies the mass of objects in an event using the formula
 * \f[
 * M_{\text{corrected},i} = M_{\text{raw},i} \times \frac{p_{T,\text{corrected},i}}{p_{T,\text{raw},i}}
 * \f]
 * for each object of an object collection in the event. The correction is applied 
 * element-wise to the mass vector and is needed as part of for example energy scale 
 * corrections that were before to the transverse momenta.
 *
 * @param df input dataframe
 * @param corrected_mass name of the output column storing the corrected masses
 * @param raw_mass name of the column containing raw object masses
 * @param raw_pt name of the column containing raw object transverse momenta
 * @param corrected_pt name of the column containing corrected transverse momenta
 *
 * @return a dataframe with a new column
 */
ROOT::RDF::RNode MassCorrectionWithPt(ROOT::RDF::RNode df,
                                      const std::string &corrected_mass,
                                      const std::string &raw_mass,
                                      const std::string &raw_pt,
                                      const std::string &corrected_pt) {
    auto mass_correction_lambda =
        [](const ROOT::RVec<float> &masses,
           const ROOT::RVec<float> &pts,
           const ROOT::RVec<float> &corrected_pts) {
                ROOT::RVec<float> corrected_masses(masses.size());
                for (int i = 0; i < masses.size(); i++) {
                    corrected_masses[i] = masses.at(i) *
                                          corrected_pts.at(i) /
                                          pts.at(i);
                }
                return corrected_masses;
            };
    auto df1 = df.Define(corrected_mass, mass_correction_lambda,
                         {raw_mass, raw_pt, corrected_pt});
    return df1;
}

/**
 * @brief This function checks for the presence of same-flavor opposite-sign 
 * lepton pairs with a given deltaR separation (`delta_r_cut`). If such a pair 
 * is found, a veto flag is set and can be later used to remove event with such
 * lepton pairs.
 *
 * @param df input dataframe
 * @param outputname name of the output flag column
 * @param lepton_pt name of the column containing lepton transverse momenta
 * @param lepton_eta name of the column containing lepton pseudorapidities
 * @param lepton_phi name of the column containing lepton azimuthal angles
 * @param lepton_mass name of the column containing lepton masses
 * @param lepton_charge name of the column containing lepton charges
 * @param lepton_mask name of the column containing the mask of selected lepton 
 * objects
 * @param delta_r_cut maximal deltaR distance allowed between lepton objects to 
 * not being vetoed
 *
 * @return a dataframe with a new flag column
 */
ROOT::RDF::RNode VetoLeptonPairs(ROOT::RDF::RNode df, 
                                 const std::string &outputname, 
                                 const std::string &lepton_pt, 
                                 const std::string &lepton_eta,
                                 const std::string &lepton_phi, 
                                 const std::string &lepton_mass,
                                 const std::string &lepton_charge,
                                 const std::string &lepton_mask, 
                                 const float delta_r_cut) {
    auto pair_finder_lambda = 
        [delta_r_cut](const ROOT::RVec<float> &pts,
                     const ROOT::RVec<float> &etas,
                     const ROOT::RVec<float> &phis,
                     const ROOT::RVec<float> &masses,
                     const ROOT::RVec<int> &charges,
                     const ROOT::RVec<int> &mask) {
        const auto selected_indices = ROOT::VecOps::Nonzero(mask);
        for (auto index_1 = selected_indices.begin();
             index_1 != selected_indices.end(); index_1++) {
            for (auto index_2 = index_1 + 1; 
                 index_2 != selected_indices.end(); index_2++) {
                if (charges.at(*index_1) != charges.at(*index_2)) {
                    auto p4_1 = ROOT::Math::PtEtaPhiMVector(
                        pts.at(*index_1), etas.at(*index_1),
                        phis.at(*index_1), masses.at(*index_1));
                    auto p4_2 = ROOT::Math::PtEtaPhiMVector(
                        pts.at(*index_2), etas.at(*index_2),
                        phis.at(*index_2), masses.at(*index_2));
                    if (ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2) >= delta_r_cut)
                        return true;
                }
            }
        }
        return false;
    };
    auto df1 = df.Define(outputname, pair_finder_lambda,
                         {lepton_pt, lepton_eta, lepton_phi, lepton_mass, 
                          lepton_charge, lepton_mask});
    return df1;
}
} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */
