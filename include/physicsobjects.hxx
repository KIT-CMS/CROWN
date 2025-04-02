#ifndef GUARD_PHYSICSOBJECTS_H
#define GUARD_PHYSICSOBJECTS_H

#include "../include/utility/utility.hxx"

namespace physicsobject {

/**
 * @brief This function defines a mask for objects that satisfy a minimum threshold 
 * requirement. The mask is created by comparing the values in the specified 
 * column with the given threshold, marking elements as `1` if they pass 
 * the cut and `0` otherwise.
 *
 * @tparam T type of the threshold and input quantity (e.g. `float`, `int`)
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the object column in the NanoAOD for which the 
 * cut should be applied, expected to be of type `ROOT::RVec<T>`
 * @param threshold minimum threshold value of type `T`
 * 
 * @return a dataframe containing the new mask as a column
 */
template <typename T> 
inline ROOT::RDF::RNode CutMin(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &quantity, const T &threshold) {
    return df.Define(outputname, 
        [threshold](const ROOT::RVec<T> &values) {
            ROOT::RVec<int> mask = values >= threshold;
            return mask;
        }, 
        {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy a minimum threshold 
 * requirement. The mask is created by comparing the absolute values in the specified 
 * column with the given threshold, marking elements as `1` if they pass 
 * the cut and `0` otherwise.
 *
 * @tparam T type of the threshold and input quantity (e.g. `float`, `int`)
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the object column in the NanoAOD for which the 
 * cut should be applied, expected to be of type `ROOT::RVec<T>`
 * @param threshold minimum threshold value of type `T`
 * 
 * @return a dataframe containing the new mask as a column
 */
template <typename T> 
inline ROOT::RDF::RNode CutAbsMin(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &quantity, const T &threshold) {
    return df.Define(outputname, 
        [threshold](const ROOT::RVec<T> &values) {
            ROOT::RVec<int> mask = abs(values) >= threshold;
            return mask;
        }, 
        {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy a maximum threshold 
 * requirement. The mask is created by comparing the values in the specified 
 * column with the given threshold, marking elements as `1` if they pass 
 * the cut and `0` otherwise.
 *
 * @tparam T type of the threshold and input quantity (e.g. `float`, `int`)
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the object column in the NanoAOD for which the 
 * cut should be applied, expected to be of type `ROOT::RVec<T>`
 * @param threshold maximum threshold value of type `T`
 * 
 * @return a dataframe containing the new mask as a column
 */
template <typename T> 
inline ROOT::RDF::RNode CutMax(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &quantity, const T &threshold) {
    return df.Define(outputname, 
        [threshold](const ROOT::RVec<T> &values) {
            ROOT::RVec<int> mask = values < threshold;
            return mask;
        }, 
        {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy a maximum threshold 
 * requirement. The mask is created by comparing the absolute values in the specified 
 * column with the given threshold, marking elements as `1` if they pass 
 * the cut and `0` otherwise.
 *
 * @tparam T type of the threshold and input quantity (e.g. `float`, `int`)
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the object column in the NanoAOD for which the 
 * cut should be applied, expected to be of type `ROOT::RVec<T>`
 * @param threshold maximum threshold value of type `T`
 * 
 * @return a dataframe containing the new mask as a column
 */
template <typename T> 
inline ROOT::RDF::RNode CutAbsMax(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &quantity, const T &threshold) {
    return df.Define(outputname, 
        [threshold](const ROOT::RVec<T> &values) {
            ROOT::RVec<int> mask = abs(values) < threshold;
            return mask;
        }, 
        {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy an exact threshold 
 * requirement. The mask is created by comparing the values in the specified 
 * column with the given threshold, marking elements as `1` if they pass 
 * the cut and `0` otherwise.
 *
 * @tparam T type of the threshold and input quantity (e.g. `float`, `int`)
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the object column in the NanoAOD for which the 
 * cut should be applied, expected to be of type `ROOT::RVec<T>`
 * @param threshold exact threshold value of type `T`
 * 
 * @return a dataframe containing the new mask as a column
 */
template <typename T> 
inline ROOT::RDF::RNode CutEqual(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &quantity, const T &threshold) {
    return df.Define(outputname, 
        [threshold](const ROOT::RVec<T> &values) {
            ROOT::RVec<int> mask = values == threshold;
            return mask;
        }, 
        {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy an exact threshold 
 * requirement. The mask is created by comparing the absolute values in the specified 
 * column with the given threshold, marking elements as `1` if they pass 
 * the cut and `0` otherwise.
 *
 * @tparam T type of the threshold and input quantity (e.g. `float`, `int`)
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the object column in the NanoAOD for which the 
 * cut should be applied, expected to be of type `ROOT::RVec<T>`
 * @param threshold exact threshold value of type `T`
 * 
 * @return a dataframe containing the new mask as a column
 */
template <typename T> 
inline ROOT::RDF::RNode CutAbsEqual(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &quantity, const T &threshold) {
    return df.Define(outputname, 
        [threshold](const ROOT::RVec<T> &values) {
            ROOT::RVec<int> mask = abs(values) == threshold;
            return mask;
        }, 
        {quantity});
}
/**
 * @brief This function takes multiple masks and applies a logical operation 
 * ("any", "all", or "none") elemet-wise to generate a combined mask. The 
 * function ensures that elements are correctly merged based on the given mode:
 *
 * - "any": The resulting mask contains true values if any element of 
 * the input masks is true (element-wise)
 * - "all": The resulting mask contains true values if all elements of 
 * the input masks are true (element-wise)
 * - "none": The resulting mask contains true values if no element of 
 * the input masks is true (element-wise)
 *
 * @tparam Args variadic template parameter pack representing mask columns plus mode
 * @param df input dataframe
 * @param outputname name of the output column containing the combined mask
 * @param args parameter pack of column names that contain the considered masks of 
 * type `ROOT::RVec<int>`, with the last argument being the mode (`"any"`, `"all"`, 
 * or `"none"`)
 *
 * @return a dataframe containing the new mask as a column
 *
 * @note The masks must have the same size, as element-wise operations are performed.
 */
template <class... Args>
inline ROOT::RDF::RNode CombineMasks(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const Args ...args) {
    
    auto argTuple = std::make_tuple(args...);
    auto mode = utility::extractLastArgument(argTuple);

    std::vector<std::string> MaskList{args...};
    MaskList.pop_back();  
    const auto nMasks = sizeof...(Args) - 1;

    auto lambda = [mode](const ROOT::RVec<ROOT::RVec<int>> &masks) {
        if (mode == std::string("any")) {
            ROOT::RVec<int> result(masks[0].size(), 0);
            for (auto &mask : masks) {
                result += mask;
            }
            result = ROOT::VecOps::Map(result, [](int x) { return x != 0; });
            return result;
        }
        else if (mode == std::string("all")) {
            ROOT::RVec<int> result(masks[0].size(), 1);
            for (auto &mask : masks) {
                result *= mask;
            }
            return result;
        }
        else if (mode == std::string("none")) {
            ROOT::RVec<int> result(masks[0].size(), 0);
            for (auto &mask : masks) {
                result += mask;
            }
            result = ROOT::VecOps::Map(result, [](int x) { return x != 0; });
            return result;
        }
        else {
            Logger::get("physicsobject::CombineMasks")->error("Mode {} is not defined!", mode);
            throw std::runtime_error("Mode is not defined!");
        }        
    };
    return df.Define(
        outputname, ROOT::RDF::PassAsVec<nMasks, ROOT::RVec<int>>(lambda),
        MaskList);
}

ROOT::RDF::RNode CutQuantityBarrelEndcap(ROOT::RDF::RNode df, 
                 const std::string &outputname, const std::string &eta, 
                 const std::string &quantity, const float &barrel_endcap_boundary, 
                 const float &lower_threshold_barrel, const float &upper_threshold_barrel, 
                 const float &lower_threshold_endcap, const float &upper_threshold_endcap);
ROOT::RDF::RNode VetoSingleObject(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &object_mask,
                            const int index);
ROOT::RDF::RNode VetoSingleObject(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &object_mask,
                            const std::string &index_vector, 
                            const int position);
ROOT::RDF::RNode Number(ROOT::RDF::RNode df, 
                        const std::string &outputname,
                        const std::string &object_mask);
ROOT::RDF::RNode NumberFlag(ROOT::RDF::RNode df, 
                            const std::string &outputname,
                            const std::string &object_mask, 
                            const int &number);
ROOT::RDF::RNode Veto(ROOT::RDF::RNode df,
                      const std::string &outputname,
                      const std::string &object_mask);
ROOT::RDF::RNode GetIndices(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &inputmaskname);
ROOT::RDF::RNode OverlapVeto(ROOT::RDF::RNode df, 
                 const std::string &outputname, const std::string &target_p4,
                 const std::string &object_pt, const std::string &object_eta, 
                 const std::string &object_phi, const std::string &object_mass, 
                 const std::string &object_mask, const float delta_r_cut);
ROOT::RDF::RNode MassCorrectionWithPt(ROOT::RDF::RNode df,
                                      const std::string &corrected_mass,
                                      const std::string &raw_mass,
                                      const std::string &raw_pt,
                                      const std::string &corrected_pt);
ROOT::RDF::RNode VetoLeptonPairs(ROOT::RDF::RNode df, 
                 const std::string &outputname, const std::string &lepton_pt, 
                 const std::string &lepton_eta, const std::string &lepton_phi, 
                 const std::string &lepton_mass, const std::string &lepton_charge, 
                 const std::string &lepton_mask, const float delta_r_cut);
} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */
