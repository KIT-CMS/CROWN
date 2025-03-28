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
 * @brief This function defines a new mask by multiplying multiple input masks 
 * element-wise. Each mask is represented as a `ROOT::RVec<int>`, and the 
 * result is computed by iterating through all masks and performing 
 * element-wise multiplication.
 *
 * @tparam Masks types of the input mask columns, expected to be `ROOT::RVec<int>`
 * @param df input dataframe
 * @param outputname name of the new column containing the combined mask
 * @param masks parameter pack of column names containing the masks to be combined
 * 
 * @return a dataframe containing the new mask as a column
 *
 * @note The masks must have the same size, as element-wise multiplication is performed.
 */
template <class... Masks>
inline ROOT::RDF::RNode CombineMasks(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const Masks &...masks) {
    auto multiplyMasks = [](const ROOT::RVec<ROOT::RVec<int>> &masks) {
        ROOT::RVec<int> result(masks[0].size(), 1);
        for (auto &mask : masks) {
            result *= mask;
        }
        return result;
    };

    std::vector<std::string> MaskList;
    utility::appendParameterPackToVector(MaskList, masks...);
    const auto nMasks = sizeof...(Masks);
    
    return df.Define(
        outputname, ROOT::RDF::PassAsVec<nMasks, ROOT::RVec<int>>(multiplyMasks),
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
                 const std::string &object_mask, const float deltaR_cut);
ROOT::RDF::RNode MassCorrectionWithPt(ROOT::RDF::RNode df,
                                      const std::string &corrected_mass,
                                      const std::string &raw_mass,
                                      const std::string &raw_pt,
                                      const std::string &corrected_pt);
ROOT::RDF::RNode VetoLeptonPairs(ROOT::RDF::RNode df, 
                 const std::string &outputname, const std::string &lepton_pt, 
                 const std::string &lepton_eta, const std::string &lepton_phi, 
                 const std::string &lepton_mass, const std::string &lepton_charge, 
                 const std::string &lepton_mask, const float deltaR_cut);
} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */
