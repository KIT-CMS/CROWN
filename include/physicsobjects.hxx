#ifndef GUARD_PHYSICSOBJECTS_H
#define GUARD_PHYSICSOBJECTS_H

#include "../include/utility/utility.hxx"

namespace physicsobject {

/**
 * @brief This function takes multiple masks and applies a logical operation
 * (`"any_of"`, `"all_of"`, or `"none_of"`) elemet-wise to generate a combined mask.
 * The function ensures that elements are correctly merged based on the given
 * mode:
 *
 * - `"any_of"`: The resulting mask contains true values if any element of
 * the input masks is true (element-wise)
 * - `"all_of"`: The resulting mask contains true values if all elements of
 * the input masks are true (element-wise)
 * - `"none_of"`: The resulting mask contains true values if no element of
 * the input masks is true (element-wise)
 *
 * @tparam Args variadic template parameter pack representing mask columns plus
 * mode
 * @param df input dataframe
 * @param outputname name of the output column containing the combined mask
 * @param args parameter pack of column names that contain the considered masks
 * of type `ROOT::RVec<int>`, with the last argument being the mode (`"any_of"`,
 * `"all_of"`, or `"none_of"`)
 *
 * @return a dataframe containing the new mask as a column
 *
 * @note The masks must have the same size, as element-wise operations are
 * performed.
 */
template <class... Args>
inline ROOT::RDF::RNode CombineMasks(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const Args... args) {

    auto argTuple = std::make_tuple(args...);
    auto mode = utility::extractLastArgument(argTuple);

    std::vector<std::string> MaskList{args...};
    MaskList.pop_back();
    const auto nMasks = sizeof...(Args) - 1;

    auto lambda = [mode](const ROOT::RVec<ROOT::RVec<int>> &masks) {
        if (mode == std::string("any_of")) {
            ROOT::RVec<int> result(masks[0].size(), 0);
            for (auto &mask : masks) {
                result += mask;
            }
            result = ROOT::VecOps::Map(result, [](int x) { return x != 0; });
            return result;
        } else if (mode == std::string("all_of")) {
            ROOT::RVec<int> result(masks[0].size(), 1);
            for (auto &mask : masks) {
                result *= mask;
            }
            return result;
        } else if (mode == std::string("none_of")) {
            ROOT::RVec<int> result(masks[0].size(), 0);
            for (auto &mask : masks) {
                result += mask;
            }
            result = ROOT::VecOps::Map(result, [](int x) { return x != 0; });
            return result;
        } else {
            Logger::get("physicsobject::CombineMasks")
                ->error("Mode {} is not defined!", mode);
            throw std::runtime_error("Mode is not defined!");
        }
    };
    return df.Define(outputname,
                     ROOT::RDF::PassAsVec<nMasks, ROOT::RVec<int>>(lambda),
                     MaskList);
}

/**
 * @brief This function defines a mask for objects that satisfy a minimum
 * threshold requirement. The mask is created by comparing the values in the
 * specified column with the given threshold, marking elements as `1` if they
 * pass the cut and `0` otherwise.
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
inline ROOT::RDF::RNode
CutMin(ROOT::RDF::RNode df, const std::string &outputname,
       const std::string &quantity, const T &threshold) {
    return df.Define(outputname,
                     [threshold](const ROOT::RVec<T> &values) {
                         ROOT::RVec<int> mask = values >= threshold;
                         return mask;
                     },
                     {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy a minimum
 * threshold requirement. The mask is created by comparing the absolute values
 * in the specified column with the given threshold, marking elements as `1` if
 * they pass the cut and `0` otherwise.
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
inline ROOT::RDF::RNode
CutAbsMin(ROOT::RDF::RNode df, const std::string &outputname,
          const std::string &quantity, const T &threshold) {
    return df.Define(outputname,
                     [threshold](const ROOT::RVec<T> &values) {
                         ROOT::RVec<int> mask = abs(values) >= threshold;
                         return mask;
                     },
                     {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy a maximum
 * threshold requirement. The mask is created by comparing the values in the
 * specified column with the given threshold, marking elements as `1` if they
 * pass the cut and `0` otherwise.
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
inline ROOT::RDF::RNode
CutMax(ROOT::RDF::RNode df, const std::string &outputname,
       const std::string &quantity, const T &threshold) {
    return df.Define(outputname,
                     [threshold](const ROOT::RVec<T> &values) {
                         ROOT::RVec<int> mask = values < threshold;
                         return mask;
                     },
                     {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy a maximum
 * threshold requirement. The mask is created by comparing the absolute values
 * in the specified column with the given threshold, marking elements as `1` if
 * they pass the cut and `0` otherwise.
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
inline ROOT::RDF::RNode
CutAbsMax(ROOT::RDF::RNode df, const std::string &outputname,
          const std::string &quantity, const T &threshold) {
    return df.Define(outputname,
                     [threshold](const ROOT::RVec<T> &values) {
                         ROOT::RVec<int> mask = abs(values) < threshold;
                         return mask;
                     },
                     {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy an exact
 * threshold requirement. The mask is created by comparing the values in the
 * specified column with the given threshold, marking elements as `1` if they
 * pass the cut and `0` otherwise.
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
inline ROOT::RDF::RNode
CutEqual(ROOT::RDF::RNode df, const std::string &outputname,
         const std::string &quantity, const T &threshold) {
    return df.Define(outputname,
                     [threshold](const ROOT::RVec<T> &values) {
                         ROOT::RVec<int> mask = values == threshold;
                         return mask;
                     },
                     {quantity});
}

/**
 * @brief This function defines a mask for objects that satisfy an exact
 * threshold requirement. The mask is created by comparing the absolute values
 * in the specified column with the given threshold, marking elements as `1` if
 * they pass the cut and `0` otherwise.
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
inline ROOT::RDF::RNode
CutAbsEqual(ROOT::RDF::RNode df, const std::string &outputname,
            const std::string &quantity, const T &threshold) {
    return df.Define(outputname,
                     [threshold](const ROOT::RVec<T> &values) {
                         ROOT::RVec<int> mask = abs(values) == threshold;
                         return mask;
                     },
                     {quantity});
}

/**
 * @brief This function checks whether a specified bit (given by `threshold`) is
 * set in each element of the input column and creates a new mask column with
 * values of `1` if the bit is set, and `0` otherwise. If `threshold` is `0` or
 * negative, all values in the mask are set to `1`.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the input column containing the bitmasks of
 * type `ROOT::RVec<UChar_t>`
 * @param threshold bit position to check in each bitmask
 *
 * @return a dataframe containing the new mask as a column
 */
inline ROOT::RDF::RNode CutBitmask(ROOT::RDF::RNode df,
                                   const std::string &outputname,
                                   const std::string &quantity,
                                   const int &threshold) {
    return df.Define(outputname,
                     [threshold](const ROOT::RVec<UChar_t> &values) {
                         ROOT::RVec<int> mask;

                         for (auto const value : values) {
                             if (threshold > 0)
                                 mask.push_back(std::min(
                                     1, int(value & 1 << (threshold - 1))));
                             else
                                 mask.push_back(int(1));
                         }
                         return mask;
                     },
                     {quantity});
}

/**
 * @brief This function compares each element in the input column against a
 * provided `selection` list and returns a mask column where each entry is `1`
 * if the value is in the `selection` list and `0` otherwise.
 *
 * @tparam T type of values in the input column and the selection list
 * @param df input dataframe
 * @param outputname name of the new column containing the selected object mask
 * @param quantity name of the input column containing values to be checked of
 * type `ROOT::RVec<T>`
 * @param selection a vector containing the selection of values of type `T`
 *
 * @return a dataframe containing the new mask as a column
 */
template <typename T>
inline ROOT::RDF::RNode
CutQuantity(ROOT::RDF::RNode df, const std::string &outputname,
            const std::string &quantity, const std::vector<T> &selection) {
    return df.Define(
        outputname,
        [selection](const ROOT::RVec<T> &values) {
            ROOT::RVec<int> mask;
            for (auto const value : values) {
                mask.push_back(int(std::find(selection.begin(), selection.end(),
                                             value) != selection.end()));
            }
            return mask;
        },
        {quantity});
}

/**
 * @brief This function returns the length of the vector containers in columns
 * that contain `ROOT::VecOps::RVec<T>` objects. The function is templated
 * with `T`, which must be set to the underlying type of the vector object.
 *
 * @param df input dataframe
 * @param outputname name of the output column storing the object count
 * @param vector_quantity name of the quantity stored in vector structures
 *
 * @return a dataframe with a new column
 */
template <typename T>
ROOT::RDF::RNode Size(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &vector_quantity) {
    return df.Define(outputname,
                     [](const ROOT::RVec<T> &quantity) {
                         int length = quantity.size();
                         return length;
                     },
                     {vector_quantity});
}

ROOT::RDF::RNode CutQuantityBarrelEndcap(
    ROOT::RDF::RNode df, const std::string &outputname, const std::string &eta,
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
ROOT::RDF::RNode Count(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &object_mask);
ROOT::RDF::RNode CountFlag(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &object_mask, const int &number);
ROOT::RDF::RNode Veto(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &object_mask);
ROOT::RDF::RNode
OverlapVeto(ROOT::RDF::RNode df, const std::string &outputname,
            const std::string &target_p4, const std::string &object_pt,
            const std::string &object_eta, const std::string &object_phi,
            const std::string &object_mass, const std::string &object_mask,
            const float min_delta_r);
ROOT::RDF::RNode
LeptonPairVeto(ROOT::RDF::RNode df, const std::string &outputname,
               const std::string &lepton_pt, const std::string &lepton_eta,
               const std::string &lepton_phi, const std::string &lepton_mass,
               const std::string &lepton_charge, const std::string &lepton_mask,
               const float max_delta_r);
ROOT::RDF::RNode GetIndices(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &inputmaskname);
ROOT::RDF::RNode OrderByPt(ROOT::RDF::RNode df, const std::string &outputname,
                           const std::string &object_pt,
                           const std::string &object_mask);
ROOT::RDF::RNode MassCorrectionWithPt(ROOT::RDF::RNode df,
                                      const std::string &outputname,
                                      const std::string &raw_mass,
                                      const std::string &raw_pt,
                                      const std::string &corrected_pt);
} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */
