#ifndef GUARDEVENT_H
#define GUARDEVENT_H

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

namespace event {
ROOT::RDF::RNode
GoldenJSONFilter(ROOT::RDF::RNode df,
                 correctionManager::CorrectionManager &correctionManager,
                 const std::string &filtername, const std::string &run,
                 const std::string &luminosity, const std::string &json_path);

/**
 * @brief Function to apply a flag filter to the dataframe. The input flag should 
 * be a flag in the NanoAOD e.g. the noise filters recommended by the CMS JetMET 
 * group (https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2).
 * 
 * @param df input dataframe
 * @param filtername name of the filter, used in the dataframe report
 * @param flagname name of the filter flag
 *
 * @return a filtered dataframe
 */
inline ROOT::RDF::RNode FilterFlag(ROOT::RDF::RNode df,
                                   const std::string &filtername,
                                   const std::string &flagname) {
    return df.Filter([](const bool flag) { return flag; }, {flagname},
                     filtername);
}

/**
 * @brief This function filters events, where none of the input `flags` are true.
 * This can be used to filter events which do not pass an underlying requirement 
 * in any systematic variation.
 *
 * @param df input dataframe
 * @param filtername name of the filter, used in the dataframe report
 * @param flags parameter pack of column names that contain the considered
 * flags of type bool
 *
 * @return a filtered dataframe
 */
template <class... Flags>
inline ROOT::RDF::RNode FilterFlagsAny(ROOT::RDF::RNode df,
                                       const std::string &filtername,
                                       const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Filter(
        utility::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return Any(flags); }),
        FlagList, filtername);
}

/**
 * @brief This function filters events, where at least one of the input `flags` 
 * is false. This is used to filter events which do not pass an underlying 
 * requirement in all systematic variation.
 *
 * @param df input dataframe
 * @param filtername name of the filter, used in the dataframe report
 * @param flags parameter pack of column names that contain the considered
 * flags of type bool
 *
 * @return a filtered dataframe
 */
template <class... Flags>
inline ROOT::RDF::RNode FilterFlagsAll(ROOT::RDF::RNode df,
                                       const std::string &filtername,
                                       const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Filter(
        utility::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return All(flags); }),
        FlagList, filtername);
}

/**
 * @brief Function to apply a `selection` on a `quantity` of type `T` (e.g. `int`,
 * `float`, `double`, `bool`). The type has to be defined when calling this function
 * `event::FilterQuantity<T>`. If the `quantity` value is in the given `selection`
 * value vector (with values of type `T`), an event is kept. If not it is filtered
 * out.
 *
 * @param df input dataframe
 * @param filtername name of the filter, used in the dataframe report
 * @param quantity quantity for the selection
 * @param selection values of the quantity that are accepted
 *
 * @return a filtered dataframe
 */
template <typename T>
inline ROOT::RDF::RNode
FilterQuantity(ROOT::RDF::RNode df, const std::string &filtername,
               const std::string &quantity, const std::vector<T> &selection) {
    return df.Filter(
        [selection](const T quantity) {
            return std::find(selection.begin(), selection.end(), quantity) !=
                   selection.end();
        },
        {quantity}, filtername);
}

/**
 * @brief Function to write out a quantity from a vector quantity based
 * on an index taken from a index vector.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the quantity value
 * @param column name of the column containing the vector quantity
 * @param index_vector name of the column containing the index vector
 * @param position index of the position in the index vector
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode
GetQuantity(ROOT::RDF::RNode df, const std::string &outputname,
            const std::string &column, const std::string &index_vector,
            const int &position) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &idx_vec, const ROOT::RVec<T> &col) {
            T out = default_value<T>();

            try {
                const int index = idx_vec.at(position);
                out = col.at(index, default_value<T>());
            } catch (const std::out_of_range &e) {
                Logger::get("GetQuantity")
                    ->debug("Index not found, returning dummy value!");
            }
            return out;
        },
        {index_vector, column});
}

/**
 * @brief Function to write out a quantity from a vector quantity based
 * on an index.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the quantity value
 * @param column name of the column containing the vector quantity
 * @param index index in the vector quantity
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode
GetQuantity(ROOT::RDF::RNode df, const std::string &outputname,
            const std::string &column, const int &index) {
    return df.Define(outputname,
                     [index](const ROOT::RVec<T> &col) {
                         return col.at(index, default_value<T>());
                     },
                     {column});
}

/**
 * @brief Function to sum all elements of the column with name `quantity`
 * containing `ROOT::VecOps::RVec<T>` objects.
 *
 * This function is a template implementation, i.e., call SumVectorQuantity<T>
 * if the column `quantity` contains `ROOT::VecOps::RVec<T>` objects.
 *
 * Elements of the `ROOT::VecOps::RVec<T>`, which should enter the sum, can be
 * selected with index lists from the column `collection_index` as
 * `ROOT::VecOps::RVec<int>` objects per entry.
 *
 * Internally, `ROOT::VecOps::Sum` is used to calculate the sum. A custom zero
 * element, which is a second optional argument of `ROOT::VecOps::Sum`, can be
 * passed to this function by setting the parameter `zero`. Its default value is
 * `T(0)`. For instance, when dealing with `ROOT::Math::PtEtaPhiMVector`
 * objects, the `zero` parameter must be set to `ROOT::Math::PtEtaPhiMVector(0.,
 * 0., 0., 0)` in order to enable summation with this function.
 *
 * @param df input dataframe
 * @param outputname name of the output column
 * @param quantity column name of the vector quantity which is summed per event
 * @param collection_index column name for index lists of the elements to be
 * summed
 * @param zero zero element passed as the second argument to the
 * `ROOT::VecOps::Sum` function
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode
SumVectorQuantity(ROOT::RDF::RNode df, const std::string &outputname,
                  const std::string &quantity,
                  const std::string &collection_index, const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity,
                                const ROOT::RVec<int> &collection_index) {
        Logger::get("SumVectorQuantity")
            ->debug("sum values {} at indices {}", quantity, collection_index);
        T sum = ROOT::VecOps::Sum(
            ROOT::VecOps::Take(quantity, collection_index), zero);
        Logger::get("SumVectorQuantity")->debug("sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity, collection_index});
}

/**
 * @brief Function to sum all elements of the column with name `quantity`
 * containing `ROOT::VecOps::RVec<T>` objects.
 *
 * This function is a template implementation, i.e., call SumVectorQuantity<T>
 * if the column `quantity` contains `ROOT::VecOps::RVec<T>` objects.
 *
 * Internally, `ROOT::VecOps::Sum` is used to calculate the sum. A custom zero
 * element, which is a second optional argument of `ROOT::VecOps::Sum`, can be
 * passed to this function by setting the parameter `zero`. Its default value is
 * `T(0)`.
 *
 * @param df input dataframe
 * @param outputname name of the output column
 * @param quantity column name of the vector quantity which is summed per event
 * @param zero zero element passed as the second argument to the
 * `ROOT::VecOps::Sum` function
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode
SumVectorQuantity(ROOT::RDF::RNode df, const std::string &outputname,
                  const std::string &quantity, const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity) {
        Logger::get("SumVectorQuantity")->debug("sum values {}", quantity);
        T sum = ROOT::VecOps::Sum(quantity, zero);
        Logger::get("SumVectorQuantity")->debug("sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity});
}

/**
 * @brief Helper function to recursively define columns for each entry of a
 * vector quantity.
 *
 * @param df input dataframe
 * @param outputnames vector of names for the new columns
 * @param name name of the vector quantity
 * @param idx index of the current recursion loop, should not be set outside
 * this function
 *
 * @return a lambda function to be used in RDF Define and a new column for each
 * vector entry
 */
template <typename T>
inline ROOT::RDF::RNode
UnrollVectorQuantity(ROOT::RDF::RNode df, const std::string &name,
                     const std::vector<std::string> &names,
                     const size_t &idx = 0) {
    if (idx >= names.size()) {
        return df;
    }
    auto df1 = df.Define(
        names.at(idx),
        [idx](const std::vector<T> &quantities) { return quantities.at(idx); },
        {name});
    return UnrollVectorQuantity<T>(df1, name, names, idx + 1);
}

/**
 * @brief This function defines a flag that is true if at least one of the input
 * `flags` is true.
 *
 * @param df input dataframe
 * @param outputflag name of the new column
 * @param flags parameter pack of column names that contain the considered
 * flags of type bool
 *
 * @return a dataframe with the new column
 */
template <class... Flags>
inline ROOT::RDF::RNode CombineFlagsAny(ROOT::RDF::RNode df,
                                        const std::string &outputflag,
                                        const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Define(
        outputflag,
        utility::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return Any(flags); }),
        FlagList);
}

/**
 * @brief This function defines a flag that is true if all of the input
 * `flags` are true.
 *
 * @param df input dataframe
 * @param outputflag name of the new column
 * @param flags parameter pack of column names that contain the considered
 * flags of type bool
 *
 * @return a dataframe with the new column
 */
template <class... Flags>
inline ROOT::RDF::RNode CombineFlagsAll(ROOT::RDF::RNode df,
                                        const std::string &outputflag,
                                        const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Define(
        outputflag,
        utility::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return All(flags); }),
        FlagList);
}
} // namespace event

#endif /* GUARDEVENT_H */