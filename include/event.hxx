#ifndef GUARD_EVENT_H
#define GUARD_EVENT_H

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <TRandom3.h>
#include <type_traits>

namespace event {

/**
 * @brief This function combines multiple boolean flags into a single boolean
 * value based on the selected mode ("any_of", "all_of", or "none_of"). The mode
 * determines how the flags are evaluated:
 * - `"any_of"`: Returns `true` if at least one of the flags is `true`
 * - `"all_of"`: Returns `true` if all flags are `true`
 * - `"none_of"`: Returns `true` if none of the flags are `true`
 *
 * @tparam Args variadic template parameter pack representing the flag columns
 * plus mode
 * @param df input dataframe
 * @param outputname name of the output column containing the combined flag
 * @param args parameter pack of column names that contain the considered flags
 * of type `bool`, with the last argument being the mode (`"any_of"`,
 * `"all_of"`, or
 * `"none_of"`)
 *
 * @return a dataframe with a new column
 *
 * @note The mode (`"any_of"`, `"all_of"`, or `"none_of"`) is extracted as the
 * last argument in the `args` parameter pack, and the rest of the arguments are
 * treated as individual flag columns.
 */
template <typename... Args>
inline auto CombineFlags(ROOT::RDF::RNode df, const std::string &outputname,
                         Args... args) {
    auto argTuple = std::make_tuple(args...);
    auto mode = utility::extractLastArgument(argTuple);

    std::vector<std::string> FlagList{args...};
    FlagList.pop_back();
    const auto nFlags = sizeof...(Args) - 1;

    using namespace ROOT::VecOps;
    return df.Define(
        outputname,
        utility::PassAsVec<nFlags, bool>([mode](const ROOT::RVec<bool> &flags) {
            if (mode == std::string("any_of")) {
                return Any(flags);
            } else if (mode == std::string("all_of")) {
                return All(flags);
            } else if (mode == std::string("none_of")) {
                return !Any(flags);
            } else {
                Logger::get("event::CombineFlags")
                    ->error("Mode {} is not defined!", mode);
                throw std::runtime_error("Mode is not defined!");
            }
        }),
        FlagList);
}

namespace quantity {

/**
 * @brief This function creates a new column in the dataframe with the specified
 * `outputname`, copying the values from an existing `quantity` column. The
 * original column remains unchanged.
 *
 * @tparam T type of the input quantity values
 * @param df input dataframe
 * @param outputname name of the new column
 * @param quantity name of the existing column to copy values from
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Rename(ROOT::RDF::RNode df,
                               const std::string &outputname,
                               const std::string &quantity) {
    return df.Define(outputname, [](const T &q) { return q; }, {quantity});
}

/**
 * @brief This function adds a new column to the dataframe, assigning it a
 * constant value for all entries.
 *
 * @tparam T type of the value to be assigned
 * @param df input dataframe
 * @param outputname name of the new column
 * @param value constant value to be assigned to the new column
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Define(ROOT::RDF::RNode df,
                               const std::string &outputname, T const &value) {
    return df.Define(outputname, [value]() { return value; }, {});
}

/**
 * @brief This function defines a new column in the dataframe, where each
 * element is a randomly generated number. The random values are generated using
 * `TRandom3`, seeded with a user-specified value and uniformly distributed in
 * the range [0,1]. The number of generated values matches the size of the input
 * column vector.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the generated random
 * vector
 * @param quantity name of the input column whose size determines the length of
 * the random vector
 * @param seed seed value for the random number generator, if not set the answer 
 * to everything is used as default `42`
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode
GenerateRandomVector(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &quantity, const int seed = 42) {
    return df.Define(
        outputname,
        [rndm = TRandom3(seed)](const ROOT::RVec<T> &quantity) mutable {
            ROOT::RVec<float> rndm_array(quantity.size());
            rndm.RndmArray(quantity.size(), rndm_array.data());
            return rndm_array;
        },
        {quantity});
}

ROOT::RDF::RNode
GenerateSeed(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &lumi,
    const std::string &run,
    const std::string &event,
    const UInt_t &master_seed = 42
);

/**
 * @brief This function creates a new column in the dataframe by applying
 * element-wise negation to an existing `quantity` column.
 *
 * @tparam T type of the input quantity values
 * @param df input dataframe
 * @param outputname name of the new column
 * @param quantity name of the existing column to be negated
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Negate(ROOT::RDF::RNode df,
                               const std::string &outputname,
                               const std::string &quantity) {
    auto negative = [](const T &quantity) { return -quantity; };
    return df.Define(outputname, negative, {quantity});
}

/**
 * @brief This function extracts values from the given quantity at the indices
 * specified in a collection index. The order of the output values reflects the
 * order of the indices. The function uses `ROOT::VecOps::Take` internally,
 * leading to the following behavior:
 * 
 * ```C++
 * auto values = ROOT::RVec<float>({0.1, 0.2, 0.3, 0.4});
 * auto index = ROOT::RVec<int>({2, 3, 1});
 * auto result = ROOT::VecOps::Take(values, index);
 * result
 * // (ROOT::VecOps::RVec<float>) {0.3, 0.4, 0.2}
 * ```
 * 
 * The column `index_vector` must contain the indices for which values
 * should be extracted, and the `quantity` column must contain the values
 * of the quantity.
 * 
 * Note that `T` is the type of the values stored in the `RVec` containers in
 * the `quantity` column, e.g., if the column has type `RVec<float>`, you
 * must use `T = float`.
 *
 * @tparam T underlying type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the extracted value
 * @param quantity name of the column from which the value is retrieved
 * @param index_vector index list for values to be extracted
 *
 * @return a dataframe with the new column
 *
 * @note If the index is out of range, a default value of type `T` is returned.
 */
template <typename T>
inline ROOT::RDF::RNode Take(
    ROOT::RDF::RNode df,
    const std::string &outputname,
    const std::string &quantity,
    const std::string &index_vector
) {
    auto take = [] (
        const ROOT::RVec<T> &quantity,
        const ROOT::RVec<int> &index_vector
    ) {
        Logger::get("event::quantity::Take")
            ->debug("Taking quantity {} at indices {}", quantity, index_vector);
        auto result = ROOT::VecOps::Take(quantity, index_vector);
        Logger::get("event::quantity::Take")
            ->debug("Result {}", result);

        return result;
    };

    return df.Define(
        outputname,
        take,
        {quantity, index_vector}
    );
}

/**
 * @brief This function extracts a value from the given column at a specified
 * index. If the index is out of range, a default value of type `T` is returned.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the extracted value
 * @param quantity name of the column from which the value is retrieved
 * @param index fixed index position used to extract the value
 *
 * @return a dataframe with the new column
 *
 * @note If the index is out of range, a default value of type `T` is returned.
 */
template <typename T>
inline ROOT::RDF::RNode Get(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &quantity, const int &index) {
    return df.Define(outputname,
                     [index](const ROOT::RVec<T> &quantity) {
                        T result = default_value<T>();

                        try {
                            result = quantity.at(index);
                        } catch (const std::out_of_range &e) {
                            Logger::get("event::quantity::Get")
                                ->debug(
                                    "Index not found, returning dummy value!");
                        }
                        // the static_cast is used because some types of nanoAOD 
                        // branches changed from Int_t to UChar_t for nanoAOD 
                        // versions > 9
                        if constexpr (std::is_same<T, UChar_t>::value || std::is_same<T, Short_t>::value) {
                            int cast_result = static_cast<int>(result);
                            Logger::get("event::quantity::Get")
                                ->debug("Returning UChar_t/Short_t quantity as int: {}",
                                        cast_result);
                            return cast_result;
                        } else {
                            return result;
                        }
                     },
                     {quantity});
}

/**
 * @brief This function extracts a value from the given column based on an index
 * stored in another column. If the index is out of range, a default value
 * is returned.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the extracted value
 * @param quantity name of the column from which the value is retrieved
 * @param index_vector name of the column containing index values
 * @param position position within the index vector used to retrieve the index
 *
 * @return a dataframe with the new column
 *
 * @note If the index is out of range, a default value of type T is returned.
 */
template <typename T>
inline ROOT::RDF::RNode Get(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &quantity,
                            const std::string &index_vector,
                            const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<T> &quantity,
                                const ROOT::RVec<int> &indices) {
                        T result = default_value<T>();

                        try {
                            const int index = indices.at(position);
                            result = quantity.at(index);
                        } catch (const std::out_of_range &e) {
                            Logger::get("event::quantity::Get")
                                ->debug(
                                    "Index not found, returning dummy value!");
                        }
                        // the static_cast is used because some types of nanoAOD 
                        // branches changed from Int_t to UChar_t for nanoAOD 
                        // versions > 9
                        if constexpr (std::is_same<T, UChar_t>::value || std::is_same<T, Short_t>::value) {
                            int cast_result = static_cast<int>(result);
                            Logger::get("event::quantity::Get")
                                ->debug("Returning UChar_t/Short_t quantity as int: {}",
                                        cast_result);
                            return cast_result;
                        } else {
                            return result;
                        }
                    },
                    {quantity, index_vector});
}

/**
 * @brief Get the value of a generator-level quantity associated with a
 * reconstructed object.
 *
 * The function operates as follows:
 *
 * - The reconstructed object under consideration is chosen by providing its
 *   `"position"` within an `"index_list"`.
 *
 * - The `"genobj_index"` column of the reconstructed collection must contain
 *   the index of the matched object in the associated generator-level
 *   collection.
 *
 * - The aforementioned `"genobj_index"` is used to access the generator-level
 *   object. From this object, the value of the `"gen_quantity"` column is
 *   taken.
 *
 * If the generator-level object cannot be accessed, the function returns a
 * default value.
 *
 * _Example:_ Let the column `"good_jet_indices"` contain the indices of
 * selected AK4 jets. For the `Jet` collection, the column `"Jet_genJetIdx"`
 * contains the index of the matched generator-level jet in the `"GenJet"`
 * collection. To define the generator-level p<sub>T</sub> of the leading
 * reconstructed AK4 jet, one needs to call:
 *
 * ```c++
 * event::quantity::GetFromGenObject(
 *     df,
 *     "jet_gen_pt_1",
 *     "Jet_genJetIdx",
 *     "GenJet_pt",
 *     "good_jet_indices",
 *     0
 *  )
 * ```
 *
 * @param df the input dataframe
 * @param output_name the name of the produced column
 * @param genobj_idx column that holds the index of the associated generator-level object
 * @param gen_quantity name of the quantity of the generator-level object
 * @param index_vector name of the column containing index values
 * @param position position within the index vector used to retrieve the index 
 *
 * @return a dataframe with the new column
 */
template <typename T>
ROOT::RDF::RNode GetFromGenObject(
    ROOT::RDF::RNode df,
    const std::string &output_name,
    const std::string &genobj_idx,
    const std::string &gen_quantity,
    const std::string &index_vector,
    const int &position
) {

    auto get_gen_quantity = [genobj_idx, gen_quantity, index_vector, position] (
        const ROOT::RVec<Short_t> &genobj_idx_val,
        const ROOT::RVec<T> &gen_quantity_val,
        const ROOT::RVec<int> &index_vector_val,
    ) {
        // Log column names and respective values in debug mode
        Logger::get("event::quantity::GetFromGenObject")->debug(
            "Trying to access value of matched generator-level object"
        );
        Logger::get("event::quantity::GetFromGenObject")->debug(
            "    genobj_idx '{}': {}"
            genobj_idx,
            genobj_idx_val
        );
        Logger::get("event::quantity::GetFromGenObject")->debug(
            "    gen_quantity '{}': {}"
            gen_quantity,
            gen_quantity_val
        );
        Logger::get("event::quantity::GetFromGenObject")->debug(
            "    index_vector '{}': {}"
            index_vector,
            index_vector_val
        );

        // Define the result with the default value, which is returned when
        // accessing the entries in the vectors fails
        T result = default_value<T>();

        // Get the index to access in the fatjet list
        if (position >= 0 && position < index_vector_val.size()) {
            auto index = index_vector_val.at(position);
                Logger::get("event::quantity::GetFromGenObject")->debug(
                    "    Reco object index {} at position {}",
                    index,
                    position
                );
            if (index >= 0 && index < genjet_idx_val.size()) {
                auto gen_index = genjet_idx_val.at(index);
                Logger::get("event::quantity::GetFromGenObject")->debug(
                    "    Gen object index {}",
                    gen_index
                );
                result = gen_quantity_val.at(gen_index);
                Logger::get("event::quantity::GetFromGenObject")->debug(
                    "    Retrieved quantity value {}",
                    result
                );
            }
        } else {
            Logger::get("event::quantity::GetFromGenObject")->debug(
                "Index not found, returning dummy value!"
            );
        }

        return result;
    };

    return df.Define(
        output_name,
        get_gen_quantity,
        {genjet_idx, gen_quantity, index_vector}
    );
}

/**
 * @brief This function computes the sum of the elements in the `quantity`
 * column for each event. If no elements are selected, a default value (provided
 * by `zero`) is used as the sum for that event.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the summed values
 * @param quantity name of the column containing the vector of values to be
 * summed
 * @param zero default value to use in `ROOT::VecOps::Sum` (default is `T(0)`)
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Sum(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &quantity, const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity) {
        Logger::get("event::quantity::Sum")->debug("sum values {}", quantity);
        T sum = ROOT::VecOps::Sum(quantity, zero);
        Logger::get("event::quantity::Sum")->debug("sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity});
}

/**
 * @brief This function computes the sum of the elements in the `quantity`
 * column, selected by the indices from the `indices` column. The sum is
 * computed per event, and a default value (provided by `zero`) is used if no
 * elements are selected.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the summed values
 * @param quantity name of the column containing the vector of values to be
 * summed
 * @param index_vector name of the column containing the indices used to select
 * values from `quantity`
 * @param zero default value to use in `ROOT::VecOps::Sum` (default is `T(0)`)
 *
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Sum(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &quantity,
                            const std::string &index_vector, const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity,
                                const ROOT::RVec<int> &indices) {
        Logger::get("event::quantity::Sum")
            ->debug("sum values {} at indices {}", quantity, indices);
        T sum = ROOT::VecOps::Sum(ROOT::VecOps::Take(quantity, indices), zero);
        Logger::get("event::quantity::Sum")->debug("sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity, index_vector});
}

/**
 * @brief This function calculates the scalar sum of an arbitrary set of quantities
 * of type `float`.
 *
 * @tparam Quantities variadic template parameter pack representing the quantity columns
 * @param df input dataframe
 * @param outputname name of the output column containing the scalar sum
 * @param quantities parameter pack of column names that contain the considered quantities
 *
 * @return a dataframe with a new column
 */
template <typename... Quantities>
inline ROOT::RDF::RNode 
ScalarSum(ROOT::RDF::RNode df, const std::string &outputname,
          Quantities... quantities) {
    auto argTuple = std::make_tuple(quantities...);
    std::vector<std::string> QuantityList{quantities...};
    const auto nQuantities = sizeof...(Quantities);

    using namespace ROOT::VecOps;
    return df.Define(
        outputname,
        utility::PassAsVec<nQuantities, float>([](const ROOT::RVec<float> &quantities) {
            for (const auto &quantity : quantities) {
                if (quantity < 0.0) {
                    Logger::get("event::quantity::ScalarSum")
                        ->debug("Negative quantity found, returning default value!");
                    return default_float;
                }
            }
            const auto sum = Sum(quantities, float(0.0));
            return sum;
        }),
        QuantityList);
}

/**
 * @brief This function recursively unrolls a vector (`std::vector<T>`) from the
 * `quantity` column into individual columns in the dataframe. Each element of
 * the vector is stored in a separate column with names provided in the
 * `outputnames` vector. The function works recursively to define a new column
 * for each element in the vector.
 *
 * @warning The length of the quantity vector has to be the same for each event.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputnames a vector of names for the new columns where the individual
 * elements of the vector will be stored
 * @param quantity name of the column containing the vector of values to unroll
 * @param index index of the current element to unroll (defaults to 0).
 *
 * @return a dataframe with the new columns containing each individual element
 * of the vector from the `quantity` column
 *
 * @note The function is recursive and will create one column for each element
 * of the vector in `quantity`. If `outputnames` has fewer entries than the
 * number of elements in the vector, the function will stop at the end of
 * `outputnames`. The `index` should not be set outside this function.
 */
template <typename T>
inline ROOT::RDF::RNode
Unroll(ROOT::RDF::RNode df, const std::vector<std::string> &outputnames,
       const std::string &quantity, const size_t &index = 0) {
    if (index >= outputnames.size()) {
        return df;
    }
    auto df1 = df.Define(
        outputnames.at(index),
        [index](const std::vector<T> &quantities) { return quantities.at(index); },
        {quantity});
    return Unroll<T>(df1, outputnames, quantity, index + 1);
}

} // end namespace quantity

namespace filter {

/**
 * @brief This function applies a filter to the input dataframe based on a
 * boolean flag column. It returns only the rows where the flag value is `true`.
 *
 * Use case examples are the noise filters recommended by the CMS JetMET
 * group
 * (https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2).
 *
 * @param df input dataframe
 * @param filtername name of the filter to be applied (used in the dataframe
 * report)
 * @param flagname name of the boolean flag column to use for filtering
 *
 * @return a filtered dataframe
 */
inline ROOT::RDF::RNode Flag(ROOT::RDF::RNode df, const std::string &filtername,
                             const std::string &flagname) {
    return df.Filter([](const bool flag) { return flag; }, {flagname},
                     filtername);
}

/**
 * @brief This function filters the rows of the input dataframe by evaluating
 * multiple boolean flags according to a specified mode. The filtering mode
 * can be "any_of", "all_of", or "none_of":
 * - `"any_of"`: Keeps the rows where at least one flag is `true`
 * - `"all_of"`: Keeps the rows where all flags are `true`
 * - `"none_of"`: Keeps the rows where none of the flags are `true`
 *
 * @tparam Args variadic template parameter pack representing the flag columns
 * plus mode
 * @param df input dataframe
 * @param filtername name of the filter to be applied (used in the dataframe
 * report)
 * @param args parameter pack of column names that contain the considered flags
 * of type `bool`, with the last argument being the mode (`"any_of"`,
 * `"all_of"`, or
 * `"none_of"`)
 *
 * @return a filtered dataframe
 *
 * @note The last argument must be the mode, while the preceding arguments are
 * the boolean flag columns to be evaluated.
 */
template <typename... Args>
inline auto Flags(ROOT::RDF::RNode df, const std::string &filtername,
                  Args... args) {
    auto argTuple = std::make_tuple(args...);
    auto mode = utility::extractLastArgument(argTuple);

    std::vector<std::string> FlagList{args...};
    FlagList.pop_back();
    const auto nFlags = sizeof...(Args) - 1;

    using namespace ROOT::VecOps;
    return df.Filter(
        utility::PassAsVec<nFlags, bool>([mode](const ROOT::RVec<bool> &flags) {
            if (mode == std::string("any_of")) {
                return Any(flags);
            } else if (mode == std::string("all_of")) {
                return All(flags);
            } else if (mode == std::string("none_of")) {
                return !Any(flags);
            } else {
                Logger::get("event::filter::Flags")
                    ->error("Mode {} is not defined!", mode);
                throw std::runtime_error("Mode is not defined!");
            }
        }),
        FlagList, filtername);
}

/**
 * @brief This function filters the rows of the input dataframe by checking if a
 * specified `quantity` exists in the provided `selection` vector. Rows where
 * the quantity is found in the selection vector are kept, while others are
 * removed.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param filtername name of the filter to be applied (used in the dataframe
 * report)
 * @param quantity name of the quantity column in the dataframe of type `T`
 * @param selection a vector containing the selection of values of type `T` to
 * filter the quantity against
 *
 * @return a filtered dataframe
 */
template <typename T>
inline ROOT::RDF::RNode
Quantity(ROOT::RDF::RNode df, const std::string &filtername,
         const std::string &quantity, const std::vector<T> &selection) {
    return df.Filter(
        [selection](const T quantity) {
            return std::find(selection.begin(), selection.end(), quantity) !=
                   selection.end();
        },
        {quantity}, filtername);
}

ROOT::RDF::RNode
GoldenJSON(ROOT::RDF::RNode df,
           correctionManager::CorrectionManager &correctionManager,
           const std::string &filtername, const std::string &run,
           const std::string &luminosity, const std::string &json_path);
} // end namespace filter
} // end namespace event

#endif /* GUARD_EVENT_H */
