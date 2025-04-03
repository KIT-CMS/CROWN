#ifndef GUARD_EVENT_H
#define GUARD_EVENT_H

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <TRandom3.h>

namespace event {

/**
 * @brief This function combines multiple boolean flags into a single boolean value 
 * based on the selected mode ("any", "all", or "none"). The mode determines 
 * how the flags are evaluated:
 * - `"any"`: Returns `true` if at least one of the flags is `true`
 * - `"all"`: Returns `true` if all flags are `true`
 * - `"none"`: Returns `true` if none of the flags are `true`
 *
 * @tparam Args variadic template parameter pack representing the flag columns plus mode
 * @param df input dataframe
 * @param outputname name of the output column containing the combined flag
 * @param args parameter pack of column names that contain the considered flags of 
 * type `bool`, with the last argument being the mode (`"any"`, `"all"`, or `"none"`)
 * 
 * @return a dataframe with a new column
 *
 * @note The mode (`"any"`, `"all"`, or `"none"`) is extracted as the last argument 
 * in the `args` parameter pack, and the rest of the arguments are treated as 
 * individual flag columns.
 */
template <typename... Args>
inline auto CombineFlags(ROOT::RDF::RNode df,
                         const std::string &outputname,
                         Args... args) {
    auto argTuple = std::make_tuple(args...);
    auto mode = utility::extractLastArgument(argTuple);

    std::vector<std::string> FlagList{args...};
    FlagList.pop_back();  
    const auto nFlags = sizeof...(Args) - 1;

    using namespace ROOT::VecOps;
    return df.Define(
        outputname,
        utility::PassAsVec<nFlags, bool>(
            [mode](const ROOT::RVec<bool> &flags) {
                if (mode == std::string("any")) {
                    return Any(flags);
                }
                else if (mode == std::string("all")) {
                    return All(flags);
                }
                else if (mode == std::string("none")) {
                    return !Any(flags);
                }
                else {
                    Logger::get("event::CombineFlags")->error("Mode {} is not defined!", mode);
                    throw std::runtime_error("Mode is not defined!");
                }
            }
        ),
        FlagList);
}

namespace quantity {

/**
 * @brief This function creates a new column in the dataframe with the specified 
 * `outputname`, copying the values from an existing `quantity` column. The original 
 * column remains unchanged.
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
                               const std::string &outputname,
                               T const &value) {
    return df.Define(outputname, [value]() { return value; }, {});
}

/**
 * @brief This function defines a new column in the dataframe, where each element is 
 * a randomly generated number. The random values are generated using `TRandom3`, 
 * seeded with a user-specified value and uniformly distributed in the range [0,1]. 
 * The number of generated values matches the size of the input column vector.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the generated random vector
 * @param quantity name of the input column whose size determines the length of the 
 * random vector
 * @param seed seed value for the random number generator
 * 
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode GenerateRandomVector(ROOT::RDF::RNode df,
                                      const std::string &outputname,
                                      const std::string &quantity, 
                                      const int seed) {
    return df.Define(outputname, 
        [rndm = TRandom3(seed)](const ROOT::RVec<T> &quantity) mutable {        
            ROOT::RVec<float> rndm_array(quantity.size());
            rndm.RndmArray(quantity.size(), rndm_array.data());
            return rndm_array;
        }, 
        {quantity});
}

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
 * @brief This function extracts a value from the given column at a specified index.
 * If the index is out of range, a default value of type `T` is returned.
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
inline ROOT::RDF::RNode Get(ROOT::RDF::RNode df, 
                            const std::string &outputname,
                            const std::string &quantity, 
                            const int &index) {
    return df.Define(outputname,
        [index](const ROOT::RVec<T> &quantity) {
            return quantity.at(index, default_value<T>());
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
inline ROOT::RDF::RNode Get(ROOT::RDF::RNode df, 
                            const std::string &outputname,
                            const std::string &quantity, 
                            const std::string &index_vector,
                            const int &position) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<T> &quantity, const ROOT::RVec<int> &indices) {
            T result = default_value<T>();

            try {
                const int index = indices.at(position);
                result = quantity.at(index, default_value<T>());
            } catch (const std::out_of_range &e) {
                Logger::get("event::quantity::Get")
                    ->debug("Index not found, returning dummy value!");
            }
            return result;
        },
        {quantity, index_vector});
}

/**
 * @brief This function computes the sum of the elements in the `quantity` column, 
 * selected by the indices from the `indices` column. The sum is computed 
 * per event, and a default value (provided by `zero`) is used if no elements 
 * are selected.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the summed values
 * @param quantity name of the column containing the vector of values to be summed
 * @param indices name of the column containing the indices used to select values 
 * from `quantity`
 * @param zero default value to use in `ROOT::VecOps::Sum` (default is `T(0)`)
 * 
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Sum(ROOT::RDF::RNode df, 
                                  const std::string &outputname,
                                  const std::string &quantity,
                                  const std::string &indices, 
                                  const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity,
                                const ROOT::RVec<int> &indices) {
        Logger::get("event::quantity::SumVector")
            ->debug("sum values {} at indices {}", quantity, indices);
        T sum = ROOT::VecOps::Sum(
            ROOT::VecOps::Take(quantity, indices), zero);
        Logger::get("event::quantity::SumVector")->debug("sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity, indices});
}

/**
 * @brief This function computes the sum of the elements in the `quantity` column 
 * for each event. If no elements are selected, a default value (provided 
 * by `zero`) is used as the sum for that event.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputname name of the new column containing the summed values
 * @param quantity name of the column containing the vector of values to be summed
 * @param zero default value to use in `ROOT::VecOps::Sum` (default is `T(0)`)
 * 
 * @return a dataframe with the new column
 */
template <typename T>
inline ROOT::RDF::RNode Sum(ROOT::RDF::RNode df, 
                                  const std::string &outputname,
                                  const std::string &quantity, 
                                  const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity) {
        Logger::get("event::quantity::Sum")->debug("sum values {}", quantity);
        T sum = ROOT::VecOps::Sum(quantity, zero);
        Logger::get("event::quantity::Sum")->debug("sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity});
}

/**
 * @brief This function recursively unrolls a vector (`std::vector<T>`) from the 
 * `quantity` column into individual columns in the dataframe. Each element of the 
 * vector is stored in a separate column with names provided in the `outputnames` vector.
 * The function works recursively to define a new column for each element in the vector.
 * 
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param outputnames a vector of names for the new columns where the individual 
 * elements of the vector will be stored
 * @param quantity name of the column containing the vector of values to unroll
 * @param idx index of the current element to unroll (defaults to 0).
 * 
 * @return a dataframe with the new columns containing each individual element of 
 * the vector from the `quantity` column
 * 
 * @note The function is recursive and will create one column for each element of 
 * the vector in `quantity`. If `outputnames` has fewer entries than the 
 * number of elements in the vector, the function will stop at the end of 
 * `outputnames`. The `idx` should not be set outside this function.
 */
template <typename T>
inline ROOT::RDF::RNode Unroll(ROOT::RDF::RNode df, 
                               const std::vector<std::string> &outputnames,
                               const std::string &quantity,
                               const size_t &idx = 0) {
    if (idx >= outputnames.size()) {
        return df;
    }
    auto df1 = df.Define(
        outputnames.at(idx),
        [idx](const std::vector<T> &quantities) { return quantities.at(idx); },
        {quantity});
    return Unroll<T>(df1, outputnames, quantity, idx + 1);
}
} // namespace quantity

namespace filter {

/**
 * @brief This function applies a filter to the input dataframe based on a boolean 
 * flag column. It returns only the rows where the flag value is `true`. 
 *
 * Use case examples are the noise filters recommended by the CMS JetMET 
 * group (https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2).
 *
 * @param df input dataframe
 * @param filtername name of the filter to be applied (used in the dataframe report)
 * @param flagname name of the boolean flag column to use for filtering
 *
 * @return a filtered dataframe
 */
inline ROOT::RDF::RNode Flag(ROOT::RDF::RNode df,
                             const std::string &filtername,
                             const std::string &flagname) {
    return df.Filter([](const bool flag) { return flag; }, {flagname},
        filtername);
}

/**
 * @brief This function filters the rows of the input dataframe by evaluating 
 * multiple boolean flags according to a specified mode. The filtering mode 
 * can be "any", "all", or "none":
 * - `"any"`: Keeps the rows where at least one flag is `true`
 * - `"all"`: Keeps the rows where all flags are `true`
 * - `"none"`: Keeps the rows where none of the flags are `true`
 *
 * @tparam Args variadic template parameter pack representing the flag columns plus mode
 * @param df input dataframe
 * @param filtername name of the filter to be applied (used in the dataframe report)
 * @param args parameter pack of column names that contain the considered flags of 
 * type `bool`, with the last argument being the mode (`"any"`, `"all"`, or `"none"`)
 *
 * @return a filtered dataframe
 *
 * @note The last argument must be the mode, while the preceding arguments are the 
 * boolean flag columns to be evaluated.
 */
template <typename... Args>
inline auto Flags(ROOT::RDF::RNode df,
                  const std::string &filtername,
                  Args... args) {
    auto argTuple = std::make_tuple(args...);
    auto mode = utility::extractLastArgument(argTuple);

    std::vector<std::string> FlagList{args...};
    FlagList.pop_back();  
    const auto nFlags = sizeof...(Args) - 1;

    using namespace ROOT::VecOps;
    return df.Filter(
        utility::PassAsVec<nFlags, bool>(
            [mode](const ROOT::RVec<bool> &flags) {
                if (mode == std::string("any")) {
                    return Any(flags);
                }
                else if (mode == std::string("all")) {
                    return All(flags);
                }
                else if (mode == std::string("none")) {
                    return !Any(flags);
                }
                else {
                    Logger::get("FilterFlags")->error("Mode {} is not defined!", mode);
                    throw std::runtime_error("Mode is not defined!");
                }
            }
        ),
        FlagList, filtername);
}

/**
 * @brief This function filters the rows of the input dataframe by checking if a 
 * specified `quantity` exists in the provided `selection` vector. Rows where the 
 * quantity is found in the selection vector are kept, while others are removed.
 *
 * @tparam T type of the input column values
 * @param df input dataframe
 * @param filtername name of the filter to be applied (used in the dataframe report)
 * @param quantity name of the quantity column in the dataframe of type `T`
 * @param selection a vector containing the selection of values of type `T` to filter 
 * the quantity against
 *
 * @return a filtered dataframe
 */
template <typename T>
inline ROOT::RDF::RNode Quantity(ROOT::RDF::RNode df, 
                                 const std::string &filtername,
                                 const std::string &quantity, 
                                 const std::vector<T> &selection) {
    return df.Filter(
        [selection](const T quantity) {
            return std::find(selection.begin(), selection.end(), quantity) !=
                   selection.end();
        },
        {quantity}, filtername);
}

ROOT::RDF::RNode GoldenJSON(ROOT::RDF::RNode df,
           correctionManager::CorrectionManager &correctionManager,
           const std::string &filtername, const std::string &run,
           const std::string &luminosity, const std::string &json_path);
} // namespace filter
} // namespace event

#endif /* GUARD_EVENT_H */