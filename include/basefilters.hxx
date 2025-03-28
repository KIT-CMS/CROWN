#ifndef GUARDBASEFILTERS_H
#define GUARDBASEFILTERS_H

#include "../include/utility/Logger.hxx"
#include "ROOT/RVec.hxx"

namespace basefilters {

/**
 * @brief Function to apply a maximal filter requirement to a quantity
 * of type `T`. The type has to be defined when calling this function
 * `basefilters::Max<T>`. It returns true if the value is
 * smaller than the given `cut` value.
 *
 * @param cut cut value of the filter
 *
 * @return a lambda function to be used in RDF Define
 */
template <typename T> 
inline auto Max(const T &cut) {
    return [cut](const ROOT::RVec<T> &values) {
        ROOT::RVec<int> mask = values < cut;
        return mask;
    };
}

/**
 * @brief Function to apply a maximal filter requirement to a quantity
 * of type `T`. The type has to be defined when calling this function
 * `basefilters::AbsMax<T>`. It returns true if the absolute
 * value is smaller than the given `cut` value.
 *
 * @param cut cut value of the filter
 *
 * @return a lambda function to be used in RDF Define
 */
template <typename T> 
inline auto AbsMax(const T &cut) {
    return [cut](const ROOT::RVec<T> &values) {
        ROOT::RVec<int> mask = abs(values) < cut;
        return mask;
    };
}

/**
 * @brief Function to apply a minimal filter requirement to a quantity
 * of type `T`. The type has to be defined when calling this function
 * `basefilters::Min<T>`. It returns true if the value is
 * larger than the given `cut` value or equal to it.
 *
 * @param cut cut value of the filter
 *
 * @return a lambda function to be used in RDF Define
 */
template <typename T> 
inline auto Min(const T &cut) {
    return [cut](const ROOT::RVec<T> &values) {
        ROOT::RVec<int> mask = values >= cut;
        return mask;
    };
}

/**
 * @brief Function to apply a minimal filter requirement to a quantity
 * of type `T`. The type has to be defined when calling this function
 * `basefilters::AbsMin<T>`. It returns true if the absolute
 * value is larger than the given `cut` value or equal to it.
 *
 * @param cut cut value of the filter
 *
 * @return a lambda function to be used in RDF Define
 */
template <typename T> 
inline auto AbsMin(const T &cut) {
    return [cut](const ROOT::RVec<T> &values) {
        ROOT::RVec<int> mask = abs(values) >= cut;
        return mask;
    };
}

/**
 * @brief Function to apply an exact filter requirement to a quantity
 * of type `T`. The type has to be defined when calling this function
 * `basefilters::Equal<T>`. It returns true if the value is
 * equal to the given `cut` value.
 *
 * @param cut cut value of the filter
 *
 * @return a lambda function to be used in RDF Define
 */
template <typename T> 
inline auto Equal(const T &cut) {
    return [cut](const ROOT::RVec<T> &values) {
        ROOT::RVec<int> mask = values == cut;
        return mask;
    };
}

/**
 * @brief Function to apply an exact filter requirement to a quantity
 * of type `T`. The type has to be defined when calling this function
 * `basefilters::AbsEqual<T>`. It returns true if the absolute
 * value is equal to the given `cut` value.
 *
 * @param cut cut value of the filter
 *
 * @return a lambda function to be used in RDF Define
 */
template <typename T> 
inline auto AbsEqual(const T &cut) {
    return [cut](const ROOT::RVec<T> &values) {
        ROOT::RVec<int> mask = abs(values) == cut;
        return mask;
    };
}

/**
 * @brief Function to apply a filter requirement to a bitmask quantity
 * of type `Int_t` or `UChar_t`. The type has to be defined when calling this
 * function `basefilters::Bitmask<T>`. It returns true if the bitmask
 * value is larger than the given `cut` bitmask value.
 *
 * @param cut bitmask cut value of the filter
 *
 * @return a lambda function to be used in RDF Define
 */
template <typename T> 
inline auto Bitmask(const int &cut) {
    return [cut](const ROOT::RVec<T> &values) {
        ROOT::RVec<int> mask;
        if (std::is_same<T, int>::value) {
            mask = values >= cut;
        } else if (std::is_same<T, UChar_t>::value) {
            for (auto const value : values) {
                if (cut > 0)
                    mask.push_back(std::min(1, int(value & 1 << (cut - 1))));
                else
                    mask.push_back(int(1));
            }
        } else {
            Logger::get("FilterBitmask")
                ->debug("Type of {} not implemented.", values);
            mask = ROOT::RVec<int>(values.size(), int(0));
        }
        return mask;
    };
}
} // namespace basefilters

#endif /* GUARDBASEFILTERS_H */
