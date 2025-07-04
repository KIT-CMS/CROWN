#ifndef GUARD_DEFAULTS_H
#define GUARD_DEFAULTS_H

#include "ROOT/RDataFrame.hxx"
#include <Math/Vector4D.h>
#include <type_traits>

const int default_int = -10;
const float default_float = -10.0;
// casting this default UChar_t to an int results in a value of 246 
// which should still be out of range for the usual use cases 
const UChar_t default_uchar = -10; 
const bool default_bool = false;
const auto default_lorentzvector = ROOT::Math::PtEtaPhiMVector(0., 0., 0., 0.);

template <typename T> const T default_value() {
    if (std::is_same<T, int>::value) {
        return default_int;
    } else if (std::is_same<T, const int>::value) {
        return default_int;
    } else if (std::is_same<T, float>::value) {
        return default_float;
    } else if (std::is_same<T, const float>::value) {
        return default_float;
    } else if (std::is_same<T, bool>::value) {
        return default_bool;
    } else if (std::is_same<T, const bool>::value) {
        return default_bool;
    } else if (std::is_same<T, UChar_t>::value) {
        return default_uchar;
    } else if (std::is_same<T, const UChar_t>::value) {
        return default_uchar;
    }
    return static_cast<T>(default_int);
};

#endif /* GUARD_DEFAULTS_H */
