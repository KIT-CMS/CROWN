#ifndef GUARDDEFAULTS_H
#define GUARDDEFAULTS_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include <Math/Vector4D.h>
#include <type_traits>

const int default_int = -10;
const int default_pdgid = -999;
const float default_float = -10.0;
const UChar_t default_uchar = -10;
const bool default_bool = false;
const auto default_lorentzvector = ROOT::Math::PtEtaPhiMVector(
    default_float, default_float, default_float, default_float);

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
    // is there a better way to handle this?
    return static_cast<T>(default_int);
};

#endif /* GUARDDEFAULTS_H */
