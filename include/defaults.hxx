#ifndef GUARDDEFAULTS_H
#define GUARDDEFAULTS_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include <Math/Vector4D.h>

const int default_int = -10;
const int default_pdgid = -999;
const float default_float = -10.0;
const UChar_t default_uchar = -10;
const auto default_lorentzvector = ROOT::Math::PtEtaPhiMVector(
    default_float, default_float, default_float, default_float);

#endif /* GUARDDEFAULTS_H */
