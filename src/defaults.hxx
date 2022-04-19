#ifndef GUARDDEFAULTS_H
#define GUARDDEFAULTS_H

#include <Math/Vector4D.h>
#include <Rtypes.h> // UChar_t

const int default_int = -10;
const int default_pdgid = -999;
const float default_float = -10.0;
const UChar_t default_uchar = -10;
auto default_lorentzvector = ROOT::Math::PtEtaPhiMVector(
    default_float, default_float, default_float, default_float);

#endif /* GUARDDEFAULTS_H */
