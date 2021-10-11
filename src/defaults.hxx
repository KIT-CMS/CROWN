#ifndef GUARDDEFAULTS_H
#define GUARDDEFAULTS_H

int default_int = -10;
int default_pdgid = -999;
float default_float = -10.0;
UChar_t default_uchar = -10;
auto default_lorentzvector = ROOT::Math::PtEtaPhiMVector(
    default_float, default_float, default_float, default_float);

#endif /* GUARDDEFAULTS_H */
