#ifndef GUARD_FATJETS_H
#define GUARD_FATJETS_H

#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <algorithm>

namespace physicsobject {
namespace fatjet {
namespace quantity {

/**
 * @brief This function calculates a discriminator score from two ParticleNet 
 * tagger outputs (a signal X and the QCD background). The signal can e.g. be 
 * \f$X\rightarrow bb\f$ or \f$X\rightarrow cc\f$ etc. 
 * 
 * The score is computed using the formula: 
 * \f[
 * Score = \frac{P(X)}{P(X) + P(QCD)}
 * \f]
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the XvsQCD score
 * @param pNet_X_decay name of the column containing the ParticleNet score for 
 * the signal process (\f$X\rightarrow ...\f$)
 * @param pNet_QCD name of the column containing the ParticleNet score for the 
 * QCD background
 * @param fatjet_collection name of the column containing a collection (vector) 
 * of good fatjet indices
 * @param position position of the fatjet in the collection (vector) that should 
 * be used for the score calculation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
ParticleNet_XvsQCD(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pNet_X_decay, const std::string &pNet_QCD,
                     const std::string &fatjet_collection, const int &position) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<float> &X_tagger,
                                const ROOT::RVec<float> &QCD_tagger,
                                const ROOT::RVec<int> &fatjets) {
                        float X_decay = default_float;
                        float QCD = default_float;
                        float X_vs_QCD = default_float;
                        if (position >= 0) {
                            const int index = fatjets.at(position);
                            if (index >= 0) {
                                X_decay = X_tagger.at(index);
                                QCD = QCD_tagger.at(index);
                                X_vs_QCD = X_decay / (X_decay + QCD);
                            }
                        }
                        return X_vs_QCD;
                     },
                     {pNet_X_decay, pNet_QCD, fatjet_collection});
}

/**
 * @brief This function calculates the ratio of two N-subjettiness variables 
 * of one fatjet, typically used for discriminating jet substructure:
 * \f[
 * \tau_{ratio} = \frac{\tau_N}{\tau_{N-1}}
 * \f]
 * For example, \f$\tau_{21} = \tau_2 / \tau_1\f$ for 2-prong vs. 1-prong decay 
 * identification. The resulting ratio is always between 0 and 1. If it is close 
 * to 0, the jet is more likely to be a N-prong jet, while a value close to 1 
 * indicates a (N-1)-prong jet.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the N-subjettiness ratio
 * @param tau_N name of the column containing the N-subjettiness variable for 
 * the numerator with number `N`
 * @param tau_Nm1 name of the column containing the N-subjettiness variable for 
 * the denominator with number `N-1`
 * @param fatjet_collection name of the column containing a collection (vector) 
 * of good fatjet indices
 * @param position position of the fatjet in the collection (vector) that should 
 * be used for the ratio calculation
 *
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
NsubjettinessRatio(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tau_N, const std::string &tau_Nm1,
                    const std::string &fatjet_collection, const int &position) {
    return df.Define(outputname,
                    [position](const ROOT::RVec<float> &tau_N,
                               const ROOT::RVec<float> &tau_Nm1,
                               const ROOT::RVec<int> &fatjets) {
                        float nsubjet_N = default_float;
                        float nsubjet_Nm1 = default_float;
                        float ratio = default_float;
                        if (position >= 0) { 
                            const int index = fatjets.at(position);
                            if (index >= 0) {
                                nsubjet_N = tau_N.at(index);
                                nsubjet_Nm1 = tau_Nm1.at(index);
                                ratio = nsubjet_N / nsubjet_Nm1;
                            }
                        }
                        return ratio;
                    },
                    {tau_N, tau_Nm1, fatjet_collection});
}
} // end namespace quantity
} // end namespace fatjet
} // end namespace physicsobject
#endif /* GUARD_FATJETS_H */