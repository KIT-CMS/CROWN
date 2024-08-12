#ifndef GUARDFINALSTATE_H
#define GUARDFINALSTATE_H

#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include <typeinfo>


namespace finalstate{

ROOT::RDF::RNode n_taus(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &tau_pt,
                                 const std::string &tau_eta,
                                 const std::string &tau_dz,
                                 const std::string &tau_vs_jet,
                                 const std::string &tau_vs_ele,
                                 const std::string &tau_vs_mu,
                                 const std::string &tau_dm,
                                 const float &TauPtThreshold,
                                 const float &TauEtaThreshold,
                                 const float &TauDZThreshold,
                                 const float &TauVsJetThreshold,
                                 const float &TauVsMuThreshold,
                                 const float &TauVsEleThreshold,
                                 const std::vector<int> &SelectedDecayModes) {


                        // Lambda function to count taus passing the cuts
                        auto count_taus_lambda = [=](const ROOT::RVec<float> &pt,
                                                    const ROOT::RVec<float> &eta,
                                                    const ROOT::RVec<float> &dz,
                                                    const ROOT::RVec<unsigned char> &vs_jet,
                                                    const ROOT::RVec<unsigned char> &vs_muon,
                                                    const ROOT::RVec<unsigned char> &vs_ele,
                                                    const ROOT::RVec<unsigned char> &dm) {
                            int counts = 0;
                            Logger::get("finalstate::n_taus")->debug("Before starting counting we have {} counts", counts);
                            for (unsigned int i = 0; i < pt.size(); ++i) {
                                if (pt[i] > TauPtThreshold &&
                                    fabs(eta[i]) < TauEtaThreshold &&
                                    fabs(dz[i]) < TauDZThreshold &&
                                    vs_jet[i] > TauVsJetThreshold &&
                                    vs_muon[i] > TauVsMuThreshold &&
                                    vs_ele[i] > TauVsEleThreshold &&
                                    std::find(SelectedDecayModes.begin(), SelectedDecayModes.end(), dm[i]) != SelectedDecayModes.end()) {
                                    ++counts;
                                }
                            }
                            return counts;

                        
                        };
                        
                        // Apply the lambda function and define a new column with the counts
                        auto df1 =
                         df.Define(outputname, count_taus_lambda,
                          {tau_pt, tau_eta, tau_dz, tau_vs_jet, tau_vs_mu, tau_vs_ele, tau_dm});


                        return df1;

                                 }


ROOT::RDF::RNode n_muons(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &muon_pt,
                                 const std::string &muon_eta,
                                 const float &MuPtThreshold,
                                 const float &MuEtaThreshold) {


                        // Lambda function to count taus passing the cuts
                        auto count_muons_lambda = [=](const ROOT::RVec<float> &pt,
                                                    const ROOT::RVec<float> &eta) {
                            int counts = 0;
                            
                            for (unsigned int i = 0; i < pt.size(); ++i) {
                                if (pt[i] > MuPtThreshold &&
                                    fabs(eta[i]) < MuEtaThreshold ) {
                                    ++counts;
                                }
                            }
                            return counts;

                        
                        };
                        
                        // Apply the lambda function and define a new column with the counts
                        auto df1 =
                         df.Define(outputname, count_muons_lambda,
                          {muon_pt, muon_eta});


                        return df1;

                                 }


ROOT::RDF::RNode n_eles(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &ele_pt,
                                 const std::string &ele_eta,
                                 const float &ElePtThreshold,
                                 const float &EleEtaThreshold) {


                        // Lambda function to count taus passing the cuts
                        auto count_eles_lambda = [=](const ROOT::RVec<float> &pt,
                                                    const ROOT::RVec<float> &eta) {
                            int counts = 0;
                            
                            for (unsigned int i = 0; i < pt.size(); ++i) {
                                if (pt[i] > EleEtaThreshold &&
                                    fabs(eta[i]) < EleEtaThreshold ) {
                                    ++counts;
                                }
                            }
                            return counts;

                        
                        };
                        
                        // Apply the lambda function and define a new column with the counts
                        auto df1 =
                         df.Define(outputname, count_eles_lambda,
                          {ele_pt, ele_eta});


                        return df1;

                                 }

    
} // end namespace finalstate
#endif /* GUARDFINALSTATE_H */

