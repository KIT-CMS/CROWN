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
namespace mutau{

// ROOT::RDF::RNode single_mu_in_fatjet_mutau(ROOT::RDF::RNode df,
//                                              const std::string &outputname,
//                                              const std::string &fatjet_p4,
//                                              const std::string &muon_pt,
//                                              const std::string &muon_eta,
//                                              const std::string &muon_phi,
//                                              const std::string &muon_mass,
//                                              const float &pt_threshold,
//                                              const float &eta_threshold) {

//     auto single_muon_in_fatjet = [=]
//                                  (const ROOT::Math::PtEtaPhiMVector &fatjet_p4,
//                                   const ROOT::RVec<float> &muon_pts,
//                                   const ROOT::RVec<float> &muon_etas,
//                                   const ROOT::RVec<float> &muon_phis,
//                                   const ROOT::RVec<float> &muon_masses) {
//         // Check if there is exactly one muon in the event

//         float mutau_flag = 0;

//         if (muon_pts.size() != 1) {
//             Logger::get("fatjet::trigger_single_mu_in_fatjet")
//                 ->info("Muon count not equal to 1: found {} muons", muon_pts.size());
//             // return 0;
//         }

//         // Access the single muon's properties
//         float muon_pt = muon_pts[0];
//         float muon_eta = muon_etas[0];

//         // Apply pt, eta, and isolation cuts
//         if (muon_pt < pt_threshold || std::abs(muon_eta) > eta_threshold ) {
//             Logger::get("fatjet::trigger_single_mu_in_fatjet")
//                 ->info("Single muon does not pass thresholds. pt: {}, eta: {}", 
//                         muon_pt, muon_eta);
//             // return 0;
//         }else{

//         // Construct the muon 4-vector
//         ROOT::Math::PtEtaPhiMVector muon_p4(muon_pt, muon_eta, muon_phis[0], muon_masses[0]);

//         // Calculate deltaR between the single muon and the fatjet
//         float delta_r = ROOT::Math::VectorUtil::DeltaR(muon_p4, fatjet_p4);

//         if (delta_r < 0.8) {
//             Logger::get("fatjet::trigger_single_mu_in_fatjet")
//                 ->info("Single muon found in fatjet with pt {} and eta {}", muon_pt, muon_eta);
//             // return 1;
//             mutau_flag = 1;
//         } else {
//             Logger::get("fatjet::trigger_single_mu_in_fatjet")
//                 ->info("Single muon NOT in fatjet. Muon pt: {}, eta: {}", muon_pt, muon_eta);
//             // return 0;
//         }
//             }
//         return mutau_flag;
//     };

//     auto df1 = df.Define(outputname, single_muon_in_fatjet,
//                          {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass});
//     return df1;
// }

ROOT::RDF::RNode single_mu_in_fatjet_mutau(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const float &pt_threshold,
                                             const float &eta_threshold) {

    auto single_muon_in_fatjet = [=]
                                 (const ROOT::Math::PtEtaPhiMVector &fatjet_p4,
                                  const ROOT::RVec<float> &muon_pts,
                                  const ROOT::RVec<float> &muon_etas,
                                  const ROOT::RVec<float> &muon_phis,
                                  const ROOT::RVec<float> &muon_masses) {

        float mutau_flag = 0;
        // Check if there is exactly one muon in the event
        if (muon_pts.size() != 1) {
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count NOT equal to 1: found {} muons", muon_pts.size());
            
        }
        else{

        // Access the single muon's properties
        float muon_pt1 = muon_pts[0];
        float muon_eta1 = muon_etas[0];

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
            ->debug("Muon count equal to 1: found {} muons", muon_pts.size());

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count equal to 1  pt: {}, eta: {}", muon_pt1, muon_eta1);

        if (muon_pt1 > pt_threshold && std::abs(muon_eta1) < eta_threshold){

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon passed pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);

            // Construct the muon 4-vector
            ROOT::Math::PtEtaPhiMVector muon_p4(muon_pt1, muon_eta1, muon_phis[0], muon_masses[0]);

            // Calculate deltaR between the single muon and the fatjet
            float delta_r = ROOT::Math::VectorUtil::DeltaR(muon_p4, fatjet_p4);

            if (delta_r < 0.8){
                mutau_flag = 1;
            }

            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Delta R {} between fatjet and muon with pt: {}, eta: {}", delta_r, muon_pt1, muon_eta1);

        }else{
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon didn't pass pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);
        }

        }
        
        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("MuTau final state flag: {}", mutau_flag);

        return mutau_flag;
    };

    auto df1 = df.Define(outputname, single_muon_in_fatjet,
                         {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass});
    return df1;
}


ROOT::RDF::RNode single_mu_in_fatjet_mutau_deltaR(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const float &pt_threshold,
                                             const float &eta_threshold) {

    auto single_muon_in_fatjet_deltaR = [=]
                                 (const ROOT::Math::PtEtaPhiMVector &fatjet_p4,
                                  const ROOT::RVec<float> &muon_pts,
                                  const ROOT::RVec<float> &muon_etas,
                                  const ROOT::RVec<float> &muon_phis,
                                  const ROOT::RVec<float> &muon_masses) {

        float mutau_deltaR = 27.0;
        // Check if there is exactly one muon in the event
        if (muon_pts.size() != 1) {
            Logger::get("fatjet::single_mu_in_fatjet_mutau_deltaR")
                ->debug("Muon count NOT equal to 1: found {} muons", muon_pts.size());
            
        }
        else{

        // Access the single muon's properties
        float muon_pt1 = muon_pts[0];
        float muon_eta1 = muon_etas[0];

        Logger::get("fatjet::single_mu_in_fatjet_mutau_deltaR")
            ->debug("Muon count equal to 1: found {} muons", muon_pts.size());

        Logger::get("fatjet::single_mu_in_fatjet_mutau_deltaR")
                ->debug("Muon count equal to 1  pt: {}, eta: {}", muon_pt1, muon_eta1);

        if (muon_pt1 > pt_threshold && std::abs(muon_eta1) < eta_threshold){

        Logger::get("fatjet::single_mu_in_fatjet_mutau_deltaR")
                ->debug("Muon passed pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);

            // Construct the muon 4-vector
            ROOT::Math::PtEtaPhiMVector muon_p4(muon_pt1, muon_eta1, muon_phis[0], muon_masses[0]);

            // Calculate deltaR between the single muon and the fatjet
            float delta_r = ROOT::Math::VectorUtil::DeltaR(muon_p4, fatjet_p4);
            mutau_deltaR = delta_r;

        }
        
        }

        Logger::get("fatjet::single_mu_in_fatjet_mutau_deltaR")
                ->debug("MuTau final state flag: {}", mutau_deltaR);

        return mutau_deltaR;
    };

    auto df1 = df.Define(outputname, single_muon_in_fatjet_deltaR,
                         {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass});
    return df1;

}

ROOT::RDF::RNode single_mu_in_fatjet_mutau_mu_pT(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const float &pt_threshold,
                                             const float &eta_threshold) {

    auto single_muon_in_fatjet = [=]
                                 (const ROOT::Math::PtEtaPhiMVector &fatjet_p4,
                                  const ROOT::RVec<float> &muon_pts,
                                  const ROOT::RVec<float> &muon_etas,
                                  const ROOT::RVec<float> &muon_phis,
                                  const ROOT::RVec<float> &muon_masses) {

        float mu_pt_flag = -10;
        // Check if there is exactly one muon in the event
        if (muon_pts.size() != 1) {
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count NOT equal to 1: found {} muons", muon_pts.size());
            
        }
        else{

        // Access the single muon's properties
        float muon_pt1 = muon_pts[0];
        float muon_eta1 = muon_etas[0];

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
            ->debug("Muon count equal to 1: found {} muons", muon_pts.size());

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count equal to 1  pt: {}, eta: {}", muon_pt1, muon_eta1);

        if (muon_pt1 > pt_threshold && std::abs(muon_eta1) < eta_threshold){

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon passed pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);

            // Construct the muon 4-vector
            ROOT::Math::PtEtaPhiMVector muon_p4(muon_pt1, muon_eta1, muon_phis[0], muon_masses[0]);

            // Calculate deltaR between the single muon and the fatjet
            float delta_r = ROOT::Math::VectorUtil::DeltaR(muon_p4, fatjet_p4);

            if (delta_r < 0.8){
                mu_pt_flag = muon_pt1;
            }

            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Delta R {} between fatjet and muon with pt: {}, eta: {}", delta_r, muon_pt1, muon_eta1);

        }else{
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon didn't pass pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);
        }

        }
        
        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Mu pT: {}", mu_pt_flag);

        return mu_pt_flag;
    };

    auto df1 = df.Define(outputname, single_muon_in_fatjet,
                         {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass});
    return df1;
}


ROOT::RDF::RNode single_mu_in_fatjet_mutau_mu_eta(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const float &pt_threshold,
                                             const float &eta_threshold) {

    auto single_muon_in_fatjet = [=]
                                 (const ROOT::Math::PtEtaPhiMVector &fatjet_p4,
                                  const ROOT::RVec<float> &muon_pts,
                                  const ROOT::RVec<float> &muon_etas,
                                  const ROOT::RVec<float> &muon_phis,
                                  const ROOT::RVec<float> &muon_masses) {

        float mu_eta_flag = -10;
        // Check if there is exactly one muon in the event
        if (muon_pts.size() != 1) {
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count NOT equal to 1: found {} muons", muon_pts.size());
            
        }
        else{

        // Access the single muon's properties
        float muon_pt1 = muon_pts[0];
        float muon_eta1 = muon_etas[0];

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
            ->debug("Muon count equal to 1: found {} muons", muon_pts.size());

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count equal to 1  pt: {}, eta: {}", muon_pt1, muon_eta1);

        if (muon_pt1 > pt_threshold && std::abs(muon_eta1) < eta_threshold){

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon passed pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);

            // Construct the muon 4-vector
            ROOT::Math::PtEtaPhiMVector muon_p4(muon_pt1, muon_eta1, muon_phis[0], muon_masses[0]);

            // Calculate deltaR between the single muon and the fatjet
            float delta_r = ROOT::Math::VectorUtil::DeltaR(muon_p4, fatjet_p4);

            if (delta_r < 0.8){
                mu_eta_flag = muon_eta1;
            }

            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Delta R {} between fatjet and muon with pt: {}, eta: {}", delta_r, muon_pt1, muon_eta1);

        }else{
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon didn't pass pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);
        }

        }
        
        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Mu eta: {}", mu_eta_flag);

        return mu_eta_flag;
    };

    auto df1 = df.Define(outputname, single_muon_in_fatjet,
                         {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass});
    return df1;
}

ROOT::RDF::RNode single_mu_in_fatjet_mutau_mu_iso(ROOT::RDF::RNode df,
                                             const std::string &outputname,
                                             const std::string &fatjet_p4,
                                             const std::string &muon_pt,
                                             const std::string &muon_eta,
                                             const std::string &muon_phi,
                                             const std::string &muon_mass,
                                             const std::string &muon_iso,
                                             const float &pt_threshold,
                                             const float &eta_threshold) {

    auto single_muon_in_fatjet = [=]
                                 (const ROOT::Math::PtEtaPhiMVector &fatjet_p4,
                                  const ROOT::RVec<float> &muon_pts,
                                  const ROOT::RVec<float> &muon_etas,
                                  const ROOT::RVec<float> &muon_phis,
                                  const ROOT::RVec<float> &muon_masses,
                                  const ROOT::RVec<float> &muon_isos) {

        float mu_iso_flag = -10;
        // Check if there is exactly one muon in the event
        if (muon_pts.size() != 1) {
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count NOT equal to 1: found {} muons", muon_pts.size());
            
        }
        else{

        // Access the single muon's properties
        float muon_pt1 = muon_pts[0];
        float muon_eta1 = muon_etas[0];
        float muon_iso1 = muon_isos[0];

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
            ->debug("Muon count equal to 1: found {} muons", muon_pts.size());

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon count equal to 1  pt: {}, eta: {}", muon_pt1, muon_eta1);

        if (muon_pt1 > pt_threshold && std::abs(muon_eta1) < eta_threshold){

        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon passed pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);

            // Construct the muon 4-vector
            ROOT::Math::PtEtaPhiMVector muon_p4(muon_pt1, muon_eta1, muon_phis[0], muon_masses[0]);

            // Calculate deltaR between the single muon and the fatjet
            float delta_r = ROOT::Math::VectorUtil::DeltaR(muon_p4, fatjet_p4);

            if (delta_r < 0.8){
                mu_iso_flag = muon_iso1;
            }

            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Delta R {} between fatjet and muon with pt: {}, eta: {}", delta_r, muon_pt1, muon_eta1);

        }else{
            Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Muon didn't pass pT-eta threshold  pt: {}, eta: {}", muon_pt1, muon_eta1);
        }

        }
        
        Logger::get("fatjet::trigger_single_mu_in_fatjet")
                ->debug("Mu eta: {}", mu_iso_flag);

        return mu_iso_flag;
    };

    auto df1 = df.Define(outputname, single_muon_in_fatjet,
                         {fatjet_p4, muon_pt, muon_eta, muon_phi, muon_mass, muon_iso});
    return df1;
}

} // end of the mutau namespace
    
} // end namespace finalstate
#endif /* GUARDFINALSTATE_H */

