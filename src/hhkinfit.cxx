#ifndef GUARDHHKINFIT_H
#define GUARDHHKINFIT_H
/// The namespace that contains the HHKinFit function.
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "../include/HHKinFit/YHKinFitMaster.hxx"

#include "TH1F.h"
#include <math.h>

namespace hhkinfit {
/**
 * @brief Function to get a specific value from a vector
 *
 * @param idx position of the wanted quantity in a vector
 * @returns a specific result from the results vector
 */
auto single_output(const int &idx) {
    return [idx](const ROOT::RVec<float> &result) { return result[idx]; };
};
/**
 * @brief Function to run a kinematic fit of a X -> YH di-Higgs system with a
 * bb+tautau final state. Code for calculation based on
 * https://github.com/janekbechtel/HHKinFit
 *
 * @param df the input dataframe
 * @param outputname_1 name of the output column for the convergence status of
 * the fit
 * @param outputname_2 name of the output column for the estimated mass of the X
 * resonance
 * @param outputname_3 name of the output column for the estimated mass of the Y
 * resonance
 * @param outputname_4 name of the output column for the estimated mass of the
 * SM Higgs boson
 * @param outputname_5 name of the output column for the ch2 value
 * @param outputname_6 name of the output column for the probability of chi2
 * value
 * @param tau_pt_1 name of the column containing the pt of the first tau in the
 * tau pair
 * @param tau_eta_1 name of the column containing the eta of the first tau in
 * the tau pair
 * @param tau_phi_1 name of the column containing the phi of the first tau in
 * the tau pair
 * @param tau_mass_1 name of the column containing the mass of the first tau in
 * the tau pair
 * @param tau_pt_2 name of the column containing the pt of the second tau in the
 * tau pair
 * @param tau_eta_2 name of the column containing the eta of the second tau in
 * the tau pair
 * @param tau_phi_2 name of the column containing the phi of the second tau in
 * the tau pair
 * @param tau_mass_2 name of the column containing the mass of the second tau in
 * the tau pair
 * @param b_pt_1 name of the column containing the pt of the first b-jet in the
 * bb pair
 * @param b_eta_1 name of the column containing the eta of the first b-jet in
 * the bb pair
 * @param b_phi_1 name of the column containing the phi of the first b-jet in
 * the bb pair
 * @param b_mass_1 name of the column containing the mass of the first b-jet in
 * the bb pair
 * @param b_reso_1 name of the column containing the pt resolution of the first
 * b-jet in the bb pair
 * @param b_pt_2 name of the column containing the pt of the second b-jet in the
 * bb pair
 * @param b_eta_2 name of the column containing the eta of the second b-jet in
 * the bb pair
 * @param b_phi_2 name of the column containing the phi of the second b-jet in
 * the bb pair
 * @param b_mass_2 name of the column containing the mass of the second b-jet in
 * the bb pair
 * @param b_reso_2 name of the column containing the pt resolution of the second
 * b-jet in the bb pair
 * @param met name of the column containing the met pt
 * @param met_phi name of the column containing the met phi
 * @param met_cov00 name of the column containing the met covariance xx
 * @param met_cov01 name of the column containing the met covariance xy
 * @param met_cov10 name of the column containing the met covariance yx
 * @param met_cov11 name of the column containing the met covariance yy
 * @param YDecay name of the Y resonace decay, either "YToTauTau" or "YToBB"
 * @returns a dataframe with all outputs of the kinematic fit
 */
ROOT::RDF::RNode
YHKinFit(ROOT::RDF::RNode df, const std::string &outputname_1,
         const std::string &outputname_2, const std::string &outputname_3,
         const std::string &outputname_4, const std::string &outputname_5,
         const std::string &outputname_6,
         const std::string &tau_pt_1, const std::string &tau_eta_1,
         const std::string &tau_phi_1, const std::string &tau_mass_1,
         const std::string &tau_pt_2, const std::string &tau_eta_2,
         const std::string &tau_phi_2, const std::string &tau_mass_2,
         const std::string &b_pt_1, const std::string &b_eta_1,
         const std::string &b_phi_1, const std::string &b_mass_1,
         const std::string &b_reso_1, const std::string &b_pt_2,
         const std::string &b_eta_2, const std::string &b_phi_2,
         const std::string &b_mass_2, const std::string &b_reso_2,
         const std::string &met, const std::string &met_phi,
         const std::string &met_cov00, const std::string &met_cov01,
         const std::string &met_cov10, const std::string &met_cov11,
         const std::string &YDecay) {
    Logger::get("YHKinFit" + YDecay)
        ->debug("Fitting bbtautau system to get estimation for X mass.");

    auto kin_fit = [YDecay](const float &tau_pt_1, const float &tau_eta_1,
                            const float &tau_phi_1, const float &tau_mass_1,
                            const float &tau_pt_2, const float &tau_eta_2,
                            const float &tau_phi_2, const float &tau_mass_2,
                            const float &b_pt_1, const float &b_eta_1,
                            const float &b_phi_1, const float &b_mass_1,
                            const float &b_reso_1, const float &b_pt_2,
                            const float &b_eta_2, const float &b_phi_2,
                            const float &b_mass_2, const float &b_reso_2,
                            const float &met, const float &met_phi,
                            const float &met_cov00, const float &met_cov01,
                            const float &met_cov10, const float &met_cov11) {
        auto kinfit_mX = -10.;
        auto kinfit_mY = -10.;
        auto kinfit_mh = -10.;
        auto kinfit_chi2 = 999.;
        auto kinfit_prob = 0.;
        auto kinfit_convergence = -1.;

        if ((tau_pt_1 > 0.) && (tau_pt_2 > 0.) && (b_pt_1 > 0.) &&
            (b_pt_2 > 0.)) {
            ROOT::Math::PtEtaPhiEVector tau_1 = (ROOT::Math::PtEtaPhiEVector) ROOT::Math::PtEtaPhiMVector(
                tau_pt_1, tau_eta_1, tau_phi_1, tau_mass_1);
            ROOT::Math::PtEtaPhiEVector tau_2 = (ROOT::Math::PtEtaPhiEVector) ROOT::Math::PtEtaPhiMVector(
                tau_pt_2, tau_eta_2, tau_phi_2, tau_mass_2);
            ROOT::Math::PtEtaPhiEVector b_1 =
                (ROOT::Math::PtEtaPhiEVector) ROOT::Math::PtEtaPhiMVector(b_pt_1, b_eta_1, b_phi_1, b_mass_1);
            ROOT::Math::PtEtaPhiEVector b_2 =
                (ROOT::Math::PtEtaPhiEVector) ROOT::Math::PtEtaPhiMVector(b_pt_2, b_eta_2, b_phi_2, b_mass_2);

            ROOT::Math::PtEtaPhiEVector met_LV =
                (ROOT::Math::PtEtaPhiEVector) ROOT::Math::PtEtaPhiMVector(met, 0., met_phi, 0.);
            TMatrixD met_cov(2, 2);
            met_cov[0][0] = met_cov00;
            met_cov[1][0] = met_cov10;
            met_cov[0][1] = met_cov01;
            met_cov[1][1] = met_cov11;

            std::vector<int> hypo_mh = {125};
            // std::vector<int> hypo_mY =
            // {5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,125,130,140,150,180,210,240,270,300,330,360,390,420,450,480,510,540,570,600,630,660,690,720,750,780,810,840,870,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000};
            std::vector<int> hypo_mY = {
                50,   60,   70,   80,   90,   95,   100,  125,
                150,  250,  300,  400,  500,  600,  700,  800,
                900,  1000, 1100, 1200, 1300, 1400, 1600, 1800,
                2000, 2200, 2400, 2500, 2600, 2800, 3000, 3500};
            bool YToTauTau = false;
            if (YDecay == "YToTauTau") {
                YToTauTau = true;
            }

            YHKinFitMaster kinFits =
                YHKinFitMaster(b_1, b_reso_1, b_2, b_reso_2, tau_1, tau_2,
                               met_LV, met_cov, YToTauTau);
            kinFits.addMhHypothesis(hypo_mh);
            kinFits.addMYHypothesis(hypo_mY);

            kinFits.doFullFit();

            std::pair<int, int> bestHypo = kinFits.getBestHypoFullFit();
            Logger::get("YHKinFit" + YDecay)
                ->debug("best hypothesis: tautau {}, bb {}", bestHypo.first,
                        bestHypo.second);

            if (bestHypo.second > 0) {
                std::map<std::pair<int, int>, double> fit_results_chi2 =
                    kinFits.getChi2FullFit();
                std::map<std::pair<int, int>, double> fit_results_fitprob =
                    kinFits.getFitProbFullFit();
                std::map<std::pair<int, int>, double> fit_results_mX =
                    kinFits.getMXFullFit();
                std::map<std::pair<int, int>, double> fit_results_mY =
                    kinFits.getMYFullFit();
                std::map<std::pair<int, int>, double> fit_results_mh =
                    kinFits.getMhFullFit();
                std::map<std::pair<int, int>, int> fit_convergence =
                    kinFits.getConvergenceFullFit();

                kinfit_convergence = fit_convergence.at(bestHypo);
                kinfit_mX = fit_results_mX.at(bestHypo);
                kinfit_mY = fit_results_mY.at(bestHypo);
                kinfit_mh = fit_results_mh.at(bestHypo);
                kinfit_chi2 = fit_results_chi2.at(bestHypo);
                kinfit_prob = fit_results_fitprob.at(bestHypo);
            }
            Logger::get("YHKinFit" + YDecay)
                ->debug("kinfit_convergence: {}", kinfit_convergence);
            Logger::get("YHKinFit" + YDecay)->debug("kinfit_mX: {}", kinfit_mX);
            Logger::get("YHKinFit" + YDecay)->debug("kinfit_mY: {}", kinfit_mY);
            Logger::get("YHKinFit" + YDecay)->debug("kinfit_mh: {}", kinfit_mh);
            Logger::get("YHKinFit" + YDecay)
                ->debug("kinfit_chi2: {}", kinfit_chi2);
            Logger::get("YHKinFit" + YDecay)
                ->debug("kinfit_prob: {}", kinfit_prob);
        }

        ROOT::RVec<float> result = {
            (float)kinfit_convergence, (float)kinfit_mX,
            (float)kinfit_mY,          (float)kinfit_mh,
            (float)kinfit_chi2,        (float)kinfit_prob};
        return result;
    };

    std::string variation = "";
    if (outputname_1.find("__") != std::string::npos) {
        size_t pos = outputname_1.find("__");
        variation = outputname_1.substr(pos);;
    }
    
    std::string result_vec_name = "HYKinFit_vector_" + YDecay + "_resolved" + variation;
    if (outputname_1.find("boosted") != std::string::npos) {
        result_vec_name = "HYKinFit_vector_" + YDecay + "_boosted" + variation;
    }

    auto df1 = df.Define(
        result_vec_name, kin_fit,
        {tau_pt_1,  tau_eta_1,  tau_phi_1, tau_mass_1, tau_pt_2,  tau_eta_2,
         tau_phi_2, tau_mass_2, b_pt_1,    b_eta_1,    b_phi_1,   b_mass_1,
         b_reso_1,  b_pt_2,     b_eta_2,   b_phi_2,    b_mass_2,  b_reso_2,
         met,       met_phi,    met_cov00, met_cov01,  met_cov10, met_cov11});

    auto df2 =
        df1.Define(outputname_1, hhkinfit::single_output(0), {result_vec_name});
    auto df3 =
        df2.Define(outputname_2, hhkinfit::single_output(1), {result_vec_name});
    auto df4 =
        df3.Define(outputname_3, hhkinfit::single_output(2), {result_vec_name});
    auto df5 =
        df4.Define(outputname_4, hhkinfit::single_output(3), {result_vec_name});
    auto df6 =
        df5.Define(outputname_5, hhkinfit::single_output(4), {result_vec_name});
    auto df7 =
        df6.Define(outputname_6, hhkinfit::single_output(5), {result_vec_name});

    return df7;
}
/**
 * @brief Function to compare the chi2 results of kinematic fits for two
 * different decays X -> Y(tautau)H(bb) and X -> Y(bb)H(tautau)
 *
 * @param outputname_1 name of the output column for the convergence status of
 * the fit
 * @param outputname_2 name of the output column for the estimated mass of the X
 * resonance
 * @param outputname_3 name of the output column for the estimated mass of the Y
 * resonance
 * @param outputname_4 name of the output column for the estimated mass of the
 * SM Higgs boson
 * @param outputname_5 name of the output column for the ch2 value
 * @param outputname_6 name of the output column for the probability of chi2
 * value
 * @param kinfit_convergence_YToBB name of the column containing the convergence
 * status of the fit for X -> Y(bb)H(tautau)
 * @param kinfit_mX_YToBB name of the column containing the estimated mass of
 * the X resonance for X -> Y(bb)H(tautau)
 * @param kinfit_mY_YToBB name of the column containing the estimated mass of
 * the Y resonance for X -> Y(bb)H(tautau)
 * @param kinfit_mh_YToBB name of the column containing the estimated mass of
 * the SM Higgs boson for X -> Y(bb)H(tautau)
 * @param kinfit_chi2_YToBB name of the column containing the ch2 value for X ->
 * Y(bb)H(tautau)
 * @param kinfit_prob_YToBB name of the column containing the probability of
 * chi2 value for X -> Y(bb)H(tautau)
 * @param kinfit_convergence_YToTauTau name of the column containing the
 * convergence status of the fit for X -> Y(tautau)H(bb)
 * @param kinfit_mX_YToTauTau name of the column containing the estimated mass
 * of the X resonance for X -> Y(tautau)H(bb)
 * @param kinfit_mY_YToTauTau name of the column containing the estimated mass
 * of the Y resonance for X -> Y(tautau)H(bb)
 * @param kinfit_mh_YToTauTau name of the column containing the estimated mass
 * of the SM Higgs boson for X -> Y(tautau)H(bb)
 * @param kinfit_chi2_YToTauTau name of the column containing the ch2 value for
 * X -> Y(tautau)H(bb)
 * @param kinfit_prob_YToTauTau name of the column containing the probability of
 * chi2 value for X -> Y(tautau)H(bb)
 * @returns a dataframe with all outputs of the best kinematic fit between the
 * two considered decays
 */
ROOT::RDF::RNode BestYHKinFit(
    ROOT::RDF::RNode df, const std::string &outputname_1,
    const std::string &outputname_2, const std::string &outputname_3,
    const std::string &outputname_4, const std::string &outputname_5,
    const std::string &outputname_6,
    const std::string &kinfit_convergence_YToBB,
    const std::string &kinfit_mX_YToBB, const std::string &kinfit_mY_YToBB,
    const std::string &kinfit_mh_YToBB, const std::string &kinfit_chi2_YToBB,
    const std::string &kinfit_prob_YToBB,
    const std::string &kinfit_convergence_YToTauTau,
    const std::string &kinfit_mX_YToTauTau,
    const std::string &kinfit_mY_YToTauTau,
    const std::string &kinfit_mh_YToTauTau,
    const std::string &kinfit_chi2_YToTauTau,
    const std::string &kinfit_prob_YToTauTau) {
    Logger::get("BestYHKinFit")
        ->debug("Decide on the best YHKinFit between the Y(tautau)H(bb) and "
                "Y(bb)H(tautau) cases.");

    auto best_kin_fit =
        [](const float &kinfit_convergence_YToBB, const float &kinfit_mX_YToBB,
           const float &kinfit_mY_YToBB, const float &kinfit_mh_YToBB,
           const float &kinfit_chi2_YToBB, const float &kinfit_prob_YToBB,
           const float &kinfit_convergence_YToTauTau,
           const float &kinfit_mX_YToTauTau, const float &kinfit_mY_YToTauTau,
           const float &kinfit_mh_YToTauTau, const float &kinfit_chi2_YToTauTau,
           const float &kinfit_prob_YToTauTau) {
            auto kinfit_mX = -10.;
            auto kinfit_mY = -10.;
            auto kinfit_mh = -10.;
            auto kinfit_chi2 = 999.;
            auto kinfit_prob = 0.;
            auto kinfit_convergence = -1.;

            if ((kinfit_mX_YToBB > 0.) || (kinfit_mX_YToTauTau > 0.)) {
                if (kinfit_chi2_YToBB < kinfit_chi2_YToTauTau) {
                    ROOT::RVec<float> result = {
                        kinfit_convergence_YToBB, kinfit_mX_YToBB,
                        kinfit_mY_YToBB,          kinfit_mh_YToBB,
                        kinfit_chi2_YToBB,        kinfit_prob_YToBB};
                    return result;
                } else if (kinfit_chi2_YToBB >= kinfit_chi2_YToTauTau) {
                    ROOT::RVec<float> result = {
                        kinfit_convergence_YToTauTau, kinfit_mX_YToTauTau,
                        kinfit_mY_YToTauTau,          kinfit_mh_YToTauTau,
                        kinfit_chi2_YToTauTau,        kinfit_prob_YToTauTau};
                    return result;
                }
            }

            ROOT::RVec<float> result = {
                (float)kinfit_convergence, (float)kinfit_mX,
                (float)kinfit_mY,          (float)kinfit_mh,
                (float)kinfit_chi2,        (float)kinfit_prob};
            return result;
        };

    std::string variation = "";
    if (outputname_1.find("__") != std::string::npos) {
        size_t pos = outputname_1.find("__");
        variation = outputname_1.substr(pos);;
    }
    
    std::string result_vec_name = "HYKinFit_vector_resolved" + variation;
    if (outputname_1.find("boosted") != std::string::npos) {
        result_vec_name = "HYKinFit_vector_boosted" + variation;
    }

    auto df1 = df.Define(
        result_vec_name, best_kin_fit,
        {kinfit_convergence_YToBB, kinfit_mX_YToBB, kinfit_mY_YToBB,
         kinfit_mh_YToBB, kinfit_chi2_YToBB, kinfit_prob_YToBB,
         kinfit_convergence_YToTauTau, kinfit_mX_YToTauTau, kinfit_mY_YToTauTau,
         kinfit_mh_YToTauTau, kinfit_chi2_YToTauTau, kinfit_prob_YToTauTau});

    auto df2 =
        df1.Define(outputname_1, hhkinfit::single_output(0), {result_vec_name});
    auto df3 =
        df2.Define(outputname_2, hhkinfit::single_output(1), {result_vec_name});
    auto df4 =
        df3.Define(outputname_3, hhkinfit::single_output(2), {result_vec_name});
    auto df5 =
        df4.Define(outputname_4, hhkinfit::single_output(3), {result_vec_name});
    auto df6 =
        df5.Define(outputname_5, hhkinfit::single_output(4), {result_vec_name});
    auto df7 =
        df6.Define(outputname_6, hhkinfit::single_output(5), {result_vec_name});

    return df7;
}

} // namespace hhkinfit
#endif /* GUARDHHKINFIT_H */