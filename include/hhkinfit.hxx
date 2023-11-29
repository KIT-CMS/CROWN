#ifndef GUARDHHKINFIT_H
#define GUARDHHKINFIT_H

namespace hhkinfit {
auto single_output(const int &idx);

ROOT::RDF::RNode
YHKinFit(ROOT::RDF::RNode df, const std::string &outputname_1, const std::string &outputname_2,
            const std::string &outputname_3, const std::string &outputname_4,
            const std::string &outputname_5, const std::string &outputname_6,
            const std::string &outputname_7, const std::string &outputname_8,
            const std::string &outputname_9,
            const std::string &tau_pt_1, const std::string &tau_eta_1,
            const std::string &tau_phi_1, const std::string &tau_mass_1,
            const std::string &tau_pt_2, const std::string &tau_eta_2,
            const std::string &tau_phi_2, const std::string &tau_mass_2,
            const std::string &b_pt_1, const std::string &b_eta_1,
            const std::string &b_phi_1, const std::string &b_mass_1, const std::string &b_reso_1,
            const std::string &b_pt_2, const std::string &b_eta_2,
            const std::string &b_phi_2, const std::string &b_mass_2, const std::string &b_reso_2,
            const std::string &met, const std::string &met_phi,
            const std::string &met_cov00, const std::string &met_cov01,
            const std::string &met_cov10, const std::string &met_cov11,
            const std::string &YDecay);
ROOT::RDF::RNode
BestYHKinFit(ROOT::RDF::RNode df, const std::string &outputname_1, const std::string &outputname_2,
            const std::string &outputname_3, const std::string &outputname_4,
            const std::string &outputname_5, const std::string &outputname_6,
            const std::string &outputname_7, const std::string &outputname_8,
            const std::string &outputname_9, const std::string &kinfit_convergence_YToBB,
            const std::string &kinfit_mX_YToBB, const std::string &kinfit_mY_YToBB,
            const std::string &kinfit_mh_YToBB, const std::string &kinfit_chi2_YToBB,
            const std::string &kinfit_prob_YToBB, const std::string &kinfit_pull1_YToBB,
            const std::string &kinfit_pull2_YToBB, const std::string &kinfit_pullBalance_YToBB,
            const std::string &kinfit_convergence_YToTauTau,
            const std::string &kinfit_mX_YToTauTau, const std::string &kinfit_mY_YToTauTau,
            const std::string &kinfit_mh_YToTauTau, const std::string &kinfit_chi2_YToTauTau,
            const std::string &kinfit_prob_YToTauTau, const std::string &kinfit_pull1_YToTauTau,
            const std::string &kinfit_pull2_YToTauTau, const std::string &kinfit_pullBalance_YToTauTau);
} // namespace hhkinfit
#endif /* GUARDHHKINFIT_H */