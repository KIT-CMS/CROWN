#ifndef YHKinFitMaster_H
#define YHKinFitMaster_H

#include "ROOT/RVec.hxx"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TMatrixD.h>
#include <map>

class YHKinFitMaster {
  public:
    YHKinFitMaster(ROOT::Math::PtEtaPhiEVector bjet1, float bjet_reso_1,
                   ROOT::Math::PtEtaPhiEVector bjet2, float bjet_reso_2,
                   ROOT::Math::PtEtaPhiEVector tauvis1,
                   ROOT::Math::PtEtaPhiEVector tauvis2,
                   ROOT::Math::PtEtaPhiEVector met, TMatrixD met_cov,
                   bool Ytautau);

    void doFullFit();
    void Fit(int mh, int mY);

    // Hypotheses
    void addMhHypothesis(std::vector<int> v);
    void addMYHypothesis(std::vector<int> v);

    // Covariance matrix
    TMatrixD CalcCov(ROOT::Math::PtEtaPhiEVector p4, double dE);

    // Resolution
    double CalcBjetResolution(ROOT::Math::PtEtaPhiEVector p4, double res);
    // double GetBjetResolution(double eta, double et);

    double ConstrainEnergy(ROOT::Math::PtEtaPhiEVector p4_mother,
                           ROOT::Math::PtEtaPhiEVector p4_1,
                           ROOT::Math::PtEtaPhiEVector p4_2, int mHypo);

    double Chi2_V4(ROOT::Math::PtEtaPhiEVector p4_reco,
                   ROOT::Math::PtEtaPhiEVector p4_fit, double res);
    double Chi2_Balance(ROOT::Math::PtEtaPhiEVector p4_X_fit);

    double GetPullBalanceX(ROOT::Math::PtEtaPhiEVector p4_X_fit);
    double GetPullBalanceY(ROOT::Math::PtEtaPhiEVector p4_X_fit);

    // Getters for fit results
    std::map<std::pair<int, int>, double> getChi2FullFit() {
        return m_fullFitResultChi2;
    }
    std::map<std::pair<int, int>, double> getChi2B1FullFit() {
        return m_fullFitResultChi2BJet1;
    }
    std::map<std::pair<int, int>, double> getChi2B2FullFit() {
        return m_fullFitResultChi2BJet2;
    }
    std::map<std::pair<int, int>, double> getChi2BalanceFullFit() {
        return m_fullFitResultChi2Balance;
    }
    std::map<std::pair<int, int>, double> getFitProbFullFit() {
        return m_fullFitResultFitProb;
    }
    std::map<std::pair<int, int>, double> getMXFullFit() {
        return m_fullFitResultMX;
    }
    std::map<std::pair<int, int>, double> getMYFullFit() {
        return m_fullFitResultMY;
    }
    std::map<std::pair<int, int>, double> getMhFullFit() {
        return m_fullFitResultMh;
    }
    std::map<std::pair<int, int>, double> getPullBalanceFullFitX() {
        return m_fullFitPullBalanceX;
    }
    std::map<std::pair<int, int>, double> getPullBalanceFullFitY() {
        return m_fullFitPullBalanceY;
    }
    std::map<std::pair<int, int>, int> getConvergenceFullFit() {
        return m_fullFitConvergence;
    }

    double getBestChi2FullFit() { return m_bestChi2FullFit; }
    double getBestMXFullFit() { return m_bestMXFullFit; }
    std::pair<int, int> getBestHypoFullFit() { return m_bestHypoFullFit; }

  private:
    // input vectors
    ROOT::Math::PtEtaPhiEVector m_bjet1;
    ROOT::Math::PtEtaPhiEVector m_bjet2;
    float m_bjet_reso1;
    float m_bjet_reso2;
    ROOT::Math::PtEtaPhiEVector m_tauvis1;
    ROOT::Math::PtEtaPhiEVector m_tauvis2;
    bool m_Ytautau;

    ROOT::Math::PtEtaPhiEVector m_MET;
    TMatrixD m_MET_COV;
    TMatrixD m_covRecoil;
    TMatrixD m_V_inv;

    // hypotheses
    std::vector<int> m_mh;
    std::vector<int> m_mY;

    double mHttHypo;
    double mHbbHypo;

    double m_chi2;
    double m_chi2_b1;
    double m_chi2_b2;
    double m_chi2_balance;
    ROOT::Math::PtEtaPhiEVector p4_X_fit;
    ROOT::Math::PtEtaPhiEVector p4_Y_fit;
    ROOT::Math::PtEtaPhiEVector p4_h_fit;

    int m_convergence;
    double m_fitted_mX;
    double m_fitted_mY;
    double m_fitted_mh;

    // full event fit
    std::map<std::pair<int, int>, double> m_fullFitResultChi2;
    std::map<std::pair<int, int>, double> m_fullFitResultChi2BJet1;
    std::map<std::pair<int, int>, double> m_fullFitResultChi2BJet2;
    std::map<std::pair<int, int>, double> m_fullFitResultChi2Balance;
    std::map<std::pair<int, int>, double> m_fullFitResultFitProb;
    std::map<std::pair<int, int>, double> m_fullFitResultMX;
    std::map<std::pair<int, int>, double> m_fullFitResultMY;
    std::map<std::pair<int, int>, double> m_fullFitResultMh;
    std::map<std::pair<int, int>, double> m_fullFitPullBalanceX;
    std::map<std::pair<int, int>, double> m_fullFitPullBalanceY;
    std::map<std::pair<int, int>, int> m_fullFitConvergence;

    double m_bestChi2FullFit;
    double m_bestMXFullFit;
    double m_bestMYFullFit;
    double m_bestMhFullFit;
    std::pair<int, int> m_bestHypoFullFit;
};

#endif