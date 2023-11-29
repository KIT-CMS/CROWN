#include "../../include/HHKinFit/YHKinFitMaster.hxx"
#include "../../include/HHKinFit/PSFit.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include "TMatrixDEigen.h"
#include "../include/utility/Logger.hxx"


YHKinFitMaster::YHKinFitMaster(ROOT::Math::PtEtaPhiMVector bjet1, float bjet_reso_1, ROOT::Math::PtEtaPhiMVector bjet2, float bjet_reso_2, ROOT::Math::PtEtaPhiMVector tauvis1, ROOT::Math::PtEtaPhiMVector tauvis2, ROOT::Math::PtEtaPhiMVector met, TMatrixD met_cov, bool Ytautau):
    m_bjet1(bjet1),
    m_bjet2(bjet2),
    m_bjet_reso1(bjet_reso_1),
    m_bjet_reso2(bjet_reso_2),
    m_tauvis1(tauvis1),
    m_tauvis2(tauvis2),
    m_Ytautau(Ytautau),
    m_MET(met),
    m_MET_COV(met_cov),

    m_covRecoil(2, 2),
    m_mh(std::vector<int>()),
    m_mY(std::vector<int>()),
    mHttHypo(-10.),
    mHbbHypo(-10.),
    p4_X_fit(ROOT::Math::PtEtaPhiEVector(-10, -10, -10, -10)),
    p4_Y_fit(ROOT::Math::PtEtaPhiEVector(-10, -10, -10, -10)),
    p4_h_fit(ROOT::Math::PtEtaPhiEVector(-10, -10, -10, -10)),

    m_convergence(0),
    m_fitted_mX(-10),
    m_fitted_mY(-10),

    m_fullFitResultChi2(std::map< std::pair<int, int> , double>()),
    m_fullFitResultMX(std::map< std::pair<int, int> , double>()),
    m_bestChi2FullFit(999),
    m_bestMXFullFit(-1),
    m_bestHypoFullFit(std::pair<int, int>(-1,-1))
{}

void YHKinFitMaster::doFullFit()
{
    //loop over all hypotheses
    for (auto &mh: m_mh) {
        for (auto &mY: m_mY) {
            m_convergence = 0;
            Fit(mh, mY);

            double prob = TMath::Prob(m_chi2, 2);
            std::pair< int, int > hypo_full (mHttHypo, mHbbHypo);
            std::pair< std::pair<int, int>, double > entry_chi2_full (hypo_full, m_chi2);
            std::pair< std::pair<int, int>, double > entry_chi2_bjet1 (hypo_full, m_chi2_b1);
            std::pair< std::pair<int, int>, double > entry_chi2_bjet2 (hypo_full, m_chi2_b2);
            std::pair< std::pair<int, int>, double > entry_chi2_balance (hypo_full, m_chi2_balance);
            std::pair< std::pair<int, int>, double > entry_fitprob_full (hypo_full, prob);
            std::pair< std::pair<int, int>, double > entry_mX_full (hypo_full, m_fitted_mX);
            std::pair< std::pair<int, int>, double > entry_mY_full (hypo_full, m_fitted_mY);
            std::pair< std::pair<int, int>, double > entry_mh_full (hypo_full, m_fitted_mh);
            std::pair< std::pair<int, int>, double > entry_pullb1_full (hypo_full, sqrt(m_chi2_b1));
            std::pair< std::pair<int, int>, double > entry_pullb2_full (hypo_full, sqrt(m_chi2_b2));
            std::pair< std::pair<int, int>, double > entry_pullbalance_full (hypo_full, sqrt(m_chi2_balance));
            std::pair< std::pair<int, int>, double > entry_pullbalance_fullX (hypo_full, GetPullBalanceX(p4_X_fit));
            std::pair< std::pair<int, int>, double > entry_pullbalance_fullY (hypo_full, GetPullBalanceY(p4_X_fit));
            std::pair< std::pair<int, int>, int > entry_convergence_full (hypo_full, m_convergence);

            m_fullFitResultChi2.insert(entry_chi2_full);
            m_fullFitResultChi2BJet1.insert(entry_chi2_bjet1);
            m_fullFitResultChi2BJet2.insert(entry_chi2_bjet2);
            m_fullFitResultChi2Balance.insert(entry_chi2_balance);
            m_fullFitResultFitProb.insert(entry_fitprob_full);
            m_fullFitResultMX.insert(entry_mX_full);
            m_fullFitResultMY.insert(entry_mY_full);
            m_fullFitResultMh.insert(entry_mh_full);
            m_fullFitPullB1.insert(entry_pullb1_full);
            m_fullFitPullB2.insert(entry_pullb2_full);
            m_fullFitPullBalance.insert(entry_pullbalance_full);
            m_fullFitPullBalanceX.insert(entry_pullbalance_fullX);
            m_fullFitPullBalanceY.insert(entry_pullbalance_fullY);
            m_fullFitConvergence.insert(entry_convergence_full);

            if (m_chi2 < m_bestChi2FullFit) {
                m_bestChi2FullFit = m_chi2;
                m_bestMXFullFit = m_fitted_mX;
                m_bestMYFullFit = m_fitted_mY;
                m_bestMhFullFit = m_fitted_mh;
                m_bestHypoFullFit = hypo_full;
            }
        }
    }
}

void YHKinFitMaster::Fit(int mh, int mY)
{
    //  ----------  for PSfit ----------
    const int np = 2;
    double a[np];
    double a_start[np];
    double a_limit[np][2];
    double a_precision[np];
    double daN[np];
    double h[np];
    double chi2_iter[1], a_Memory[np][5], g[np], H[np * np], H_inv[np * np];
    bool noNewtonShifts = false;

    int iter = 0;             //  number of iterations
    int method = 1;           //  initial fit method, see PSfitter()
    int mode = 1;             //  mode =1 for start of a new fit by PSfitter()
    //  --------------------------------
    
    // double bjet1_reso = GetBjetResolution(m_bjet1.Eta(), m_bjet1.Et());
    // double bjet2_reso = GetBjetResolution(m_bjet2.Eta(), m_bjet2.Et());
    double bjet1_reso = CalcBjetResolution(m_bjet1, m_bjet_reso1);
    double bjet2_reso = CalcBjetResolution(m_bjet2, m_bjet_reso2);

    double bjet1UpperLimit = m_bjet1.E() + 5.0 * bjet1_reso;
    double bjet1LowerLimit = m_bjet1.E() - 5.0 * bjet1_reso;
    double bjet2UpperLimit = m_bjet2.E() + 5.0 * bjet2_reso;
    double bjet2LowerLimit = m_bjet2.E() - 5.0 * bjet2_reso;

    double tau1LowerLimit = m_tauvis1.E();
    double tau2LowerLimit = m_tauvis2.E();
    
    if (m_Ytautau) {
        mHttHypo = mY;
        mHbbHypo = 125.0;
    }
    else {
        mHttHypo = 125.0;
        mHbbHypo = mY;
    }
    
    double mtau  = 1.777;

    //Calculate Recoil CovMatrix
    TMatrixD Cov_MET = m_MET_COV;
    TMatrixD Cov_b1 = CalcCov(m_bjet1, bjet1_reso);
    TMatrixD Cov_b2 = CalcCov(m_bjet2, bjet2_reso);

    //Tau resolution assumed to be exact
    TMatrixD tau_cov(2, 2);
    tau_cov[0][0] = 0.;
    tau_cov[1][1] = 0.;
    tau_cov[0][1] = 0.;
    tau_cov[1][0] = 0.;
    TMatrixD Cov_tauvis1 = tau_cov;
    TMatrixD Cov_tauvis2 = tau_cov;

    m_covRecoil = Cov_MET - (Cov_b1 + Cov_b2 + Cov_tauvis1 + Cov_tauvis2);
    
    TMatrixDEigen eigenmatrix = TMatrixDEigen(m_covRecoil);
    if (eigenmatrix.GetEigenValues()(0,0)<0 || eigenmatrix.GetEigenValues()(1,1)<0) {
        m_covRecoil[0][0] = 100.;
        m_covRecoil[1][1] = 100.;
        m_covRecoil[1][0] = 0.;
        m_covRecoil[0][1] = 0.;
    }

    if (m_tauvis1.E() < mtau) {
        m_tauvis1.SetM(mtau);
    }
    if (m_tauvis2.E() < mtau) {
        m_tauvis2.SetM(mtau);
    }

    ROOT::Math::PtEtaPhiMVector Htt_vis = m_tauvis1 + m_tauvis2;
    ROOT::Math::PtEtaPhiMVector Hbb = m_bjet1 + m_bjet2;

    //Compute upper limit of E(tau1) by having set E(tau2)=E(tau2)_min and compute E(tau1)
    double tau1UpperLimit = ConstrainEnergy(Htt_vis, m_tauvis2, m_tauvis1, mHttHypo);
    double tau2UpperLimit = ConstrainEnergy(Htt_vis, m_tauvis1, m_tauvis2, mHttHypo);

    if (tau1UpperLimit < m_tauvis1.E()) {
        m_convergence=-1;
        m_chi2=9999;
        m_chi2_b1=9999;
        m_chi2_b2=9999;
        m_chi2_balance=9999;
        m_fitted_mX=-1;
        return;
    }
    if (tau2UpperLimit < m_tauvis2.E()) {
        m_convergence=-1;
        m_chi2=9999;
        m_chi2_b1=9999;
        m_chi2_b2=9999;
        m_chi2_balance=9999;
        m_fitted_mX=-1;
        return;
    }

    if (tau1UpperLimit > 13000.) {
        tau1UpperLimit = 13000.;
    }
    if (tau2UpperLimit > 13000.) {
        tau2UpperLimit = 13000.;
    }
    
    // fill initial tau fit parameters
    a_start[0] = m_bjet1.E();          // energy of first b-jet
    a_precision[0] = a_start[0]*0.002;             // precision for fit
    a_start[1] = tau1LowerLimit;       // energy of first tau
    a_precision[1] = 0.1;                          // precision for fit
    
    // fill initial step width
    h[0] = 0.5 * bjet1_reso;
    h[1] = 0.1 * tau1LowerLimit;                

    daN[0] = 1.0;                  // initial search direction in E_b-E_tau diagonal
    daN[1] = 1.0; 

    // fit range
    if (bjet1LowerLimit<0.01) {
        a_limit[0][0] = 0.01;
    } 
    else {
        a_limit[0][0] = bjet1LowerLimit; 
    }
    a_limit[0][1] = bjet1UpperLimit;

    // adjusting lower and upper energy limits for bjet 1
    ROOT::Math::PtEtaPhiEVector tmp_bjet_2 = ROOT::Math::PtEtaPhiEVector(m_bjet2.Pt(), m_bjet2.Eta(), m_bjet2.Phi(), m_bjet2.E());
    tmp_bjet_2.SetPxPyPzE(tmp_bjet_2.Px()*bjet2LowerLimit/tmp_bjet_2.E(), tmp_bjet_2.Py()*bjet2LowerLimit/tmp_bjet_2.E(), tmp_bjet_2.Pz()*bjet2LowerLimit/tmp_bjet_2.E(), bjet2LowerLimit); //Set b2 to lower b2 limit

    double max_E_b1 = ConstrainEnergy(Hbb, ROOT::Math::PtEtaPhiMVector(tmp_bjet_2.Pt(), tmp_bjet_2.Eta(), tmp_bjet_2.Phi(), tmp_bjet_2.M()), m_bjet1, mHbbHypo);       //Calculate b1 with b2 fixed
    if(a_limit[0][1] > max_E_b1){
        a_limit[0][1] = max_E_b1;                                   //Modify b1 limit if b2 limit is tighter
    }

    tmp_bjet_2 = ROOT::Math::PtEtaPhiEVector(m_bjet2.Pt(), m_bjet2.Eta(), m_bjet2.Phi(), m_bjet2.E());
    tmp_bjet_2.SetPxPyPzE(tmp_bjet_2.Px()*bjet2UpperLimit/tmp_bjet_2.E(), tmp_bjet_2.Py()*bjet2UpperLimit/tmp_bjet_2.E(), tmp_bjet_2.Pz()*bjet2UpperLimit/tmp_bjet_2.E(), bjet2UpperLimit); //Set b2 to lower b2 limit

    double min_E_b1 = ConstrainEnergy(Hbb, ROOT::Math::PtEtaPhiMVector(tmp_bjet_2.Pt(), tmp_bjet_2.Eta(), tmp_bjet_2.Phi(), tmp_bjet_2.M()), m_bjet1, mHbbHypo);       //Calculate b1 with b2 fixed
    if(a_limit[0][0] < min_E_b1){
        a_limit[0][0] = min_E_b1;                                   //Modify b1 limit if b2 limit is tighter
    }

    //Fill initial bjet fit parameters
    if((a_limit[0][1] - a_limit[0][0]) > 0.5*bjet1_reso) {
        a_start[0] = a_limit[0][0] + (a_limit[0][1] - a_limit[0][0]) / 2.0;         // energy of first b-jet
    }
    else {
        m_convergence=-2;
        m_chi2=9999;
        m_chi2_b1=9999;
        m_chi2_b2=9999;
        m_chi2_balance=9999;
        m_fitted_mX=-1;
        return;
    }

    a_limit[1][0] = tau1LowerLimit;              // tau: minimum is visible tau1 energy
    a_limit[1][1] = tau1UpperLimit;              //      maximum as computed above

    // tau: check initial values against fit range
    if (a_start[1] - h[1] < a_limit[1][0]) {
        a_start[1] = a_limit[1][0] + h[1];
    }
    else if (a_start[1] + h[1] > a_limit[1][1]) {
        a_start[1] = a_limit[1][1] - h[1];
    }

    for (int ip = 0; ip < np; ip++) {
        a[ip] = a_start[ip];

        a_Memory[ip][0] = -999.0;
        a_Memory[ip][1] = -995.0;
        a_Memory[ip][2] = -990.0;
        a_Memory[ip][3] = -985.0;
        a_Memory[ip][4] = -980.0;
    }

    static const int nloopmax = 1000;

    ROOT::Math::PtEtaPhiEVector bjet_1_fit = ROOT::Math::PtEtaPhiEVector(m_bjet1.Pt(), m_bjet1.Eta(), m_bjet1.Phi(), m_bjet1.E());
    ROOT::Math::PtEtaPhiEVector bjet_2_fit = ROOT::Math::PtEtaPhiEVector(m_bjet2.Pt(), m_bjet2.Eta(), m_bjet2.Phi(), m_bjet2.E());
    ROOT::Math::PtEtaPhiEVector tau_1_fit = ROOT::Math::PtEtaPhiEVector(m_tauvis1.Pt(), m_tauvis1.Eta(), m_tauvis1.Phi(), m_tauvis1.E());
    ROOT::Math::PtEtaPhiEVector tau_2_fit = ROOT::Math::PtEtaPhiEVector(m_tauvis2.Pt(), m_tauvis2.Eta(), m_tauvis2.Phi(), m_tauvis2.E());
    p4_X_fit = bjet_1_fit + bjet_2_fit + tau_1_fit + tau_2_fit;
    double new_E_b2, new_E_tau2;
    
    for (int iloop = 0; iloop < nloopmax; iloop++) { // FIT loop
        // Logger::get("YHKinFit")->debug("kinfit_convergence: {}, energy b {}, energy tau {}", m_convergence, new_E_b2, new_E_tau2); 
        bjet_1_fit.SetPxPyPzE(bjet_1_fit.Px()*a[0]/bjet_1_fit.E(), bjet_1_fit.Py()*a[0]/bjet_1_fit.E(), bjet_1_fit.Pz()*a[0]/bjet_1_fit.E(), a[0]);

        tau_1_fit.SetPxPyPzE(tau_1_fit.Px()*a[1]/tau_1_fit.E(), tau_1_fit.Py()*a[1]/tau_1_fit.E(), tau_1_fit.Pz()*a[1]/tau_1_fit.E(), a[1]);
        
        new_E_b2 = ConstrainEnergy(Hbb, ROOT::Math::PtEtaPhiMVector(bjet_1_fit.Pt(), bjet_1_fit.Eta(), bjet_1_fit.Phi(), bjet_1_fit.M()), m_bjet2, mHbbHypo);
        bjet_2_fit.SetPxPyPzE(bjet_2_fit.Px()*new_E_b2/bjet_2_fit.E(), bjet_2_fit.Py()*new_E_b2/bjet_2_fit.E(), bjet_2_fit.Pz()*new_E_b2/bjet_2_fit.E(), new_E_b2);
        
        new_E_tau2 = ConstrainEnergy(Htt_vis, ROOT::Math::PtEtaPhiMVector(tau_1_fit.Pt(), tau_1_fit.Eta(), tau_1_fit.Phi(), tau_1_fit.M()), m_tauvis2, mHttHypo);
        if (new_E_tau2 > tau2UpperLimit) {
            new_E_tau2 = tau2UpperLimit;
        }
        tau_2_fit.SetPxPyPzE(tau_2_fit.Px()*new_E_tau2/tau_2_fit.E(), tau_2_fit.Py()*new_E_tau2/tau_2_fit.E(), tau_2_fit.Pz()*new_E_tau2/tau_2_fit.E(), new_E_tau2);

        m_chi2_b1 = Chi2_V4(ROOT::Math::PtEtaPhiEVector(m_bjet1.Pt(), m_bjet1.Eta(), m_bjet1.Phi(), m_bjet1.E()), bjet_1_fit, m_bjet_reso1);
        m_chi2_b2 = Chi2_V4(ROOT::Math::PtEtaPhiEVector(m_bjet2.Pt(), m_bjet2.Eta(), m_bjet2.Phi(), m_bjet2.E()), bjet_2_fit, m_bjet_reso2);

        if (m_Ytautau) {
            p4_Y_fit = tau_1_fit + tau_2_fit;
            p4_h_fit = bjet_1_fit + bjet_2_fit;
        }
        else {
            p4_Y_fit = bjet_1_fit + bjet_2_fit;
            p4_h_fit = tau_1_fit + tau_2_fit;
        }
        
        p4_X_fit = bjet_1_fit + bjet_2_fit + tau_1_fit + tau_2_fit;
        m_chi2_balance = Chi2_Balance(p4_X_fit);
        m_chi2 = m_chi2_b1 + m_chi2_b2 + m_chi2_balance; // chi2 calculation

        if (m_convergence != 0) {
            break;
        }
        m_convergence = PSFit::PSfitter(iloop, iter, method, mode, noNewtonShifts,
                            np, a, a_start, a_limit, a_precision,
                            daN, h, a_Memory, m_chi2, chi2_iter, g, H, H_inv);
    }

    m_fitted_mX = p4_X_fit.M();
    m_fitted_mY = p4_Y_fit.M();
    m_fitted_mh = p4_h_fit.M();

    if(m_convergence != 0 && m_convergence != 5) {
        if(a[0] < (a_limit[0][0] + 2*a_precision[0]) ) {
            // std::cout << "Convergence at lower bjet limit!" << std::endl;
            m_convergence = 3;
        }
        if(a[0] > (a_limit[0][1] - 2*a_precision[0]) ) {
            // std::cout << "Convergence at upper bjet limit!" << std::endl;
            m_convergence = 3;
        }
        if(a[1] < (a_limit[1][0] + 2*a_precision[1]) ) {
            if(m_convergence == 3) {
                m_convergence = 5;
            }  
            else{
                // std::cout << "Convergence at lower tau limit!" << std::endl;
                m_convergence = 4;
            }
        }
        if(a[1] > (a_limit[1][1] - 2*a_precision[1]) ) {
            if(m_convergence == 3) {
                m_convergence = 5;
            }
            else{
                // std::cout << "Convergence at upper tau limit!" << std::endl;
                m_convergence = 4;
            }
        }
    }
}

void YHKinFitMaster::addMhHypothesis(std::vector<int> v)
{
    m_mh = v;
}

void YHKinFitMaster::addMYHypothesis(std::vector<int> v)
{
    m_mY = v;
}

TMatrixD YHKinFitMaster::CalcCov(ROOT::Math::PtEtaPhiMVector p4, double dE)
{
    //NOTE: only dE is used for Cov calculation!
    double dp = p4.E() / p4.P() * dE; // error propagation p=sqrt(e^2-m^2)
    double dpt = sin(p4.Theta()) * dp;

    TMatrixD cov(2, 2);
    cov[0][0] = pow(cos(p4.Phi()) * dpt, 2);
    cov[1][1] = pow(sin(p4.Phi()) * dpt, 2);
    cov[0][1] = sin(p4.Phi()) * cos(p4.Phi()) *dpt*dpt;
    cov[1][0] = sin(p4.Phi()) * cos(p4.Phi()) *dpt*dpt;

    return cov;
}

double YHKinFitMaster::CalcBjetResolution(ROOT::Math::PtEtaPhiMVector p4, double res){
    double pt_res = p4.Pt()*res;
    double dE = pt_res * p4.P() / sin(p4.Theta()) / p4.E();
    return dE;
}

double YHKinFitMaster::GetBjetResolution(double eta, double et){
    double det=0;
    double de=10;

    if(0.000<=abs(eta) && abs(eta)<0.087){
        det = et * (sqrt(0.0686*0.0686 + (1.03/sqrt(et))*(1.03/sqrt(et)) + (1.68/et)*(1.68/et)));
        de = 1.0/sin(2 * atan(exp(-(0.000+0.087)/2))) * det;
    }

    if(0.087<=abs(eta) && abs(eta)<0.174){
        det = et * (sqrt(0.0737*0.0737 + (1.01/sqrt(et))*(1.01/sqrt(et)) + (1.74/et)*(1.74/et)));
        de = 1.0/sin(2 * atan(exp(-(0.087+0.174)/2))) * det;
    }

    if(0.174<=abs(eta) && abs(eta)<0.261){
        det = et * (sqrt(0.0657*0.0657 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (5.16e-06/et)*(5.16e-06/et)));
        de = 1.0/sin(2 * atan(exp(-(0.174+0.261)/2))) * det;
    }

    if(0.261<=abs(eta) && abs(eta)<0.348){
        det = et * (sqrt(0.062*0.062 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (0.000134/et)*(0.000134/et)));
        de = 1.0/sin(2 * atan(exp(-(0.261+0.348)/2))) * det;
    }

    if(0.348<=abs(eta) && abs(eta)<0.435){
        det = et * (sqrt(0.0605*0.0605 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (1.84e-07/et)*(1.84e-07/et)));
        de = 1.0/sin(2 * atan(exp(-(0.348+0.435)/2))) * det;
    }

    if(0.435<=abs(eta) && abs(eta)<0.522){
        det = et * (sqrt(0.059*0.059 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (9.06e-09/et)*(9.06e-09/et)));
        de = 1.0/sin(2 * atan(exp(-(0.435+0.522)/2))) * det;
    }

    if(0.522<=abs(eta) && abs(eta)<0.609){
        det = et * (sqrt(0.0577*0.0577 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (5.46e-06/et)*(5.46e-06/et)));
        de = 1.0/sin(2 * atan(exp(-(0.522+0.609)/2))) * det;
    }

    if(0.609<=abs(eta) && abs(eta)<0.696){
        det = et * (sqrt(0.0525*0.0525 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (4.05e-05/et)*(4.05e-05/et)));
        de = 1.0/sin(2 * atan(exp(-(0.609+0.696)/2))) * det;
    }

    if(0.696<=abs(eta) && abs(eta)<0.783){
        det = et * (sqrt(0.0582*0.0582 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.17e-05/et)*(1.17e-05/et)));
        de = 1.0/sin(2 * atan(exp(-(0.696+0.783)/2))) * det;
    }

    if(0.783<=abs(eta) && abs(eta)<0.870){
        det = et * (sqrt(0.0649*0.0649 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (7.85e-06/et)*(7.85e-06/et)));
        de = 1.0/sin(2 * atan(exp(-(0.783+0.870)/2))) * det;
    }

    if(0.870<=abs(eta) && abs(eta)<0.957){
        det = et * (sqrt(0.0654*0.0654 + (1.1/sqrt(et))*(1.1/sqrt(et)) + (1.09e-07/et)*(1.09e-07/et)));
        de = 1.0/sin(2 * atan(exp(-(0.870+0.957)/2))) * det;
    }

    if(0.957<=abs(eta) && abs(eta)<1.044){
        det = et * (sqrt(0.0669*0.0669 + (1.11/sqrt(et))*(1.11/sqrt(et)) + (1.87e-06/et)*(1.87e-06/et)));
        de = 1.0/sin(2 * atan(exp(-(0.957+1.044)/2))) * det;
    }

    if(1.044<=abs(eta) && abs(eta)<1.131){
        det = et * (sqrt(0.0643*0.0643 + (1.15/sqrt(et))*(1.15/sqrt(et)) + (2.76e-05/et)*(2.76e-05/et)));
        de = 1.0/sin(2 * atan(exp(-(1.044+1.131)/2))) * det;
    }

    if(1.131<=abs(eta) && abs(eta)<1.218){
        det = et * (sqrt(0.0645*0.0645 + (1.16/sqrt(et))*(1.16/sqrt(et)) + (1.04e-06/et)*(1.04e-06/et)));
        de = 1.0/sin(2 * atan(exp(-(1.131+1.218)/2))) * det;
    }

    if(1.218<=abs(eta) && abs(eta)<1.305){
        det = et * (sqrt(0.0637*0.0637 + (1.19/sqrt(et))*(1.19/sqrt(et)) + (1.08e-07/et)*(1.08e-07/et)));
        de = 1.0/sin(2 * atan(exp(-(1.218+1.305)/2))) * det;
    }

    if(1.305<=abs(eta) && abs(eta)<1.392){
        det = et * (sqrt(0.0695*0.0695 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5.75e-06/et)*(5.75e-06/et)));
        de = 1.0/sin(2 * atan(exp(-(1.305+1.392)/2))) * det;
    }

    if(1.392<=abs(eta) && abs(eta)<1.479){
        det = et * (sqrt(0.0748*0.0748 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (5.15e-08/et)*(5.15e-08/et)));
        de = 1.0/sin(2 * atan(exp(-(1.392+1.479)/2))) * det;
    }

    if(1.479<=abs(eta) && abs(eta)<1.566){
        det = et * (sqrt(0.0624*0.0624 + (1.23/sqrt(et))*(1.23/sqrt(et)) + (2.28e-05/et)*(2.28e-05/et)));
        de = 1.0/sin(2 * atan(exp(-(1.479+1.566)/2))) * det;
    }

    if(1.566<=abs(eta) && abs(eta)<1.653){
        det = et * (sqrt(0.0283*0.0283 + (1.25/sqrt(et))*(1.25/sqrt(et)) + (4.79e-07/et)*(4.79e-07/et)));
        de = 1.0/sin(2 * atan(exp(-(1.566+1.653)/2))) * det;
    }

    if(1.653<=abs(eta) && abs(eta)<1.740){
        det = et * (sqrt(0.0316*0.0316 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5e-05/et)*(5e-05/et)));
        de = 1.0/sin(2 * atan(exp(-(1.653+1.740)/2))) * det;
    }

    if(1.740<=abs(eta) && abs(eta)<1.830){
        det = et * (sqrt(2.29e-07*2.29e-07 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (1.71e-05/et)*(1.71e-05/et)));
        de = 1.0/sin(2 * atan(exp(-(1.740+1.830)/2))) * det;
    }

    if(1.830<=abs(eta) && abs(eta)<1.930){
        det = et * (sqrt(5.18e-09*5.18e-09 + (1.14/sqrt(et))*(1.14/sqrt(et)) + (1.7/et)*(1.7/et)));
        de = 1.0/sin(2 * atan(exp(-(1.830+1.930)/2))) * det;
    }

    if(1.930<=abs(eta) && abs(eta)<2.043){
        det = et * (sqrt(2.17e-07*2.17e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (2.08/et)*(2.08/et)));
        de = 1.0/sin(2 * atan(exp(-(1.930+2.043)/2))) * det;
    }

    if(2.043<=abs(eta) && abs(eta)<2.172){
        det = et * (sqrt(3.65e-07*3.65e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.63/et)*(1.63/et)));
        de = 1.0/sin(2 * atan(exp(-(2.043+2.172)/2))) * det;
    }

    if(2.172<=abs(eta) && abs(eta)<2.322){
        det = et * (sqrt(2.02e-07*2.02e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.68/et)*(1.68/et)));
        de = 1.0/sin(2 * atan(exp(-(2.172+2.322)/2))) * det;
    }

    if(2.322<=abs(eta) && abs(eta)<2.500){
        det = et * (sqrt(5.27e-07*5.27e-07 + (1.12/sqrt(et))*(1.12/sqrt(et)) + (1.78/et)*(1.78/et)));
        de = 1.0/sin(2 * atan(exp(-(2.322+2.500)/2))) * det;
    }

    return de;
}

double YHKinFitMaster::ConstrainEnergy(ROOT::Math::PtEtaPhiMVector p4_mother, ROOT::Math::PtEtaPhiMVector p4_1, ROOT::Math::PtEtaPhiMVector p4_2, int mHypo)
{
    ROOT::Math::PtEtaPhiEVector new_p4_1 = ROOT::Math::PtEtaPhiEVector(p4_1.Pt(), p4_1.Eta(), p4_1.Phi(), p4_1.E());
    ROOT::Math::PtEtaPhiEVector new_p4_2 = ROOT::Math::PtEtaPhiEVector(p4_2.Pt(), p4_2.Eta(), p4_2.Phi(), p4_2.E());
    ROOT::Math::PtEtaPhiEVector new_p4_mother = ROOT::Math::PtEtaPhiEVector(p4_mother.Pt(), p4_mother.Eta(), p4_mother.Phi(), p4_mother.E());
    double M_truth = mHypo;
    double M_reco = new_p4_mother.M();
    double E_1 = new_p4_1.E();
    double M_1 = new_p4_1.M();
    double E_2 = new_p4_2.E();
    double M_2 = new_p4_2.M();

    double beta_2 = -10.;
    double gamma2_2 = -10.;
    double c = -10.;
    double E_part1 = -10.;
    double E_part2 = -10.;
    double new_E_2 = -10.;

    int loopCount = 0;

    while (abs(M_reco - M_truth) > 0.0001) {
        loopCount++;
     
        if (loopCount>=10) {
            m_convergence=-3;
            m_chi2=9999;
            m_chi2_b1=9999;
            m_chi2_b2=9999;
            m_chi2_balance=9999;
            m_fitted_mX=-1;
            return new_E_2;
        }
        
        if ( (M_2 < (1.e-3*E_2)) || (M_2 < 0.) ) { // massless case; if mass is negative or much smaller than the energy
            new_E_2 = E_2 * (M_truth / M_reco) * (M_truth / M_reco);
            return new_E_2;
        }
        else {
            beta_2 = sqrt(E_2 * E_2 - M_2 * M_2) / E_2;
            gamma2_2 = 1 / (1 - beta_2 * beta_2);

            c = (M_reco * M_reco - M_1 * M_1 - M_2 * M_2) / (2. * E_1 * E_2);
            E_part1 = E_1 * c * gamma2_2;
            E_part2 = (M_truth * M_truth - M_1 * M_1) / (gamma2_2 * E_1 * E_1 * c * c);
            
            new_E_2 = E_part1 * (-1. + sqrt(1. + E_part2));
        }
       
        new_p4_2.SetPxPyPzE(new_p4_2.Px()*new_E_2/E_2, new_p4_2.Py()*new_E_2/E_2, new_p4_2.Pz()*new_E_2/E_2, new_E_2);
        
        new_p4_mother = new_p4_2 + new_p4_1;

        M_reco = new_p4_mother.M();
        E_2 = new_p4_2.E();
        M_2 = new_p4_2.M();
    }

    return new_E_2;
}

double YHKinFitMaster::Chi2_V4(ROOT::Math::PtEtaPhiEVector p4_reco, ROOT::Math::PtEtaPhiEVector p4_fit, double res)
{
    double chi2_E = 0;
    // double dE_fit = GetBjetResolution(p4_fit.Eta(), p4_fit.Et());
    double dE_fit = CalcBjetResolution(ROOT::Math::PtEtaPhiMVector(p4_fit.Pt(), p4_fit.Eta(), p4_fit.Phi(), p4_fit.M()), res);

    if (dE_fit > 0.) {
        chi2_E = ((p4_reco.E() - p4_fit.E()) / dE_fit) * ((p4_reco.E() - p4_fit.E()) / dE_fit);
    }

    return chi2_E;
}

double YHKinFitMaster::Chi2_Balance(ROOT::Math::PtEtaPhiEVector p4_X_fit)
{
    //reco objects
    ROOT::Math::PtEtaPhiMVector p4_X_reco = m_MET + m_tauvis1 + m_tauvis2 + m_bjet1 + m_bjet2;

    double res_px, res_py;
    res_px = p4_X_fit.Px() - p4_X_reco.Px();    // residuum in Pt_H
    res_py = p4_X_fit.Py() - p4_X_reco.Py();    // residuum in Pt_H
    
    double Vxx = m_covRecoil[0][0];
    double Vyy = m_covRecoil[1][1];
    double Vxy = m_covRecoil[0][1];

    double det, Vxx_inv, Vyy_inv, Vxy_inv;
    det = Vxx * Vyy - Vxy * Vxy;
    Vxx_inv = Vyy / det;
    Vyy_inv = Vxx / det;
    Vxy_inv = -Vxy / det;

    TMatrixD V_inv(2, 2);
    V_inv[0][0] = Vxx_inv;
    V_inv[1][0] = Vxy_inv;
    V_inv[0][1] = Vxy_inv;
    V_inv[1][1] = Vyy_inv;

    double chi2 = res_px*(V_inv[0][0]*res_px + V_inv[0][1]*res_py) + res_py*(V_inv[1][0]*res_px + V_inv[1][1]*res_py); // chi2 = res_transponiert * Vinv * res    

    return chi2;
}

double YHKinFitMaster::GetPullBalanceX(ROOT::Math::PtEtaPhiEVector p4_X_fit)
{
    ROOT::Math::PtEtaPhiMVector p4_X_reco = m_MET + m_tauvis1 + m_tauvis2 + m_bjet1 + m_bjet2;
    double dE_fit = sqrt(m_covRecoil[0][0]);
    double pull = (p4_X_fit.Px()- p4_X_reco.Px()) / dE_fit;
    return pull;
}

double YHKinFitMaster::GetPullBalanceY(ROOT::Math::PtEtaPhiEVector p4_X_fit)
{
    ROOT::Math::PtEtaPhiMVector p4_X_reco = m_MET + m_tauvis1 + m_tauvis2 + m_bjet1 + m_bjet2;
    double dE_fit = sqrt(m_covRecoil[1][1]);
    double pull = (p4_X_fit.Py()- p4_X_reco.Py()) / dE_fit;
    return pull;
}