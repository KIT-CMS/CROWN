#ifndef GUARDGGH_HTXS_H
#define GUARDGGH_HTXS_H

#include <vector>

typedef std::vector<double> NumV;

NumV qcd_ggF_uncert_wg1(int Njets30, double pTH, int STXS);
NumV qcd_ggF_uncert_stxs(int Njets30, double pTH, int STXS);
NumV qcd_ggF_uncert_2017(int Njets30, double pTH, int STXS);
NumV qcd_ggF_uncert_jve(int Njets30, double pTH, int STXS);

NumV qcd_ggF_uncertSF_wg1(int Njets30, double pTH, int STXS_Stage1,
                          double Nsigma = 1.0);
NumV qcd_ggF_uncertSF_stxs(int Njets30, double pTH, int STXS_Stage1,
                           double Nsigma = 1.0);
NumV qcd_ggF_uncertSF_2017(int Njets30, double pTH, int STXS_Stage1,
                           double Nsigma = 1.0);
NumV qcd_ggF_uncertSF_jve(int Njets30, double pTH, int STXS_Stage1,
                          double Nsigma = 1.0);

// Cross sections of ggF with =0, =1, and >=2 jets
static double g_sig0 = 30.117, g_sig1 = 12.928, g_sig_ge2 = 5.475,
              g_sig_ge1 = g_sig1 + g_sig_ge2, g_sig_tot = g_sig0 + g_sig_ge1,
              g_sig_vbfTopo = 0.630, g_sig_ge2noVBF = g_sig_ge2 - g_sig_vbfTopo,
              g_sig_ge1noVBF = g_sig_ge1 - g_sig_vbfTopo;

NumV blptw(int Njets30);
double vbf_2j(int STXS);
double vbf_3j(int STXS);
double interpol(double x, double x1, double y1, double x2, double y2);
double qm_t(double pT);
double pT120(double pT, int Njets30);
double pT60(double pT, int Njets30);
NumV jetBinUnc(int Njets30, int STXS);
NumV unc2sf(const NumV &unc, double Nsigma);

#endif /* GUARDGGH_HTXS_H */