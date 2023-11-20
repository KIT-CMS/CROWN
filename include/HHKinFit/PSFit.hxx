#ifndef PSFIT_H_
#define PSFIT_H_

#include <Rtypes.h>

class PSFit {
    public:
        PSFit();
        static int PSfitter(int iloop, int &iter, int &method, int &mode,  
                bool &noNewtonShifts, int np, double a[], double a_start[], double a_limit[][2], 
                double a_precision[], double daN[], double h[], double a_Memory[][5],
                double chi2, double chi2_iter[], double g[], double H[], double H_inv[]);

        static double PSVnorm(double x[], int n);

        static void PSLineLimit(int np, double a_start[], double daN[], double a_limit[][2], double x_limit[]);

        static double PSLineSearch(int &mode, double hh, double x_limit[], 
                    double eps_x, double eps_f, double x[4], double f[], double chi2);

        static void PSNewtonLimitShift(int sign, int np, double a[], double a_limit[][2], double a_precision[],
                    double daN[], double h[], double g[], double H[]);

        static int PSderivative(int icall, int np, double a[], double h[], double chi2, double chi2_iter[], double g[], double H[]);

        static int PSderivative1(int icall, double a[], double h[], double chi2, double g[], double H[]);

        static double PSNewtonAnalyzer(int np, double a[], double a_limit[][2], double a_precision[],
                        double daN[], double h[], double g[], double H[], double H_inv[], double chi2, bool noNewtonShifts);

        static double PSMinverse(double H[], double H_inv[], int p);

        static double PSMCholesky(double M[], double R[],  int n);

        static double PSMRTrianInvert2(double R[], int n);

        static double PSMmultiplyMRRT(double A[], int n1, int n2);
};

#endif