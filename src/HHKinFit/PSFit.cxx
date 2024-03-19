#include "../../include/HHKinFit/PSFit.hxx"

#include <TMarker.h>
#include <TPolyMarker.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

// Collection of mathematical tools
//-----------------------------
//  PSfitter            Fit tool
//  PSderivative        Tool for numerical calc. of derivative and Hesse matrix
//  PSderivative1       1-dim case of PSderivative,  see documentation there
//  PSMinverse       invert a symmetric matrix
//  PSMCholesky      Cholesky decomposition of  M = n*n  symmetric matrix
//  PSMRTrianInvert2 invert a triangular R - matrix
//  PSMmultiplyMRRT  multiply triangular matrix R with transpose,

int PSFit::PSfitter(int iloop, int &iter, int &method, int &mode,
                    bool &noNewtonShifts, int np, double a[], double a_start[],
                    double a_limit[][2], double a_precision[], double daN[],
                    double h[], double a_Memory[][5], double chi2,
                    double chi2_iter[], double g[], double H[],
                    double H_inv[]) {
    // generic fitter using Newton method and Line Search within parameter
    // limits iter:     set iter=0 at start of a new fit,
    //           iter<0 during construction of numerical derivatives
    //           iter>0 after each iteration
    // method:   preferred initial fit method:
    //           =1 line search in direction of daN[],
    //           =2 Newton in all dimensions
    // mode:     mode within method;
    // np:         number of fit parameters
    // a[np]:      fit parameters, fill start values at beginning when iloop=0,
    //             contains actual fit parameters for chi2 calculation,
    //             will contain final result
    // a_start[np]: start values for a[], will not be changed
    // a_limit[np]: upper and lower limits for a[]
    // a_precision[np]:  precision requested for a[], also minimum value for h[]
    //             and minimum distance of a[] to alimit[] (needed to calculate
    //             derivatives)
    // daN[]:      step vector from iteration to iteration,
    //             may contain initial search direction if method= line search
    // h[]:        for Newton method: initial step width for calculation of num.
    // derivatives chi2:       must contain function value at input value of a
    // to be minimized chi2_iter[]  Function (chi2) with nominal a[], i.e. no
    // shift of a[] g[np]       derivative vector  g[np] H[np*np]    Hesse
    // matrix,   g[] and H[] are also used as
    //                             intermediate storage of chi2
    // H_inv[np*np] Inverse of Hesse matrix

    static int icall_Newton, iter_Memory;
    static double chi2_Memory;
    static double x[4], f[4];
    static double xx, x_limit[2];
    static double x_h, daN_abs;
    static double eps_x = 0.1, eps_f = 0.1;
    int convergence;

    int itemp, ready;
    double temp, d;

    convergence = 0;
    if (iloop == 0) { // start of a new fit
        iter_Memory = 0;
        iter = 0;
        chi2_Memory = 123456.;
        for (int ip = 0; ip < np; ip++) { // check limits
            a_limit[ip][0] = std::max(a_limit[ip][0], -pow(10.0, 10.0));
            a_limit[ip][1] = std::min(a_limit[ip][1], pow(10.0, 10.0));
            if (a_limit[ip][0] >= a_limit[ip][1]) {
                a_limit[ip][0] = -pow(10.0, 10.0);
                a_limit[ip][1] = pow(10.0, 10.0);
            }
            if (a_start[ip] < a_limit[ip][0] ||
                a_start[ip] > a_limit[ip][1]) { // NOT WORKING astart has no
                                                // effect once the loop started!
                a_start[ip] =
                    0.5 * (a_limit[ip][0] +
                           a_limit[ip][1]); // set to mean in between limits
                return convergence;
            }
        }
    }

    if (method == 3) {
        if (std::abs(chi2_Memory - chi2) < eps_f) { // test convergence
            convergence = 1;
            iter_Memory = iter_Memory + 1;
            iter = iter_Memory;
            chi2_Memory = chi2;
        } else {
            // std::cout << "No Convergence! Back to LineSearch!" << std::endl;
            method = 1; // switch to line search
            mode = 1;
            for (int ip = 0; ip < np; ip++) {
                a_start[ip] = a[ip];
            }
        }
    }

    iter = -1; // only positive if a new iteration step (LineSearch, ..) is
               // finished
    if (method == 1) { // ----------------------------- line search
        // std::cout << "In Linesearch!" << std::endl;
        if (mode == 1) {                // start of a new line search
            daN_abs = PSVnorm(daN, np); // check initial search direction
            //      x_h     = PSVnorm(h, np) ;        // initial step width for
            //      line search
            x_h = 1.0 / PSVnorm(daN, np); // initial step width for line search

            if (daN_abs == 0.)
                return 5; // Minimum found at both limits

            else {
                PSLineLimit(np, a_start, daN, a_limit, x_limit);
            } // limits for line search
        }

        double eps_x_LS = eps_x;

        for (int ip = 0; ip < np; ip++) {
            if (a_precision[ip] * 0.5 / std::abs(daN[ip]) < eps_x_LS) {
                eps_x_LS = a_precision[ip] * 0.5 / std::abs(daN[ip]);
            }
        }

        xx = PSLineSearch(mode, x_h, x_limit, eps_x_LS, eps_f, x, f, chi2);
        for (int ip = 0; ip < np; ip++) {
            a[ip] = a_start[ip] + xx * daN[ip];
        }

        if (mode <= 0) { // Minimum found
            // std::cout << "Minimum found between:" << std::endl;
            // std::cout << "Eb1 = " << a_start[0] + x[1] * daN[0] << " ETau1 =
            // " << a_start[1] + x[1] * daN[1] << std::endl; std::cout << "and"
            // << std::endl; std::cout << "Eb1 = " << a_start[0] + x[3] * daN[0]
            // << " ETau1 = " << a_start[1] + x[3] * daN[1] << std::endl;
            // std::cout << "Finished Linesearch! Test Convergence" <<
            // std::endl;

            bool didConverge = true;
            for (int ip = 0; ip < np;
                 ip++) { // check for progress w.r.t. previous iteration
                if (std::abs(a[ip] - a_Memory[ip][0]) > a_precision[ip]) {
                    // std::cout << "No Convergence as Parameter " << np << " is
                    // " << std::abs(a[ip] - a_Memory[ip][0]) << " away from
                    // previous iteration." << std::endl;
                    didConverge = false;
                    break;
                }
            }
            if (didConverge) {
                convergence = 2;
                iter_Memory = iter_Memory + 1;
                iter = iter_Memory;
                return convergence;
            }

            for (int previous_iter = 1; previous_iter <= 4; previous_iter++) {
                bool we_have_been_here_before = true;
                for (int ip = 0; ip < np; ip++) {
                    if (std::abs(a_Memory[ip][previous_iter] - a[ip]) >
                        a_precision[ip] * 0.01) {
                        we_have_been_here_before = false;
                        // std::cout << "We have not been here before!" <<
                        // std::endl;
                        break;
                    }
                }
                if (we_have_been_here_before) {
                    // std::cout << "Pattern recognized! We have been here " <<
                    // previous_iter+1 << " iterations ago! Disable
                    // NewtonShifts!" << std::endl;
                    noNewtonShifts = true;
                }
            }

            for (int ip = 0; ip < np; ip++) {
                a_Memory[ip][4] = a_Memory[ip][3];
                a_Memory[ip][3] = a_Memory[ip][2];
                a_Memory[ip][2] = a_Memory[ip][1];
                a_Memory[ip][1] = a_Memory[ip][0];
                a_Memory[ip][0] = a[ip];
            }

            // std::cout << "No Convergence! Switch to Newton!" << std::endl;
            method = 2;
            icall_Newton = -1; // switch to Newton Method
            iter_Memory = iter_Memory + 1;
            iter = iter_Memory;
            chi2_Memory = chi2;
        }
    }

    if (method == 2) { // Newton method
        if (icall_Newton == -1) {
            PSNewtonLimitShift(1, np, a, a_limit, a_precision, daN, h, g, H);
            icall_Newton = 0;
            return convergence;
        }
        // std::cout << "Newton Method! Calc. derivative!" << std::endl;

        if (np > 1) {
            ready = PSderivative(icall_Newton, np, a, h, chi2, chi2_iter, g, H);
        } else {
            ready = PSderivative1(icall_Newton, a, h, chi2, g, H);
        }

        icall_Newton = icall_Newton + 1;
        d = 99.;

        if (ready == 1) {
            // std::cout << "Newton Method! Got derivative! Get Minimum!" <<
            // std::endl;

            PSNewtonLimitShift(-1, np, a, a_limit, a_precision, daN, h, g, H);
            d = PSNewtonAnalyzer(np, a, a_limit, a_precision, daN, h, g, H,
                                 H_inv, chi2, noNewtonShifts);

            if (d < eps_x) { // test convergence in next run by checking if
                             // delta_Chi2 is small enough
                method = 3;
                // std::cout << "Small Shift! Check for Convergence (small Chi2)
                // in next Loop!" << std::endl;
            } else {
                // std::cout << "No Minimum! Back to Linesearch!" << std::endl;
                method = 1;
                mode = 1; // switch to line search
                for (int ip = 0; ip < np; ip++) {
                    a_start[ip] = a[ip];
                }
            }
        }
    }
    return convergence;
}

double PSFit::PSVnorm(double x[], int n) {
    double xnorm = 0.;
    for (int i = 0; i < n; i++) {
        xnorm = xnorm + pow(x[i], 2.0);
    }
    return sqrt(xnorm);
}

void PSFit::PSLineLimit(int np, double a_start[], double daN[],
                        double a_limit[][2], double x_limit[]) {
    // calculate x_limit for line search along direction daN
    double temp0, temp1, temp;
    x_limit[0] = -pow(10.0, 10.0);
    x_limit[1] = pow(10.0, 10.0);

    for (int ip = 0; ip < np; ip++) { // limits
        if (daN[ip] != 0.) {
            temp0 = (a_limit[ip][0] - a_start[ip]) / daN[ip];
            temp1 = (a_limit[ip][1] - a_start[ip]) / daN[ip];

            if (temp0 > temp1) {
                temp = temp0;
                temp0 = temp1;
                temp1 = temp;
            }

            x_limit[0] = std::max(x_limit[0], temp0);
            x_limit[1] = std::min(x_limit[1], temp1);
        }
    }
}

double PSFit::PSLineSearch(
    int &mode, double hh, double x_limit[], double eps_x, double eps_f,
    double x[4], double f[],
    double chi2) { // 1-dim Line-Search, Method from Blobel textbook p. 252
    static double xt, ft;
    double d31, d32, d21;
    double g, H;
    double tau = 0.618034;
    double tau_fac = tau - tau * tau;
    double close = 0.004;
    double h;
    double temp;
    int i_limit;

    if (mode == 1) {
        // std::cout << "In Linesearch Mode 1!" << std::endl;
        // set up for "+x" search direction
        h = std::min(std::abs(hh),
                     0.25 * (x_limit[1] - x_limit[0])); // limit h to x_limit
        x[1] = 0.;
        f[1] = chi2;
        x[2] = x[1] + h;
        if (x[2] > x_limit[1] - eps_x) {
            x[2] = x[1] - h;
        } // too close to limit: invert direction
        x[3] = 0.;
        f[2] = 0.;
        f[3] = 0.;

        mode = 2;
        return x[2];
    }

    else if (mode == 2) {
        // std::cout << "In Linesearch Mode 2!" << std::endl;
        f[2] = chi2;
        x[3] = x[2] + x[2] - x[1]; // step in same direction as previous step
        if (x[3] > x_limit[1] - eps_x) {
            x[3] = x[1] + x[1] - x[2];
        } // too close to limit: invert direction

        mode = 3;
        return x[3];
    }

    else if (mode == 3) {
        // std::cout << "In Linesearch Mode 3!" << std::endl;
        f[3] = chi2;
        // due to limits, the order x1->x2->x3 might be screwed up: repair this
        // choose f1 to be the worst point
        if (f[1] < f[2]) { // Switch points 1 and 2
            temp = f[2];
            f[2] = f[1];
            f[1] = temp;
            temp = x[2];
            x[2] = x[1];
            x[1] = temp;
        }                  // order f1>f2
        if (f[1] < f[3]) { // Switch points 1 and 3
            temp = f[3];
            f[3] = f[1];
            f[1] = temp;
            temp = x[3];
            x[3] = x[1];
            x[1] = temp;
        } // order f1>f3
        if ((x[1] > x[2] && x[1] < x[3]) ||
            (x[1] < x[2] &&
             x[1] > x[3])) { // Local maximum between x[2] and x[3]! Continue in
                             // direction of steepest descent.
            // std::cout << "Local Maximum!" << std::endl;
            if ((f[1] - f[2]) > (f[1] - f[3])) { // Continue in direction of
                                                 // x[2]
                temp = x[3];
                x[3] = x[2] + x[2] - x[1];
                if ((x[3] > x_limit[1] - eps_x) ||
                    (x[3] < x_limit[0] +
                                eps_x)) { // To close to limit. Other direction.
                    x[3] = temp;
                    x[2] = x[3];
                    f[2] = f[3];
                    x[3] = x[2] + x[2] - x[1];
                    return x[3];
                }
                return x[3];
            } else { // Continue in direction of x[3]
                temp = x[2];
                double temp_f = f[2];
                x[2] = x[3];
                f[2] = f[3];
                x[3] = x[2] + x[2] - x[1];
                if ((x[3] > x_limit[1] - eps_x) ||
                    (x[3] < x_limit[0] +
                                eps_x)) { // To close to limit. Other direction.
                    x[2] = temp;
                    f[2] = temp_f;
                    x[3] = x[2] + x[2] - x[1];
                    return x[3];
                }
                return x[3];
            }
        }
        // sort x1 -> x2 -> x3 in one direction
        else if (std::abs(x[1] - x[2]) >
                 std::abs(x[1] - x[3])) { // exchange x2 <-> x3
            temp = f[3];
            f[3] = f[2];
            f[2] = temp;
            temp = x[3];
            x[3] = x[2];
            x[2] = temp;
        }

        // check if minimum found
        if (f[2] < f[3]) { // minimum found
            xt = 0.5 *
                 (x[3] +
                  x[2]); // simple step, and switch to normal "+x" direction
            if (x[1] > x[2]) {
                temp = f[1];
                f[1] = f[3];
                f[3] = temp;
                temp = x[1];
                x[1] = x[3];
                x[3] = temp;
            }
            mode = 5;
        } else { // no minimum found, check if x3 is close to limit
            xt = x[3] + x[3] - x[2];
            i_limit = -1; // check limits
            if (xt < x_limit[0]) {
                i_limit = 0;
            } else if (xt > x_limit[1]) {
                i_limit = 1;
            }
            if (i_limit > -1) { // limit xt to x_limit0 <= xt <= x_limit1
                xt = x_limit[i_limit] - tau_fac * (x_limit[i_limit] - x[3]);
            }
            mode = 4;
            x[3] = xt;
        }
        return xt;
    }

    else if (mode == 4) { //  find intervall around minimum
        // std::cout << "Mode 4!" << std::endl;
        f[3] = chi2;
        if (f[3] <=
            f[2]) { //  no minimum found:  shift x1 <- x2 <- x3 and get next x3
            f[1] = f[2];
            f[2] = f[3];
            x[1] = x[2];
            x[2] = x[3];
            x[3] = x[2] + (1 + tau) * (x[2] - x[1]);
            f[3] = 0; // Goldener Schnitt
            // check limits
            i_limit = -1; // check limits
            if (x[3] < x_limit[0]) {
                i_limit = 0;
            } else if (x[3] > x_limit[1]) {
                i_limit = 1;
            }
            if (i_limit > -1) { // limit x3 to x_limit0 <= x3 <= x_limit1
                if (fabs(x[2] - x_limit[i_limit]) <
                    eps_x) { // Minimum found at limit : finish !
                    // std::cout << "Minimum found near limit!" << std::endl;
                    f[1] = f[2];
                    f[2] = f[3];
                    x[1] = x[2];
                    x[2] = x_limit[i_limit];
                    x[3] = x_limit[i_limit];
                    mode = -1;
                    return x_limit[i_limit];
                }
                x[3] = x_limit[i_limit] - tau_fac * (x_limit[i_limit] - x[2]);
                // std::cout << "Searching close to limit!" << std::endl;
                return x[3];
            }
            return x[3];
        } // case when not at limit and no minimum found

        else { // search intervall found
            //      cout << "search intervall found " << endl;
            mode = 5;
            xt =
                0.5 * (x[3] + x[2]); // ?  step not needed but easier to code...
            // switch to normal "+x" direction
            if (x[1] > x[2]) {
                temp = f[1];
                f[1] = f[3];
                f[3] = temp;
                temp = x[1];
                x[1] = x[3];
                x[3] = temp;
            }
        }
        return xt;
    }

    else if (mode == 5) { // ------------ find minimum within search range
                          // --------------
        // std::cout << "Mode 5!" << std::endl;
        ft = chi2;
        double delta_Chi2 = std::abs(ft - f[2]);

        if (ft < f[2]) { // better point found: narrow intervall
            if (xt < x[2]) {
                f[3] = f[2];
                x[3] = x[2];
                f[2] = ft;
                x[2] = xt;
            } else {
                f[1] = f[2];
                x[1] = x[2];
                f[2] = ft;
                x[2] = xt;
            }
        } else { // worse point
            if (xt < x[2]) {
                f[1] = ft;
                x[1] = xt;
            } else {
                f[3] = ft;
                x[3] = xt;
            }
        }

        // std::cout << "x[1]: " << x[1] <<  "  f[1]: " << f[1] << std::endl;
        // std::cout << "x[2]: " << x[2] <<  "  f[2]: " << f[2] << std::endl;
        // std::cout << "x[3]: " << x[3] <<  "  f[3]: " << f[3] << std::endl;

        if ((std::abs(x[1] - x[3]) < eps_x && delta_Chi2 < eps_f) ||
            std::abs(x[1] - x[3]) <
                0.00001) { // CONVERGENCE reached   Requirement of small
                           // delta_Chi2 difficult in very steep regions
            mode = 0;
            return x[2];
        }

        d21 = x[2] - x[1];
        d31 = x[3] - x[1];
        d32 = x[3] - x[2];
        g = 1. / d31 *
            ((f[3] - f[2]) * d21 / d32 +
             (f[2] - f[1]) * d32 / d21); // derivative
        H = 2. / d31 * ((f[3] - f[2]) / d32 - (f[2] - f[1]) / d21); // Hesse
        xt = x[2] - (g / H) * 1.02; // Newton Method

        if (std::abs(x[2] - xt) <
            close * d31) { // safety for numerical precision
            if (d21 < d32) {
                xt = x[2] + tau_fac * d31;
            } // and numerical precision
            else {
                xt = x[2] - tau_fac * d31;
            }
        }
        return xt;
    }
    return -1; // default return should never be reached
}

void PSFit::PSNewtonLimitShift(int sign, int np, double a[],
                               double a_limit[][2], double a_precision[],
                               double daN[], double h[], double g[],
                               double H[]) {
    // ------- for Newton Method: if close to limit shift central point by less
    // than precision
    if (sign > 0) {
        for (int ip = 0; ip < np; ip++) { // check limits
            daN[ip] = 0.;
            if (a[ip] < a_limit[ip][0] + a_precision[ip]) {
                daN[ip] = a_precision[ip] * 0.9;
                a[ip] = a_limit[ip][0] + daN[ip];
                h[ip] = a_precision[ip] * 0.8;
            } else if (a[ip] < a_limit[ip][0] + h[ip]) {
                h[ip] = std::abs(a_limit[ip][0] - a[ip]) * 0.99;
            }
            if (a[ip] > a_limit[ip][1] - a_precision[ip]) {
                daN[ip] = -a_precision[ip] * 0.9;
                a[ip] = a_limit[ip][1] + daN[ip];
                h[ip] = a_precision[ip] * 0.8;
            } else if (a[ip] > a_limit[ip][1] - h[ip]) {
                h[ip] = std::abs(a_limit[ip][1] - a[ip]) * 0.99;
            }
        }
    } else if (sign < 0) { //  ------------------    shift back to limit
        for (int ip = 0; ip < np; ip++) {
            daN[ip] = -daN[ip];
            a[ip] = a[ip] + daN[ip];
            for (int jp = 0; jp < np; jp++) {
                g[ip] =
                    g[ip] + H[ip * np + jp] * daN[jp]; // also shift derivative
            }
        }
    }
}

int PSFit::PSderivative(int icall, int np, double a[], double h[], double chi2,
                        double chi2_iter[], double g[], double H[]) {
    // Tool for numerical calculation of derivative and Hesse matrix
    //  depending on icall, the components of a[] are shifted up and down.
    //  In the following call the corresponding function (chi2) value is stored
    //  and the next shift to a[] is applied.
    //  After all shifts are done the derivative and Hesse matrix are calculated
    //  and the step size is recalculated (reduced).

    //  Int_t icall  ;                 //     current call number
    //  Int_t np  ;                    //     number of fit parameters
    //  Double_t a[np], Double_t h[np] ; //  fit parameters and step width
    //                                 //  (will be updated during the fit)
    //  Double_t chi2                   //  Function (chi2) value for current
    //  a[] Double_t chi2iter               //  Function (chi2) with nominal a[]
    //                                 //     (i.e. no shift of a[])
    //  Double_t g[4] ;                 //  derivative vector  g[np]
    //  Double_t H[4*4] ;               //  Hesse matrix       H[np*np]
    //                                 //  g[] and H[] are also used as
    //                                 //    intermediate storage of chi2
    //     icall = 0             for nominal a[]
    //     icall = 1...np        for + shift of a
    //     icall = np+1...2np    for - shift of a
    //     icall = 2np+1 ... (np*np+3np+2)/2        for ++ shift of a[0]
    //     icall = (np*np+3np+2)/2 ... np*np+np+1   for -- shift of a[0]

    //  Double_t chi2,  d ;             // chi2 value and stop criteria for fit
    if (np == 1) {
        // std::cout << "WARNING! Using PSderivative for 1-dim Case! Use
        // PSderivative1 instead!" << std::endl;
    }

    int ready = -1;
    int nstep = 1 + 2 * np + np * (np - 1);
    int iter;  // number of finished iteration
    int icalc; // current number of chi2 calculations
    //     within current iteration
    int ia;         // index of fit parameter to be changed
    int ia_i, ia_j; // same for Hesse_ij

    int shift_nom = 0; //     icall = 0             for nominal a[] //0
    int shift_p1 = 1;  //     etc see above                          //1
    int shift_p2 = shift_p1 + np - 1;              // 1
    int shift_m1 = shift_p2 + 1;                   // 2
    int shift_m2 = shift_m1 + np - 1;              // 2
    int shift_pp1 = shift_m2 + 1;                  // 3
    int shift_pp2 = shift_m2 + np * (np - 1) / 2;  // 3
    int shift_mm1 = shift_pp2 + 1;                 // 4
    int shift_mm2 = shift_pp2 + np * (np - 1) / 2; // 3

    iter = icall / nstep;
    icalc = icall - iter * nstep;

    int ia_old = -1, ia1_old = -1, ia_i_old = -1, ia_j_old = -1;
    int ia_new = -1, ia1_new = -1, ia_i_new = -1, ia_j_new = -1;
    double sign_old = 0., sign_new = 0.;

    if (icalc == shift_nom) {
        ia_old = -1;
        ia_new = 0;
        sign_old = 0.;
        sign_new = +1.;
    } else if (icalc < shift_p2) {
        ia_old = icalc - shift_p1;
        ia_new = ia_old + 1;
        sign_old = -1.;
        sign_new = +1.;
    } else if (icalc == shift_p2) {
        ia_old = np - 1;
        ia_new = 0;
        sign_old = -1.;
        sign_new = -1.;
    } else if (icalc < shift_m2) {
        ia_old = icalc - shift_m1;
        ia_new = ia_old + 1;
        sign_old = +1.;
        sign_new = -1.;

    } else if (icalc == shift_m2) {
        ia_old = np - 1;
        sign_old = +1.;
        ia_i_new = 0;
        ia_j_new = 1;
        sign_new = +1.;
    } else if (icalc < shift_pp2) {
        ia = icalc - shift_m2 - 1;
        int sum = 0;
        int test = -1;
        for (ia_i = 0; test < 0; ia_i++) {
            sum = sum + np - 1 - ia_i;
            if (ia < sum) {
                ia_i = ia_i - 1;
                ia_j = ia - sum + np;
                test = 1;
            };
        };

        ia_i_old = ia_i;
        ia_j_old = ia_j;
        sign_old = -1.;
        if (ia_j < np - 1) {
            ia_i_new = ia_i;
            ia_j_new = ia_j + 1;
            sign_new = +1.;
        } else {
            ia_i_new = ia_i + 1;
            ia_j_new = ia_i + 2;
            sign_new = +1.;
        };
    } else if (icalc == shift_pp2) {
        ia_i_old = np - 2;
        ia_j_old = np - 1;
        sign_old = -1.;
        ia_i_new = 0;
        ia_j_new = 1;
        sign_new = -1.;

    } else if (icalc < shift_mm2) {
        ia = icalc - shift_pp2 - 1;
        int sum = 0;
        int test = -1;
        for (ia_i = 0; test < 0; ia_i++) {
            sum = sum + np - 1 - ia_i;
            if (ia < sum) {
                ia_i = ia_i - 1;
                ia_j = ia - sum + np;
                test = 1;
            };
        };

        ia_i_old = ia_i;
        ia_j_old = ia_j;
        sign_old = +1.;
        if (ia_j < np - 1) {
            ia_i_new = ia_i;
            ia_j_new = ia_j + 1;
            sign_new = -1.;
        } else {
            ia_i_new = ia_i + 1;
            ia_j_new = ia_i + 2;
            sign_new = -1.;
        };
    } else if (icalc == shift_mm2) {
        ia_i_old = np - 2;
        ia_j_old = np - 1;
        sign_old = +1.;
    } else {
        // std::cout << "ERROR \n" << std::endl;
    };

    if (ia_old >= 0) {
        a[ia_old] = a[ia_old] + sign_old * h[ia_old];
    };
    if (ia_new >= 0) {
        a[ia_new] = a[ia_new] + sign_new * h[ia_new];
    };
    if (ia_i_old >= 0) {
        a[ia_i_old] = a[ia_i_old] + sign_old * h[ia_i_old];
    };
    if (ia_i_new >= 0) {
        a[ia_i_new] = a[ia_i_new] + sign_new * h[ia_i_new];
    };
    if (ia_j_old >= 0) {
        a[ia_j_old] = a[ia_j_old] + sign_old * h[ia_j_old];
    };
    if (ia_j_new >= 0) {
        a[ia_j_new] = a[ia_j_new] + sign_new * h[ia_j_new];
    };

    if (icalc == shift_nom) {
        chi2_iter[0] = chi2;
        for (int ii = 0; ii < np; ii++) {
            g[ii] = 0.; // build up derivative
            for (int jj = 0; jj < np; jj++) {
                H[ii * np + jj] = -2. * chi2; // build up Hesse
            };
        };
    };

    if (ia_old >= 0) {
        g[ia_old] = g[ia_old] - sign_old * chi2; // build up derivative
        H[ia_old * np + ia_old] =
            H[ia_old * np + ia_old] + chi2; // build up Hesse diagonal
    }
    if (ia_i_old >= 0) {
        H[ia_i_old * np + ia_j_old] =
            H[ia_i_old * np + ia_j_old] + chi2; // Hesse (non-diagonal)
    }
    if (icalc == shift_mm2) {
        for (int ii = 0; ii < np; ii++) {
            g[ii] = g[ii] / (2. * h[ii]); // final derivative
            H[ii * np + ii] =
                H[ii * np + ii] / pow(h[ii], 2.0); // final Hesse diagonal
        };
        for (int ii = 0; ii < np - 1; ii++) { // Hesse non-diagonal
            for (int jj = ii + 1; jj < np; jj++) {
                H[ii * np + jj] = H[ii * np + jj] -
                                  pow(h[ii], 2.0) * H[ii * np + ii] -
                                  pow(h[jj], 2.0) * H[jj * np + jj];
            };
        };
        for (int ii = 0; ii < np - 1; ii++) { // final Hesse non-diagonal
            for (int jj = ii + 1; jj < np; jj++) {
                H[ii * np + jj] = H[ii * np + jj] / (2. * h[ii] * h[jj]);
                H[jj * np + ii] = H[ii * np + jj];
            };
        };
        ready = 1;
    };
    return ready;
}

int PSFit::PSderivative1(int icall, double a[], double h[], double F,
                         double g[], double H[]) {
    // 1-dim case of PSderivative,  see documentation there
    // g,H will contain derivative and second derivative of F(a)
    // ready = 1 if derivatives are ready

    int nstep = 3;
    int iter;  // number of finished iterations
    int icalc; // current number of function calculations
               //     within current iteration
    int ready = -1;
    iter = icall / nstep;
    icalc = icall - iter * nstep;
    if (icalc == 0) {
        a[0] = a[0] + h[0];
    } else if (icalc == 1) {
        a[0] = a[0] - 2. * h[0];
    } else if (icalc == 2) {
        a[0] = a[0] + h[0];
    }

    if (icalc == 0) {
        g[0] = 0.;      // build up derivative
        H[0] = -2. * F; // build up Hesse
    } else if (icalc == 1) {
        g[0] = g[0] + F; // build up derivative
        H[0] = H[0] + F; // build up Hesse
    } else if (icalc == 2) {
        g[0] = g[0] - F;            // build up derivative
        H[0] = H[0] + F;            // build up Hesse
        g[0] = g[0] / (2. * h[0]);  // at last step: divide by step width
        H[0] = H[0] / pow(h[0], 2); // at last step: divide by step width
        ready = 1;
    }
    return ready;
}

double PSFit::PSNewtonAnalyzer(int np, double a[], double a_limit[][2],
                               double a_precision[], double daN[], double h[],
                               double g[], double H[], double H_inv[],
                               double chi2, bool noNewtonShifts) {
    double d;    // convergence test value
    int iNewton; // =1 if Newton Method can be used
    double g_norm, h_norm;
    double x_limit[2]; // straight line distance from a[] to alimits[][2]

    d = 100.;
    for (int ip = 0; ip < np; ip++) {
        daN[ip] = 0.;
    }

    //  --------- test validity of Hesse for Newton method --------
    iNewton = 1;
    if (np == 2) {
        if (H[0] * H[3] <= H[1] * H[2]) {
            iNewton = 0;
            // std::cout << "PSNewton Analyzer check 1 (det<0)" << std::endl;
        }
    }                                 // check only works for np=2
    for (int ip = 0; ip < np; ip++) { // check for positive Diagonal of Hesse
        if (H[ip * np + ip] < 0.) {
            iNewton = 0;
            // std::cout << "PSNewton Analyzer check 2 (negative diagonal
            // entries)" << std::endl;
            break;
        }
        if (H[ip * np + ip] * h[ip] < 0.07 * std::abs(g[ip])) {
            iNewton = 0;
            // std::cout << "PSNewton Analyzer check 3" << std::endl;
            break;
        }
    }
    if (iNewton == 1) { // Newton Step width
        double temp = PSMinverse(H, H_inv, np);
        for (int ii = 0; ii < np; ii++) {
            for (int jj = 0; jj < np; jj++) {
                daN[ii] = daN[ii] - H_inv[ii * np + jj] * g[jj];
            }
        }
        d = 0.;
        for (int ip = 0; ip < np; ip++) {
            d = d - g[ip] * daN[ip];
        } // convergence test value

        if (d < 0.) {
            iNewton = 0;
        } // daN wrong since not in direction of -g
    }

    // ------------- Newton not ok ------------------------
    if (iNewton != 1 || noNewtonShifts) {
        // std::cout << "PSNewton Analyzer ===== WARNING ========== H is not
        // positive definit" << std::endl;
        if (noNewtonShifts) {
            // std::cout << "Newton Shifts are disabled!" << std::endl;
        }

        g_norm = PSVnorm(g, np); // continue in direction of derivative
        h_norm = PSVnorm(h, np); // continue with previous step width
        for (int ip = 0; ip < np; ip++) {
            daN[ip] = -g[ip] * h_norm / g_norm;
            // std::cout << "daN[" << np << "] set to " << daN[ip] << std::endl;
        }
        for (int ip = 0; ip < np; ip++) { // check limits
            if (a[ip] < a_limit[ip][0] + a_precision[ip]) {
                h[ip] = a_precision[ip];
                daN[ip] = std::max(daN[ip], 0.);
                // std::cout << "At lower limit for parameter " << ip <<
                // std::endl;
            } else if (a[ip] > a_limit[ip][1] - a_precision[ip]) {
                h[ip] = a_precision[ip];
                daN[ip] = std::min(daN[ip], 0.);
                // std::cout << "At upper limit for parameter " << ip <<
                // std::endl;
            }
        }
        // std::cout << "daN set to: daN[0]: " << daN[0] << " daN[1]: " <<
        // daN[1] << std::endl;
        d = 110.;
        return d;
    }

    // ------------ Newton step seems to be fine --------------
    PSLineLimit(np, a, daN, a_limit, x_limit); // get distance to limit

    // std::cout << "daN shift for Newton in daN[0]: " << daN[0] << " daN[1]: "
    // << daN[1] << std::endl;

    double x_max =
        std::min(1., x_limit[1]); //  Newton step would lead out of limits

    for (int ip = 0; ip < np; ip++) {    // apply Newton step
        a[ip] = a[ip] + daN[ip] * x_max; // update parameter value
        h[ip] =
            100. * sqrt(0.000001 *
                        std::abs(chi2 / H[ip * np + ip])); // update step width
    }

    // std::cout << "Newton shifted to: a[0] = " << a[0] << " a[1] = " << a[1]
    // << std::endl;

    for (int ip = 0; ip < np; ip++) { // check limits
        if (a[ip] < a_limit[ip][0] + a_precision[ip]) {
            h[ip] = a_precision[ip];
            daN[ip] = std::max(daN[ip], 0.);
            // std::cout << "At lower limit for parameter " << ip << std::endl;
        } else if (a[ip] > a_limit[ip][1] - a_precision[ip]) {
            h[ip] = a_precision[ip];
            daN[ip] = std::min(daN[ip], 0.);
            // std::cout << "At upper limit for parameter " << ip << std::endl;
        }
    }

    d = 0.;
    for (int ip = 0; ip < np; ip++) { // check limits
        d = d -
            g[ip] * daN[ip]; // convergence test value (only when not at limit)
    }
    if (d < 0.) {
        d = 111.;
    }

    double d2 = 0.;
    for (int ip = 0; ip < np; ip++) {
        for (int jp = 0; jp < np; jp++) {
            d2 = d2 + daN[ip] * H[ip * np + jp] * daN[jp];
        }
    }

    // std::cout << "daN set to: dan[0]: " << daN[0] << " daN[1]: " << daN[1] <<
    // std::endl;
    return d;
}

double PSFit::PSMinverse(double H[], double H_inv[],
                         int p) { // invert a symmetric matrix H[p,p]
    // for 2*2 and 3*3 inversion also works for a non-symmetric matrix
    double temp;
    if (p == 2) {
        double det = H[0] * H[3] - H[1] * H[2];
        H_inv[0] = H[3] / det;
        H_inv[3] = H[0] / det;
        H_inv[1] = -H[1] / det;
        H_inv[2] = -H[2] / det;
    } else if (p == 3) {
        H_inv[0 * p + 0] =
            H[1 * p + 1] * H[2 * p + 2] - H[1 * p + 2] * H[2 * p + 1];
        H_inv[0 * p + 1] =
            H[1 * p + 2] * H[2 * p + 0] - H[1 * p + 0] * H[2 * p + 2];
        H_inv[0 * p + 2] =
            H[1 * p + 0] * H[2 * p + 1] - H[1 * p + 1] * H[2 * p + 0];
        H_inv[1 * p + 0] =
            H[2 * p + 1] * H[0 * p + 2] - H[2 * p + 2] * H[0 * p + 1];
        H_inv[1 * p + 1] =
            H[2 * p + 2] * H[0 * p + 0] - H[2 * p + 0] * H[0 * p + 2];
        H_inv[1 * p + 2] =
            H[2 * p + 0] * H[0 * p + 1] - H[2 * p + 1] * H[0 * p + 0];
        H_inv[2 * p + 0] =
            H[0 * p + 1] * H[1 * p + 2] - H[0 * p + 2] * H[1 * p + 1];
        H_inv[2 * p + 1] =
            H[0 * p + 2] * H[1 * p + 0] - H[0 * p + 0] * H[1 * p + 2];
        H_inv[2 * p + 2] =
            H[0 * p + 0] * H[1 * p + 1] - H[0 * p + 1] * H[1 * p + 0];
        double det = H[0 * p + 0] * H_inv[0 * p + 0] +
                     H[1 * p + 0] * H_inv[0 * p + 1] +
                     H[2 * p + 0] * H_inv[0 * p + 2];

        for (int ii = 0; ii < p; ii++) {
            for (int jj = 0; jj < p; jj++) {
                H_inv[ii * p + jj] = H_inv[ii * p + jj] / det;
            }
        }
    } else {
        temp = PSMCholesky(H, H_inv, p);
        temp = PSMRTrianInvert2(H_inv, p);
        temp = PSMmultiplyMRRT(H_inv, p, p);
    };
    return 0.;
}

double PSFit::PSMCholesky(
    double M[], double R[],
    int n) { // Cholesky decomposition of  M = n*n  symmetric matrix
    // M = R^T * R
    for (int ii = 0; ii < n; ii++) {
        R[ii * n + ii] = M[ii * n + ii];
        for (int ll = 0; ll < ii; ll++) {
            R[ii * n + ii] = R[ii * n + ii] - pow(R[ll * n + ii], 2.0);
        };
        if (R[ii * n + ii] <= 0.) {
            return -1.;
        }
        R[ii * n + ii] = sqrt(R[ii * n + ii]);

        for (int jj = ii + 1; jj < n; jj++) {
            R[ii * n + jj] = 0.;
            for (int kk = 0; kk < ii; kk++) {
                R[ii * n + jj] =
                    R[ii * n + jj] + R[kk * n + ii] * R[kk * n + jj];
            };
            R[ii * n + jj] =
                1. / R[ii * n + ii] * (M[ii * n + jj] - R[ii * n + jj]);
        }
        for (int jj = 0; jj < ii; jj++) {
            R[ii * n + jj] = 0.;
        }
    }
    return 0.;
}

double PSFit::PSMRTrianInvert2(double R[],
                               int n) { // invert a triangular R - matrix (n*n)
    //  store result in position of input (R)
    for (int ii = 0; ii < n; ii++) {
        for (int jj = ii + 1; jj < n; jj++) {
            R[jj * n + ii] = R[ii * n + jj];
            R[ii * n + jj] = 0.;
        }
        for (int jj = ii + 1; jj < n; jj++) {
            R[ii * n + jj] =
                R[ii * n + jj] - 1. / R[ii * n + ii] * R[jj * n + ii];
            for (int kk = ii + 1; kk < jj; kk++) {
                R[ii * n + jj] =
                    R[ii * n + jj] - R[ii * n + kk] * R[kk * n + jj];
            }
            R[ii * n + jj] = R[ii * n + jj] / R[jj * n + jj];
        }
        R[ii * n + ii] = 1. / R[ii * n + ii];
        for (int jj = ii + 1; jj < n; jj++) {
            R[jj * n + ii] = 0.;
        }
    }
    return 0.;
}

double PSFit::PSMmultiplyMRRT(
    double A[], int n1,
    int n2) { // multiply triangular matrix R=A  with transpose,
    //  store result in input, B = R * R^T,  A = B,
    //  Aij= Bik * Bjk,  A(n1,n1)   B(n1,n2)
    for (int ii = 0; ii < n1; ii++) {
        for (int jj = ii + 1; jj < n1; jj++) {
            for (int kk = jj; kk < n2; kk++) {
                A[jj * n1 + ii] =
                    A[jj * n1 + ii] + A[ii * n1 + kk] * A[jj * n2 + kk];
            }
        }
        A[ii * n1 + ii] = pow(A[ii * n1 + ii], 2.0);
        for (int ll = ii + 1; ll < n2; ll++) {
            A[ii * n1 + ii] = A[ii * n1 + ii] + pow(A[ii * n1 + ll], 2);
        }
        for (int jj = ii + 1; jj < n1; jj++) {
            A[ii * n1 + jj] = A[jj * n1 + ii];
        }
    }
    return 0;
}