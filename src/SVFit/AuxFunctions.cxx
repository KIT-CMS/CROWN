#include "../../include/SVFit/AuxFunctions.hxx"

#include <TF1.h>
#include <TFitResult.h>
#include <TMath.h>

namespace fastmtt {

double roundToNdigits(double x, int n) {
    double tmp = TMath::Power(10., n);
    if (x != 0.) {
        tmp /= TMath::Power(10., TMath::Floor(TMath::Log10(TMath::Abs(x))));
    }
    double x_rounded = TMath::Nint(x * tmp) / tmp;
    // std::cout << "<roundToNdigits>: x = " << x << ", x_rounded = " <<
    // x_rounded << std::endl;
    return x_rounded;
}

TGraphErrors *makeGraph(const std::string &graphName,
                        const std::vector<GraphPoint> &graphPoints) {
    // std::cout << "<makeGraph>:" << std::endl;
    size_t numPoints = graphPoints.size();
    // std::cout << " numPoints = " << numPoints << std::endl;
    TGraphErrors *graph = new TGraphErrors(numPoints);
    graph->SetName(graphName.data());
    for (size_t iPoint = 0; iPoint < numPoints; ++iPoint) {
        const GraphPoint &graphPoint = graphPoints[iPoint];
        graph->SetPoint(iPoint, graphPoint.x_, graphPoint.y_);
        graph->SetPointError(iPoint, graphPoint.xErr_, graphPoint.yErr_);
    }
    return graph;
}

void extractResult(TGraphErrors *graph, double &mass, double &massErr,
                   double &Lmax, int verbosity) {
    // determine range of mTest values that are within ~2 sigma interval within
    // maximum of likelihood function
    double x_Lmax = 0.;
    double y_Lmax = 0.;
    double idxPoint_Lmax = -1;
    for (int iPoint = 0; iPoint < graph->GetN(); ++iPoint) {
        double x, y;
        graph->GetPoint(iPoint, x, y);
        if (y > y_Lmax) {
            x_Lmax = x;
            y_Lmax = y;
            idxPoint_Lmax = iPoint;
        }
    }

    double xMin = 1.e+6;
    double xMax = 0.;
    for (int iPoint = 0; iPoint < graph->GetN(); ++iPoint) {
        double x, y;
        graph->GetPoint(iPoint, x, y);
        if (x < xMin)
            xMin = x;
        if (x > xMax)
            xMax = x;
    }

    // fit log-likelihood function within ~2 sigma interval within maximum
    // with parabola
    std::vector<GraphPoint> graphPoints_forFit;
    double xMin_fit = 1.e+6;
    double xMax_fit = 0.;
    for (int iPoint = 0; iPoint < graph->GetN(); ++iPoint) {
        double x, y;
        graph->GetPoint(iPoint, x, y);
        double xErr = graph->GetErrorX(iPoint);
        double yErr = graph->GetErrorY(iPoint);
        // std::cout << "point #" << iPoint << ": x = " << x << " +/- " << xErr
        // << ", y = " << y << " +/- " << yErr << std::endl;
        if (y > (1.e-1 * y_Lmax) && TMath::Abs(iPoint - idxPoint_Lmax) <= 5) {
            GraphPoint graphPoint;
            graphPoint.x_ = x;
            graphPoint.xErr_ = xErr;
            if ((x - xErr) < xMin_fit)
                xMin_fit = x - xErr;
            if ((x + xErr) > xMax_fit)
                xMax_fit = x + xErr;
            graphPoint.y_ = -TMath::Log(y);
            graphPoint.yErr_ = yErr / y;
            graphPoints_forFit.push_back(graphPoint);
        }
    }

    TGraphErrors *likelihoodGraph_forFit =
        makeGraph("svFitLikelihoodGraph_forFit", graphPoints_forFit);
    int numPoints = likelihoodGraph_forFit->GetN();
    bool useFit = false;
    if (numPoints >= 3) {
        TF1 *fitFunction =
            new TF1("fitFunction", "TMath::Power((x - [0])/[1], 2.) + [2]",
                    xMin_fit, xMax_fit);
        fitFunction->SetParameter(0, x_Lmax);
        fitFunction->SetParameter(1, 0.20 * x_Lmax);
        fitFunction->SetParameter(2, -TMath::Log(y_Lmax));

        std::string fitOptions = "NSQ";
        // if ( !verbosity ) fitOptions.append("Q");
        TFitResultPtr fitResult =
            likelihoodGraph_forFit->Fit(fitFunction, fitOptions.data());
        if (fitResult.Get()) {
            if (verbosity >= 1) {
                std::cout << "fitting graph of p versus M(test) in range "
                          << xMin_fit << ".." << xMax_fit
                          << ", result:" << std::endl;
                std::cout << " parameter #0 = " << fitFunction->GetParameter(0)
                          << " +/- " << fitFunction->GetParError(0)
                          << std::endl;
                std::cout << " parameter #1 = " << fitFunction->GetParameter(1)
                          << " +/- " << fitFunction->GetParError(1)
                          << std::endl;
                std::cout << " parameter #2 = " << fitFunction->GetParameter(2)
                          << " +/- " << fitFunction->GetParError(2)
                          << std::endl;
                std::cout << "chi^2 = " << fitResult->Chi2() << std::endl;
            }
            if (fitResult->Chi2() < (10. * numPoints) &&
                fitFunction->GetParameter(0) > xMin &&
                fitFunction->GetParameter(0) < xMax &&
                TMath::Abs(fitFunction->GetParameter(0) - x_Lmax) <
                    (0.10 * x_Lmax)) {
                mass = fitFunction->GetParameter(0);
                massErr = TMath::Sqrt(square(fitFunction->GetParameter(1)) +
                                      square(fitFunction->GetParError(0)));
                Lmax = TMath::Exp(-fitFunction->GetParameter(2));
                // std::cout << "fit: mass = " << mass << " +/- " << massErr <<
                // " (Lmax = " << Lmax << ")" << std::endl;
                useFit = true;
            }
        } else {
            std::cerr << "Warning in <extractResult>: Fit did not converge !!"
                      << std::endl;
        }
        delete fitFunction;
    }
    if (!useFit) {
        mass = x_Lmax;
        massErr = TMath::Sqrt(0.5 * (square(x_Lmax - xMin_fit) +
                                     square(xMax_fit - x_Lmax))) /
                  TMath::Sqrt(2. * TMath::Log(10.));
        Lmax = y_Lmax;
        // std::cout << "graph: mass = " << mass << " +/- " << massErr << "
        // (Lmax = " << Lmax << ")" << std::endl;
    }

    delete likelihoodGraph_forFit;
}

Vector normalize(const Vector &p) {
    double p_x = p.x();
    double p_y = p.y();
    double p_z = p.z();
    double mag2 = square(p_x) + square(p_y) + square(p_z);
    if (mag2 <= 0.)
        return p;
    double mag = TMath::Sqrt(mag2);
    return Vector(p_x / mag, p_y / mag, p_z / mag);
}

double compScalarProduct(const Vector &p1, const Vector &p2) {
    return (p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z());
}

Vector compCrossProduct(const Vector &p1, const Vector &p2) {
    double p3_x = p1.y() * p2.z() - p1.z() * p2.y();
    double p3_y = p1.z() * p2.x() - p1.x() * p2.z();
    double p3_z = p1.x() * p2.y() - p1.y() * p2.x();
    return Vector(p3_x, p3_y, p3_z);
}

double compCosThetaNuNu(double visEn, double visP, double visMass2,
                        double nunuEn, double nunuP, double nunuMass2) {
    double cosThetaNuNu =
        (visEn * nunuEn - 0.5 * (tauLeptonMass2 - (visMass2 + nunuMass2))) /
        (visP * nunuP);
    return cosThetaNuNu;
}

double compPSfactor_tauToLepDecay(double x, double visEn, double visP,
                                  double visMass, double nunuEn, double nunuP,
                                  double nunuMass) {
    // std::cout << "<compPSfactor_tauToLepDecay>:" << std::endl;
    // std::cout << " x = " << x << std::endl;
    // std::cout << " visEn = " << visEn << std::endl;
    // std::cout << " visP = " << visP << std::endl;
    // std::cout << " visMass = " << visMass << std::endl;
    // std::cout << " nunuEn = " << nunuEn << std::endl;
    // std::cout << " nunuP = " << nunuP << std::endl;
    // std::cout << " nunuMass = " << nunuMass << std::endl;
    double visMass2 = square(visMass);
    double nunuMass2 = square(nunuMass);
    if (x >= (visMass2 / tauLeptonMass2) && x <= 1. &&
        nunuMass2 < ((1. - x) * tauLeptonMass2)) { // physical solution
        double tauEn_rf =
            (tauLeptonMass2 + nunuMass2 - visMass2) / (2. * nunuMass);
        double visEn_rf = tauEn_rf - nunuMass;
        if (!(tauEn_rf >= tauLeptonMass && visEn_rf >= visMass))
            return 0.;
        double I =
            nunuMass2 *
            (2. * tauEn_rf * visEn_rf -
             (2. / 3.) * TMath::Sqrt((square(tauEn_rf) - tauLeptonMass2) *
                                     (square(visEn_rf) - visMass2)));
#ifdef XSECTION_NORMALIZATION
        I *= GFfactor;
#endif
        double cosThetaNuNu =
            compCosThetaNuNu(visEn, visP, visMass2, nunuEn, nunuP, nunuMass2);
        if (!(cosThetaNuNu >= (-1. + epsilon) && cosThetaNuNu <= +1.))
            return 0.;
        double PSfactor =
            (visEn + nunuEn) * I /
            (8. * visP * square(x) *
             TMath::Sqrt(square(visP) + square(nunuP) +
                         2. * visP * nunuP * cosThetaNuNu + tauLeptonMass2));
//-------------------------------------------------------------------------
// CV: fudge factor to reproduce literature value for cross-section times
// branching fraction
#ifdef XSECTION_NORMALIZATION
        PSfactor *= 2.;
#endif
        //-------------------------------------------------------------------------
        return PSfactor;
    } else {
        return 0.;
    }
}

double compPSfactor_tauToHadDecay(double x, double visEn, double visP,
                                  double visMass, double nuEn, double nuP) {
    // std::cout << "<compPSfactor_tauToHadDecay>:" << std::endl;
    // std::cout << " x = " << x << std::endl;
    // std::cout << " visEn = " << visEn << std::endl;
    // std::cout << " visP = " << visP << std::endl;
    // std::cout << " visMass = " << visMass << std::endl;
    // std::cout << " nuEn = " << nuEn << std::endl;
    // std::cout << " nuP = " << nuP << std::endl;
    double visMass2 = square(visMass);
    if (x >= (visMass2 / tauLeptonMass2) && x <= 1.) { // physical solution
        double cosThetaNu =
            compCosThetaNuNu(visEn, visP, visMass2, nuEn, nuP, 0.);
        // std::cout << "cosThetaNu = " << cosThetaNu << std::endl;
        if (!(cosThetaNu >= (-1. + epsilon) && cosThetaNu <= +1.))
            return 0.;
        double PSfactor =
            (visEn + nuEn) /
            (8. * visP * square(x) *
             TMath::Sqrt(square(visP) + square(nuP) +
                         2. * visP * nuP * cosThetaNu + tauLeptonMass2));
        PSfactor *= 1.0 / (tauLeptonMass2 - visMass2);
//-------------------------------------------------------------------------
// CV: multiply by constant matrix element,
//     chosen such that the branching fraction of the tau to decay into hadrons
//     is reproduced
// const double M2
// = 16.*TMath::Pi()*cube(tauLeptonMass)*GammaTauToHad/(tauLeptonMass2 -
// visMass2); Remove multiplication as it add to execution time, and does not
// alter the result.
#ifdef XSECTION_NORMALIZATION
        PSfactor *= M2;
#endif
        //-------------------------------------------------------------------------
        return PSfactor;
    } else {
        return 0.;
    }
}

void integrationParameters::reset() {
    idx_X_ = -1;
    idx_phi_ = -1;
    idx_VisPtShift_ = -1;
    idx_mNuNu_ = -1;
}

} // namespace fastmtt
