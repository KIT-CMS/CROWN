#include <algorithm>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#include "TMath.h"
#include "TMatrixD.h"
#include "TVector2.h"

#include "Math/BasicMinimizer.h"
#include "TF1.h"

#include "../../include/SVFit/FastMTT.hxx"
#include "../../include/SVFit/MeasuredTauLepton.hxx"

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::Likelihood() {

    covMET.ResizeTo(2, 2);

    compnentsBitWord.reset();
    enableComponent(fastMTT::MASS);
    enableComponent(fastMTT::MET);
    /// experimental components. Disabled by default.
    disableComponent(fastMTT::PX);
    disableComponent(fastMTT::PY);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::~Likelihood() {}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
template <unsigned int max_order>
std::array<double, max_order + 1> getPowTable(double x) {
    std::array<double, max_order + 1> powerTable{};
    powerTable[0] = 1.;
    for (unsigned int i = 1; i <= max_order; ++i) {
        powerTable[i] = powerTable[i - 1] * x;
    }
    return powerTable;
};
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setLeptonInputs(const LorentzVector &aLeg1P4,
                                 const LorentzVector &aLeg2P4,
                                 int aLeg1DecayType, int aLeg2DecayType,
                                 int aLeg1DecayMode, int aLeg2DecayMode) {
    leg1P4 = aLeg1P4;
    leg2P4 = aLeg2P4;
    using PowTable = std::array<double, 5u>;
    auto getPowTableLam = [](double x) {
        PowTable powerTable{};
        powerTable[0] = 1.;
        for (unsigned int i = 1; i <= 4; ++i) {
            powerTable[i] = powerTable[i - 1] * x;
        }
        return powerTable;
    };
    allpTpows = std::array<std::array<PowTable, 3u>, 2u>{
        {{{getPowTableLam(leg1P4.Px()), getPowTableLam(leg1P4.Py()),
           getPowTableLam(leg1P4.Pz())}},
         {{getPowTableLam(leg2P4.Px()), getPowTableLam(leg2P4.Py()),
           getPowTableLam(leg2P4.Pz())}}}};

    mVis = (leg1P4 + leg2P4).M();

    mVisLeg1 = leg1P4.M();
    mVisLeg2 = leg2P4.M();

    if (aLeg1DecayType == fastmtt::MeasuredTauLepton::kTauToHadDecay &&
        mVisLeg1 > 1.5) {
        mVisLeg1 = 0.3;
    }
    if (aLeg2DecayType == fastmtt::MeasuredTauLepton::kTauToHadDecay &&
        mVisLeg2 > 1.5) {
        mVisLeg2 = 0.3;
    }

    const double &mTau = fastmtt::tauLeptonMass;
    mVis1OverTauSquare = std::pow(mVisLeg1 / mTau, 2);
    mVis2OverTauSquare = std::pow(mVisLeg2 / mTau, 2);

    leg1DecayType = aLeg1DecayType;
    leg2DecayType = aLeg2DecayType;

    leg1DecayMode = aLeg1DecayMode;
    leg2DecayMode = aLeg2DecayMode;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setMETInputs(const LorentzVector &aMET,
                              const TMatrixD &aCovMET) {
    recoMET = aMET;
    covMET = aCovMET;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setParameters(const std::vector<double> &aPars) {

    parameters = aPars;
    if (parameters.size() < 2)
        parameters = std::vector<double>{6, 1.0 / 1.15};
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::enableComponent(fastMTT::likelihoodComponent aCompIndex) {

    compnentsBitWord.set(aCompIndex);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::disableComponent(fastMTT::likelihoodComponent aCompIndex) {

    compnentsBitWord.reset(aCompIndex);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::massLikelihood(const double &m) const {

    double coeff1 = parameters[0];
    double coeff2 = parameters[1];
    double mScaled = m * coeff2;

    if (mScaled < mVis)
        return 0.0;

    const double mVS2 = std::pow(mVis / mScaled, 2);

    double x1Min = std::min(1.0, mVis1OverTauSquare);
    double x2Min = std::max(mVis2OverTauSquare, mVS2);
    double x2Max = std::min(1.0, mVS2 / x1Min);

    if (x2Max < x2Min)
        return 0.0;

    double jacobiFactor = 2.0 * std::pow(mVis, 2) * std::pow(mScaled, -coeff1);
    double x2IntegralTerm = log(x2Max) - log(x2Min);

    double value = x2IntegralTerm;
    if (leg1DecayType != fastmtt::MeasuredTauLepton::kTauToHadDecay) {

        double mNuNuIntegralTermLeg1 =
            mVS2 * (std::pow(x2Max, -1) - std::pow(x2Min, -1));
        value += mNuNuIntegralTermLeg1;
    }
    if (leg2DecayType != fastmtt::MeasuredTauLepton::kTauToHadDecay) {
        double mNuNuIntegralTermLeg2 = mVS2 * x2IntegralTerm - (x2Max - x2Min);
        value += mNuNuIntegralTermLeg2;
    }

    /// The E9 factor to get values around 1.0
    value *= 1E9 * jacobiFactor;

    return value;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::ptLikelihood(const double &pTTauTau, int type) const {

    /// Protection against numerical singularity in phase space volume.
    if (std::abs(pTTauTau) < 0.5)
        return 0.0;

    const auto pT1pow = allpTpows[0][type];
    const auto pT2pow = allpTpows[1][type];
    const auto pT1 = pT1pow[1];
    const auto pT2 = pT2pow[1];
    const auto pTTauTauPow = getPowTable<5>(pTTauTau);

    Double_t x1Min = std::min(1.0, mVis1OverTauSquare);
    Double_t x2Min = std::min(1.0, mVis2OverTauSquare);

    Double_t x2Max = 1.0;
    Double_t x1Max = 1.0;

    Double_t a_x2 = x1Min * pT2 / (x1Min * pTTauTau - pT1);
    Double_t b_x2 = x1Max * pT2 / (x1Max * pTTauTau - pT1);

    const bool is_x1_vs_x2_falling = (-pT2 * pT1) < 0;
    const double x1_singularity = pT1 / pTTauTau;
    const bool x2_vs_x1_hasSingularity =
        x1_singularity > 0.0 && x1_singularity < 1.0;
    if (x2_vs_x1_hasSingularity && x1_singularity < x1Min)
        return 0.0;

    if (is_x1_vs_x2_falling) {
        x2Min = std::max(x2Min, b_x2);
        x2Max = std::min(x2Max, a_x2);
        if (x2_vs_x1_hasSingularity && x2Max < 0)
            x2Max = 1.0;
    } else {
        x2Min = std::max(x2Min, a_x2);
        x2Max = std::min(x2Max, b_x2);
        if (x2_vs_x1_hasSingularity && x2Max < 0)
            x2Max = 1.0;
    }

    if (x2Min < 0)
        x2Min = 0.0;
    if (x2Min > x2Max)
        return 0.0;

    Double_t mNuNuIntegral = 0.0;
    Double_t x2 = std::min(1.0, x2Max);

    const Double_t term1 = pT2 - pTTauTau * x2;
    const Double_t log_term1 = log(std::abs(term1));
    const Double_t term1Square = std::pow(term1, 2);

    Double_t integralMax =
        (pT1 * (pTTauTau * x2 + pT2pow[2] / term1 + 2 * pT2 * log_term1)) /
        pTTauTauPow[3];
    if (leg1DecayType != fastmtt::MeasuredTauLepton::kTauToHadDecay) {
        mNuNuIntegral =
            -pT1pow[2] *
            (2 * pTTauTau * x2 +
             (pT2pow[2] * (5 * pT2 - 6 * pTTauTau * x2)) / term1Square +
             6 * pT2 * log_term1) /
            (2 * pTTauTauPow[4]);
    }
    if (leg2DecayType != fastmtt::MeasuredTauLepton::kTauToHadDecay) {
        mNuNuIntegral += -pT1 / (2 * pTTauTauPow[5]) *
                         (2 * pT2 * pTTauTau * (-3 * pT1 + 2 * pTTauTau) * x2 +
                          pTTauTauPow[2] * (-pT1 + pTTauTau) * pow(x2, 2) +
                          (pT2pow[4] * pT1) / term1Square +
                          (2 * pT2pow[3] * (-4 * pT1 + pTTauTau)) / term1 +
                          6 * pT2pow[2] * (-2 * pT1 + pTTauTau) * log_term1);
    }
    integralMax += mNuNuIntegral;

    x2 = x2Min;
    const Double_t term2 = pT2 - pTTauTau * x2;
    const Double_t log_term2 = log(std::abs(term2));
    const Double_t term2Square = std::pow(term2, 2);

    Double_t integralMin =
        (pT1 * (pTTauTau * x2 + pT2pow[2] / term2 + 2 * pT2 * log_term2)) /
        pTTauTauPow[3];
    if (leg1DecayType != fastmtt::MeasuredTauLepton::kTauToHadDecay) {
        mNuNuIntegral =
            -pT1pow[2] *
            (2 * pTTauTau * x2 +
             (pT2pow[2] * (5 * pT2 - 6 * pTTauTau * x2)) / term2Square +
             6 * pT2 * log_term2) /
            (2 * pTTauTauPow[4]);
    }
    if (leg2DecayType != fastmtt::MeasuredTauLepton::kTauToHadDecay) {
        mNuNuIntegral += -pT1 / (2 * pTTauTauPow[5]) *
                         (2 * pT2 * pTTauTau * (-3 * pT1 + 2 * pTTauTau) * x2 +
                          pTTauTauPow[2] * (-pT1 + pTTauTau) * pow(x2, 2) +
                          (pT2pow[4] * pT1) / term2Square +
                          (2 * pT2pow[3] * (-4 * pT1 + pTTauTau)) / term2 +
                          6 * pT2pow[2] * (-2 * pT1 + pTTauTau) * log_term2);
    }

    integralMin += mNuNuIntegral;

    Double_t value = integralMax - integralMin;

    /// The 1E4 factor to get values around 1.0
    value *= 1E4;

    return std::abs(value);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::metTF(const LorentzVector &metP4, const LorentzVector &nuP4,
                         const TMatrixD &covMET) const {

    const double aMETx = metP4.X();
    const double aMETy = metP4.Y();

    double invCovMETxx = covMET(1, 1);
    double invCovMETxy = -covMET(0, 1);
    double invCovMETyx = -covMET(1, 0);
    double invCovMETyy = covMET(0, 0);
    double covDet = invCovMETxx * invCovMETyy - invCovMETxy * invCovMETyx;

    if (std::abs(covDet) < 1E-10) {
        std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!"
                  << "METx: " << aMETy << " METy: " << aMETy << std::endl;
        return 0;
    }
    double const_MET = 1. / (2. * M_PI * TMath::Sqrt(covDet));
    double residualX = aMETx - (nuP4.X());
    double residualY = aMETy - (nuP4.Y());

    double pull2 =
        residualX * (invCovMETxx * residualX + invCovMETxy * residualY) +
        residualY * (invCovMETyx * residualX + invCovMETyy * residualY);
    pull2 /= covDet;

    return const_MET * TMath::Exp(-0.5 * pull2);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double Likelihood::value(const double *x) const {

    double x1Min = std::min(1.0, mVis1OverTauSquare);
    double x2Min = std::min(1.0, mVis2OverTauSquare);

    if (x[0] < x1Min || x[1] < x2Min)
        return 0.0;

    const auto testP4 = leg1P4 * (1.0 / x[0]) + leg2P4 * (1.0 / x[1]);
    const auto testMET = testP4 - leg1P4 - leg2P4;

    double value = -1.0;
    if (compnentsBitWord.test(fastMTT::MET))
        value *= metTF(recoMET, testMET, covMET);
    if (compnentsBitWord.test(fastMTT::MASS))
        value *= massLikelihood(testP4.M());
    if (compnentsBitWord.test(fastMTT::PX))
        value *= ptLikelihood(testP4.Px(), 0);
    if (compnentsBitWord.test(fastMTT::PY))
        value *= ptLikelihood(testP4.Py(), 1);
    return value;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
FastMTT::FastMTT() {

    minimizerName = "Minuit2";
    minimizerAlgorithm = "Migrad";
    initialize();
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
FastMTT::~FastMTT() { delete minimizer; }
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::initialize() {

    minimizer =
        ROOT::Math::Factory::CreateMinimizer(minimizerName, minimizerAlgorithm);
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(100000);
    minimizer->SetTolerance(0.01);

    std::vector<std::string> varNames = {"x1", "x2"};
    nVariables = varNames.size();
    std::vector<double> initialValues(nVariables, 0.5);
    std::vector<double> stepSizes(nVariables, 0.01);
    minimumPosition = initialValues;
    minimumValue = 999.0;

    for (unsigned int iVar = 0; iVar < nVariables; ++iVar) {
        minimizer->SetVariable(iVar, varNames[iVar].c_str(),
                               initialValues[iVar], stepSizes[iVar]);
    }

    std::vector<double> shapeParams = {6, 1.0 / 1.15};
    setLikelihoodParams(shapeParams);
    likelihoodFunctor =
        new ROOT::Math::Functor(&myLikelihood, &Likelihood::value, nVariables);
    minimizer->SetFunction(*likelihoodFunctor);

    verbosity = 0;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::setLikelihoodParams(const std::vector<double> &aPars) {

    myLikelihood.setParameters(aPars);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::enableComponent(fastMTT::likelihoodComponent aCompIndex) {

    myLikelihood.enableComponent(aCompIndex);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::disableComponent(fastMTT::likelihoodComponent aCompIndex) {

    myLikelihood.disableComponent(aCompIndex);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
bool FastMTT::compareLeptons(
    const fastmtt::MeasuredTauLepton &measuredTauLepton1,
    const fastmtt::MeasuredTauLepton &measuredTauLepton2) {

    using namespace fastmtt;

    if ((measuredTauLepton1.type() == MeasuredTauLepton::kTauToElecDecay ||
         measuredTauLepton1.type() == MeasuredTauLepton::kTauToMuDecay) &&
        measuredTauLepton2.type() == MeasuredTauLepton::kTauToHadDecay)
        return true;
    if ((measuredTauLepton2.type() == MeasuredTauLepton::kTauToElecDecay ||
         measuredTauLepton2.type() == MeasuredTauLepton::kTauToMuDecay) &&
        measuredTauLepton1.type() == MeasuredTauLepton::kTauToHadDecay)
        return false;
    return (measuredTauLepton1.pt() > measuredTauLepton2.pt());
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::run(
    const std::vector<fastmtt::MeasuredTauLepton> &measuredTauLeptons,
    const double &measuredMETx, const double &measuredMETy,
    const TMatrixD &covMET) {

    bestP4 = LorentzVector();
    ////////////////////////////////////////////

    if (measuredTauLeptons.size() != 2) {
        std::cout << "Number of MeasuredTauLepton is "
                  << measuredTauLeptons.size()
                  << " a user shouls pass exactly two leptons." << std::endl;
        return;
    }

    std::vector<fastmtt::MeasuredTauLepton> sortedMeasuredTauLeptons =
        measuredTauLeptons;
    std::sort(sortedMeasuredTauLeptons.begin(), sortedMeasuredTauLeptons.end(),
              compareLeptons);

    double metLength =
        sqrt(std::pow(measuredMETx, 2) + std::pow(measuredMETy, 2));
    LorentzVector aMET =
        LorentzVector(measuredMETx, measuredMETy, 0, metLength);

    const fastmtt::MeasuredTauLepton &aLepton1 = measuredTauLeptons[0];
    const fastmtt::MeasuredTauLepton &aLepton2 = measuredTauLeptons[1];

    myLikelihood.setLeptonInputs(aLepton1.p4(), aLepton2.p4(), aLepton1.type(),
                                 aLepton2.type(), aLepton1.decayMode(),
                                 aLepton2.decayMode());
    myLikelihood.setMETInputs(aMET, covMET);

    scan();
    // minimize();

    tau1P4 = aLepton1.p4() * (1.0 / minimumPosition[0]);
    tau2P4 = aLepton2.p4() * (1.0 / minimumPosition[1]);
    bestP4 = tau1P4 + tau2P4;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::minimize() {

    clock.Reset();
    clock.Start("minimize");

    minimizer->SetVariableLimits(0, 0.01, 1.0);
    minimizer->SetVariableLimits(1, 0.01, 1.0);

    minimizer->SetVariable(0, "x1", 0.5, 0.1);
    minimizer->SetVariable(1, "x2", 0.5, 0.1);

    minimizer->Minimize();

    const double *theMinimum = minimizer->X();
    minimumPosition[0] = theMinimum[0];
    minimumPosition[1] = theMinimum[1];
    minimumValue = minimizer->MinValue();

    if (true || minimizer->Status() != 0) {
        std::cout << " minimizer "
                  << " Status: " << minimizer->Status()
                  << " nCalls: " << minimizer->NCalls()
                  << " nIterations: " << minimizer->NIterations()
                  << " x1Max: " << theMinimum[0] << " x2Max: " << theMinimum[1]
                  << " max LLH: " << minimizer->MinValue()
                  << " m: " << bestP4.M() << std::endl;
    }
    clock.Stop("minimize");
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::scan() {

    clock.Reset();
    clock.Start("scan");

    double lh = 0.0;
    double bestLH = 0.0;

    double x[2] = {0.5, 0.5};

    double theMinimum[2] = {0.75, 0.75};
    const int nGridPoints = 100;
    const double gridFactor = 1. / nGridPoints;
    int nCalls = 0;
    for (int iX2 = 1; iX2 < nGridPoints; ++iX2) {
        x[1] = iX2 * gridFactor;
        for (int iX1 = 1; iX1 < nGridPoints; ++iX1) {
            x[0] = iX1 * gridFactor;
            lh = myLikelihood.value(x);
            ++nCalls;
            if (lh < bestLH) {
                bestLH = lh;
                theMinimum[0] = x[0];
                theMinimum[1] = x[1];
            }
        }

        minimumPosition[0] = theMinimum[0];
        minimumPosition[1] = theMinimum[1];
        minimumValue = bestLH;
    }
    clock.Stop("scan");
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getCpuTime(const std::string &method) {

    return clock.GetCpuTime(method.c_str());
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getRealTime(const std::string &method) {

    return clock.GetRealTime(method.c_str());
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
std::tuple<double, double> FastMTT::getBestX() const {

    return std::make_tuple(minimumPosition[0], minimumPosition[1]);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getBestLikelihood() const { return minimumValue; }
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getLikelihoodForX(double *x) const {

    return myLikelihood.value(x);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
