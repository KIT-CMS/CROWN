#ifndef FastMTT_FastMTT_H
#define FastMTT_FastMTT_H

#include "Math/LorentzVector.h"
#include "TBenchmark.h"
#include "TMatrixD.h"
#include <Math/Vector4D.h>
#include <bitset>
#include <string>
#include <tuple>
#include <vector>

namespace fastmtt {
class MeasuredTauLepton;
}

namespace ROOT {
namespace Math {
class Minimizer;
class Functor;
} // namespace Math
} // namespace ROOT
class TVector2;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;

namespace fastMTT {
double likelihoodFunc(double *x, double *par);

enum likelihoodComponent { MET, MASS, PX, PY, ENERGY, IP };
} // namespace fastMTT

class Likelihood {

  public:
    Likelihood();

    ~Likelihood();

    double value(const double *x) const;

    void setLeptonInputs(const LorentzVector &aLeg1P4,
                         const LorentzVector &aLeg2P4, int aLeg1DecayType,
                         int aLeg2DecayType, int aLeg1DecayMode,
                         int aLeg2DecayMode);

    void setMETInputs(const LorentzVector &aMET, const TMatrixD &aCovMET);

    void setParameters(const std::vector<double> &parameters);

    void enableComponent(fastMTT::likelihoodComponent aCompIndex);

    void disableComponent(fastMTT::likelihoodComponent aCompIndex);

    double massLikelihood(const double &m) const;

    double ptLikelihood(const double &pTTauTau, int type) const;

    double metTF(const LorentzVector &metP4, const LorentzVector &nuP4,
                 const TMatrixD &covMET) const;

  private:
    LorentzVector leg1P4, leg2P4;
    LorentzVector recoMET;

    TMatrixD covMET;

    double mVis, mVisLeg1, mVisLeg2;

    int leg1DecayType, leg2DecayType;
    int leg1DecayMode, leg2DecayMode;

    std::vector<double> parameters;

    /// Bit word coding enabled likelihood components
    std::bitset<128> componentsBitWord;

    // precomputed values used to reduce the number of redundant calculations
    double mVis1OverTauSquare;
    double mVis2OverTauSquare;
    using PowTable = std::array<double, 5u>; // first powers of a number
    // array with dimensions 2x3x5 used to store the powers of the coefficients
    // of two vectors
    std::array<std::array<PowTable, 3u>, 2u> allpTpows;
};
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
class FastMTT {

  public:
    FastMTT();

    ~FastMTT();

    void initialize();

    /// Run fastMTT algorithm for given input
    ROOT::Math::PtEtaPhiMVector
    run(const std::vector<fastmtt::MeasuredTauLepton> &, const double &,
        const double &, const TMatrixD &);

  private:
    static bool
    compareLeptons(const fastmtt::MeasuredTauLepton &measuredTauLepton1,
                   const fastmtt::MeasuredTauLepton &measuredTauLepton2);

    /// Run a scan over x1 and x2 [0,1] rectangle for given inputs.
    /// Results are stored in internal variables accesed by
    /// relevant get methods.
    void scan(Likelihood &myLikelihood);

    /// Minimum location
    std::vector<double> minimumPosition;

    /// Mimimum value
    double minimumValue;

    /// Dimension of minimalization space
    unsigned int nVariables;

    /// Names of variables to be minimized
    std::vector<std::string> varNames;

    /// Values of variables to be minimized

    std::vector<double> variables;

    /// Step sizes for each minimized variable
    std::vector<double> stepSizes;

    LorentzVector tau1P4, tau2P4, bestP4;

    int verbosity;
};

#endif
