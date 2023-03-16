#include "../../include/SVFit/MeasuredTauLepton.hxx"
#include "../../include/SVFit/AuxFunctions.hxx"

#include <TMath.h>

using namespace fastmtt;

MeasuredTauLepton::MeasuredTauLepton()
    : type_(kUndefinedDecayType), pt_(0.), eta_(0.), phi_(0.), mass_(0.),
      decayMode_(-1) {
    initialize();
}

MeasuredTauLepton::MeasuredTauLepton(int type, double pt, double eta,
                                     double phi, double mass, int decayMode)
    : type_(type), pt_(pt), eta_(eta), phi_(phi), mass_(mass),
      decayMode_(decayMode) {
    double minVisMass = electronMass;
    double maxVisMass = tauLeptonMass;
    if (type_ == kTauToElecDecay) {
        minVisMass = electronMass;
        maxVisMass = minVisMass;
    } else if (type_ == kTauToMuDecay) {
        minVisMass = muonMass;
        maxVisMass = minVisMass;
    } else if (type_ == kTauToHadDecay) {
        if (decayMode_ == -1) {
            minVisMass = chargedPionMass;
            maxVisMass = 1.5;
        } else if (decayMode_ == 0) {
            minVisMass = chargedPionMass;
            maxVisMass = minVisMass;
        } else {
            minVisMass = 0.3;
            maxVisMass = 1.5;
        }
    }
    preciseVisMass_ = mass_;
    if (preciseVisMass_ < minVisMass)
        preciseVisMass_ = minVisMass;
    if (preciseVisMass_ > maxVisMass)
        preciseVisMass_ = maxVisMass;
    initialize();
}

MeasuredTauLepton::MeasuredTauLepton(const MeasuredTauLepton &measuredTauLepton)
    : type_(measuredTauLepton.type()), pt_(measuredTauLepton.pt()),
      eta_(measuredTauLepton.eta()), phi_(measuredTauLepton.phi()),
      mass_(measuredTauLepton.mass()),
      decayMode_(measuredTauLepton.decayMode()) {
    preciseVisMass_ = measuredTauLepton.mass();

    initialize();
}

MeasuredTauLepton::~MeasuredTauLepton() {}

int MeasuredTauLepton::type() const { return type_; }

double MeasuredTauLepton::pt() const { return pt_; }
double MeasuredTauLepton::eta() const { return eta_; }
double MeasuredTauLepton::phi() const { return phi_; }
double MeasuredTauLepton::mass() const { return preciseVisMass_; }

double MeasuredTauLepton::energy() const { return energy_; }
double MeasuredTauLepton::px() const { return px_; }
double MeasuredTauLepton::py() const { return py_; }
double MeasuredTauLepton::pz() const { return pz_; }

double MeasuredTauLepton::p() const { return p_; }

int MeasuredTauLepton::decayMode() const { return decayMode_; }

LorentzVector MeasuredTauLepton::p4() const { return p4_; }

Vector MeasuredTauLepton::p3() const { return p3_; }

double MeasuredTauLepton::cosPhi_sinTheta() const { return cosPhi_sinTheta_; }
double MeasuredTauLepton::sinPhi_sinTheta() const { return sinPhi_sinTheta_; }
double MeasuredTauLepton::cosTheta() const { return cosTheta_; }

void MeasuredTauLepton::roundToNdigits(unsigned int nDigis) {
    pt_ = fastmtt::roundToNdigits(pt_, nDigis);
    eta_ = fastmtt::roundToNdigits(eta_, nDigis);
    phi_ = fastmtt::roundToNdigits(phi_, nDigis);
    mass_ = fastmtt::roundToNdigits(mass_, nDigis);
    initialize();
}

void MeasuredTauLepton::initialize() {
    // CV: relations between pT and p, energy taken from
    // http://en.wikipedia.org/wiki/Pseudorapidity
    p_ = pt_ * TMath::CosH(eta_);
    px_ = pt_ * TMath::Cos(phi_);
    py_ = pt_ * TMath::Sin(phi_);
    pz_ = pt_ * TMath::SinH(eta_);
    energy_ = TMath::Sqrt(p_ * p_ + preciseVisMass_ * preciseVisMass_);
    p4_ = LorentzVector(px_, py_, pz_, energy_);
    p3_ = Vector(px_, py_, pz_);
    double theta = p4_.theta();
    cosPhi_sinTheta_ = TMath::Cos(phi_) * TMath::Sin(theta);
    sinPhi_sinTheta_ = TMath::Sin(phi_) * TMath::Sin(theta);
    cosTheta_ = TMath::Cos(theta);
}
