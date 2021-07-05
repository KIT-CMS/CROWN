#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "RecoilCorrections/RecoilCorrector.cxx"
#include "basefunctions.hxx"
#include "bitset"
#include "utility/Logger.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>

typedef std::bitset<20> IntBits;

namespace met {
/**
 * @brief Function to get pT from the met Lorentz vector
 *
 * @param df the input dataframe
 * @param metvector the name of the column containing the met lorentz vector
 * @param outputname name of the new column containing the met pT value
 * @return a new dataframe containing the met pT value
 */
auto metPt(auto df, const std::string &metvector,
           const std::string &outputname) {
    auto GetMetPt = [](const ROOT::Math::PtEtaPhiEVector &met) {
        return met.Pt();
    };
    return df.Define(outputname, GetMetPt, {metvector});
}
/**
 * @brief Function to get Phi from the met Lorentz vector
 *
 * @param df the input dataframe
 * @param metvector the name of the column containing the met lorentz vector
 * @param outputname name of the new column containing the met phi value
 * @return a new dataframe containing the met phi value
 */
auto metPhi(auto df, const std::string &metvector,
            const std::string &outputname) {
    auto GetMetPhi = [](const ROOT::Math::PtEtaPhiEVector met) {
        return met.Phi();
    };
    return df.Define(outputname, GetMetPhi, {metvector});
}
/**
 * @brief Function used to propagate lepton corrections to the met. If the
 energy of a lepton is corrected (via some scale factor) or due to a shift, this
 change in energy has to be propagated to the met vector, and the met vector has
 to be adapted accordingly. The met is recalculated via
 @code
  Recalculate MET with corrected lepton energies :
  MetX_corrected = MetX + Px - Px_corrected
  MetY_corrected = MetY + Py - Py_corrected
  MET_corrected = sqrt(MetX_corrected * MetX_corrected + MetY_corrected *
 MetY_corrected)
 @endcode
 * @param df the input dataframe
 * @param met the uncorrected met lorentz vector
 * @param p4_1_uncorrected the uncorrected lorentz vector of the first lepton
 * @param p4_2_uncorrected the uncorrected lorentz vector of the second lepton
 * @param p4_1 the corrected lorentz vector of the first lepton
 * @param p4_2 the corrected lorentz vector of the second lepton
 * @param outputname name of the column containing the corrected met lorentz
 * @param apply_propagation if bool is set, the propagation is applied, if not,
 the outputcolumn contains the original met value vector
 * @return a new df containing the corrected met lorentz vector
 */
auto propagateLeptonsToMET(auto df, const std::string &met,
                           const std::string &p4_1_uncorrected,
                           const std::string &p4_2_uncorrected,
                           const std::string &p4_1, const std::string &p4_2,
                           const std::string &outputname,
                           bool apply_propagation) {
    auto scaleMet = [](const ROOT::Math::PtEtaPhiEVector &met,
                       const ROOT::Math::PtEtaPhiMVector &uncorrected_object,
                       const ROOT::Math::PtEtaPhiMVector &corrected_object) {
        // We propagate the lepton corrections to the MET by scaling the x and y
        // component of the MET according to the correction of the lepton
        // Recalculate MET with corrected lepton energies :
        // MetX_corrected = MetX + Px - Px_corrected
        // MetY_corrected = MetY + Py - Py_corrected
        // MET_corrected = sqrt(MetX_corrected * MetX_corrected + MetY_corrected
        // * MetY_corrected)
        float corr_x = uncorrected_object.Px() - corrected_object.Px();
        float corr_y = uncorrected_object.Py() - corrected_object.Py();
        float MetX = met.Px() - corr_x;
        float MetY = met.Py() - corr_y;
        Logger::get("propagateLeptonsToMET")->debug("corr_x {}", corr_x);
        Logger::get("propagateLeptonsToMET")->debug("corr_y {}", corr_y);
        Logger::get("propagateLeptonsToMET")->debug("MetX {}", MetX);
        Logger::get("propagateLeptonsToMET")->debug("MetY {}", MetY);
        ROOT::Math::PtEtaPhiEVector corrected_met;
        corrected_met.SetPxPyPzE(MetX, MetY, 0,
                                 std::sqrt(MetX * MetX + MetY * MetY));
        Logger::get("propagateLeptonsToMET")
            ->debug("corrected_object pt - {}", corrected_object.Pt());
        Logger::get("propagateLeptonsToMET")
            ->debug("uncorrected_object pt - {}", uncorrected_object.Pt());
        Logger::get("propagateLeptonsToMET")->debug("old met {}", met.Pt());
        Logger::get("propagateLeptonsToMET")
            ->debug("corrected met {}", corrected_met.Pt());
        return corrected_met;
    };
    if (apply_propagation) {
        // first correct for the first lepton, store the met in an
        // intermediate column
        Logger::get("propagateLeptonsToMET")
            ->debug("Setting up correction for first lepton {}", p4_1);
        auto df1 = df.Define(outputname + "_intermidiate", scaleMet,
                             {met, p4_1_uncorrected, p4_1});
        // after the second lepton correction, the correct output column is used
        Logger::get("propagateLeptonsToMET")
            ->debug("Setting up correction for second lepton {}", p4_2);
        return df1.Define(
            outputname, scaleMet,
            {outputname + "_intermidiate", p4_2_uncorrected, p4_2});
    } else {
        // if we do not apply the propagation, just rename the met column to the
        // new outputname and dont change anything else
        return basefunctions::rename<ROOT::Math::PtEtaPhiEVector>(df, met,
                                                                  outputname);
    }
}
/**
 * @brief Function used to propagate jet corrections to the met. If the
 energy of a jet is corrected, this
 change in energy has to be propagated to the met vector, and the met vector has
 to be adapted accordingly. The met is recalculated via
 @code
  Recalculate MET with corrected jet energies :
  MetX_corrected = MetX + Px - Px_corrected
  MetY_corrected = MetY + Py - Py_corrected
  MET_corrected = sqrt(MetX_corrected * MetX_corrected + MetY_corrected *
 MetY_corrected)
 @endcode
 The correction is done for all valid jets.
 * @param df the input dataframe
 * @param met the uncorrected met lorentz vector
 * @param jet_pt_corrected pt of the corrected jet
 * @param jet_eta_corrected eta of the corrected jet
 * @param jet_phi_corrected phi of the corrected jet
 * @param jet_mass_corrected mass of the corrected jet
 * @param jet_pt pt of the uncorrected jet
 * @param jet_eta eta of the uncorrected jet
 * @param jet_phi phi of the uncorrected jet
 * @param jet_mass mass of the uncorrected jet
 * @param outputname name of the column containing the corrected met lorentz
 vector
  * @param apply_propagation if bool is set, the propagation is applied, if not,
 the outputcolumn contains the original met value
 * @param min_jet_pt minimal pt, the corrected jet has to have, in order for the
 met propagation to be applied
 * @return a new df containing the corrected met lorentz vector
 */
auto propagateJetsToMET(auto df, const std::string &met,
                        const std::string &jet_pt_corrected,
                        const std::string &jet_eta_corrected,
                        const std::string &jet_phi_corrected,
                        const std::string &jet_mass_corrected,
                        const std::string &jet_pt, const std::string &jet_eta,
                        const std::string &jet_phi, const std::string &jet_mass,
                        const std::string &outputname, bool apply_propagation,
                        float min_jet_pt) {
    // propagate jet corrections to met, since we can have an arbitrary
    // amount of jets, this has to be done per event
    auto scaleMet = [min_jet_pt](const ROOT::Math::PtEtaPhiEVector &met,
                                 const ROOT::RVec<float> &jet_pt_corrected,
                                 const ROOT::RVec<float> &jet_eta_corrected,
                                 const ROOT::RVec<float> &jet_phi_corrected,
                                 const ROOT::RVec<float> &jet_mass_corrected,
                                 const ROOT::RVec<float> &jet_pt,
                                 const ROOT::RVec<float> &jet_eta,
                                 const ROOT::RVec<float> &jet_phi,
                                 const ROOT::RVec<float> &jet_mass) {
        ROOT::Math::PtEtaPhiEVector corrected_met;
        ROOT::Math::PtEtaPhiMVector uncorrected_jet;
        ROOT::Math::PtEtaPhiMVector corrected_jet;
        float corr_x = 0.0;
        float corr_y = 0.0;
        // now loop though all jets in the event
        for (std::size_t index = 0; index < jet_pt.size(); ++index) {
            // only propagate jets above the given pt threshold
            if (jet_pt_corrected.at(index) > min_jet_pt) {
                // construct the uncorrected and the corrected lorentz
                // vectors
                uncorrected_jet = ROOT::Math::PtEtaPhiMVector(
                    jet_pt_corrected.at(index), jet_eta_corrected.at(index),
                    jet_phi_corrected.at(index), jet_mass_corrected.at(index));
                corrected_jet = ROOT::Math::PtEtaPhiMVector(
                    jet_pt.at(index), jet_eta.at(index), jet_phi.at(index),
                    jet_mass.at(index));
                // update the correction factors that are applied to the met
                corr_x += uncorrected_jet.Px() - corrected_jet.Px();
                corr_y += uncorrected_jet.Py() - corrected_jet.Py();
            }
        }
        float MetX = met.Px() - corr_x;
        float MetY = met.Py() - corr_y;
        Logger::get("propagateJetsToMET")->debug("corr_x {}", corr_x);
        Logger::get("propagateJetsToMET")->debug("corr_y {}", corr_y);
        Logger::get("propagateJetsToMET")->debug("MetX {}", MetX);
        Logger::get("propagateJetsToMET")->debug("MetY {}", MetY);
        corrected_met.SetPxPyPzE(MetX, MetY, 0,
                                 std::sqrt(MetX * MetX + MetY * MetY));
        Logger::get("propagateJetsToMET")->debug("old met {}", met.Pt());
        Logger::get("propagateJetsToMET")
            ->debug("corrected met {}", corrected_met.Pt());
        return corrected_met;
    };
    if (apply_propagation) {
        return df.Define(outputname, scaleMet,
                         {met, jet_pt_corrected, jet_eta_corrected,
                          jet_phi_corrected, jet_mass_corrected, jet_pt,
                          jet_eta, jet_phi, jet_mass});
    } else {
        // if we do not apply the propagation, just rename the met column to the
        // new outputname and dont change anything else
        return basefunctions::rename<ROOT::Math::PtEtaPhiEVector>(df, met,
                                                                  outputname);
    }
}
/**
 * @brief
 *
 * @param df
 * @param met
 * @param genparticle_pt
 * @param genparticle_eta
 * @param genparticle_phi
 * @param genparticle_mass
 * @param genparticle_id
 * @param genparticle_status
 * @param jet_pt
 * @param outputname
 * @param inputfile
 * @param applyRecoilCorrections
 * @return auto
 */
auto applyRecoilCorrections(
    auto df, const std::string &met, const std::string &genparticle_pt,
    const std::string &genparticle_eta, const std::string &genparticle_phi,
    const std::string &genparticle_mass, const std::string &genparticle_id,
    const std::string &genparticle_status, const std::string &jet_pt,
    const std::string &outputname, const std::string &inputfile,
    bool applyRecoilCorrections) {
    if (applyRecoilCorrections) {
        Logger::get("applyRecoilCorrections")
            ->debug("Will run recoil corrections");
        const auto corrector = new RecoilCorrector(inputfile);
        auto RecoilCorrections = [corrector](
                                     ROOT::Math::PtEtaPhiEVector &met,
                                     const ROOT::RVec<float> &genparticle_pt,
                                     const ROOT::RVec<float> &genparticle_eta,
                                     const ROOT::RVec<float> &genparticle_phi,
                                     const ROOT::RVec<float> &genparticle_mass,
                                     const ROOT::RVec<int> &genparticle_id,
                                     const ROOT::RVec<int> &genparticle_status,
                                     const ROOT::RVec<float> &jet_pt) {
            // TODO is this the correct number of jets ?
            auto jets_above30 =
                Filter(jet_pt, [](float pt) { return pt > 30; });
            int nJets30 = jets_above30.size() + 1;
            float genPx = 0.; // generator Z(W) px
            float genPy = 0.; // generator Z(W) py
            float visPx = 0.; // visible (generator) Z(W) px
            float visPy = 0.; // visible (generator) Z(W) py
            float MetX = met.Px();
            float MetY = met.Py();
            float correctedMetX = 0.;
            float correctedMetY = 0.;
            ROOT::Math::PtEtaPhiMVector genparticle;
            ROOT::Math::PtEtaPhiEVector corrected_met;
            // now loop though all genparticles in the event
            for (std::size_t index = 0; index < genparticle_id.size();
                 ++index) {
                // consider a genparticle,
                // 1. if it is a lepton and fromHardProcess --> bit 8 from
                // status
                // 2. if it is isDirectHardProcessTauDecayProduct --> bit 10
                // in status
                Logger::get("applyRecoilCorrections")
                    ->debug("Checking particle {} ", genparticle_id.at(index));
                if ((abs(genparticle_id.at(index)) >= 11 &&
                     abs(genparticle_id.at(index)) <= 16) ||
                    (IntBits(genparticle_status.at(index)).test(10))) {
                    Logger::get("applyRecoilCorrections")
                        ->debug("Adding to gen p*");
                    genparticle = ROOT::Math::PtEtaPhiMVector(
                        genparticle_pt.at(index), genparticle_eta.at(index),
                        genparticle_phi.at(index), genparticle_mass.at(index));
                    genPx += genparticle.Px();
                    genPy += genparticle.Py();
                    // if the genparticle is no neutrino, we add the x and y
                    // component to the visible generator component as well
                    if (abs(genparticle_id.at(index)) != 12 &&
                        abs(genparticle_id.at(index)) != 14 &&
                        abs(genparticle_id.at(index)) != 16) {
                        Logger::get("applyRecoilCorrections")
                            ->debug("Adding to vis p*");
                        visPx += genparticle.Px();
                        visPy += genparticle.Py();
                    }
                }
            };
            Logger::get("applyRecoilCorrections")->debug("Corrector Inputs");
            Logger::get("applyRecoilCorrections")
                ->debug("nJets30 {} ", nJets30);
            Logger::get("applyRecoilCorrections")->debug("genPx {} ", genPx);
            Logger::get("applyRecoilCorrections")->debug("genPy {} ", genPy);
            Logger::get("applyRecoilCorrections")->debug("visPx {} ", visPx);
            Logger::get("applyRecoilCorrections")->debug("visPy {} ", visPy);
            Logger::get("applyRecoilCorrections")->debug("MetX {} ", MetX);
            Logger::get("applyRecoilCorrections")->debug("MetY {} ", MetY);
            corrector->CorrectWithHist(MetX, MetY, genPx, genPy, visPx, visPy,
                                       nJets30, correctedMetX, correctedMetY);
            Logger::get("applyRecoilCorrections")
                ->debug("correctedMetX {} ", correctedMetX);
            Logger::get("applyRecoilCorrections")
                ->debug("correctedMetY {} ", correctedMetY);
            corrected_met.SetPxPyPzE(correctedMetX, correctedMetY, 0,
                                     std::sqrt(correctedMetX * correctedMetX +
                                               correctedMetY * correctedMetY));
            Logger::get("applyRecoilCorrections")
                ->debug("old met {} ", met.Pt());
            Logger::get("applyRecoilCorrections")
                ->debug("corrected met {} ", corrected_met.Pt());
            return corrected_met;
        };
        return df.Define(outputname, RecoilCorrections,
                         {met, genparticle_pt, genparticle_eta, genparticle_phi,
                          genparticle_mass, genparticle_id, genparticle_status,
                          jet_pt});
    } else {
        // if we do not apply the recoil corrections, just rename the met
        // column to the new outputname and dont change anything else
        return basefunctions::rename<ROOT::Math::PtEtaPhiEVector>(df, met,
                                                                  outputname);
    }
}
} // end namespace met
