#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "RecoilCorrections/MetSystematics.cxx"
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
 * @brief function used to calculate the GenBosonVector and the
visibleGenBosonVector for an event.

The meaning of the genparticle statusflag codes is listed in the table below.

 Meaning                             | Value | Bit (value used in the config)
-------------------------------------|-------|-------
 isPrompt                            | 1     | 0
 isDecayedLeptonHadron               | 2     | 1
 isTauDecayProduct                   | 4     | 2
 isPromptTauDecayProduct             | 8     | 3
 isDirectTauDecayProduct             | 16    | 4
 isDirectPromptTauDecayProduct       | 32    | 5
 isDirectHadronDecayProduct          | 64    | 6
 isHardProcess                       | 128   | 7
 fromHardProcess                     | 256   | 8
 isHardProcessTauDecayProduct        | 512   | 9
 isDirectHardProcessTauDecayProduct  | 1024  | 10
 fromHardProcessBeforeFSR            | 2048  | 11
 isFirstCopy                         | 4096  | 12
 isLastCopy                          | 8192  | 13
 isLastCopyBeforeFSR                 | 16384 | 14

 *
 * @param df the input dataframe
 * @param genparticle_pt genparticle pt
 * @param genparticle_eta genparticle eta
 * @param genparticle_phi genparticle phi
 * @param genparticle_mass genparticle mass
 * @param genparticle_id genparticle PDG ID
 * @param genparticle_status genparticle status
 * @param genparticle_statusflag genparticle statusflag bit (see above)
 * @param outputname name of the new column containing the corrected met
 * @return a new dataframe containing a pair of lorentz vectors, first is the
GenBosonVector, second is the visibleGenBosonVector
 */
auto calculateGenBosonVector(auto df, const std::string &genparticle_pt,
                             const std::string &genparticle_eta,
                             const std::string &genparticle_phi,
                             const std::string &genparticle_mass,
                             const std::string &genparticle_id,
                             const std::string &genparticle_status,
                             const std::string &genparticle_statusflag,
                             const std::string outputname) {
    auto calculateGenBosonVector =
        [](const ROOT::RVec<float> &genparticle_pt,
           const ROOT::RVec<float> &genparticle_eta,
           const ROOT::RVec<float> &genparticle_phi,
           const ROOT::RVec<float> &genparticle_mass,
           const ROOT::RVec<int> &genparticle_id,
           const ROOT::RVec<int> &genparticle_status,
           const ROOT::RVec<int> &genparticle_statusflag) {
            ROOT::Math::PtEtaPhiMVector genparticle;
            float genPx = 0.; // generator Z(W) px
            float genPy = 0.; // generator Z(W) py
            float visPx = 0.; // visible (generator) Z(W) px
            float visPy = 0.; // visible (generator) Z(W) py
            // now loop though all genparticles in the event
            for (std::size_t index = 0; index < genparticle_id.size();
                 ++index) {
                // consider a genparticle,
                // 1. if it is a lepton and fromHardProcessFinalState --> bit 8
                // from statusflag and 1 from status
                // 2. if it is isDirectHardProcessTauDecayProduct --> bit 10
                // in statusflag
                Logger::get("getGenMet")
                    ->debug("Checking particle {} ", genparticle_id.at(index));
                if ((abs(genparticle_id.at(index)) >= 11 &&
                     abs(genparticle_id.at(index)) <= 16 &&
                     (IntBits(genparticle_status.at(index)).test(8)) &&
                     genparticle_status.at(index) == 1) ||
                    (IntBits(genparticle_status.at(index)).test(10))) {
                    Logger::get("getGenMet")->debug("Adding to gen p*");
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
                        Logger::get("getGenMet")->debug("Adding to vis p*");
                        visPx += genparticle.Px();
                        visPy += genparticle.Py();
                    }
                }
            }
            ROOT::Math::PtEtaPhiMVector genBoson;
            ROOT::Math::PtEtaPhiMVector visgenBoson;

            genBoson.SetPxPyPzE(genPx, genPy, 0,
                                std::sqrt(genPx * genPx + genPy * genPy));
            visgenBoson.SetPxPyPzE(visPx, visPy, 0,
                                   std::sqrt(visPx * visPx + visPy * visPy));
            std::pair<ROOT::Math::PtEtaPhiMVector, ROOT::Math::PtEtaPhiMVector>
                metpair = {genBoson, visgenBoson};
            return metpair;
        };
    return df.Define(outputname, calculateGenBosonVector,
                     {genparticle_pt, genparticle_eta, genparticle_phi,
                      genparticle_mass, genparticle_id, genparticle_status,
                      genparticle_statusflag});
}
/**
 * @brief Function used to propagate lepton corrections to the met. If the
 energy of a lepton is corrected (via some scale factor) or due to a shift,
 this change in energy has to be propagated to the met vector, and the met
 vector has to be adapted accordingly. The met is recalculated via
 @code
  Recalculate Met with corrected lepton energies :
  MetX_corrected = MetX + Px - Px_corrected
  MetY_corrected = MetY + Py - Py_corrected
  Met_corrected = sqrt(MetX_corrected * MetX_corrected + MetY_corrected *
 MetY_corrected)
 @endcode
 * @param df the input dataframe
 * @param met the uncorrected met lorentz vector
 * @param p4_1_uncorrected the uncorrected lorentz vector of the first
 lepton
 * @param p4_2_uncorrected the uncorrected lorentz vector of the second
 lepton
 * @param p4_1 the corrected lorentz vector of the first lepton
 * @param p4_2 the corrected lorentz vector of the second lepton
 * @param outputname name of the column containing the corrected met lorentz
 * @param apply_propagation if bool is set, the propagation is applied, if
 not, the outputcolumn contains the original met value vector
 * @return a new df containing the corrected met lorentz vector
 */
auto propagateLeptonsToMet(auto df, const std::string &met,
                           const std::string &p4_1_uncorrected,
                           const std::string &p4_2_uncorrected,
                           const std::string &p4_1, const std::string &p4_2,
                           const std::string &outputname,
                           bool apply_propagation) {
    auto scaleMet = [](const ROOT::Math::PtEtaPhiMVector &met,
                       const ROOT::Math::PtEtaPhiMVector &uncorrected_object,
                       const ROOT::Math::PtEtaPhiMVector &corrected_object) {
        // We propagate the lepton corrections to the Met by scaling the x
        // and y component of the Met according to the correction of the
        // lepton Recalculate Met with corrected lepton energies :
        // MetX_corrected = MetX + Px - Px_corrected
        // MetY_corrected = MetY + Py - Py_corrected
        // Met_corrected = sqrt(MetX_corrected * MetX_corrected +
        // MetY_corrected
        // * MetY_corrected)
        float corr_x = uncorrected_object.Px() - corrected_object.Px();
        float corr_y = uncorrected_object.Py() - corrected_object.Py();
        float MetX = met.Px() + corr_x;
        float MetY = met.Py() + corr_y;
        Logger::get("propagateLeptonsToMet")->debug("corr_x {}", corr_x);
        Logger::get("propagateLeptonsToMet")->debug("corr_y {}", corr_y);
        Logger::get("propagateLeptonsToMet")->debug("MetX {}", MetX);
        Logger::get("propagateLeptonsToMet")->debug("MetY {}", MetY);
        ROOT::Math::PtEtaPhiMVector corrected_met;
        corrected_met.SetPxPyPzE(MetX, MetY, 0,
                                 std::sqrt(MetX * MetX + MetY * MetY));
        Logger::get("propagateLeptonsToMet")
            ->debug("corrected_object pt - {}", corrected_object.Pt());
        Logger::get("propagateLeptonsToMet")
            ->debug("uncorrected_object pt - {}", uncorrected_object.Pt());
        Logger::get("propagateLeptonsToMet")->debug("old met {}", met.Pt());
        Logger::get("propagateLeptonsToMet")
            ->debug("corrected met {}", corrected_met.Pt());
        return corrected_met;
    };
    if (apply_propagation) {
        // first correct for the first lepton, store the met in an
        // intermediate column
        Logger::get("propagateLeptonsToMet")
            ->debug("Setting up correction for first lepton {}", p4_1);
        auto df1 = df.Define(outputname + "_intermediate", scaleMet,
                             {met, p4_1_uncorrected, p4_1});
        // after the second lepton correction, the correct output column is
        // used
        Logger::get("propagateLeptonsToMet")
            ->debug("Setting up correction for second lepton {}", p4_2);
        return df1.Define(
            outputname, scaleMet,
            {outputname + "_intermediate", p4_2_uncorrected, p4_2});
    } else {
        // if we do not apply the propagation, just rename the met column to
        // the new outputname and dont change anything else
        return basefunctions::rename<ROOT::Math::PtEtaPhiMVector>(df, met,
                                                                  outputname);
    }
}
/**
 * @brief Function used to propagate jet corrections to the met. If the
 energy of a jet is corrected, this
 change in energy has to be propagated to the met vector, and the met vector
 has to be adapted accordingly. The met is recalculated via
 @code
  Recalculate Met with corrected jet energies :
  MetX_corrected = MetX + Px - Px_corrected
  MetY_corrected = MetY + Py - Py_corrected
  Met_corrected = sqrt(MetX_corrected * MetX_corrected + MetY_corrected *
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
  * @param apply_propagation if bool is set, the propagation is applied, if
 not, the outputcolumn contains the original met value
 * @param min_jet_pt minimal pt, the corrected jet has to have, in order for
 the met propagation to be applied
 * @return a new df containing the corrected met lorentz vector
 */
auto propagateJetsToMet(auto df, const std::string &met,
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
    auto scaleMet = [min_jet_pt](const ROOT::Math::PtEtaPhiMVector &met,
                                 const ROOT::RVec<float> &jet_pt_corrected,
                                 const ROOT::RVec<float> &jet_eta_corrected,
                                 const ROOT::RVec<float> &jet_phi_corrected,
                                 const ROOT::RVec<float> &jet_mass_corrected,
                                 const ROOT::RVec<float> &jet_pt,
                                 const ROOT::RVec<float> &jet_eta,
                                 const ROOT::RVec<float> &jet_phi,
                                 const ROOT::RVec<float> &jet_mass) {
        ROOT::Math::PtEtaPhiMVector corrected_met;
        ROOT::Math::PtEtaPhiMVector uncorrected_jet;
        ROOT::Math::PtEtaPhiMVector corrected_jet;
        float corr_x = 0.0;
        float corr_y = 0.0;
        // now loop through all jets in the event
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
        float MetX = met.Px() + corr_x;
        float MetY = met.Py() + corr_y;
        Logger::get("propagateJetsToMet")->debug("corr_x {}", corr_x);
        Logger::get("propagateJetsToMet")->debug("corr_y {}", corr_y);
        Logger::get("propagateJetsToMet")->debug("MetX {}", MetX);
        Logger::get("propagateJetsToMet")->debug("MetY {}", MetY);
        corrected_met.SetPxPyPzE(MetX, MetY, 0,
                                 std::sqrt(MetX * MetX + MetY * MetY));
        Logger::get("propagateJetsToMet")->debug("old met {}", met.Pt());
        Logger::get("propagateJetsToMet")
            ->debug("corrected met {}", corrected_met.Pt());
        return corrected_met;
    };
    if (apply_propagation) {
        return df.Define(outputname, scaleMet,
                         {met, jet_pt_corrected, jet_eta_corrected,
                          jet_phi_corrected, jet_mass_corrected, jet_pt,
                          jet_eta, jet_phi, jet_mass});
    } else {
        // if we do not apply the propagation, just rename the met column to
        // the new outputname and dont change anything else
        return basefunctions::rename<ROOT::Math::PtEtaPhiMVector>(df, met,
                                                                  outputname);
    }
}
/**
 * @brief function used to apply Recoil corrections on a given sample. For
more information on recoil corrections, check [this
repo](https://github.com/KIT-CMS/RecoilCorrections/). The code for the
application of these corrections can be found in `src/RecoilCorrections`.
 *
 * @param df the input dataframe
 * @param met the input met
 * @param genboson the genboson vector, this is a std::pair where the first
value is the genboson vector and the seconds value is the visible genboson
vector
 * @param jet_pt the pt of all jets
 * @param outputname name of the new column containing the corrected met
 * @param recoilfile path to the recoil corrections file
 * @param systematicsfile path to the systematics corrections file
 * @param applyRecoilCorrections if bool is set, the recoil correction is
applied, if not, the outputcolumn contains the original met value
 * @param resolution bool - if set, resolution corrections are applied
 * @param response bool - if set, response corrections are applied
 * @param shiftUp bool - if set, the up shift is applied
 * @param shiftDown bool - if set, the down shift is applied
 * @param isWjets bool - if set, the number of jets is incrased by one (this
is only needed for WJets samples)
 * @return a new dataframe containing the new met column
 */
auto applyRecoilCorrections(
    auto df, const std::string &met, const std::string &genmet,
    const std::string &jet_pt, const std::string &outputname,
    const std::string &recoilfile, const std::string &systematicsfile,
    bool applyRecoilCorrections, bool resolution, bool response, bool shiftUp,
    bool shiftDown, bool isWjets) {
    if (applyRecoilCorrections) {
        Logger::get("RecoilCorrections")->debug("Will run recoil corrections");
        const auto corrector = new RecoilCorrector(recoilfile);
        const auto systematics = new MetSystematic(systematicsfile);
        auto shiftType = MetSystematic::SysShift::Nominal;
        if (shiftUp) {
            shiftType = MetSystematic::SysShift::Up;
        } else if (shiftDown) {
            shiftType = MetSystematic::SysShift::Down;
        }
        auto sysType = MetSystematic::SysType::None;
        if (response) {
            sysType = MetSystematic::SysType::Response;
        } else if (resolution) {
            sysType = MetSystematic::SysType::Resolution;
        }
        auto RecoilCorrections = [sysType, systematics, shiftType, corrector,
                                  isWjets](
                                     ROOT::Math::PtEtaPhiMVector &met,
                                     std::pair<ROOT::Math::PtEtaPhiMVector,
                                               ROOT::Math::PtEtaPhiMVector>
                                         &genboson,
                                     const ROOT::RVec<float> &jet_pt) {
            // TODO is this the correct number of jets ?
            auto jets_above30 =
                Filter(jet_pt, [](float pt) { return pt > 30; });
            int nJets30 = jets_above30.size();
            if (isWjets) {
                nJets30 = nJets30 + 1;
            }
            float MetX = met.Px();
            float MetY = met.Py();
            float correctedMetX = 0.;
            float correctedMetY = 0.;
            ROOT::Math::PtEtaPhiMVector genparticle;
            ROOT::Math::PtEtaPhiMVector corrected_met;
            float genPx = genboson.first.Px();  // generator Z(W) px
            float genPy = genboson.first.Py();  // generator Z(W) py
            float visPx = genboson.second.Px(); // visible (generator) Z(W) px
            float visPy = genboson.second.Py(); // visible (generator) Z(W) py
            Logger::get("RecoilCorrections")->debug("Corrector Inputs");
            Logger::get("RecoilCorrections")->debug("nJets30 {} ", nJets30);
            Logger::get("RecoilCorrections")->debug("genPx {} ", genPx);
            Logger::get("RecoilCorrections")->debug("genPy {} ", genPy);
            Logger::get("RecoilCorrections")->debug("visPx {} ", visPx);
            Logger::get("RecoilCorrections")->debug("visPy {} ", visPy);
            Logger::get("RecoilCorrections")->debug("MetX {} ", MetX);
            Logger::get("RecoilCorrections")->debug("MetY {} ", MetY);
            Logger::get("RecoilCorrections")
                ->debug("correctedMetX {} ", correctedMetX);
            Logger::get("RecoilCorrections")
                ->debug("correctedMetY {} ", correctedMetY);
            Logger::get("RecoilCorrections")->debug("old met {} ", met.Pt());
            corrector->CorrectWithHist(MetX, MetY, genPx, genPy, visPx, visPy,
                                       nJets30, correctedMetX, correctedMetY);
            // only apply shifts if the correpsonding variables are set
            if (sysType != MetSystematic::SysType::None &&
                shiftType != MetSystematic::SysShift::Nominal) {
                Logger::get("RecoilCorrections")
                    ->debug(" apply systematics {} {}", sysType, shiftType);
                systematics->ApplyMetSystematic(
                    correctedMetX, correctedMetY, genPx, genPy, visPx, visPy,
                    nJets30, sysType, shiftType, correctedMetX, correctedMetY);
            }
            corrected_met.SetPxPyPzE(correctedMetX, correctedMetY, 0,
                                     std::sqrt(correctedMetX * correctedMetX +
                                               correctedMetY * correctedMetY));
            Logger::get("RecoilCorrections")
                ->debug("shifted and corrected met {} ", corrected_met.Pt());

            return corrected_met;
        };
        return df.Define(outputname, RecoilCorrections, {met, genmet, jet_pt});
    } else {
        // if we do not apply the recoil corrections, just rename the met
        // column to the new outputname and dont change anything else
        return basefunctions::rename<ROOT::Math::PtEtaPhiMVector>(df, met,
                                                                  outputname);
    }
}
} // end namespace met
