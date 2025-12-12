#ifndef GUARD_MET_H
#define GUARD_MET_H

#include "../include/RecoilCorrections/MetSystematics.hxx"
#include "../include/RecoilCorrections/RecoilCorrector.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "../include/event.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "bitset"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include "correction.h"

typedef std::bitset<20> IntBits;

namespace met {

/**
 * @brief This function is used to apply recoil corrections on a given sample.
 * It is only recommented to apply it to single boson samples, like Z, W or
 * Higgs boson. The corrections are provided by the HLepRare group. For more
 * information on how the recoil corrections were calculated and have to be used,
 * check out [this presentation](https://indico.cern.ch/event/1583951/contributions/6751916/attachments/3159171/5612627/HLepRare_25.10.22.pdf). 
 * and the HLepRare group [documentation](https://cms-higgs-leprare.docs.cern.ch/htt-common/V_recoil/).
 *
 * @note This function is intended to be used in Run3 where recoil corrections are
 * provided in json files to be read by correctionlib. For Run2 recoil correction
 * see the overloaded version of this function.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * recoil corrections file
 * @param outputname name of the new column containing the corrected MET Lorentz
 * vector
 * @param p4_met initial/uncorrected MET Lorentz vector
 * @param p4_gen_boson name of the column containing the generator-level boson
 * Lorentz vector
 * @param p4_vis_gen_boson name of the column containing the visible part of
 * the generator-level boson Lorentz vector
 * @param n_jets name of the column containing the number of good jets in an event 
 * @param corr_file path to the json file with the recoil corrections
 * @param corr_name name of the recoil correction, this is the first part of the
 * correction string in the json file (e.g. "Recoil_correction")
 * @param method method to be used to apply the corrections, possible options are
 * "Rescaling", "QuantileMapHist" and "Uncertainty" (second part of the correction string)
 * @param order order of the used DY samples: "LO" for madgraph, "NLO" for amc@nlo,
 * "NNLO" for powheg
 * @param variation name of the variation that should be evaluated, options are
 * "nom", "RespUp", "RespDown", "ResolUp", "ResolDown". This is only used if
 * `method` is set to "Uncertainty".
 * @param apply_correction if bool is set to true, the recoil correction is
 * applied, if not, the output column contains the original MET vector
 *
 * @return a new dataframe containing the new MET column
 */
ROOT::RDF::RNode RecoilCorrection(
    ROOT::RDF::RNode df, correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname,
    const std::string &p4_met, const std::string &p4_gen_boson,
    const std::string &p4_vis_gen_boson, const std::string &n_jets,
    const std::string &corr_file, const std::string &corr_name,
    const std::string &method, const std::string &order,
    const std::string &variation, bool apply_correction) {
    if (apply_correction) {
        Logger::get("met::RecoilCorrection")->debug("Will run recoil corrections with correctionlib");
        
        // Rescaling method is the only recommended method for outside
        // -150 GeV < UPara, UPerp < 150 GeV
        auto RecoilCorr = correction_manager.loadCorrection(corr_file, corr_name + "_" + method);

        auto Correction = [RecoilCorr, order, method, variation](
                            ROOT::Math::PtEtaPhiMVector &met,
                            ROOT::Math::PtEtaPhiMVector &gen_boson,
                            ROOT::Math::PtEtaPhiMVector &vis_gen_boson,
                            const int &n_jets) {
            // For Run3 jets with pT > 30 GeV and |eta| < 2.5 or pT > 50 GeV outside the tracker region
            // need to be considered (avoid jet horn region)
            // type of n_jets needs to be a float for the correctionlib evaluation
            float nJets = static_cast<float>(n_jets);
            float genPt = gen_boson.Pt();  // generator Z(W) pT
            ROOT::Math::PtEtaPhiMVector met_new;

            Logger::get("met::RecoilCorrection")->debug("Corrector Inputs");
            Logger::get("met::RecoilCorrection")->debug("N jets {} ", nJets);
            Logger::get("met::RecoilCorrection")->debug("gen. boson Px {} ", gen_boson.Px());
            Logger::get("met::RecoilCorrection")->debug("gen. boson Py {} ", gen_boson.Py());
            Logger::get("met::RecoilCorrection")->debug("vis. gen. boson Px {} ", vis_gen_boson.Px());
            Logger::get("met::RecoilCorrection")->debug("vis. gen. boson Py {} ", vis_gen_boson.Py());
            Logger::get("met::RecoilCorrection")->debug("old MET X {} ", met.Px());
            Logger::get("met::RecoilCorrection")->debug("old MET Y {} ", met.Py());
            Logger::get("met::RecoilCorrection")->debug("old MET {} ", met.Pt());

            if (std::isnan(met.Px()) || std::isnan(met.Py())) {
                Logger::get("met::RecoilCorrection")->debug("NaN detected in MET, returning uncorrected MET");
                return met;
            } 
            else {
                if (method == "Rescaling" || method == "QuantileMapHist") {
                    ROOT::Math::PtEtaPhiMVector U = met + vis_gen_boson - gen_boson;
                    float dPhi_U = U.Phi() - gen_boson.Phi();
                    float Upara = U.Pt() * std::cos(dPhi_U);
                    float Uperp = U.Pt() * std::sin(dPhi_U);

                    float Upara_new = RecoilCorr->evaluate({order, nJets, genPt, std::string("Upara"), Upara});
                    float Uperp_new = RecoilCorr->evaluate({order, nJets, genPt, std::string("Uperp"), Uperp});
                    
                    float Upt_new = std::sqrt(Upara_new*Upara_new + Uperp_new*Uperp_new);
                    float Uphi_new = std::atan2(Uperp_new, Upara_new) + gen_boson.Phi();
                    ROOT::Math::PtEtaPhiMVector U_new = ROOT::Math::PtEtaPhiMVector(Upt_new, 0., Uphi_new, 0.);
                    ROOT::Math::PtEtaPhiMVector met_new = U_new - vis_gen_boson + gen_boson;
                }
                else if (method == "Uncertainty") {
                    if (std::set<std::string>{"RespUp", "RespDown", "ResolUp", "ResolDown"}.count(variation)) {
                        ROOT::Math::PtEtaPhiMVector H = - met - vis_gen_boson;
                        float dPhi_H = H.Phi() - gen_boson.Phi();
                        float Hpara = H.Pt() * std::cos(dPhi_H);
                        float Hperp = H.Pt() * std::sin(dPhi_H);

                        float Hpara_new = RecoilCorr->evaluate({order, nJets, genPt, std::string("Hpara"), Hpara, variation});
                        float Hperp_new = RecoilCorr->evaluate({order, nJets, genPt, std::string("Hperp"), Hperp, variation});

                        float Hpt_new = std::sqrt(Hpara_new*Hpara_new + Hperp_new*Hperp_new);
                        float Hphi_new = std::atan2(Hperp_new, Hpara_new) + gen_boson.Phi();
                        if (Hphi_new > M_PI) Hphi_new -= 2*M_PI;
                        if (Hphi_new < -M_PI) Hphi_new += 2*M_PI;
                        ROOT::Math::PtEtaPhiMVector H_new = ROOT::Math::PtEtaPhiMVector(Hpt_new, 0., Hphi_new, 0.);
                        ROOT::Math::PtEtaPhiMVector met_new = - H_new - vis_gen_boson;
                    }
                    else {
                        Logger::get("met::RecoilCorrection")
                            ->error("Variation {} not known. Choose either 'RespUp', 'RespDown', 'ResolUp' or 'ResolDown'", variation);
                        throw std::runtime_error("Invalid variation for Recoil corrections");
                    }
                }
                else if (method == "QuantileMapFit") {
                    Logger::get("met::RecoilCorrection")
                        ->debug("QuantileMapFit method not yet implemented, returning uncorrected MET");
                    return met;
                }
                else {
                    Logger::get("met::RecoilCorrection")
                        ->error("Method {} not known. Choose either 'Rescaling', 'QuantileMapHist' or 'Uncertainty'", method);
                    throw std::runtime_error("Invalid method for Recoil corrections");
                }
                
                Logger::get("met::RecoilCorrection")->debug("corrected MET X {} ", met_new.Px());
                Logger::get("met::RecoilCorrection")->debug("corrected MET Y {} ", met_new.Py());
                Logger::get("met::RecoilCorrection")->debug("corrected MET {} ", met_new.Pt());
                return met_new;
            }
        };
        return df.Define(outputname, Correction,
                         {p4_met, p4_gen_boson, p4_vis_gen_boson, n_jets});
    } else {
        // if we do not apply the recoil corrections, just rename the met
        // column to the new outputname and dont change anything else
        return event::quantity::Rename<ROOT::Math::PtEtaPhiMVector>(df, outputname, p4_met);
    }
}

/**
 * @brief This function applies Recoil corrections on a given sample. For
 * more information on recoil corrections, check [this
 * repo](https://github.com/KIT-CMS/RecoilCorrections/). The code for the
 * application of these corrections can be found in `src/RecoilCorrections`.
 * This correction is applicable for processes involving W and Z bosons.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected MET
 * @param p4_met initial/uncorrected MET Lorentz vector
 * @param p4_gen_boson name of the column containing the generator-level boson
 * Lorentz vector
 * @param p4_vis_gen_boson name of the column containing the visible part of
 * the generator-level boson Lorentz vector
 * @param jet_pt name of the column containing the \f$p_T\f$ of all jets in
 * an event
 * @param corr_file path to the recoil corrections file
 * @param syst_file path to the systematics corrections file
 * @param apply_correction if bool is set to true, the recoil correction is
 * applied, if not, the output column contains the original MET vector
 * @param resolution bool - if set, resolution corrections are applied
 * @param response bool - if set, response corrections are applied
 * @param shift_up bool - if set, the up shift is applied
 * @param shift_down bool - if set, the down shift is applied
 * @param is_Wjets bool - if set, the number of jets is incrased by one (this
 * is only needed for W+jets samples)
 *
 * @return a dataframe containing a new column with the new MET vector
 *
 * @warning The corrections are derived by the H(tautau) group for the legacy
 * analyses of Run2. They are not up to standard with the correctionlib json
 * format, therefore, to be used with caution.
 */
ROOT::RDF::RNode RecoilCorrection(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &p4_met, const std::string &p4_gen_boson,
    const std::string &p4_vis_gen_boson, const std::string &jet_pt,
    const std::string &corr_file, const std::string &syst_file,
    const bool apply_correction, const bool resolution, const bool response,
    const bool shift_up, const bool shift_down, const bool is_Wjets) {
    if (apply_correction) {
        Logger::get("met::RecoilCorrection")->debug("Will run recoil corrections");
        const auto corrector = new RecoilCorrector(corr_file);
        const auto systematics = new MetSystematic(syst_file);
        auto shiftType = MetSystematic::SysShift::Nominal;
        if (shift_up) {
            shiftType = MetSystematic::SysShift::Up;
        } else if (shift_down) {
            shiftType = MetSystematic::SysShift::Down;
        }
        auto sysType = MetSystematic::SysType::None;
        if (response) {
            sysType = MetSystematic::SysType::Response;
        } else if (resolution) {
            sysType = MetSystematic::SysType::Resolution;
        }
        auto RecoilCorrections = [sysType, systematics, shiftType, corrector,
                                  is_Wjets](
                                     ROOT::Math::PtEtaPhiMVector &met,
                                     ROOT::Math::PtEtaPhiMVector &gen_boson,
                                     ROOT::Math::PtEtaPhiMVector &vis_gen_boson,
                                     const ROOT::RVec<float> &jet_pt) {
            // TODO is this the correct number of jets ?
            auto jets_above30 =
                Filter(jet_pt, [](float pt) { return pt > 30; });
            int nJets30 = jets_above30.size();
            if (is_Wjets) {
                nJets30 = nJets30 + 1;
            }
            float MetX = met.Px();
            float MetY = met.Py();
            float correctedMetX = 0.;
            float correctedMetY = 0.;
            ROOT::Math::PtEtaPhiMVector corrected_met;
            float genPx = gen_boson.Px();  // generator Z(W) px
            float genPy = gen_boson.Py();  // generator Z(W) py
            float visPx = vis_gen_boson.Px(); // visible (generator) Z(W) px
            float visPy = vis_gen_boson.Py(); // visible (generator) Z(W) py
            Logger::get("met::RecoilCorrection")->debug("Corrector Inputs");
            Logger::get("met::RecoilCorrection")->debug("nJets30 {} ", nJets30);
            Logger::get("met::RecoilCorrection")->debug("genPx {} ", genPx);
            Logger::get("met::RecoilCorrection")->debug("genPy {} ", genPy);
            Logger::get("met::RecoilCorrection")->debug("visPx {} ", visPx);
            Logger::get("met::RecoilCorrection")->debug("visPy {} ", visPy);
            Logger::get("met::RecoilCorrection")->debug("MetX {} ", MetX);
            Logger::get("met::RecoilCorrection")->debug("MetY {} ", MetY);
            Logger::get("met::RecoilCorrection")->debug("old met {} ", met.Pt());
            corrector->CorrectWithHist(MetX, MetY, genPx, genPy, visPx, visPy,
                                       nJets30, correctedMetX, correctedMetY);
            Logger::get("met::RecoilCorrection")
                ->debug("correctedMetX {} ", correctedMetX);
            Logger::get("met::RecoilCorrection")
                ->debug("correctedMetY {} ", correctedMetY);
            // only apply shifts if the correpsonding variables are set
            if (sysType != MetSystematic::SysType::None &&
                shiftType != MetSystematic::SysShift::Nominal) {
                Logger::get("met::RecoilCorrection")
                    ->debug(" apply systematics {} {}", (int)sysType,
                            (int)shiftType);
                systematics->ApplyMetSystematic(
                    correctedMetX, correctedMetY, genPx, genPy, visPx, visPy,
                    nJets30, sysType, shiftType, correctedMetX, correctedMetY);
            }
            corrected_met.SetPxPyPzE(correctedMetX, correctedMetY, 0,
                                     std::sqrt(correctedMetX * correctedMetX +
                                               correctedMetY * correctedMetY));
            Logger::get("met::RecoilCorrection")
                ->debug("shifted and corrected met {} ", corrected_met.Pt());

            return corrected_met;
        };
        return df.Define(outputname, RecoilCorrections,
                         {p4_met, p4_gen_boson, p4_vis_gen_boson, jet_pt});
    } else {
        // if we do not apply the recoil corrections, just rename the met
        // column to the new outputname and dont change anything else
        return event::quantity::Rename<ROOT::Math::PtEtaPhiMVector>(df, outputname, p4_met);
    }
}

/**
 * @brief This function applies MET \f$\phi\f$ corrections as provided by JME POG
 * for Run2 in correction JSON files.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected MET Lorentz
 * vector
 * @param p4_met initial/uncorrected MET Lorentz vector
 * @param n_pv name of the column containing the number of primary vertices
 * in the event
 * @param run name of the column containing the run number
 * @param corr_file path to the file containing the correction
 * @param corr_name name of the correction to be applied, e.g. "metphicorr_pfmet_mc",
 * possible options are "puppimet" instead of "pfmet" or "data" instead of "mc"
 *
 * @return a dataframe with the new column containing the corrected MET Lorentz
 * vector
 *
 * @note This function is only valid for Run2 data and MC. For Run3 the overloaded
 * version of this function should be used.
 */
ROOT::RDF::RNode
METPhiCorrection(ROOT::RDF::RNode df, const std::string &outputname,
       const std::string &p4_met, const std::string &n_pv,
       const std::string &run, const std::string &corr_file,
       const std::string &corr_name) {

    auto evaluator_met_pt =
        correction::CorrectionSet::from_file(corr_file)->at("pt_" + corr_name);
    auto evaluator_met_phi =
        correction::CorrectionSet::from_file(corr_file)->at("phi_" + corr_name);

    auto xyCorrection = [evaluator_met_pt, evaluator_met_phi](
                      const ROOT::Math::PtEtaPhiMVector &met,
                      const int npv, const UInt_t run) {
        Logger::get("met::METPhiCorrection")
            ->debug("before: pt {} phi {}", met.Pt(), met.Phi());

        double corr_met_pt = evaluator_met_pt->evaluate(
            {met.Pt(), met.Phi(), float(npv), float(run)});
        double corr_met_phi = evaluator_met_phi->evaluate(
            {met.Pt(), met.Phi(), float(npv), float(run)});

        double corr_met_X = corr_met_pt * cos(corr_met_phi);
        double corr_met_Y = corr_met_pt * sin(corr_met_phi);
        ROOT::Math::PtEtaPhiMVector corr_met;
        corr_met.SetPxPyPzE(
            corr_met_X, corr_met_Y, 0,
            std::sqrt(corr_met_X * corr_met_X + corr_met_Y * corr_met_Y));

        Logger::get("met::METPhiCorrection")
            ->debug("after: pt {} phi {}", corr_met.Pt(), corr_met.Phi());

        return corr_met;
    };
    return df.Define(outputname, xyCorrection, {p4_met, n_pv, run});
}

/**
 * @brief This function applies MET \f$\phi\f$ corrections as provided by JME POG
 * for Run3 in correction JSON files.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected MET Lorentz
 * vector
 * @param p4_met initial/uncorrected MET Lorentz vector
 * @param n_pv name of the column containing the number of primary vertices
 * in the event
 * @param corr_file path to the file containing the correction
 * @param corr_name name of the correction to be applied, e.g. "met_xy_corrections"
 * @param met_type type of the MET, possible options are "PuppiMET" and "MET" (for
 * PF)
 * @param era data-taking period, possible options are e.g. "2022", "2022EE",
 * "2023", "2023BPix"
 * @param is_mc boolean indicating whether the sample is MC or data
 * @param stat_variation name of the statistical variation, possible options
 * are "nom", "xup", "xdn", "yup", "ydn"
 * @param pileup_variation name of the pileup variation, possible options
 * are "nom", "pu_up", "pu_dn"
 *
 * @return a dataframe with the new column containing the corrected MET Lorentz
 * vector
 *
 * @note This function is only valid for Run3 data and MC. For Run2 the overloaded
 * version of this function should be used.
 */
ROOT::RDF::RNode
METPhiCorrection(ROOT::RDF::RNode df, const std::string &outputname,
       const std::string &p4_met, const std::string &n_pv,
       const std::string &corr_file, const std::string &corr_name,
       const std::string &met_type, const std::string &era,
       const bool is_mc, const std::string &stat_variation,
       const std::string &pileup_variation) {

    auto evaluator =
        correction::CorrectionSet::from_file(corr_file)->at(corr_name);

    const std::string data_mc_key = is_mc ? "MC" : "DATA";
    const std::string pt_key = stat_variation != "nom" ?
                               "pt_stat_" + stat_variation : "pt";
    const std::string phi_key = stat_variation != "nom" ?
                                "phi_stat_" + stat_variation : "phi";

    auto xyCorrection = [evaluator, met_type, era, data_mc_key,
                         pt_key, phi_key, pileup_variation](
                      const ROOT::Math::PtEtaPhiMVector &met,
                      const int npv) {
        Logger::get("met::METPhiCorrection")
            ->debug("before: pt {} phi {}", met.Pt(), met.Phi());

        double corr_met_pt = evaluator->evaluate(
            {pt_key, met_type, era, data_mc_key, pileup_variation,
             met.Pt(), met.Phi(), float(npv)});
        double corr_met_phi = evaluator->evaluate(
            {phi_key, met_type, era, data_mc_key, pileup_variation,
             met.Pt(), met.Phi(), float(npv)});

        double corr_met_X = corr_met_pt * cos(corr_met_phi);
        double corr_met_Y = corr_met_pt * sin(corr_met_phi);
        ROOT::Math::PtEtaPhiMVector corr_met;
        corr_met.SetPxPyPzE(
            corr_met_X, corr_met_Y, 0,
            std::sqrt(corr_met_X * corr_met_X + corr_met_Y * corr_met_Y));

        Logger::get("met::METPhiCorrection")
            ->debug("after: pt {} phi {}", corr_met.Pt(), corr_met.Phi());

        return corr_met;
    };
    return df.Define(outputname, xyCorrection, {p4_met, n_pv});
}
} // end namespace met

namespace physicsobject {

/**
 * @brief This function propagates object corrections to the MET based on
 * vectors of physics objects. The objects can be e.g. a collection/vector of jets.
 * If the energy of an object is corrected/changed (e.g. via some scale factor) or
 * due to a shift, this change in energy has to be propagated to the MET vector, and
 * the MET vector has to be adapted accordingly. The MET is recalculated via
 *
 * \f[
 *  E_{T,miss,x}^{\text{corrected}} = E_{T,miss,x} + p_{x,\text{object}}
 *        - p_{x,\text{object}}^{\text{corrected}} \\
 *  E_{T,miss,y}^{\text{corrected}} = E_{T,miss,y} + p_{y,\text{object}}
 *        - p_{y,\text{object}}^{\text{corrected}} \\
 *  E_{T,miss}^{\text{corrected}} = \sqrt{E_{T,miss,x}^{\text{corrected}} * E_{T,miss,x}^{\text{corrected}}
 *        + E_{T,miss,y}^{\text{corrected}} * E_{T,miss,y}^{\text{corrected}}}
 * \f]
 *
 * The correction is done for all jets above a certain \f$p_T\f$ threshold.
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected MET Lorentz
 * vector
 * @param p4_met initial/uncorrected MET Lorentz vector
 * @param pt_corrected name of the column containing \f$p_T\f$ vector of the
 * corrected objects
 * @param eta_corrected name of the column containing \f$\eta\f$ vector of the
 * corrected objects
 * @param phi_corrected name of the column containing \f$\phi\f$ vector of the
 * corrected objects
 * @param mass_corrected name of the column containing mass vector of the
 * corrected objects
 * @param pt name of the column containing \f$p_T\f$ vector of the uncorrected
 * objects
 * @param eta name of the column containing \f$\eta\f$ vector of the uncorrected
 * objects
 * @param phi name of the column containing \f$\phi\f$ vector of the uncorrected
 * objects
 * @param mass name of the column containing mass vector of the uncorrected objects
 * @param apply_propagation boolean indicating whether the propagation should be
 * applied or just the original MET vector should be returned
 * @param min_pt minimal \f$p_T\f$, the corrected object has to have in order for
 * the MET propagation to be applied
 *
 * @return a dataframe with the new column containing the corrected MET Lorentz
 * vector
 */
ROOT::RDF::RNode PropagateToMET(
    ROOT::RDF::RNode df, const std::string &outputname, const std::string &p4_met,
    const std::string &pt_corrected, const std::string &eta_corrected,
    const std::string &phi_corrected, const std::string &mass_corrected,
    const std::string &pt, const std::string &eta,
    const std::string &phi, const std::string &mass,
    bool apply_propagation, float min_pt) {
    // propagate objects corrections to MET, since we can have an arbitrary
    // amount of objects, this has to be done per event
    auto scaleMet = [min_pt](const ROOT::Math::PtEtaPhiMVector &met,
                                 const ROOT::RVec<float> &pts_corrected,
                                 const ROOT::RVec<float> &etas_corrected,
                                 const ROOT::RVec<float> &phis_corrected,
                                 const ROOT::RVec<float> &masses_corrected,
                                 const ROOT::RVec<float> &pts,
                                 const ROOT::RVec<float> &etas,
                                 const ROOT::RVec<float> &phis,
                                 const ROOT::RVec<float> &masses) {
        ROOT::Math::PtEtaPhiMVector corrected_met;
        ROOT::Math::PtEtaPhiMVector uncorrected_object;
        ROOT::Math::PtEtaPhiMVector corrected_object;
        float corr_x = 0.0;
        float corr_y = 0.0;
        // now loop through all objects in the event
        for (std::size_t index = 0; index < pts.size(); ++index) {
            // only propagate objects above the given pt threshold
            if (pts_corrected.at(index) > min_pt) {
                // construct the uncorrected and the corrected Lorentz
                // vectors
                corrected_object = ROOT::Math::PtEtaPhiMVector(
                    pts_corrected.at(index), etas_corrected.at(index),
                    phis_corrected.at(index), masses_corrected.at(index));
                uncorrected_object = ROOT::Math::PtEtaPhiMVector(
                    pts.at(index), etas.at(index), phis.at(index),
                    masses.at(index));
                // update the correction factors that are applied to the MET
                corr_x += uncorrected_object.Px() - corrected_object.Px();
                corr_y += uncorrected_object.Py() - corrected_object.Py();
            }
        }
        float MetX = met.Px() + corr_x;
        float MetY = met.Py() + corr_y;
        Logger::get("physicsobject::PropagateToMET")
            ->debug("corr_x {}, corr_y {}", corr_x, corr_y);
        Logger::get("physicsobject::PropagateToMET")
            ->debug("MetX {}, MetY {}", MetX, MetY);
        corrected_met.SetPxPyPzE(MetX, MetY, 0,
                                 std::sqrt(MetX * MetX + MetY * MetY));
        Logger::get("physicsobject::PropagateToMET")->debug("old met {}", met.Pt());
        Logger::get("physicsobject::PropagateToMET")
            ->debug("corrected met {}", corrected_met.Pt());
        return corrected_met;
    };
    if (apply_propagation) {
        return df.Define(outputname, scaleMet,
                         {p4_met, pt_corrected, eta_corrected,
                          phi_corrected, mass_corrected, pt,
                          eta, phi, mass});
    } else {
        // if we do not apply the propagation, just rename the met column to
        // the new outputname and dont change anything else
        return event::quantity::Rename<ROOT::Math::PtEtaPhiMVector>(df, outputname, p4_met);
    }
}
} // end namespace physicsobject
#endif /* GUARD_MET_H */
