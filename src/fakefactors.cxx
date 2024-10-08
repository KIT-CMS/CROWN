#ifndef GUARDFAKEFACTORS_H
#define GUARDFAKEFACTORS_H
/// The namespace that contains the fake factor function.
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "correction.h"

namespace fakefactors {
/**
 * @brief Function to calculate raw fake factors without corrections with
 * correctionlib for the semileptonic channels
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param nbtags number of good b-tagged jets in the event
 * @param qcd_variation name of the QCD FF uncertainty variation or nominal
 * @param wjets_variation name of the Wjets FF uncertainty variation or nominal
 * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
 * @param fraction_variation name of the process fraction uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
raw_fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &tau_pt, const std::string &njets,
                        const std::string &lep_mt, const std::string &nbtags,
                        const std::string &qcd_variation, const std::string &wjets_variation,
                        const std::string &ttbar_variation, const std::string &fraction_variation,
                        const std::string &ff_file) {
    Logger::get("RawFakeFactor")
        ->debug("Setting up functions for raw fake factor (without "
                "corrections) evaluation with correctionlib");
    Logger::get("RawFakeFactor")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("RawFakeFactor")->debug("Wjets variation - Name {}", wjets_variation);
    Logger::get("RawFakeFactor")->debug("ttbar variation - Name {}", ttbar_variation);
    Logger::get("RawFakeFactor")->debug("Fraction variation - Name {}", fraction_variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");
    auto calc_fake_factor = [qcd_variation, wjets_variation, ttbar_variation, fraction_variation, 
                             qcd, wjets, ttbar, fractions](const float &pt_2, const int &njets,
                                        const float &mt_1, const int &nbtag) {
        float ff = 0.;
        if (pt_2 >= 0.) {
            Logger::get("RawFakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, qcd_variation});
            Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff = wjets->evaluate({pt_2, (float)njets, wjets_variation});
            Logger::get("RawFakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, ttbar_variation});
            Logger::get("RawFakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("RawFakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("RawFakeFactor")->debug("N b-jets - value {}", nbtag);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)nbtag, fraction_variation});
            Logger::get("RawFakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)nbtag, fraction_variation});
            Logger::get("RawFakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)nbtag, fraction_variation});
            Logger::get("RawFakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);

            ff = qcd_frac * std::max(qcd_ff, (float)0.) + wjets_frac * std::max(wjets_ff, (float)0.) +
                 ttbar_frac * std::max(ttbar_ff, (float)0.);
        }

        Logger::get("RawFakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {tau_pt, njets, lep_mt, nbtags});
    return df1;
}
/**
 * @brief Function to calculate raw fake factors without corrections with
 * correctionlib for the NMSSM Di-Higgs analysis for the full hadronic channel
 *
 * @param df the dataframe to add the quantity to
 * @param outputname name of the output column for the fake factor
 * @param tau_idx index of the tau, leading/subleading
 * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
 * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param m_vis visible di-tau mass of the tau pair
 * @param nbtag number of good b-tagged jets in the event
 * @param qcd_variation name of the QCD FF uncertainty variation or nominal
 * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
 * @param fraction_variation name of the process fraction uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
raw_fakefactor_nmssm_tt(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &tau_idx, const std::string &tau_pt_1,
                        const std::string &tau_pt_2, const std::string &njets, 
                        const std::string &m_vis, const std::string &nbtag,
                        const std::string &qcd_variation, const std::string &ttbar_variation, 
                        const std::string &fraction_variation, const std::string &ff_file) {

    Logger::get("RawFakeFactor")
        ->debug("Setting up functions for raw fake factor (without "
                "corrections) evaluation with correctionlib");
    Logger::get("RawFakeFactor")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("RawFakeFactor")->debug("ttbar variation - Name {}", ttbar_variation);
    Logger::get("RawFakeFactor")->debug("Fraction variation - Name {}", fraction_variation);

    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto qcd_subleading = correction::CorrectionSet::from_file(ff_file)->at(
        "QCD_subleading_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto ttbar_subleading = correction::CorrectionSet::from_file(ff_file)->at(
        "ttbar_subleading_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");
    auto fractions_subleading =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions_subleading");
    
    auto calc_fake_factor = [tau_idx, qcd_variation, ttbar_variation, fraction_variation, qcd, qcd_subleading, ttbar, ttbar_subleading, fractions, fractions_subleading](
                                const float &pt_1, const float &pt_2,
                                const int &njets, const float &m_vis, const int &nbtag) {
        float ff = 0.;
        if (pt_2 >= 0.) {
            Logger::get("RawFakeFactor")
                ->debug("Leading Tau pt - value {}", pt_1);
            Logger::get("RawFakeFactor")
                ->debug("Subleading Tau pt - value {}", pt_2);
            Logger::get("RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = -1.;
            float ttbar_ff = -1.;
            float qcd_frac = -1.;
            float wjets_frac = -1.;
            float ttbar_frac = -1.;
            if (tau_idx == 0) {
                qcd_ff = qcd->evaluate({pt_1, (float)njets, qcd_variation});
                Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
                ttbar_ff = ttbar->evaluate({pt_1, (float)njets, ttbar_variation});
                Logger::get("RawFakeFactor")->debug("ttbar - value {}", ttbar_ff);
                qcd_frac = fractions->evaluate({"QCD", m_vis, (float)nbtag, fraction_variation});
                Logger::get("RawFakeFactor")->debug("QCD - fraction {}", qcd_frac);
                wjets_frac = fractions->evaluate({"Wjets", m_vis, (float)nbtag, fraction_variation});
                Logger::get("RawFakeFactor")->debug("Wjets - fraction {}", wjets_frac);
                ttbar_frac = fractions->evaluate({"ttbar", m_vis, (float)nbtag, fraction_variation});
                Logger::get("RawFakeFactor")->debug("ttbar - fraction {}", ttbar_frac);
                
                ff = (qcd_frac + wjets_frac) * std::max(qcd_ff, (float)0.) + ttbar_frac * std::max(ttbar_ff, (float)0.);
            } else if (tau_idx == 1) {
                qcd_ff =
                    qcd_subleading->evaluate({pt_2, (float)njets, qcd_variation});
                Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
                ttbar_ff = ttbar_subleading->evaluate({pt_1, (float)njets, ttbar_variation});
                Logger::get("RawFakeFactor")->debug("ttbar - value {}", ttbar_ff);
                qcd_frac = fractions_subleading->evaluate({"QCD", m_vis, (float)nbtag, fraction_variation});
                Logger::get("RawFakeFactor")->debug("QCD - fraction {}", qcd_frac);
                wjets_frac = fractions_subleading->evaluate({"Wjets", m_vis, (float)nbtag, fraction_variation});
                Logger::get("RawFakeFactor")->debug("Wjets - fraction {}", wjets_frac);
                ttbar_frac = fractions_subleading->evaluate({"ttbar", m_vis, (float)nbtag, fraction_variation});
                Logger::get("RawFakeFactor")->debug("ttbar - fraction {}", ttbar_frac);

                ff = (qcd_frac + wjets_frac) * std::max(qcd_ff, (float)0.) + ttbar_frac * std::max(ttbar_ff, (float)0.);
            }
        }

        Logger::get("RawFakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 =
        df.Define(outputname, calc_fake_factor, {tau_pt_1, tau_pt_2, njets, m_vis, nbtag});
    return df1;
}
/**
 * @brief Function to calculate fake factors with correctionlib for the
 * semileptonic channels
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param nbtags number of good b-tagged jets in the event
 * @param lep_pt pt of the leptonic tau in the tau pair
 * @param tau_mass mass of the hadronic tau in the tau pair
 * @param m_vis visible mass of the tau pair
 * @param qcd_variation name of the QCD FF uncertainty variation or nominal
 * @param wjets_variation name of the Wjets FF uncertainty variation or nominal
 * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
 * @param fraction_variation name of the process fraction uncertainty variation or nominal
 * @param qcd_corr_leppt_variation name of the QCD lepton pt correction uncertainty variation or nominal
 * @param qcd_corr_taumass_variation name of the QCD tau mass correction uncertainty variation or nominal
 * @param qcd_corr_drsr_variation name of the QCD DR to SR correction uncertainty variation or nominal
 * @param wjets_corr_leppt_variation name of the Wjets lepton pt correction uncertainty variation or nominal
 * @param wjets_corr_taumass_variation name of the Wjets tau mass correction uncertainty variation or nominal
 * @param wjets_corr_drsr_variation name of the Wjets DR to SR correction uncertainty variation or nominal
 * @param ttbar_corr_leppt_variation name of the ttbar lepton pt correction uncertainty variation or nominal
 * @param ttbar_corr_taumass_variation name of the ttbar tau mass correction uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tau_pt, const std::string &njets,
                    const std::string &lep_mt, const std::string &nbtags,
                    const std::string &lep_pt, const std::string &tau_mass,
                    const std::string &m_vis,
                    const std::string &qcd_variation, const std::string &wjets_variation,
                    const std::string &ttbar_variation, const std::string &fraction_variation,
                    const std::string &qcd_corr_leppt_variation, const std::string &qcd_corr_taumass_variation,
                    const std::string &qcd_corr_drsr_variation, const std::string &wjets_corr_leppt_variation,
                    const std::string &wjets_corr_taumass_variation, const std::string &wjets_corr_drsr_variation, 
                    const std::string &ttbar_corr_leppt_variation, const std::string &ttbar_corr_taumass_variation,
                    const std::string &ff_file, const std::string &ff_corr_file) {

    Logger::get("FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("FakeFactor")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("FakeFactor")->debug("Wjets variation - Name {}", wjets_variation);
    Logger::get("FakeFactor")->debug("ttbar variation - Name {}", ttbar_variation);
    Logger::get("FakeFactor")->debug("Fraction variation - Name {}", fraction_variation);
    Logger::get("FakeFactor")->debug("QCD lep pt corr variation - Name {}", qcd_corr_leppt_variation);
    Logger::get("FakeFactor")->debug("QCD tau mass corr variation - Name {}", qcd_corr_taumass_variation);
    Logger::get("FakeFactor")->debug("QCD DRSR corr variation - Name {}", qcd_corr_drsr_variation);
    Logger::get("FakeFactor")->debug("Wjets lep pt corr variation - Name {}", wjets_corr_leppt_variation);
    Logger::get("FakeFactor")->debug("Wjets tau mass corr variation - Name {}", wjets_corr_taumass_variation);
    Logger::get("FakeFactor")->debug("Wjets DRSR corr variation - Name {}", wjets_corr_drsr_variation);
    Logger::get("FakeFactor")->debug("ttbar lep pt corr variation - Name {}", ttbar_corr_leppt_variation);
    Logger::get("FakeFactor")->debug("ttbar tau mass corr variation - Name {}", ttbar_corr_taumass_variation);

    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto qcd_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_leading_lep_pt_correction");
    auto qcd_tau_mass_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_subleading_lep_mass_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto wjets_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_leading_lep_pt_correction");
    auto wjets_tau_mass_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_subleading_lep_mass_correction");
    auto wjets_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                           ->at("Wjets_DR_SR_correction");
    auto ttbar_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_leading_lep_pt_correction");
    auto ttbar_tau_mass_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_subleading_lep_mass_correction");
    auto calc_fake_factor = [qcd_variation, wjets_variation, ttbar_variation, fraction_variation, 
                             qcd_corr_leppt_variation, qcd_corr_taumass_variation, qcd_corr_drsr_variation, 
                             wjets_corr_leppt_variation, wjets_corr_taumass_variation, wjets_corr_drsr_variation, 
                             ttbar_corr_leppt_variation, ttbar_corr_taumass_variation,
                             qcd, wjets, ttbar, fractions, qcd_lep_pt_closure, qcd_tau_mass_closure, qcd_DR_SR,
                             wjets_lep_pt_closure, wjets_tau_mass_closure, wjets_DR_SR, ttbar_lep_pt_closure, ttbar_tau_mass_closure](
                                const float &pt_2, const int &njets,
                                const float &mt_1, const int &nbtag,
                                const float &pt_1, const float &mass_2,
                                const float &m_vis) {
        float ff = 0.;
        if (pt_2 >= 0.) {
            Logger::get("FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, qcd_variation});
            Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff = wjets->evaluate({pt_2, (float)njets, wjets_variation});
            Logger::get("FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, ttbar_variation});
            Logger::get("FakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("FakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("FakeFactor")->debug("N b-jets - value {}", nbtag);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)nbtag, fraction_variation});
            Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)nbtag, fraction_variation});
            Logger::get("FakeFactor")->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)nbtag, fraction_variation});
            Logger::get("FakeFactor")->debug("ttbar - fraction {}", ttbar_frac);

            Logger::get("FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("FakeFactor")->debug("Tau mass - value {}", mass_2);
            Logger::get("FakeFactor")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({pt_1, qcd_corr_leppt_variation});
            Logger::get("FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_tau_mass_corr =
                qcd_tau_mass_closure->evaluate({mass_2, qcd_corr_taumass_variation});
            Logger::get("FakeFactor")
                ->debug("QCD - tau mass correction {}", qcd_tau_mass_corr);
            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, qcd_corr_drsr_variation});
            Logger::get("FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({pt_1, wjets_corr_leppt_variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_tau_mass_corr =
                wjets_tau_mass_closure->evaluate({mass_2, wjets_corr_taumass_variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - tau mass correction {}", wjets_tau_mass_corr);
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, wjets_corr_drsr_variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({pt_1, ttbar_corr_leppt_variation});
            Logger::get("FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);
            float ttbar_tau_mass_corr =
                ttbar_tau_mass_closure->evaluate({mass_2, ttbar_corr_taumass_variation});
            Logger::get("FakeFactor")
                ->debug("ttbar - tau mass correction {}", ttbar_tau_mass_corr);

            ff = qcd_frac * std::max(qcd_ff, (float)0.) * qcd_lep_pt_corr * qcd_tau_mass_corr *
                     qcd_DR_SR_corr +
                 wjets_frac * std::max(wjets_ff, (float)0.) * wjets_lep_pt_corr * wjets_tau_mass_corr * wjets_DR_SR_corr +
                 ttbar_frac * std::max(ttbar_ff, (float)0.) * ttbar_lep_pt_corr * ttbar_tau_mass_corr;
        }

        Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(
        outputname, calc_fake_factor,
        {tau_pt, njets, lep_mt, nbtags, lep_pt, tau_mass, m_vis});
    return df1;
}

/**
 * @brief Function to calculate fake factors with correctionlib for the NMSSM
 * Di-Higgs boosted analysis for the semileptonic channel
 *
 * @param df the dataframe to add the quantity to
 * @param outputname name of the output column for the fake factor
 * @param boosted_tau_pt pt of the hadronic tau in the boosted tau pair
 * @param njets number of good jets in the event
 * @param boosted_lep_mt transverse mass of the leptonic tau in the boosted tau
 * pair
 * @param nbtags number of good b-tagged jets in the event
 * @param boosted_lep_pt pt of the leptonic tau in the boosted tau pair
 * @param boosted_m_vis visible mass of the boosted tau pair
 * @param boosted_dR_ditau distance in eta-phi plane of the boosted tau pair
 * @param qcd_variation name of the QCD FF uncertainty variation or nominal
 * @param wjets_variation name of the Wjets FF uncertainty variation or nominal
 * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
 * @param fraction_variation name of the process fraction uncertainty variation or nominal
 * @param qcd_corr_leppt_variation name of the QCD lepton pt correction uncertainty variation or nominal
 * @param qcd_corr_lepmt_variation name of the QCD transverse mass (lepton+MET) correction uncertainty variation or nominal
 * @param qcd_corr_drsr_variation name of the QCD DR to SR correction uncertainty variation or nominal
 * @param wjets_corr_leppt_variation name of the Wjets lepton pt correction uncertainty variation or nominal
 * @param wjets_corr_drsr_variation name of the Wjets DR to SR correction uncertainty variation or nominal
 * @param ttbar_corr_leppt_variation name of the ttbar lepton pt correction uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode fakefactor_nmssm_boosted_lt(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &boosted_tau_pt, const std::string &njets,
    const std::string &boosted_lep_mt, const std::string &nbtags,
    const std::string &boosted_lep_pt, const std::string &boosted_m_vis,
    const std::string &boosted_dR_ditau, const std::string &qcd_variation, 
    const std::string &wjets_variation, const std::string &ttbar_variation,
    const std::string &fraction_variation, const std::string &qcd_corr_leppt_variation, 
    const std::string &qcd_corr_lepmt_variation, const std::string &qcd_corr_drsr_variation, 
    const std::string &wjets_corr_leppt_variation, const std::string &wjets_corr_drsr_variation, 
    const std::string &ttbar_corr_leppt_variation, const std::string &ff_file, 
    const std::string &ff_corr_file) {

    Logger::get("FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("FakeFactor")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("FakeFactor")->debug("Wjets variation - Name {}", wjets_variation);
    Logger::get("FakeFactor")->debug("ttbar variation - Name {}", ttbar_variation);
    Logger::get("FakeFactor")->debug("Fraction variation - Name {}", fraction_variation);
    Logger::get("FakeFactor")->debug("QCD lep pt corr variation - Name {}", qcd_corr_leppt_variation);
    Logger::get("FakeFactor")->debug("QCD lep mt corr variation - Name {}", qcd_corr_lepmt_variation);
    Logger::get("FakeFactor")->debug("QCD DRSR corr variation - Name {}", qcd_corr_drsr_variation);
    Logger::get("FakeFactor")->debug("Wjets lep pt corr variation - Name {}", wjets_corr_leppt_variation);
    Logger::get("FakeFactor")->debug("Wjets DRSR corr variation - Name {}", wjets_corr_drsr_variation);
    Logger::get("FakeFactor")->debug("ttbar lep pt corr variation - Name {}", ttbar_corr_leppt_variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto qcd_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_leading_lep_pt_correction");
    auto qcd_lep_mt_closure = correction::CorrectionSet::from_file(ff_corr_file)
                                  ->at("QCD_non_closure_lep_mt_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto wjets_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_leading_lep_pt_correction");
    auto wjets_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                           ->at("Wjets_DR_SR_correction");
    auto ttbar_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_leading_lep_pt_correction");
    // auto ttbar_m_vis_closure =
    //     correction::CorrectionSet::from_file(ff_corr_file)
    //         ->at("ttbar_non_closure_m_vis_correction");
    auto calc_fake_factor = [qcd_variation, wjets_variation, ttbar_variation, fraction_variation, 
                             qcd_corr_leppt_variation, qcd_corr_lepmt_variation, qcd_corr_drsr_variation, 
                             wjets_corr_leppt_variation, wjets_corr_drsr_variation, ttbar_corr_leppt_variation, 
                             qcd, wjets, ttbar, fractions, qcd_lep_pt_closure, qcd_lep_mt_closure, qcd_DR_SR,
                             wjets_lep_pt_closure, wjets_DR_SR, ttbar_lep_pt_closure](
                                const float &boosted_pt_2, const int &njets,
                                const float &boosted_mt_1, const int &nbtag,
                                const float &boosted_pt_1,
                                const float &boosted_m_vis,
                                const float &boosted_dR_ditau) {
        float ff = 0.;
        if (boosted_pt_2 >= 0.) {
            Logger::get("FakeFactor")->debug("Tau pt - value {}", boosted_pt_2);
            Logger::get("FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff =
                qcd->evaluate({boosted_pt_2, (float)njets, qcd_variation});
            Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff =
                wjets->evaluate({boosted_pt_2, (float)njets, wjets_variation});
            Logger::get("FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff =
                ttbar->evaluate({boosted_pt_2, (float)njets, ttbar_variation});
            Logger::get("FakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("FakeFactor")->debug("Lep mt - value {}", boosted_mt_1);
            Logger::get("FakeFactor")->debug("N b-jets - value {}", nbtag);

            float qcd_frac = fractions->evaluate(
                {"QCD", boosted_mt_1, (float)nbtag, fraction_variation});
            Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac = fractions->evaluate(
                {"Wjets", boosted_mt_1, (float)nbtag, fraction_variation});
            Logger::get("FakeFactor")->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac = fractions->evaluate(
                {"ttbar", boosted_mt_1, (float)nbtag, fraction_variation});
            Logger::get("FakeFactor")->debug("ttbar - fraction {}", ttbar_frac);

            Logger::get("FakeFactor")->debug("Lep pt - value {}", boosted_pt_1);
            Logger::get("FakeFactor")->debug("m_vis - value {}", boosted_m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({boosted_pt_1, qcd_corr_leppt_variation});
            Logger::get("FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_lep_mt_corr =
                qcd_lep_mt_closure->evaluate({boosted_mt_1, qcd_corr_lepmt_variation});
            Logger::get("FakeFactor")
                ->debug("QCD - lep mt correction {}", qcd_lep_mt_corr);
            float qcd_DR_SR_corr =
                qcd_DR_SR->evaluate({boosted_dR_ditau, qcd_corr_drsr_variation});
            Logger::get("FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({boosted_pt_1, wjets_corr_leppt_variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_DR_SR_corr =
                wjets_DR_SR->evaluate({boosted_dR_ditau, wjets_corr_drsr_variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({boosted_pt_1, ttbar_corr_leppt_variation});
            Logger::get("FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);
            // float ttbar_m_vis_corr =
            //     ttbar_m_vis_closure->evaluate({boosted_m_vis, variation});
            // Logger::get("FakeFactor")
            //     ->debug("ttbar - m_vis correction {}", ttbar_m_vis_corr);

            ff = qcd_frac * std::max(qcd_ff, (float)0.) * qcd_lep_pt_corr * qcd_lep_mt_corr * qcd_DR_SR_corr +
                 wjets_frac * std::max(wjets_ff, (float)0.) * wjets_lep_pt_corr * wjets_DR_SR_corr +
                 ttbar_frac * std::max(ttbar_ff, (float)0.) * ttbar_lep_pt_corr;
        }

        Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {boosted_tau_pt, njets, boosted_lep_mt, nbtags,
                          boosted_lep_pt, boosted_m_vis, boosted_dR_ditau});
    return df1;
}
/**
 * @brief Function to calculate fake factors with correctionlib
 *
 * @param df the dataframe to add the quantity to
 * @param outputname name of the output column for the fake factor
 * @param tau_idx index of the tau, leading/subleading
 * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
 * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param m_vis visible mass of the tau pair
 * @param nbtag number of good b-tagged jets in the event
 * @param tau_mass_1 mass of the leading hadronic tau in the tau pair
 * @param tau_mass_2 mass of the subleading hadronic tau in the tau pair
 * @param qcd_variation name of the QCD FF uncertainty variation or nominal
 * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
 * @param fraction_variation name of the process fraction uncertainty variation or nominal
 * @param qcd_corr_leppt_variation name of the QCD lepton pt correction uncertainty variation or nominal
 * @param qcd_corr_taumass_variation name of the QCD lepton mass correction uncertainty variation or nominal
 * @param qcd_corr_drsr_variation name of the QCD DR to SR correction uncertainty variation or nominal
 * @param ttbar_corr_leppt_variation name of the ttbar lepton pt correction uncertainty variation or nominal
 * @param ttbar_corr_taumass_variation name of the ttbar lepton mass correction uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
fakefactor_nmssm_tt(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &tau_idx, const std::string &tau_pt_1,
                    const std::string &tau_pt_2, const std::string &njets,
                    const std::string &m_vis, const std::string &nbtag,
                    const std::string &tau_mass_1, const std::string &tau_mass_2,
                    const std::string &qcd_variation, const std::string &ttbar_variation,
                    const std::string &fraction_variation, const std::string &qcd_corr_leppt_variation,
                    const std::string &qcd_corr_taumass_variation, const std::string &qcd_corr_drsr_variation,
                    const std::string &ttbar_corr_leppt_variation, const std::string &ttbar_corr_taumass_variation,
                    const std::string &ff_file, const std::string &ff_corr_file) {

    Logger::get("FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("FakeFactor")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("FakeFactor")->debug("ttbar variation - Name {}", ttbar_variation);
    Logger::get("FakeFactor")->debug("Fraction variation - Name {}", fraction_variation);
    Logger::get("FakeFactor")->debug("QCD lepton pt variation - Name {}", qcd_corr_leppt_variation);
    Logger::get("FakeFactor")->debug("QCD lepton mass variation - Name {}", qcd_corr_taumass_variation);
    Logger::get("FakeFactor")->debug("QCD DRSR variation - Name {}", qcd_corr_drsr_variation);
    Logger::get("FakeFactor")->debug("ttbar lepton pt variation - Name {}", ttbar_corr_leppt_variation);
    Logger::get("FakeFactor")->debug("ttbar lepton mass variation - Name {}", ttbar_corr_taumass_variation);

    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto qcd_subleading = correction::CorrectionSet::from_file(ff_file)->at(
        "QCD_subleading_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto ttbar_subleading = correction::CorrectionSet::from_file(ff_file)->at(
        "ttbar_subleading_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");
    auto fractions_subleading =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions_subleading");
    auto qcd_tau_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_subleading_lep_pt_correction");
    auto qcd_tau_mass_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_leading_lep_mass_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto ttbar_tau_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_subleading_lep_pt_correction");
    auto ttbar_tau_mass_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_leading_lep_mass_correction");
    auto qcd_tau_pt_closure_subleading =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_subleading_non_closure_leading_lep_pt_correction");
    auto qcd_tau_mass_closure_subleading =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_subleading_non_closure_subleading_lep_mass_correction");
    auto qcd_DR_SR_subleading = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_subleading_DR_SR_correction");
    auto ttbar_tau_pt_closure_subleading =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_subleading_non_closure_leading_lep_pt_correction");
    auto ttbar_tau_mass_closure_subleading =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_subleading_non_closure_subleading_lep_mass_correction");

    auto calc_fake_factor = [tau_idx, qcd_variation, ttbar_variation, fraction_variation, 
                             qcd_corr_leppt_variation, qcd_corr_taumass_variation, qcd_corr_drsr_variation, 
                             ttbar_corr_leppt_variation, ttbar_corr_taumass_variation,
                             qcd, ttbar, fractions, qcd_tau_pt_closure, qcd_tau_mass_closure,
                             qcd_DR_SR, ttbar_tau_pt_closure, ttbar_tau_mass_closure, qcd_subleading,
                             ttbar_subleading, fractions_subleading, qcd_tau_pt_closure_subleading,
                             qcd_tau_mass_closure_subleading, qcd_DR_SR_subleading, 
                             ttbar_tau_pt_closure_subleading, ttbar_tau_mass_closure_subleading](
                                const float &pt_1, const float &pt_2,
                                const int &njets, const float &m_vis, const int &nbtag,
                                const float &mass_1, const float &mass_2) {
        float ff = 0.;
        if (pt_2 >= 0.) {
            Logger::get("FakeFactor")->debug("Leading Tau pt - value {}", pt_1);
            Logger::get("FakeFactor")
                ->debug("Subleading Tau pt - value {}", pt_2);
            Logger::get("FakeFactor")->debug("m_vis - value {}", m_vis);
            Logger::get("FakeFactor")->debug("N jets - value {}", njets);
            Logger::get("FakeFactor")->debug("N btag - value {}", nbtag);
            Logger::get("FakeFactor")->debug("Leading Tau mass - value {}", mass_1);
            Logger::get("FakeFactor")
                ->debug("Subleading Tau mass - value {}", mass_2);

            float qcd_ff = -1.;
            float ttbar_ff = -1.;
            float qcd_frac = -1.;
            float wjets_frac = -1.;
            float ttbar_frac = -1.;
            float qcd_tau_pt_corr = -1.;
            float qcd_tau_mass_corr = -1.;
            float qcd_DR_SR_corr = -1.;
            float ttbar_tau_pt_corr = -1.;
            float ttbar_tau_mass_corr = -1.;
            if (tau_idx == 0) {
                qcd_ff = qcd->evaluate({pt_1, (float)njets, qcd_variation});
                Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
                ttbar_ff = ttbar->evaluate({pt_1, (float)njets, ttbar_variation});
                Logger::get("FakeFactor")->debug("ttbar - value {}", ttbar_ff);
                qcd_frac = fractions->evaluate({"QCD", m_vis, (float)nbtag, fraction_variation});
                Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);
                wjets_frac = fractions->evaluate({"Wjets", m_vis, (float)nbtag, fraction_variation});
                Logger::get("FakeFactor")->debug("Wjets - fraction {}", wjets_frac);
                ttbar_frac = fractions->evaluate({"ttbar", m_vis, (float)nbtag, fraction_variation});
                Logger::get("FakeFactor")->debug("ttbar - fraction {}", ttbar_frac);

                qcd_tau_pt_corr =
                    qcd_tau_pt_closure->evaluate({pt_2, qcd_corr_leppt_variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - lep pt correction {}", qcd_tau_pt_corr);
                qcd_tau_mass_corr =
                    qcd_tau_mass_closure->evaluate({mass_1, qcd_corr_taumass_variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - lep mass correction {}", qcd_tau_mass_corr);
                qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, qcd_corr_drsr_variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
                ttbar_tau_pt_corr =
                    ttbar_tau_pt_closure->evaluate({pt_2, ttbar_corr_leppt_variation});
                Logger::get("FakeFactor")
                    ->debug("ttbar - lep pt correction {}", ttbar_tau_pt_corr);
                ttbar_tau_mass_corr =
                    ttbar_tau_mass_closure->evaluate({mass_1, ttbar_corr_taumass_variation});
                Logger::get("FakeFactor")
                    ->debug("ttbar - lep mass correction {}", ttbar_tau_mass_corr);

                ff = (qcd_frac + wjets_frac) * std::max(qcd_ff, (float)0.) * qcd_tau_pt_corr * qcd_tau_mass_corr * qcd_DR_SR_corr + 
                      ttbar_frac * std::max(ttbar_ff, (float)0.) * ttbar_tau_pt_corr * ttbar_tau_mass_corr;
            } else if (tau_idx == 1) {
                qcd_ff = qcd_subleading->evaluate({pt_2, (float)njets, qcd_variation});
                Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
                ttbar_ff = ttbar_subleading->evaluate({pt_2, (float)njets, ttbar_variation});
                Logger::get("FakeFactor")->debug("ttbar - value {}", ttbar_ff);
                qcd_frac = fractions_subleading->evaluate({"QCD", m_vis, (float)nbtag, fraction_variation});
                Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);
                wjets_frac = fractions_subleading->evaluate({"Wjets", m_vis, (float)nbtag, fraction_variation});
                Logger::get("FakeFactor")->debug("Wjets - fraction {}", wjets_frac);
                ttbar_frac = fractions_subleading->evaluate({"ttbar", m_vis, (float)nbtag, fraction_variation});
                Logger::get("FakeFactor")->debug("ttbar - fraction {}", ttbar_frac);

                qcd_tau_pt_corr =
                    qcd_tau_pt_closure_subleading->evaluate({pt_1, qcd_corr_leppt_variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - lep pt correction {}", qcd_tau_pt_corr);
                qcd_tau_mass_corr =
                    qcd_tau_mass_closure_subleading->evaluate({mass_2, qcd_corr_taumass_variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - lep mass correction {}", qcd_tau_mass_corr);
                qcd_DR_SR_corr = qcd_DR_SR_subleading->evaluate({m_vis, qcd_corr_drsr_variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
                ttbar_tau_pt_corr =
                    ttbar_tau_pt_closure_subleading->evaluate({pt_1, ttbar_corr_leppt_variation});
                Logger::get("FakeFactor")
                    ->debug("ttbar - lep pt correction {}", ttbar_tau_pt_corr);
                ttbar_tau_mass_corr =
                    ttbar_tau_mass_closure_subleading->evaluate({mass_2, ttbar_corr_taumass_variation});
                Logger::get("FakeFactor")
                    ->debug("ttbar - lep mass correction {}", ttbar_tau_mass_corr);

                ff = (qcd_frac + wjets_frac) * std::max(qcd_ff, (float)0.) * qcd_tau_pt_corr * qcd_tau_mass_corr * qcd_DR_SR_corr + 
                      ttbar_frac * std::max(ttbar_ff, (float)0.) * ttbar_tau_pt_corr * ttbar_tau_mass_corr;
            }
        }

        Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {tau_pt_1, tau_pt_2, njets, m_vis, nbtag, tau_mass_1, tau_mass_2});
    return df1;
}
/**
 * @brief Function to calculate raw fake factors without corrections with
 * correctionlib. In difference to the NMSSM version, njets is used for the
 * fraction binning, and an additional split in deltaR is applied for wjets.
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param delta_r delta R between the two taus
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
raw_fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &tau_pt, const std::string &njets,
                     const std::string &lep_mt, const std::string &delta_r,
                     const std::string &variation, const std::string &ff_file) {
    Logger::get("SM RawFakeFactor")
        ->debug("Setting up functions for raw fake factor (without "
                "corrections) evaluation with correctionlib");
    Logger::get("SM RawFakeFactor")->debug("Variation - Name {}", variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");
    auto calc_fake_factor = [variation, qcd, wjets, ttbar, fractions](
                                const float &pt_2, const int &njets,
                                const float &mt_1, const float &delta_r) {
        float ff = 0.;
        if (pt_2 >= 0.) {
            Logger::get("SM RawFakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM RawFakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff =
                wjets->evaluate({pt_2, (float)njets, delta_r, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("ttbar - value {}", ttbar_ff);

            Logger::get("SM RawFakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("SM RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);

            ff = qcd_frac * qcd_ff + wjets_frac * wjets_ff +
                 ttbar_frac * ttbar_ff;
        }

        Logger::get("SM RawFakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {tau_pt, njets, lep_mt, delta_r});
    return df1;
}
/**
 * @brief Function to calculate fake factors with correctionlib. In difference
 * to the NMSSM version, njets is used for the fraction binning, and an
 * additional split in deltaR is applied for wjets.
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param lep_pt pt of the leptonic tau in the tau pair
 * @param lep_iso isolation of the leptonic tau in the tau pair
 * @param m_vis visible mass of the tau pair
 * @param delta_r distance in eta-phi between the two taus
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                 const std::string &tau_pt, const std::string &njets,
                 const std::string &lep_mt, const std::string &lep_pt,
                 const std::string &lep_iso, const std::string &m_vis,
                 const std::string &delta_r, const std::string &variation,
                 const std::string &ff_file, const std::string &ff_corr_file) {
    Logger::get("SM FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("SM FakeFactor")->debug("Variation - Name {}", variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto qcd_lep_pt_closure = correction::CorrectionSet::from_file(ff_corr_file)
                                  ->at("QCD_non_closure_lep_pt_correction");
    auto qcd_lep_iso_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_lep_iso_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto wjets_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_lep_pt_correction");
    auto wjets_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                           ->at("Wjets_DR_SR_correction");
    auto ttbar_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_lep_pt_correction");
    auto calc_fake_factor = [variation, qcd, wjets, ttbar, fractions,
                             qcd_lep_pt_closure, qcd_lep_iso_closure, qcd_DR_SR,
                             wjets_lep_pt_closure, wjets_DR_SR,
                             ttbar_lep_pt_closure](
                                const float &pt_2, const int &njets,
                                const float &mt_1, const float &pt_1,
                                const float &iso_1, const float &m_vis,
                                const float &delta_r) {
        float ff = 0.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff =
                wjets->evaluate({pt_2, (float)njets, delta_r, variation});
            Logger::get("SM FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("SM FakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);

            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_lep_iso_corr =
                qcd_lep_iso_closure->evaluate({iso_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep iso correction {}", qcd_lep_iso_corr);
            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);

            ff = qcd_frac * qcd_ff * qcd_lep_pt_corr * qcd_lep_iso_corr *
                     qcd_DR_SR_corr +
                 wjets_frac * wjets_ff * wjets_lep_pt_corr * wjets_DR_SR_corr +
                 ttbar_frac * ttbar_ff * ttbar_lep_pt_corr;
        }

        Logger::get("SM FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 =
        df.Define(outputname, calc_fake_factor,
                  {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis, delta_r});
    return df1;
}

} // namespace fakefactors
#endif /* GUARDFAKEFACTORS_H */