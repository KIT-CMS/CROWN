#ifndef GUARDFAKEFACTORS_H
#define GUARDFAKEFACTORS_H
/// The namespace that contains the fake factor function.
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "correction.h"

namespace fakefactors {

/// Function to calculate raw fake factors without corrections with
/// correctionlib
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the output column for the fake factor
/// \param tau_pt pt of the hadronic tau in the tau pair
/// \param njets number of good jets in the event
/// \param lep_mt transverse mass of the leptonic tau in the tau pair
/// \param nbtags number of good b-tagged jets in the event
/// \param variation name of the uncertainty variation or nominal
/// \param ff_file correctionlib json file with the fake factors
///
/// \returns a dataframe with the fake factors
ROOT::RDF::RNode
raw_fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &tau_pt, const std::string &njets,
                        const std::string &lep_mt, const std::string &nbtags,
                        const std::string &variation,
                        const std::string &ff_file) {
    auto calc_fake_factor = [variation,
                             ff_file](const float &pt_2, const int &njets,
                                      const float &mt_1, const int &nbtag) {
        Logger::get("RawFakeFactor")
            ->debug("Setting up functions for raw fake factor (without "
                    "corrections) evaluation with correctionlib");
        Logger::get("RawFakeFactor")->debug("Variation - Name {}", variation);
        auto qcd = correction::CorrectionSet::from_file(ff_file)->at(
            "QCD_fake_factors");
        auto wjets = correction::CorrectionSet::from_file(ff_file)->at(
            "Wjets_fake_factors");
        auto ttbar = correction::CorrectionSet::from_file(ff_file)->at(
            "ttbar_fake_factors");
        auto fractions = correction::CorrectionSet::from_file(ff_file)->at(
            "process_fractions");

        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("RawFakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff = wjets->evaluate({pt_2, (float)njets, variation});
            Logger::get("RawFakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("RawFakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("RawFakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("RawFakeFactor")->debug("N b-jets - value {}", nbtag);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)nbtag, variation});
            Logger::get("RawFakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)nbtag, variation});
            Logger::get("RawFakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)nbtag, variation});
            Logger::get("RawFakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);

            ff = qcd_frac * qcd_ff + wjets_frac * wjets_ff +
                 ttbar_frac * ttbar_ff;
        }

        Logger::get("RawFakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {tau_pt, njets, lep_mt, nbtags});
    return df1;
}

/// Function to calculate fake factors with correctionlib
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the output column for the fake factor
/// \param tau_pt pt of the hadronic tau in the tau pair
/// \param njets number of good jets in the event
/// \param lep_mt transverse mass of the leptonic tau in the tau pair
/// \param nbtags number of good b-tagged jets in the event
/// \param lep_pt pt of the leptonic tau in the tau pair
/// \param lep_iso isolation of the leptonic tau in the tau pair
/// \param m_vis visible mass of the tau pair
/// \param variation name of the uncertainty variation or nominal
/// \param ff_file correctionlib json file with the fake factors
/// \param ff_corr_file correctionlib json file with corrections for the fake
/// factors
///
/// \returns a dataframe with the fake factors
ROOT::RDF::RNode
fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tau_pt, const std::string &njets,
                    const std::string &lep_mt, const std::string &nbtags,
                    const std::string &lep_pt, const std::string &lep_iso,
                    const std::string &m_vis, const std::string &variation,
                    const std::string &ff_file,
                    const std::string &ff_corr_file) {
    auto calc_fake_factor = [variation, ff_file, ff_corr_file](
                                const float &pt_2, const int &njets,
                                const float &mt_1, const int &nbtag,
                                const float &pt_1, const float &iso_1,
                                const float &m_vis) {
        Logger::get("FakeFactor")
            ->debug("Setting up functions for fake factor evaluation with "
                    "correctionlib");
        Logger::get("FakeFactor")->debug("Variation - Name {}", variation);
        auto qcd = correction::CorrectionSet::from_file(ff_file)->at(
            "QCD_fake_factors");
        auto wjets = correction::CorrectionSet::from_file(ff_file)->at(
            "Wjets_fake_factors");
        auto ttbar = correction::CorrectionSet::from_file(ff_file)->at(
            "ttbar_fake_factors");
        auto fractions = correction::CorrectionSet::from_file(ff_file)->at(
            "process_fractions");

        auto qcd_lep_pt_closure =
            correction::CorrectionSet::from_file(ff_corr_file)
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

        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff = wjets->evaluate({pt_2, (float)njets, variation});
            Logger::get("FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("FakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("FakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("FakeFactor")->debug("N b-jets - value {}", nbtag);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)nbtag, variation});
            Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)nbtag, variation});
            Logger::get("FakeFactor")->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)nbtag, variation});
            Logger::get("FakeFactor")->debug("ttbar - fraction {}", ttbar_frac);

            Logger::get("FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("FakeFactor")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_lep_iso_corr =
                qcd_lep_iso_closure->evaluate({iso_1, variation});
            Logger::get("FakeFactor")
                ->debug("QCD - lep iso correction {}", qcd_lep_iso_corr);
            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, variation});
            Logger::get("FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);

            ff = qcd_frac * qcd_ff * qcd_lep_pt_corr * qcd_lep_iso_corr *
                     qcd_DR_SR_corr +
                 wjets_frac * wjets_ff * wjets_lep_pt_corr * wjets_DR_SR_corr +
                 ttbar_frac * ttbar_ff * ttbar_lep_pt_corr;
        }

        Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 =
        df.Define(outputname, calc_fake_factor,
                  {tau_pt, njets, lep_mt, nbtags, lep_pt, lep_iso, m_vis});
    return df1;
}

} // namespace fakefactors
#endif /* GUARDFAKEFACTORS_H */