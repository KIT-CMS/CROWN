// #ifndef GUARDFAKEFACTORS_H
// #define GUARDFAKEFACTORS_H
// /// The namespace that contains the fake factor function.
// #include "../include/event.hxx"
// #include "../include/utility/CorrectionManager.hxx"
// #include "../include/utility/Logger.hxx"
// #include "ROOT/RDataFrame.hxx"
// #include "correction.h"
// #include <vector>
// #include <string>
// #include <sstream>
// #include <iterator>
// #include <iostream>
// #include <random>

// namespace fakefactors {
//     /**
//      * @brief Function to join and replace characters in a vector of strings
//      * and return a single string.
//      * 
//      * @param strings Vector of strings to join.
//      * @param delimiter Delimiter to use for joining the strings.
//      * @return A single string with the joined and modified strings.
//      */
//     std::string joinAndReplace(const std::vector<std::string>& strings, const std::string& delimiter) {
//         static std::unordered_map<std::string, int> seen_strings;
//         std::ostringstream os;

//         for (const auto& str : strings) {
//             std::string modified_str = str;
//             std::replace(modified_str.begin(), modified_str.end(), '/', '_');
//             std::replace(modified_str.begin(), modified_str.end(), '.', '_');

//             os << modified_str << delimiter;
//         }

//         std::string result = os.str();
//         result = result.substr(0, result.size() - delimiter.size()); // Remove the trailing delimiter

//         if (seen_strings.find(result) != seen_strings.end()) {
//             seen_strings[result]++;
//             result += "_" + std::to_string(seen_strings[result]);
//         } else {
//             seen_strings[result] = 0;
//         }

//         return result;
//     }
//     namespace sm {
//         /**
//          * @brief Function to calculate raw fake factors without corrections with
//          * correctionlib for the semileptonic channels
//          * 
//          * @param df the input dataframe
//          * @param correctionManager the correction manager to load corrections
//          * @param outputname name of the output column for the fake factor
//          * @param pt_2 pt of the hadronic tau in the tau pair
//          * @param njets number of jets in the event
//          * @param delta_r delta R between the leptonic tau and the hadronic tau
//          * @param mt_1 transverse mass of the leptonic tau in the tau pair
//          * @param fraction_variation name of the uncertainty variation or nominal
//          * @param QCD_variation name of the uncertainty variation or nominal for QCD
//          * @param Wjets_variation name of the uncertainty variation or nominal for Wjets
//          * @param ttbar_variation name of the uncertainty variation or nominal for ttbar
//          * @param ff_file correctionlib json file with the fake factors
//          * @returns a dataframe with the fake factors
//          */
//         ROOT::RDF::RNode
//         raw_fakefactor_lt(
//             ROOT::RDF::RNode df,
//             correctionManager::CorrectionManager &correctionManager,
//             const std::string &outputname,
//             // for ff
//             const std::string &pt_2,
//             const std::string &njets,
//             //const std::string &delta_r,
//             // for fraction
//             const std::string &mt_1,
//             //
//             const std::string &fraction_variation,
//             const std::string &QCD_variation,
//             const std::string &Wjets_variation,
//             const std::string &ttbar_variation,
//             //
//             const std::string &ff_file
//         ) {
//             Logger::get("SM RawFakeFactor (lt)")->debug("Setting up functions for raw fake factor (without corrections) evaluation with correctionlib");
//             Logger::get("SM RawFakeFactor (lt)")->debug("Fraction variations: fraction={}, QCD={}, Wjets={}", fraction_variation, QCD_variation, Wjets_variation);

//             auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
//             auto wjets = correctionManager.loadCorrection(ff_file, "Wjets_fake_factors");
//             auto ttbar = correctionManager.loadCorrection(ff_file, "ttbar_fake_factors");
//             auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");

//             auto calc_fake_factor = [
//                 qcd, wjets, ttbar, fractions,
//                 QCD_variation, Wjets_variation, ttbar_variation, fraction_variation](
//                 const float &pt_2, const int &njets, const float &mt_1) {

//                 float ff = 0.0;

//                 float qcd_ff = 0.0, wjets_ff = 0.0, ttbar_ff = 0.0;
//                 float qcd_frac = 0.0, wjets_frac = 0.0, ttbar_frac = 0.0;

//                 if (pt_2 >= 0.) {
//                     Logger::get("SM RawFakeFactor (lt)")->debug("pt_tau={}, njets={}, mt={}", pt_2, njets, mt_1);

//                     qcd_ff = qcd->evaluate({pt_2, (float)njets, QCD_variation});
//                     wjets_ff = wjets->evaluate({pt_2, (float)njets, Wjets_variation});
//                     ttbar_ff = ttbar->evaluate({pt_2, (float)njets, ttbar_variation});

//                     Logger::get("SM RawFakeFactor (lt)")->debug("RawFakeFactor (lt) - QCD={}, Wjets={}", qcd_ff, wjets_ff);

//                     qcd_frac = fractions->evaluate({"QCD", mt_1, (float)njets, fraction_variation});
//                     wjets_frac = fractions->evaluate({"Wjets", mt_1, (float)njets, fraction_variation});
//                     ttbar_frac = fractions->evaluate({"ttbar", mt_1, (float)njets, fraction_variation});

//                     Logger::get("SM RawFakeFactor (lt)")->debug("Fractions: QCD={}, Wjets={}", qcd_frac, wjets_frac);

//                     ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) + 
//                         std::max(wjets_frac, (float)0.) * std::max(wjets_ff, (float)0.) +
//                         std::max(ttbar_frac, (float)0.) * std::max(ttbar_ff, (float)0.);
//                 }

//                 Logger::get("SM RawFakeFactor (lt)")->debug("Event Fake Factor {}", ff);

//                 return ff;
//             };

//             auto df1 = df.Define(outputname, calc_fake_factor, {pt_2, njets, mt_1});

//             return df1;
//         }
//         /**
//          * @brief Function to calculate fake factors with corrections with correctionlib
//          * for the semileptonic channels
//          * 
//          * @param df the input dataframe
//          * @param correctionManager the correction manager to load corrections
//          * @param outputname name of the output column for the fake factor
//          * @param pt_2 pt of the hadronic tau in the tau pair
//          * @param njets number of jets in the event
//          * @param delta_r delta R between the leptonic tau and the hadronic tau
//          * @param mt_1 transverse mass of the leptonic tau in the tau pair
//          * @param pt_1 pt of the leptonic tau in the tau pair
//          * @param decaymode_2 decay mode of the hadronic tau in the tau pair
//          * @param lep_iso isolation of the leptonic tau in the tau pair
//          * @param mass_2 mass of the hadronic tau in the tau pair
//          * @param mt_tot total transverse mass of the event
//          * @param met missing transverse energy of the event
//          * @param pt_tt visible di-tau mass of the tau pair
//          * @param fraction_variation name of the uncertainty variation or nominal
//          * @param QCD_variation name of the uncertainty variation or nominal for QCD
//          * @param Wjets_variation name of the uncertainty variation or nominal for Wjets
//          * @param ttbar_variation name of the uncertainty variation or nominal for ttbar
//          * @param QCD_DR_SR_correction_variation name of the uncertainty variation or nominal
//          * for the QCD DR to SR correction
//          * @param QCD_non_closure_correction_variation name of the uncertainty variation or
//          * nominal for the QCD non-closure correction
//          * @param Wjets_DR_SR_correction_variation name of the uncertainty variation or nominal
//          * for the Wjets DR to SR correction
//          * @param Wjets_non_closure_correction_variation name of the uncertainty variation or
//          * nominal for the Wjets non-closure correction
//          * @param ttbar_non_closure_correction_variation name of the uncertainty variation or
//          * nominal for the ttbar non-closure correction
//          * @param ff_file correctionlib json file with the fake factors
//          * @param ff_corr_file correctionlib json file with corrections for the fake factors
//          * @returns a dataframe with the fake factors
//          */
//         ROOT::RDF::RNode
//         fakefactor_lt(
//             ROOT::RDF::RNode df, 
//             correctionManager::CorrectionManager &correctionManager,
//             const std::string &outputname,
//             // for ff
//             const std::string &pt_2,
//             const std::string &njets,
//             // for fraction 
//             const std::string &mt_1,
//             // for non closure corrections
//             const std::string &delta_r,
//             const std::string &decaymode_2,
//             const std::string &mass_2,
//             const std::string &mt_tot,
//             const std::string &met,
//             // for DR SR corrections
//             const std::string &pt_tt,
//             // for corrections
//             const std::string &fraction_variation,
//             const std::string &QCD_variation,
//             const std::string &Wjets_variation,
//             const std::string &ttbar_variation,
//             //
//             const std::string &QCD_DR_SR_correction_variation,
//             const std::string &QCD_non_closure_correction_variation,
//             //
//             const std::string &Wjets_DR_SR_correction_variation,
//             const std::string &Wjets_non_closure_correction_variation,
//             //
//             const std::string &ttbar_non_closure_correction_variation,
//             //
//             const std::string &ff_file,
//             const std::string &ff_corr_file
//             ) {

//             Logger::get("SM FakeFactor (lt)")->debug("Setting up functions for fake factor evaluation with correctionlib");

//             Logger::get("SM FakeFactor (lt)")->debug("Fraction variations: fraction={}, QCD={}, Wjets={})", fraction_variation, QCD_variation, Wjets_variation);
//             Logger::get("SM FakeFactor (lt)")->debug("Correction variations: QCD_DR_SR={}, QCD_non_closure={})", QCD_DR_SR_correction_variation, QCD_non_closure_correction_variation);
//             Logger::get("SM FakeFactor (lt)")->debug("Correction variations: Wjets_DR_SR={}, Wjets_non_closure={})", Wjets_DR_SR_correction_variation, Wjets_non_closure_correction_variation);
//             Logger::get("SM FakeFactor (lt)")->debug("Correction variations: ttbar_non_closure={})", ttbar_non_closure_correction_variation);

//             auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
//             auto wjets = correctionManager.loadCorrection(ff_file, "Wjets_fake_factors");
//             auto ttbar = correctionManager.loadCorrection(ff_file, "ttbar_fake_factors");
//             auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");

//             auto qcd_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_DR_SR_correction");
//             auto qcd_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "QCD_compound_correction");

//             auto wjets_DR_SR = correctionManager.loadCorrection(ff_corr_file, "Wjets_DR_SR_correction");
//             auto wjets_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "Wjets_compound_correction");

//             auto ttbar_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "ttbar_compound_correction");

//             auto calc_fake_factor = [
//                 qcd, wjets, fractions,
//                 QCD_variation, Wjets_variation, fraction_variation,
//                 qcd_DR_SR, qcd_non_closure,
//                 QCD_DR_SR_correction_variation, QCD_non_closure_correction_variation,
//                 wjets_DR_SR, wjets_non_closure,
//                 Wjets_DR_SR_correction_variation, Wjets_non_closure_correction_variation,
//                 ttbar, ttbar_variation, ttbar_non_closure_correction_variation, ttbar_non_closure](
//                 const float &pt_2,
//                 const int &njets,
//                 const float &mt_1,
//                 const float &delta_r,
//                 const int &decaymode_2,
//                 const float &mass_2,
//                 const float &mt_tot,
//                 const float &met,
//                 const float &pt_tt) {

//                 float ff = 0.0;

//                 float qcd_ff = 0.0, wjets_ff = 0.0, ttbar_ff = 0.0;
//                 float qcd_frac = 0.0, wjets_frac = 0.0, ttbar_frac = 0.0;
//                 float qcd_DR_SR_corr = 0., qcd_non_closure_corr = 0.;
//                 float wjets_DR_SR_corr = 0., wjets_non_closure_corr = 0.;
//                 float ttbar_non_closure_corr = 0.;

//                 float qcd_correction = 0.0, wjets_correction = 0.0, ttbar_correction = 0.0;

//                 if (pt_2 >= 0.) {
//                     Logger::get("SM FakeFactor (lt)")->debug("pt_tau={}, njets={}, mt={}", pt_2, njets, mt_1);

//                     qcd_ff = qcd->evaluate({pt_2, (float)njets, QCD_variation});
//                     wjets_ff = wjets->evaluate({pt_2, (float)njets, Wjets_variation});
//                     ttbar_ff = ttbar->evaluate({pt_2, (float)njets, ttbar_variation});

//                     Logger::get("SM FakeFactor (lt)")->debug("fake factors: QCD={}, Wjets={}", qcd_ff, wjets_ff);

//                     qcd_frac = fractions->evaluate({"QCD", mt_1, (float)njets, fraction_variation});
//                     wjets_frac = fractions->evaluate({"Wjets", mt_1, (float)njets, fraction_variation});
//                     ttbar_frac = fractions->evaluate({"ttbar", mt_1, (float)njets, fraction_variation});

//                     Logger::get("SM FakeFactor (lt)")->debug("fractions: QCD={}, Wjets={}", qcd_frac, wjets_frac);

//                     qcd_DR_SR_corr = qcd_DR_SR->evaluate({pt_tt, (float)njets, QCD_DR_SR_correction_variation});
//                     //qcd_DR_SR_corr = 1.0;  // Ignoring QCD DR_SR correction for now
//                     qcd_non_closure_corr = qcd_non_closure->evaluate(
//                         {
//                             //delta_r,
//                             (float)decaymode_2,
//                             pt_1,
//                             mass_2,
//                             mt_tot,
//                             met,
//                             (float)njets,
//                             QCD_non_closure_correction_variation
//                         }
//                     );

//                     Logger::get("SM FakeFactor (lt)")->debug("QCD: DR_SR={}, non_closure={}", qcd_DR_SR_corr, qcd_non_closure_corr);

//                     wjets_DR_SR_corr = wjets_DR_SR->evaluate({pt_tt, (float)njets, Wjets_DR_SR_correction_variation});
//                     //wjets_DR_SR_corr = 1.0;  // Ignoring Wjets DR_SR correction for now
//                     wjets_non_closure_corr = wjets_non_closure->evaluate(
//                         {
//                             //delta_r,
//                             (float)decaymode_2,
//                             mass_2,
//                             mt_tot,
//                             met,
//                             (float)njets,
//                             Wjets_non_closure_correction_variation
//                         }
//                     );

//                     Logger::get("SM FakeFactor (lt)")->debug("Wjets: DR_SR={}, non_closure={}", wjets_DR_SR_corr, wjets_non_closure_corr);

//                     ttbar_non_closure_corr = ttbar_non_closure->evaluate(
//                         {
//                             //delta_r,
//                             mass_2,
//                             mt_tot,
//                             met,
//                             (float)decaymode_2,
//                             (float)njets,
//                             ttbar_non_closure_correction_variation
//                         }
//                     );

//                     Logger::get("SM FakeFactor (lt)")->debug("ttbar: non_closure={}", ttbar_non_closure_corr);

//                     qcd_correction = std::max(qcd_DR_SR_corr, (float)0.) * std::max(qcd_non_closure_corr, (float)0.);
//                     wjets_correction = std::max(wjets_DR_SR_corr, (float)0.) * std::max(wjets_non_closure_corr, (float)0.);
//                     ttbar_correction = std::max(ttbar_non_closure_corr, (float)0.);

//                     ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) * qcd_correction +
//                         std::max(wjets_frac, (float)0.) * std::max(wjets_ff, (float)0.) * wjets_correction +
//                         std::max(ttbar_frac, (float)0.) * std::max(ttbar_ff, (float)0.) * ttbar_correction;

//                 }

//                 Logger::get("SM FakeFactor (lt)")->debug("Event Fake Factor {}", ff);
                
//                 return ff;
//             };

//             auto df1 = df.Define(outputname, calc_fake_factor, {pt_2, njets, mt_1, delta_r, decaymode_2, mass_2, mt_tot, met, pt_tt});

//             return df1;
//         }
//         /**
//          * @brief Function to calculate fake factors with corrections with correctionlib
//          * for the semileptonic channels, splitting the information for further analysis
//          * 
//          * @param df the input dataframe
//          * @param correctionManager the correction manager to load corrections
//          * @param outputname name of the output column for the fake factor
//          * @param pt_2 pt of the hadronic tau in the tau pair
//          * @param njets number of jets in the event
//          * @param delta_r delta R between the leptonic tau and the hadronic tau
//          * @param mt_1 transverse mass of the leptonic tau in the tau pair
//          * @param pt_1 pt of the leptonic tau in the tau pair
//          * @param decaymode_2 decay mode of the hadronic tau in the tau pair
//          * @param lep_iso isolation of the leptonic tau in the tau pair
//          * @param mass_2 mass of the hadronic tau in the tau pair
//          * @param mt_tot total transverse mass of the event
//          * @param met missing transverse energy of the event
//          * @param pt_tt visible di-tau mass of the tau pair
//          * @param fraction_variation name of the uncertainty variation or nominal
//          * @param QCD_variation name of the uncertainty variation or nominal for QCD
//          * @param Wjets_variation name of the uncertainty variation or nominal for Wjets
//          * @param ttbar_variation name of the uncertainty variation or nominal for ttbar
//          * @param QCD_DR_SR_correction_variation name of the uncertainty variation or nominal
//          * for the QCD DR to SR correction
//          * @param QCD_non_closure_correction_variation name of the uncertainty variation or
//          * nominal for the QCD non-closure correction
//          * @param Wjets_DR_SR_correction_variation name of the uncertainty variation or nominal
//          * for the Wjets DR to SR correction
//          * @param Wjets_non_closure_correction_variation name of the uncertainty variation or
//          * nominal for the Wjets non-closure correction
//          * @param ttbar_non_closure_correction_variation name of the uncertainty variation or
//          * nominal for the ttbar non-closure correction
//          * @param ff_file correctionlib json file with the fake factors
//          * @param ff_corr_file correctionlib json file with corrections for the fake factors
//          * @returns a dataframe with the fake factors and additional split information
//          */
//         ROOT::RDF::RNode
//         fakefactor_lt_split_info(
//             ROOT::RDF::RNode df, 
//             correctionManager::CorrectionManager &correctionManager,
//             const std::vector<std::string> &outputname,
//                 // for ff
//             const std::string &pt_2,
//             const std::string &njets,
//             // for fraction 
//             const std::string &mt_1,
//             // for non closure corrections
//             const std::string &delta_r,
//             const std::string &decaymode_2,
//             const std::string &mass_2,
//             const std::string &mt_tot,
//             const std::string &met,
//             // for DR SR corrections
//             const std::string &pt_tt,
//             // for corrections
//             const std::string &fraction_variation,
//             const std::string &QCD_variation,
//             const std::string &Wjets_variation,
//             const std::string &ttbar_variation,
//             //
//             const std::string &QCD_DR_SR_correction_variation,
//             const std::string &QCD_non_closure_correction_variation,
//             //
//             const std::string &Wjets_DR_SR_correction_variation,
//             const std::string &Wjets_non_closure_correction_variation,
//             //
//             const std::string &ttbar_non_closure_correction_variation,
//             //
//             const std::string &ff_file,
//             const std::string &ff_corr_file
//             ) {

//             auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
//             auto wjets = correctionManager.loadCorrection(ff_file, "Wjets_fake_factors");
//             auto ttbar = correctionManager.loadCorrection(ff_file, "ttbar_fake_factors");
//             auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");

//             auto qcd_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_DR_SR_correction");
//             auto qcd_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "QCD_compound_correction");

//             auto wjets_DR_SR = correctionManager.loadCorrection(ff_corr_file, "Wjets_DR_SR_correction");
//             auto wjets_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "Wjets_compound_correction");

//             auto ttbar_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "ttbar_compound_correction");

//             auto calc_fake_factor = [
//                 qcd, wjets, fractions,
//                 QCD_variation, Wjets_variation, fraction_variation,
//                 qcd_DR_SR, qcd_non_closure,
//                 QCD_DR_SR_correction_variation, QCD_non_closure_correction_variation,
//                 wjets_DR_SR, wjets_non_closure,
//                 Wjets_DR_SR_correction_variation, Wjets_non_closure_correction_variation,
//                 ttbar, ttbar_variation, ttbar_non_closure_correction_variation, ttbar_non_closure](
//                 const float &pt_2,
//                 const int &njets,
//                 const float &mt_1,
//                 const float &delta_r,
//                 const int &decaymode_2,
//                 const float &mass_2,
//                 const float &mt_tot,
//                 const float &met,
//                 const float &pt_tt) {

//                 float qcd_ff = 0.0, wjets_ff = 0.0, ttbar_ff = 0.0;
//                 float qcd_frac = 0.0, wjets_frac = 0.0, ttbar_frac = 0.0;
//                 float qcd_DR_SR_corr = 0., qcd_non_closure_corr = 0.;
//                 float wjets_DR_SR_corr = 0., wjets_non_closure_corr = 0.;
//                 float ttbar_non_closure_corr = 0.;

//                 float qcd_correction = 0.0, wjets_correction = 0.0, ttbar_correction = 0.0;


//                 if (pt_2 >= 0.) {
//                     Logger::get("SM FakeFactor split (lt)")->debug("pt_tau={}, njets={}, mt={}", pt_2, njets, mt_1);

//                     qcd_ff = qcd->evaluate({pt_2, (float)njets, QCD_variation});
//                     wjets_ff = wjets->evaluate({pt_2, (float)njets, Wjets_variation});
//                     ttbar_ff = ttbar->evaluate({pt_2, (float)njets, ttbar_variation});

//                     Logger::get("SM FakeFactor split (lt)")->debug("fake factors: QCD={}, Wjets={}", qcd_ff, wjets_ff);

//                     qcd_frac = fractions->evaluate({"QCD", mt_1, (float)njets, fraction_variation});
//                     wjets_frac = fractions->evaluate({"Wjets", mt_1, (float)njets, fraction_variation});
//                     ttbar_frac = fractions->evaluate({"ttbar", mt_1, (float)njets, fraction_variation});

//                     Logger::get("SM FakeFactor split (lt)")->debug("fractions: QCD={}, Wjets={}", qcd_frac, wjets_frac);

//                     qcd_DR_SR_corr = qcd_DR_SR->evaluate({pt_tt, (float)njets, QCD_DR_SR_correction_variation});
//                     qcd_non_closure_corr = qcd_non_closure->evaluate(
//                         {
//                             //delta_r,
//                             mass_2,
//                             mt_tot, 
//                             met,
//                             (float)decaymode_2,
//                             (float)njets,
//                             QCD_non_closure_correction_variation
//                         }
//                     );

//                     Logger::get("SM FakeFactor split (lt)")->debug("QCD: DR_SR={}, non_closure={}", qcd_DR_SR_corr, qcd_non_closure_corr);

//                     wjets_DR_SR_corr = wjets_DR_SR->evaluate({pt_tt, (float)njets, Wjets_DR_SR_correction_variation});
//                     //wjets_DR_SR_corr = 1.0;
//                     wjets_non_closure_corr = wjets_non_closure->evaluate(
//                         {
//                             //delta_r,
//                             mass_2,
//                             mt_tot, 
//                             met,
//                             (float)decaymode_2,
//                             (float)njets,
//                             Wjets_non_closure_correction_variation
//                         }
//                     );

//                     Logger::get("SM FakeFactor split (lt)")->debug("Wjets: DR_SR={}, non_closure={}", wjets_DR_SR_corr, wjets_non_closure_corr);

//                     ttbar_non_closure_corr = ttbar_non_closure->evaluate(
//                         {
//                             //delta_r,
//                             mass_2,
//                             mt_tot, 
//                             met,
//                             (float)decaymode_2,
//                             (float)njets,
//                             ttbar_non_closure_correction_variation
//                         }
//                     );

//                     Logger::get("SM FakeFactor split (lt)")->debug("ttbar: non_closure={}", ttbar_non_closure_corr);

//                     qcd_correction = std::max(qcd_DR_SR_corr, (float)0.) * std::max(qcd_non_closure_corr, (float)0.);
//                     wjets_correction = std::max(wjets_DR_SR_corr, (float)0.) * std::max(wjets_non_closure_corr, (float)0.);
//                     ttbar_correction = std::max(ttbar_non_closure_corr, (float)0.);
//                 }

//                 // all of them process wise
//                 // raw_ff, factions, DR_SR, correction_wo_DR_SR, combined_correction, ff
//                 std::vector<float> result = {
//                     std::max(qcd_ff, (float)0.),
//                     std::max(wjets_ff, (float)0.),
//                     std::max(ttbar_ff, (float)0.),
//                     //
//                     std::max(qcd_frac, (float)0.),
//                     std::max(wjets_frac, (float)0.),
//                     std::max(ttbar_frac, (float)0.),
//                     //
//                     std::max(qcd_DR_SR_corr, (float)0.),
//                     std::max(wjets_DR_SR_corr, (float)0.),
//                     //
//                     std::max(qcd_non_closure_corr, (float)0.),
//                     std::max(wjets_non_closure_corr, (float)0.),
//                     std::max(ttbar_non_closure_corr, (float)0.),
//                     //
//                     std::max(qcd_correction, (float)0.),
//                     std::max(wjets_correction, (float)0.),
//                     std::max(ttbar_correction, (float)0.),
//                     //
//                     std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) * std::max(qcd_correction, (float)0.),
//                     std::max(wjets_frac, (float)0.) * std::max(wjets_ff, (float)0.) * std::max(wjets_correction, (float)0.),
//                     std::max(ttbar_frac, (float)0.) * std::max(ttbar_ff, (float)0.) * std::max(ttbar_correction, (float)0.)};
        
//                 return result;
//             };

//             std::vector<std::string> strings = {
//                 "fakefactor_lt_split_info",
//                 fraction_variation,
//                 QCD_variation,
//                 Wjets_variation,
//                 ttbar_variation,
//                 QCD_DR_SR_correction_variation,
//                 QCD_non_closure_correction_variation,
//                 Wjets_DR_SR_correction_variation,
//                 Wjets_non_closure_correction_variation,
//                 ttbar_non_closure_correction_variation,
//                 ff_file,
//                 ff_corr_file};
        
//             std::string shifted_collection_identifier =  fakefactors::joinAndReplace(strings, "_");

//             auto df1 = df.Define(shifted_collection_identifier, calc_fake_factor, {pt_2, njets, mt_1, delta_r, decaymode_2, mass_2, mt_tot, met, pt_tt});
//             auto df2 = event::quantity::Unroll<float>(df1, outputname, shifted_collection_identifier);

//             return df2;
//         }

//         /**
//         * @brief Function to calculate raw fake factors without corrections with
//         * correctionlib for the NMSSM Di-Higgs analysis for the full hadronic channel
//         *
//         * @param df the dataframe to add the quantity to
//         * @param outputname name of the output column for the fake factor
//         * @param tau_idx index of the tau, leading/subleading
//         * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
//         * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
//         * @param njets number of good jets in the event
//         * @param pt_tt visible di-tau mass of the tau pair
//         * @param qcd_variation name of the QCD FF uncertainty variation or nominal
//         * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
//         * @param fraction_variation name of the process fraction uncertainty variation or nominal
//         * @param ff_file correctionlib json file with the fake factors
//         * @returns a dataframe with the fake factors
//         */
//         ROOT::RDF::RNode
//         raw_fakefactor_tt(ROOT::RDF::RNode df, correctionManager::CorrectionManager &correctionManager,
//                                 const std::string &outputname,
//                                 const int &tau_idx,
//                                  const std::string &pt_1,
//                                 const std::string &pt_2,
//                                  const std::string &njets, 
//                                 const std::string &pt_tt,
//                                 const std::string &qcd_variation,
//                                 const std::string &fraction_variation, const std::string &ff_file) {

//             Logger::get("RawFakeFactor")->debug("Setting up functions for raw fake factor (without corrections) evaluation with correctionlib");
//             Logger::get("RawFakeFactor")->debug("QCD variation - Name {}", qcd_variation);
//             Logger::get("RawFakeFactor")->debug("Fraction variation - Name {}", fraction_variation);

//             auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
//             auto qcd_subleading = correctionManager.loadCorrection(ff_file, "QCD_subleading_fake_factors");

//             auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");
//             auto fractions_subleading = correctionManager.loadCorrection(ff_file, "process_fractions_subleading");

//             auto calc_fake_factor = [tau_idx, qcd_variation, fraction_variation, qcd, qcd_subleading, fractions, fractions_subleading](
//                                         const float &pt_1, const float &pt_2,
//                                         const int &njets, const float &pt_tt) {
                
//                 float ff = 0.;
                
//                 if (pt_2 >= 0.) {
//                     Logger::get("RawFakeFactor")
//                         ->debug("Leading Tau pt - value {}", pt_1);
//                     Logger::get("RawFakeFactor")
//                         ->debug("Subleading Tau pt - value {}", pt_2);
//                     Logger::get("RawFakeFactor")->debug("N jets - value {}", njets);

//                     float qcd_ff = 0.;
//                     float qcd_frac = 0.;
//                     if (tau_idx == 0) {
//                         qcd_ff = qcd->evaluate({pt_1, (float)njets, qcd_variation});
//                         Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
//                         qcd_frac = fractions->evaluate({"QCD", pt_tt, (float)njets, fraction_variation});
//                         Logger::get("RawFakeFactor")->debug("QCD - fraction {}", qcd_frac);
                        
//                         ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.);

//                     } else if (tau_idx == 1) {
//                         qcd_ff = qcd_subleading->evaluate({pt_2, (float)njets, qcd_variation});
//                         Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
//                         qcd_frac = fractions_subleading->evaluate({"QCD", pt_tt, (float)njets, fraction_variation});
                        
//                         ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.);
//                     }
//                 }

//                 Logger::get("RawFakeFactor")->debug("Event Fake Factor {}", ff);
//                 return ff;
//             };
//             auto df1 =
//                 df.Define(outputname, calc_fake_factor, {pt_1, pt_2, njets, pt_tt});
//             return df1;
//         }

//         /**
//         * @brief Function to calculate fake factors with correctionlib
//         *
//         * @param df the dataframe to add the quantity to
//         * @param outputname name of the output column for the fake factor
//         * @param tau_idx index of the tau, leading/subleading
//         * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
//         * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
//         * @param njets number of good jets in the event
//         * @param pt_tt visible mass of the tau pair
//         * @param tau_mass_1 mass of the leading hadronic tau in the tau pair
//         * @param tau_mass_2 mass of the subleading hadronic tau in the tau pair
//         * @param mt_tot total transverse mass of the event
//         * @param met missing transverse energy of the event
//         * @param qcd_variation name of the QCD FF uncertainty variation or nominal
//         * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
//         * @param fraction_variation name of the process fraction uncertainty variation or nominal
//         * @param qcd_non_closure_correction_variation name of the QCD lepton pt correction uncertainty variation or nominal
//         * @param qcd_corr_taumass_variation name of the QCD lepton mass correction uncertainty variation or nominal
//         * @param qcd_DR_SR_correction_variation name of the QCD DR to SR correction uncertainty variation or nominal
//         * @param ttbar_non_closure_correction_variation name of the ttbar lepton pt correction uncertainty variation or nominal
//         * @param ttbar_corr_taumass_variation name of the ttbar lepton mass correction uncertainty variation or nominal
//         * @param ff_file correctionlib json file with the fake factors
//         * @param ff_corr_file correctionlib json file with corrections for the fake
//         * factors
//         * @returns a dataframe with the fake factors
//         */
//         ROOT::RDF::RNode
//         fakefactor_tt(ROOT::RDF::RNode df, correctionManager::CorrectionManager &correctionManager,
//                             const std::string &outputname,
//                             const int &tau_idx,
//                             const std::string &tau_pt_1,
//                             const std::string &tau_pt_2,
//                             const std::string &mass_1,
//                             const std::string &mass_2,
//                             const std::string &mt_tot,
//                             const std::string &met,
//                             const std::string &decaymode_1,
//                             const std::string &decaymode_2,
//                             const std::string &njets,
//                             const std::string &pt_tt,
//                             const std::string &m_vis,
//                             const std::string &qcd_variation, 
//                             const std::string &fraction_variation, const std::string &qcd_non_closure_correction_variation,
//                             const std::string &qcd_DR_SR_correction_variation,
//                             const std::string &ff_file, const std::string &ff_corr_file) {

//             Logger::get("FakeFactor")
//                 ->debug("Setting up functions for fake factor evaluation with "
//                         "correctionlib");
//             Logger::get("FakeFactor")->debug("QCD variation - Name {}", qcd_variation);
//             Logger::get("FakeFactor")->debug("Fraction variation - Name {}", fraction_variation);
//             Logger::get("FakeFactor")->debug("QCD lepton pt variation - Name {}", qcd_non_closure_correction_variation);
//             Logger::get("FakeFactor")->debug("QCD DRSR variation - Name {}", qcd_DR_SR_correction_variation);

//             auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
//             auto qcd_subleading = correctionManager.loadCorrection(ff_file, "QCD_subleading_fake_factors");
//             auto fractions =correctionManager.loadCorrection(ff_file, "process_fractions");
//             auto fractions_subleading = correctionManager.loadCorrection(ff_file, "process_fractions_subleading");

//             auto qcd_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "QCD_compound_correction");
//             auto qcd_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_DR_SR_correction");
//             auto qcd_subleading_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "QCD_subleading_compound_correction");
//             auto qcd_subleading_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_subleading_DR_SR_correction");
            
//             auto calc_fake_factor = [tau_idx, qcd_variation, fraction_variation, 
//                                     qcd_non_closure_correction_variation, qcd_DR_SR_correction_variation, 
//                                     qcd, fractions, qcd_non_closure, qcd_DR_SR,
//                                     qcd_subleading, fractions_subleading, qcd_subleading_non_closure, qcd_subleading_DR_SR](
//                                         const float &pt_1, const float &pt_2, const float &mass_1, const float &mass_2, const float &mt_tot, const float &met,
//                                         const int &decaymode_1, const int &decaymode_2,
//                                         const int &njets, const float &pt_tt, const float &m_vis) {
                
//                 float ff = 0.;
                
//                 if (pt_2 >= 0.) {
//                     Logger::get("FakeFactor")->debug("Leading Tau pt - value {}", pt_1);
//                     Logger::get("FakeFactor")->debug("Subleading Tau pt - value {}", pt_2);
//                     Logger::get("FakeFactor")->debug("pt_tt - value {}", pt_tt);
//                     Logger::get("FakeFactor")->debug("mt_tot - value {}", mt_tot);
//                     Logger::get("FakeFactor")->debug("met - value {}", met);
//                     Logger::get("FakeFactor")->debug("N jets - value {}", njets);

//                     float qcd_ff = 0.;
//                     float qcd_frac = 0.;
//                     float qcd_non_closure_corr = 0.;
//                     float qcd_DR_SR_corr = 0.;
//                     float qcd_correction = 0.0;

//                     if (tau_idx == 0) {
//                         qcd_ff = qcd->evaluate({pt_1, (float)njets, qcd_variation});
//                         Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
//                         qcd_frac = fractions->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);

//                         qcd_non_closure_corr = qcd_non_closure->evaluate({mass_1, mt_tot, met, (float)decaymode_1, (float)njets, qcd_non_closure_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - lep pt correction {}", qcd_non_closure_correction_variation);
//                         qcd_DR_SR_corr = qcd_DR_SR->evaluate({pt_tt, (float)njets, qcd_DR_SR_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - DR to SR correction {}", qcd_DR_SR_correction_variation);

//                         qcd_correction = std::max(qcd_DR_SR_corr, (float)0.) * std::max(qcd_non_closure_corr, (float)0.);

//                         ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) * qcd_correction;

//                     } else if (tau_idx == 1) {
//                         qcd_ff = qcd_subleading->evaluate({pt_2, (float)njets, qcd_variation});
//                         Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
//                         qcd_frac = fractions_subleading->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);

//                         qcd_non_closure_corr = qcd_subleading_non_closure->evaluate({mass_2, mt_tot, met, (float)decaymode_2, (float)njets, qcd_non_closure_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - lep pt correction {}", qcd_non_closure_correction_variation);
//                         qcd_DR_SR_corr = qcd_subleading_DR_SR->evaluate({pt_tt, (float)njets, qcd_DR_SR_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - DR to SR correction {}", qcd_DR_SR_correction_variation);

//                         qcd_correction = std::max(qcd_DR_SR_corr, (float)0.) * std::max(qcd_non_closure_corr, (float)0.);

//                         ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) * qcd_correction;
//                     }
//                 }

//                 Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);
//                 return ff;
//             };
//             auto df1 = df.Define(outputname, calc_fake_factor,
//                                 {tau_pt_1, tau_pt_2, mass_1, mass_2, mt_tot, met, decaymode_1, decaymode_2, njets, pt_tt, m_vis});
//             return df1;
//         }

//         /**
//         * @brief Function to calculate fake factors with correctionlib
//         *
//         * @param df the dataframe to add the quantity to
//         * @param outputname name of the output column for the fake factor
//         * @param tau_idx index of the tau, leading/subleading
//         * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
//         * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
//         * @param njets number of good jets in the event
//         * @param pt_tt visible mass of the tau pair
//         * @param mt_tot total transverse mass of the event
//         * @param met missing transverse energy of the event
//         * @param tau_mass_1 mass of the leading hadronic tau in the tau pair
//         * @param tau_mass_2 mass of the subleading hadronic tau in the tau pair
//         * @param qcd_variation name of the QCD FF uncertainty variation or nominal
//         * @param ttbar_variation name of the ttbar FF uncertainty variation or nominal
//         * @param fraction_variation name of the process fraction uncertainty variation or nominal
//         * @param qcd_non_closure_correction_variation name of the QCD lepton pt correction uncertainty variation or nominal
//         * @param qcd_corr_taumass_variation name of the QCD lepton mass correction uncertainty variation or nominal
//         * @param qcd_DR_SR_correction_variation name of the QCD DR to SR correction uncertainty variation or nominal
//         * @param ttbar_non_closure_correction_variation name of the ttbar lepton pt correction uncertainty variation or nominal
//         * @param ttbar_corr_taumass_variation name of the ttbar lepton mass correction uncertainty variation or nominal
//         * @param ff_file correctionlib json file with the fake factors
//         * @param ff_corr_file correctionlib json file with corrections for the fake
//         * factors
//         * @returns a dataframe with the fake factors
//         */
//         ROOT::RDF::RNode
//         fakefactor_tt_split_info(ROOT::RDF::RNode df, correctionManager::CorrectionManager &correctionManager,
//                             const std::vector<std::string> &outputname,
//                             const int &tau_idx,
//                             const std::string &tau_pt_1,
//                             const std::string &tau_pt_2, 
//                             const std::string &mass_1,
//                             const std::string &mass_2,
//                             const std::string &mt_tot,
//                             const std::string &met,
//                             const std::string &decaymode_1,
//                             const std::string &decaymode_2,
//                             const std::string &njets,
//                             const std::string &pt_tt,
//                             const std::string &m_vis,
//                             const std::string &qcd_variation, 
//                             const std::string &fraction_variation, const std::string &qcd_non_closure_correction_variation,
//                             const std::string &qcd_DR_SR_correction_variation,
//                             const std::string &ff_file, const std::string &ff_corr_file) {

//             Logger::get("FakeFactor")
//                 ->debug("Setting up functions for fake factor evaluation with "
//                         "correctionlib");
//             Logger::get("FakeFactor")->debug("QCD variation - Name {}", qcd_variation);
//             Logger::get("FakeFactor")->debug("Fraction variation - Name {}", fraction_variation);
//             Logger::get("FakeFactor")->debug("QCD lepton pt variation - Name {}", qcd_non_closure_correction_variation);
//             Logger::get("FakeFactor")->debug("QCD DRSR variation - Name {}", qcd_DR_SR_correction_variation);

//             auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
//             auto qcd_subleading = correctionManager.loadCorrection(ff_file, "QCD_subleading_fake_factors");
//             auto fractions =correctionManager.loadCorrection(ff_file, "process_fractions");
//             auto fractions_subleading = correctionManager.loadCorrection(ff_file, "process_fractions_subleading");

//             auto qcd_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "QCD_compound_correction");
//             auto qcd_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_DR_SR_correction");
//             auto qcd_subleading_non_closure = correctionManager.loadCompoundCorrection(ff_corr_file, "QCD_subleading_compound_correction");
//             auto qcd_subleading_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_subleading_DR_SR_correction");
            
//             auto calc_fake_factor = [tau_idx, qcd_variation, fraction_variation, 
//                                     qcd_non_closure_correction_variation, qcd_DR_SR_correction_variation, 
//                                     qcd, fractions, qcd_non_closure, qcd_DR_SR,
//                                     qcd_subleading, fractions_subleading, qcd_subleading_non_closure, qcd_subleading_DR_SR](
//                                         const float &pt_1, const float &pt_2,
//                                         const float &mass_1, const float &mass_2, const float &mt_tot, const float &met,
//                                         const int &decaymode_1, const int &decaymode_2,
//                                         const int &njets, const float &pt_tt,  const float &m_vis) {
                
//                 float ff = 0.;
//                 float qcd_ff = 0.;
//                 float qcd_frac = 0.;
//                 float qcd_non_closure_corr = 0.;
//                 float qcd_DR_SR_corr = 0.;
//                 float qcd_correction = 0.0;
                
//                 if (pt_2 >= 0.) {
//                     Logger::get("FakeFactor")->debug("Leading Tau pt - value {}", pt_1);
//                     Logger::get("FakeFactor")->debug("Subleading Tau pt - value {}", pt_2);
//                     Logger::get("FakeFactor")->debug("pt_tt - value {}", pt_tt);
//                     Logger::get("FakeFactor")->debug("mt_tot - value {}", mt_tot);
//                     Logger::get("FakeFactor")->debug("met - value {}", met);
//                     Logger::get("FakeFactor")->debug("N jets - value {}", njets);

//                     if (tau_idx == 0) {
//                         qcd_ff = qcd->evaluate({pt_1, (float)njets, qcd_variation});
//                         Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
//                         qcd_frac = fractions->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);

//                         qcd_non_closure_corr = qcd_non_closure->evaluate({mass_1, mt_tot, met, (float)decaymode_1, (float)njets, qcd_non_closure_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - lep pt correction {}", qcd_non_closure_correction_variation);
//                         qcd_DR_SR_corr = qcd_DR_SR->evaluate({pt_tt, (float)njets, qcd_DR_SR_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - DR to SR correction {}", qcd_DR_SR_correction_variation);

//                         qcd_correction = std::max(qcd_DR_SR_corr, (float)0.) * std::max(qcd_non_closure_corr, (float)0.);

//                         ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) * qcd_correction;

//                     } else if (tau_idx == 1) {
//                         qcd_ff = qcd_subleading->evaluate({pt_2, (float)njets, qcd_variation});
//                         Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
//                         qcd_frac = fractions_subleading->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);

//                         qcd_non_closure_corr = qcd_subleading_non_closure->evaluate({mass_2, mt_tot, met, (float)decaymode_2, (float)njets, qcd_non_closure_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - lep pt correction {}", qcd_non_closure_correction_variation);
//                         qcd_DR_SR_corr = qcd_subleading_DR_SR->evaluate({pt_tt, (float)njets, qcd_DR_SR_correction_variation});
//                         Logger::get("FakeFactor")->debug("QCD - DR to SR correction {}", qcd_DR_SR_correction_variation);

//                         qcd_correction = std::max(qcd_DR_SR_corr, (float)0.) * std::max(qcd_non_closure_corr, (float)0.);

//                         ff = std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) * qcd_correction;
//                     }
//                 }

//                 Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);

//                 // all of them process wise
//                 // raw_ff, factions, DR_SR, correction_wo_DR_SR, combined_correction, ff
//                 std::vector<float> result = {
//                     std::max(qcd_ff, (float)0.),
//                     //
//                     std::max(qcd_frac, (float)0.),
//                     //
//                     std::max(qcd_DR_SR_corr, (float)0.),
//                     //
//                     std::max(qcd_non_closure_corr, (float)0.),
//                     //
//                     std::max(qcd_correction, (float)0.),
//                     //
//                     std::max(qcd_frac, (float)0.) * std::max(qcd_ff, (float)0.) * std::max(qcd_correction, (float)0.)};
        
//                 return result;
//             };

//             std::vector<std::string> strings = {
//                 "fakefactor_tt_split_info",
//                 fraction_variation,
//                 qcd_variation,
//                 qcd_DR_SR_correction_variation,
//                 qcd_non_closure_correction_variation,
//                 ff_file,
//                 ff_corr_file};
        
//             std::string shifted_collection_identifier =  fakefactors::joinAndReplace(strings, "_");

//             auto df1 = df.Define(shifted_collection_identifier, calc_fake_factor, {tau_pt_1, tau_pt_2, mass_1, mass_2, mt_tot, met, decaymode_1, decaymode_2, njets, pt_tt, m_vis});
//             auto df2 = event::quantity::Unroll<float>(df1, outputname, shifted_collection_identifier);

//             return df2;
//     }

// } // namespace sm
// } // namespace fakefactors
// #endif /* GUARDFAKEFACTORS_H */