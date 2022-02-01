#include "ROOT/RDataFrame.hxx"
#include "RooFunctor.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "basefunctions.hxx"
#include "correction.h"
#include "utility/Logger.hxx"
#include "utility/RooFunctorThreadsafe.hxx"
/// namespace used for scale factor related functions
namespace scalefactor {
namespace muon {
/**
 * @brief Function used to evaluate id scale factors from muons
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param id_output name of the id scale factor column
 * @param workspace_name path to the Rooworkspace
 * @param id_functor_name name of the function from the workspace
 * @param id_arguments arguments of the function
 * @return a new dataframe containing the new column
 */
auto id(auto &df, const std::string &pt, const std::string &eta,
        const std::string &id_output, const std::string &workspace_name,
        const std::string &id_functor_name, const std::string &id_arguments) {

    Logger::get("muonsf")->debug("Setting up functions for muon sf");
    Logger::get("muonsf")->debug("ID - Function {} // argset {}",
                                 id_functor_name, id_arguments);

    const std::shared_ptr<RooFunctorThreadsafe> id_function =
        loadFunctor(workspace_name, id_functor_name, id_arguments);
    auto df1 = basefunctions::evaluateWorkspaceFunction(df, id_output,
                                                        id_function, pt, eta);
    return df1;
}
/**
 * @brief Function used to evaluate iso scale factors from muons
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param iso muon iso
 * @param iso_output name of the iso scale factor column
 * @param workspace_name path to the Rooworkspace
 * @param iso_functor_name name of the function from the workspace
 * @param iso_arguments arguments of the function
 * @return a new dataframe containing the new column
 */
auto iso(auto &df, const std::string &pt, const std::string &eta,
         const std::string &iso, const std::string &iso_output,
         const std::string &workspace_name, const std::string &iso_functor_name,
         const std::string &iso_arguments) {

    Logger::get("muonsf")->debug("Setting up functions for muon sf");
    Logger::get("muonsf")->debug("Iso - Function {} // argset {}",
                                 iso_functor_name, iso_arguments);

    const std::shared_ptr<RooFunctorThreadsafe> iso_function =
        loadFunctor(workspace_name, iso_functor_name, iso_arguments);
    auto df1 = basefunctions::evaluateWorkspaceFunction(
        df, iso_output, iso_function, pt, eta, iso);
    return df1;
}
/**
 * @brief Function used to evaluate id scale factors from muons with
 * correctionlib
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param year_id id for the year of data taking and mc compaign
 * @param variation id for the variation of the scale factor "sf" for nominal
 * and "systup"/"systdown" for up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon id scale factor
 * @return a new dataframe containing the new column
 */
auto id_ul(auto &df, const std::string &pt, const std::string &eta,
           const std::string &year_id, const std::string &variation,
           const std::string &id_output, const std::string &sf_file,
           const std::string &sf_name) {

    Logger::get("muonIdSF")->debug("Setting up functions for muon id sf");
    Logger::get("muonIdSF")->debug("ID - Name {}", sf_name);
    auto evaluator = correction::CorrectionSet::from_file(sf_file)->at(sf_name);
    auto df1 = df.Define(
        id_output,
        [evaluator, year_id, variation](const float &pt, const float &eta) {
            Logger::get("muonIdSF")->debug("ID - pt {}, eta {}", pt, eta);
            // preventing muons with default values due to tau energy correction
            // shifts below good tau pt selection
            if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                auto sf = evaluator->evaluate(
                    {year_id, std::abs(eta), pt, variation});
                return sf;
            } else
                return (double)1.;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function used to evaluate iso scale factors from muons with
 * correctionlib
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param year_id id for the year of data taking and mc compaign
 * @param variation id for the variation of the scale factor "sf" for nominal
 * and "systup"/"systdown" the up/down variation
 * @param iso_output name of the iso scale factor column
 * @param sf_file path to the file with the muon scale factors
 * @param sf_name name of the muon iso scale factor
 * @return a new dataframe containing the new column
 */
auto iso_ul(auto &df, const std::string &pt, const std::string &eta,
            const std::string &year_id, const std::string &variation,
            const std::string &iso_output, const std::string &sf_file,
            const std::string &sf_name) {

    Logger::get("muonIsoSF")->debug("Setting up functions for muon iso sf");
    Logger::get("muonIsoSF")->debug("ISO - Name {}", sf_name);
    auto evaluator = correction::CorrectionSet::from_file(sf_file)->at(sf_name);
    auto df1 = df.Define(
        iso_output,
        [evaluator, year_id, variation](const float &pt, const float &eta) {
            Logger::get("muonIsoSF")->debug("ISO - pt {}, eta {}", pt, eta);
            // preventing muons with default values due to tau energy correction
            // shifts below good tau pt selection
            if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                auto sf = evaluator->evaluate(
                    {year_id, std::abs(eta), pt, variation});
                return sf;
            } else
                return (double)1.;
        },
        {pt, eta});
    return df1;
}
} // namespace muon
namespace tau {
/**
 * @brief Function used to evaluate id scale factors from taus with
 * correctionlib
 *
 * @param df The input dataframe
 * @param pt tau pt
 * @param decayMode decay mode of the tau
 * @param genMatch column with genmatch values (from prompt e, prompt mu,
 * tau->e, tau->mu, had. tau)
 * @param wp working point of the ID cut
 * @param variation id for the variation of the scale factor "sf" for nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_dependence "pt", "dm" or "eta" based scale factors
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param sf_name name of the tau id scale factor
 * @return a new dataframe containing the new column
 */
auto id(auto &df, const std::string &pt, const std::string &decayMode,
        const std::string &genMatch, const std::string &wp,
        const std::string &variation, const std::string &sf_dependence,
        const std::string &id_output, const std::string &sf_file,
        const std::string &sf_name, const std::vector<int> &SelectedDMs) {

    Logger::get("tauIdSF")->debug("Setting up functions for tau id sf");
    Logger::get("tauIdSF")->debug("ID - Name {}", sf_name);
    auto evaluator = correction::CorrectionSet::from_file(sf_file)->at(sf_name);
    auto idSF_calculator = [evaluator, wp, variation, sf_dependence,
                            SelectedDMs](const float &pt, const int &decayMode,
                                         const UChar_t &genMatch) {
        if (std::find(SelectedDMs.begin(), SelectedDMs.end(), decayMode) !=
            SelectedDMs.end()) {
            Logger::get("tauIdSF")->debug("ID - pt {}, genMatch {}", pt,
                                          genMatch);
            auto sf =
                evaluator->evaluate({pt, decayMode, static_cast<int>(genMatch),
                                     wp, variation, sf_dependence});
            return sf;
        } else
            return (double)1.;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {pt, decayMode, genMatch});
    return df1;
}
} // namespace tau
} // namespace scalefactor