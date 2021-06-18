#include "ROOT/RDataFrame.hxx"
#include "RooFunctor.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "basefunctions.hxx"
#include "utility/Logger.hxx"
/// namespace used for scale factor related functions
namespace scalefactor {
/**
 * @brief Function used to evaluate id scale factors from muons
 * Function used to load a
 * [`RooFunctor`](https://root.cern.ch/doc/master/classRooFunctor.html) from a
 * [`RooWorkspace`](https://root.cern.ch/doc/master/classRooWorkspace.html).
 * These can be used to store scale factors and other functions. An example,
 * how these workspaces are created can be found in [this
 * repository](https://github.com/KIT-CMS/LegacyCorrectionsWorkspace).
 *
 * @param workspace_name The path to the workspace file, from which the functor
 * should be loaded
 * @param functor_name The name of the function from the workspace to be loaded
 * @param arguments The arguments, that form the `ArgSet` of of the functor.
 * @returns A `std::shared_ptr<RooFunctor>`, which contains the functor used for
 * evaluation
 */
std::shared_ptr<RooFunctor> loadFunctor(const std::string &workspace_name,
                                        const std::string &functor_name,
                                        const std::string &arguments) {
    // first load the workspace
    auto workspacefile = TFile::Open(workspace_name.c_str(), "read");
    auto workspace = (RooWorkspace *)workspacefile->Get("w");
    workspacefile->Close();
    const std::shared_ptr<RooFunctor> functor = std::shared_ptr<RooFunctor>(
        workspace->function(functor_name.c_str())
            ->functor(workspace->argSet(arguments.c_str())));
    return functor;
}
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

    const std::shared_ptr<RooFunctor> id_function =
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

    const std::shared_ptr<RooFunctor> iso_function =
        loadFunctor(workspace_name, iso_functor_name, iso_arguments);
    auto df1 = basefunctions::evaluateWorkspaceFunction(
        df, iso_output, iso_function, pt, eta, iso);
    return df1;
}
} // namespace muon
} // namespace scalefactor