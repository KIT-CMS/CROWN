#include "ROOT/RDataFrame.hxx"
#include "RooFunctor.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "basefunctions.hxx"
#include "utility/Logger.hxx"
namespace scalefactor {
std::shared_ptr<RooFunctor> loadFunctor(const std::string &workspace_name,
                                        const std::string &functor_name,
                                        const std::string &arguments) {
    // first load the workspace
    auto workspacefile = TFile::Open(workspace_name.c_str(), "read");
    auto workspace = (RooWorkspace *)workspacefile->Get("w");
    workspacefile->Close();
    // workspace->function(id_functor_name.c_str())->Print();
    // workspace->argSet(id_arguments.c_str()).Print();
    const std::shared_ptr<RooFunctor> functor = std::shared_ptr<RooFunctor>(
        workspace->function(functor_name.c_str())
            ->functor(workspace->argSet(arguments.c_str())));
    return functor;
}
namespace muon {
auto id(auto df, const std::string &pt, const std::string &eta,
        const std::string &id_output, const std::string &workspace_name,
        const std::string &id_functor_name, const std::string &id_arguments) {

    Logger::get("muonsf")->debug("Setting up functions for muon sf");
    Logger::get("muonsf")->debug("ID - Function {} // argset {}",
                                 id_functor_name, id_arguments);

    const std::shared_ptr<RooFunctor> id_function =
        loadFunctor(workspace_name, id_functor_name, id_arguments);
    auto df1 =
            basefunctions::evaluateWorkspaceFunction(df, id_output, id_function,
                                                     pt, eta);
    return df1;
}
auto iso(auto df, const std::string &pt, const std::string &eta,
         const std::string &iso, const std::string &iso_output,
         const std::string &workspace_name, const std::string &iso_functor_name,
         const std::string &iso_arguments) {

    Logger::get("muonsf")->debug("Setting up functions for muon sf");
    Logger::get("muonsf")->debug("Iso - Function {} // argset {}",
                                 iso_functor_name, iso_arguments);

    const std::shared_ptr<RooFunctor> iso_function =
        loadFunctor(workspace_name, iso_functor_name, iso_arguments);
    auto df1 =
            basefunctions::evaluateWorkspaceFunction(
                df, iso_output, iso_function, pt, eta, iso);
    return df1;
}
} // namespace muon
} // namespace scalefactor