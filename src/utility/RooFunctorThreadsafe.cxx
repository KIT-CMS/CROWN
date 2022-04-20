#ifndef GUARDROOFUNCTORTHREADSAFE_H
#define GUARDROOFUNCTORTHREADSAFE_H

#include "RooFunctor.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include <chrono>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

/**
 * @brief Function used to load a
 * [`RooFunctor`](https://root.cern.ch/doc/master/classRooFunctor.html) from a
 * [`RooWorkspace`](https://root.cern.ch/doc/master/classRooWorkspace.html).
 * These can be used to store scale factors and other functions. An example,
 * how these workspaces are created can be found in [this
 * repository](https://github.com/KIT-CMS/LegacyCorrectionsWorkspace). This
 * version uses the threadsafe variant of loading `RooFunctors`, defined in the
 * `RooFunctorThreadsafe` class.
 *
 * @param workspace_name The path to the workspace file, from which the functor
 * should be loaded
 * @param functor_name The name of the function from the workspace to be loaded
 * @param arguments The arguments, that form the `ArgSet` of of the functor.
 * @returns A `std::shared_ptr<RooFunctorThreadsafe>`, which contains the
 * functor used for evaluation
 */

auto loadFunctor(const std::string &workspace_name,
                        const std::string &functor_name,
                        const std::string &arguments) {
    // first load the workspace
    std::unique_ptr<TFile> workspacefile{
        TFile::Open(workspace_name.c_str(), "read")};
    auto workspace = workspacefile->Get<RooWorkspace>("w");
    workspacefile->Close();
    auto func = workspace->function(functor_name.c_str());
    auto args = workspace->argSet(arguments.c_str());

    auto functor =
        std::make_shared<RooFunctorThreadsafe>(*func, args, workspace);
    return functor;
}

#endif /* GUARDROOFUNCTORTHREADSAFE_H */