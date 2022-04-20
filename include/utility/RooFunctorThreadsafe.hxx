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

class RooFunctorThreadsafe {

    // Thread-safe alternative to RooFunctor.
    // This class starts by cloning the full tree of the RooFit function and
    // creating a function for it.  If the evaluation operator is called
    // concurrently, more function clones and functors will be spawned if
    // necessary.

  public:
    RooFunctorThreadsafe(RooAbsReal const &function, RooArgSet const &args,
                         RooWorkspace *ws)
        : workspace_(ws) {
        executors_.emplace_back(function, args);
        // to avoid reallocation for executor vector which is not threadsafe
        executors_.reserve(maxNExecutors);
    }

    double operator()(const double *input) { return eval(input); }

    double eval(const double *input) { return getIdleExecutor().eval(input); }

  private:
    class Executor {

      public:
        Executor(RooAbsReal const &function, RooArgSet const &args) {
            init(function, args);
        }
        Executor(Executor const &other) {
            std::lock_guard<std::mutex> lock(other.mutex_);
            init(*other.function_, *other.args_);
        }
        Executor &operator=(Executor const &other) {
            std::lock_guard<std::mutex> lock(other.mutex_);
            init(*other.function_, *other.args_);
            return *this;
        }
        Executor(Executor &&other) = delete;
        Executor &operator=(Executor &&other) = delete;

        double eval(const double *input) {
            double result = functor_->eval(input);
            mutex_.unlock();
            return result;
        }

        bool try_lock() const { return mutex_.try_lock(); }

      private:
        void init(RooAbsReal const &function, RooArgSet const &args) {
            function_.reset(static_cast<RooAbsReal *>(function.cloneTree()));
            RooArgSet funcServers;
            function_->treeNodeServerList(&funcServers);
            args_.reset(
                static_cast<RooArgSet *>(funcServers.selectCommon(args)));
            functor_.reset(function_->functor(*args_));
        }

        std::unique_ptr<RooAbsReal> function_ = nullptr;
        std::unique_ptr<RooArgSet> args_ = nullptr;
        std::unique_ptr<RooFunctor> functor_ = nullptr;
        mutable std::mutex mutex_;
    };

    Executor &addNewExecutor() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (executors_.size() >= maxNExecutors) {
            throw std::runtime_error("Maximum number of executors reached.");
        }
        executors_.emplace_back(executors_.back());
        return executors_.back();
    }

    Executor &getIdleExecutor() {
        for (auto &executor : executors_) {
            if (executor.try_lock()) {
                return executor;
            }
        }
        auto &out = addNewExecutor();
        out.try_lock();
        return out;
    }

    mutable std::mutex mutex_;
    std::vector<Executor> executors_;
    // some of the objects indirectly used by the executors, such as
    // RooDataHists, might be owned by the RooWorkspace where they came from, so
    // we tie the workspace's lifetime to the one of the RooFunctorThreadsafe
    // object
    std::unique_ptr<RooWorkspace> workspace_;

    constexpr static int maxNExecutors = 100;
};

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

inline auto loadFunctor(const std::string &workspace_name,
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