#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "RooTrace.h"
#include "TStopwatch.h"
#include "src/htxs.hxx"
#include "src/jets.hxx"
#include "src/lorentzvectors.hxx"
#include "src/metfilter.hxx"
#include "src/pairselection.hxx"
#include "src/physicsobjects.hxx"
#include "src/quantities.hxx"
#include "src/scalefactors.hxx"
#include "src/triggers.hxx"
#include "src/utility/Logger.hxx"
#include <ROOT/RLogger.hxx>
#include <string>

static std::vector<std::string> varSet = {"run", "luminosityBlock", "event"};

int main(int argc, char *argv[]) {
    if (argc != 3) {
        Logger::get("main")->critical(
            "Require exactly two additional input arguments (the input and "
            "output paths to the ROOT files) but got {}",
            argc - 1);
        return 1;
    }

    const auto input_path = argv[1];
    Logger::get("main")->info("Input file: {}", input_path);
    const auto output_path = argv[2];
    Logger::get("main")->info("Output directory: {}", output_path);

    TStopwatch timer;
    timer.Start();
    ROOT::EnableImplicitMT(1); // for multithreading
    // ROOT logging
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(
        ROOT::Detail::RDF::RDFLogChannel(),
        ROOT::Experimental::ELogLevel::kInfo);
    RooTrace::verbose(kTRUE);

    // file logging
    ROOT::RDataFrame df0("Events", input_path);
    // for testing, we limit to 1000 events only
    // 1st stage: Good object selection
    Logger::enableFileLogging("logs/main.txt");
    Logger::setLevel(Logger::LogLevel::INFO);
    Logger::get("main")->info("Starting Setup of Dataframe");

    // auto df_final = df0;

    // {CODE_GENERATION}

    // Logger::get("main")->debug(df_final.Describe()); // <-- starting from
    // ROOT 6.25

    Logger::get("main")->info("Finished Setup");
    Logger::get("main")->info("Runtime for setup (real time: {}, CPU time: {})",
                              timer.RealTime(), timer.CpuTime());
    timer.Continue();

    Logger::get("main")->info("Starting Evaluation");
    ROOT::RDF::RSnapshotOptions dfconfig;
    dfconfig.fLazy = true;
    // {RUN_COMMANDS}
    Logger::get("main")->info("Finished Evaluation");

    const auto nruns = {NRUNS};
    if (nruns != 1) {
        Logger::get("main")->critical(
            "Analysis runs more than one event loop!");
    }

    // as a first testcase, we work on selecting good muons
    Logger::get("main")->info("Overall runtime (real time: {}, CPU time: {})",
                              timer.RealTime(), timer.CpuTime());
}
