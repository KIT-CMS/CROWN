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
#include <TFile.h>
#include <TTree.h>
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
    // Add meta-data
    const std::string outputfilename = {METADATAFILENAME};
    const std::vector<std::string> output_quanties = {OUTPUT_QUANTITIES};
    const std::vector<std::string> variations = {SYSTEMATIC_VARIATIONS};
    const std::string analysis = {ANALYSISTAG};
    const std::string era = {ERATAG};
    const std::string sample = {SAMPLETAG};
    const std::string commit_hash = {COMMITHASH};
    bool setup_clean = {CLEANSETUP};
    TFile outputfile(outputfilename.c_str(), "UPDATE");
    TTree quantities_meta = TTree("quantities", "quantities");
    for (auto quantity : output_quanties) {
        quantities_meta.Branch(quantity.c_str(), &setup_clean);
    }
    quantities_meta.Write();
    TTree variations_meta = TTree("variations", "variations");
    for (auto variation : variations) {
        variations_meta.Branch(variation.c_str(), &setup_clean);
    }
    variations_meta.Write();
    TTree conditions_meta = TTree("conditions", "conditions");
    conditions_meta.Branch(analysis.c_str(), &setup_clean);
    conditions_meta.Branch(era.c_str(), &setup_clean);
    conditions_meta.Branch(sample.c_str(), &setup_clean);
    conditions_meta.Write();
    TTree commit_meta = TTree("commit", "commit");
    commit_meta.Branch(commit_hash.c_str(), &setup_clean);
    commit_meta.Fill();
    commit_meta.Write();
    outputfile.Close();

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
