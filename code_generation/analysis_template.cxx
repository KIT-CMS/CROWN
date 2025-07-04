#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "RooTrace.h"
#include "TStopwatch.h"
#include <ROOT/RLogger.hxx>
#include "include/utility/Logger.hxx"
#include <TFile.h>
#include <TMap.h>
#include <filesystem>
#include <TObjString.h>
#include <TTree.h>
#include <TVector.h>
#include "onnxruntime_cxx_api.h"
#include <regex>
#include <string>
#include "include/utility/OnnxSessionManager.hxx"
#include "include/utility/CorrectionManager.hxx"
#include "include/genparticles.hxx"
#include "include/htxs.hxx"
#include "include/jets.hxx"
#include "include/event.hxx"
#include "include/lorentzvectors.hxx"
#include "include/met.hxx"
#include "include/ml.hxx"
#include "include/pairselection.hxx"
#include "include/tripleselection.hxx"
#include "include/embedding.hxx"
#include "include/physicsobjects.hxx"
#include "include/muons.hxx"
#include "include/electrons.hxx"
#include "include/taus.hxx"
#include "include/quantities.hxx"
#include "include/reweighting.hxx"
#include "include/topreco.hxx"
#include "include/triggers.hxx"

// {INCLUDE_ANALYSISADDONS}

// {INCLUDES}

int main(int argc, char *argv[]) {
    // {DEBUGLEVEL}
    // ROOT logging
    if (debug) {
        auto verbosity = ROOT::Experimental::RLogScopedVerbosity(
            ROOT::Detail::RDF::RDFLogChannel(),
            ROOT::Experimental::ELogLevel::kInfo);
        RooTrace::verbose(kTRUE);
        Logger::setLevel(Logger::LogLevel::DEBUG);
    } else {
        RooTrace::verbose(kFALSE);
        Logger::setLevel(Logger::LogLevel::INFO);
        gErrorIgnoreLevel = 6001; // ignore all ROOT errors
    }
    if (argc < 3) {
        Logger::get("main")->critical(
            "Require at least two arguments: a single output file and N input "
            "files \n"
            "Example:\n"
            "./analysis output.root /path/to/inputfiles/*.root");
        return 1;
    }
    // check if CROWN is run from the correct directory, if the folder "data" does not exist, exit
    if (!std::filesystem::exists("data")) {
        Logger::get("main")->critical(
            "CROWN is not run from the correct directory, "
            "data folder does not exist. Did you run CROWN from the correct directory?");
        return 1;
    }
    std::vector<std::string> input_files;
    int nevents = 0;
    Logger::get("main")->info("Checking input files");
    std::string basetree = "Events";
    for (int i = 2; i < argc; i++) {
        input_files.push_back(std::string(argv[i]));
        // Check if the input file exists and is readable, also get the number
        // of events, using a timeout of 30 seconds
        TFile *f1 = TFile::Open(argv[i], "TIMEOUT=30");
        if (!f1 || f1->IsZombie()) {
            Logger::get("main")->critical("File {} does not exist or is not "
                                          "readable",
                                          argv[i]);
            return 1;
        }
        // Get a list of all keys in the file
        TList *list = f1->GetListOfKeys();
        // Check if the Events tree exists
        if (list->FindObject("Events")) {
            TTree *t1 = (TTree *)f1->Get("Events");
            nevents += t1->GetEntries();
            Logger::get("main")->info("NanoAOD input_file {}: {} - {} Events",
                                      i - 1, argv[i], t1->GetEntries());
        } else if (list->FindObject("ntuple")) {
            TTree *t1 = (TTree *)f1->Get("ntuple");
            nevents += t1->GetEntries();
            basetree = "ntuple";
            Logger::get("main")->info("CROWN input_file {}: {} - {} Events",
                                      i - 1, argv[i], t1->GetEntries());
        } else {
            Logger::get("main")->critical("File {} does not contain a tree "
                                          "named 'Events' or 'ntuple'",
                                          argv[i]);
            return 1;
        }
    }
    const auto output_path = argv[1];
    Logger::get("main")->info("Output directory: {}", output_path);
    TStopwatch timer;
    timer.Start();
    int quantile = 10000;

    // file logging
    Logger::enableFileLogging("logs/main.txt");

    // start an onnx session manager
    OnnxSessionManager onnxSessionManager;
    // start a correction manager
    correctionManager::CorrectionManager correctionManager;

    // {MULTITHREADING}

    // initialize df
    ROOT::RDataFrame df0(basetree, input_files);
    ROOT::RDF::Experimental::AddProgressBar(df0); // add progress bar
    Logger::get("main")->info("Starting Setup of Dataframe with {} events",
                              nevents);
    std::vector<ROOT::RDF::RResultPtr<ROOT::RDF::RCutFlowReport>> cutReports;
    // setup output files
    // {OUTPUT_PATHS}
    if (nevents == 0){
        // {ZERO_EVENTS_FALLBACK}
    } else {
        // {CODE_GENERATION}
        ROOT::RDF::RSnapshotOptions dfconfig;
        dfconfig.fLazy = true;
        // {RUN_COMMANDS}
    }
    // Add meta-data
    // clang-format off
    const std::map<std::string, std::vector<std::string>> output_quantities = {OUTPUT_QUANTITIES};
    const std::map<std::string, std::vector<std::string>> variations = {SYSTEMATIC_VARIATIONS};
    std::map<std::string, std::map<std::string, std::vector<std::string>>> shift_quantities_map = {SHIFT_QUANTITIES_MAP};
    std::map<std::string, std::map<std::string, std::vector<std::string>>> quantities_shift_map = {QUANTITIES_SHIFT_MAP};
    // clang-format on
    const std::string analysis = {ANALYSISTAG};
    const std::string config = {CONFIGTAG};
    const std::string era = {ERATAG};
    const std::string sample = {SAMPLETAG};
    const std::string commit_hash = {COMMITHASH};
    bool setup_clean = {CROWN_IS_CLEAN};
    const std::string analysis_commit_hash = {ANALYSIS_COMMITHASH};
    bool analysis_setup_clean = {ANALYSIS_IS_CLEAN};
    int scope_counter = 0;
    for (auto const &output : output_quantities) {
        // output.first is the output file name
        // output.second is the list of quantities
        TFile outputfile(output.first.c_str(), "UPDATE");
        TTree quantities_meta = TTree("quantities", "quantities");
        for (auto const &quantity : output.second) {
            quantities_meta.Branch(quantity.c_str(), &setup_clean);
        }
        quantities_meta.Write();
        TTree variations_meta = TTree("variations", "variations");
        for (auto const &variation : variations.at(output.first)) {
            variations_meta.Branch(variation.c_str(), &setup_clean);
        }
        variations_meta.Write();
        TTree conditions_meta = TTree("conditions", "conditions");
        conditions_meta.Branch(analysis.c_str(), &setup_clean);
        conditions_meta.Branch(config.c_str(), &setup_clean);
        conditions_meta.Branch(era.c_str(), &setup_clean);
        conditions_meta.Branch(sample.c_str(), &setup_clean);
        conditions_meta.Write();
        TTree commit_meta = TTree("commit", "commit");
        commit_meta.Branch(commit_hash.c_str(), &setup_clean);
        commit_meta.Fill();
        commit_meta.Write();
        TTree analysis_commit_meta =
            TTree("analysis_commit", "analysis_commit");
        analysis_commit_meta.Branch(analysis_commit_hash.c_str(),
                                    &analysis_setup_clean);
        analysis_commit_meta.Fill();
        analysis_commit_meta.Write();
        if (nevents != 0){
            TH1D cutflow;
            cutflow.SetName("cutflow");
            cutflow.SetTitle("cutflow");
            // iterate through the cutflow vector and fill the histogram with the
            // .GetPass() values
            if (scope_counter >= cutReports.size()) {
                Logger::get("main")->critical(
                    "Cutflow vector is too small, this should not happen");
                return 1;
            }
            for (auto cut = cutReports[scope_counter].begin();
                cut != cutReports[scope_counter].end(); cut++) {
                cutflow.SetBinContent(
                    std::distance(cutReports[scope_counter].begin(), cut) + 1,
                    cut->GetPass());
                cutflow.GetXaxis()->SetBinLabel(
                    std::distance(cutReports[scope_counter].begin(), cut) + 1,
                    cut->GetName().c_str());
            }
            // store it in the output file
            cutflow.Write();
        }
        outputfile.Close();
        TFile *fout = TFile::Open(output.first.c_str(), "UPDATE");
        Logger::get("main")->info("Writing quantities map to {}", output.first);
        fout->WriteObject(&shift_quantities_map.at(output.first),
                          "shift_quantities_map");
        fout->WriteObject(&quantities_shift_map.at(output.first),
                          "quantities_shift_map");
        fout->Close();
        scope_counter++;
    }

    Logger::get("main")->info("Finished Evaluation");
    Logger::get("main")->info(
        "Overall runtime (real time: {0:.2f}, CPU time: {1:.2f})",
        timer.RealTime(), timer.CpuTime());
}
