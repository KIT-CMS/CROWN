#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TStopwatch.h"
#include "lorentzvectors.hxx"
#include "metfilter.hxx"
#include "pairselection.hxx"
#include "physicsobjects.hxx"
#include "utility/Logger.hxx"
#include <ROOT/RLogger.hxx>

static std::vector<std::string> varSet = {"run", "luminosityBlock", "event"};

int main() {
  TStopwatch timer;
  timer.Start();

  // ROOT logging
  auto verbosity = ROOT::Experimental::RLogScopedVerbosity(
      ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);

  // file logging
  ROOT::RDataFrame df0(
      "Events",
      "/work/sbrommer/ntuple_prototype/files/DYJetsToLL_SUmmer19_UL18.root");
  // for testing, we limit to 1000 events only
  // 1st stage: Good object selection
  Logger::enableFileLogging("logs/main.txt");
  Logger::setLevel(Logger::LogLevel::INFO);
  Logger::get("main")->info("Starting Setup of Dataframe");

  // auto df_final = df0;

  // {CODE_GENERATION}

  auto cutReport = df_final.Report();

  // Logger::get("main")->debug(df_final.Describe()); // <-- starting from
  // ROOT 6.25

  Logger::get("main")->info("Finished Setup");
  Logger::get("main")->info("Runtime for setup (real time: {}, CPU time: {})",
                            timer.RealTime(), timer.CpuTime());
  timer.Continue();

  Logger::get("main")->info("Starting Evaluation");
  df_final.Snapshot("ntuple", "test.root", varSet);
  cutReport->Print();
  Logger::get("main")->info("Finished Evaluation");

  const auto nruns = df_final.GetNRuns();
  if (nruns != 1) {
    Logger::get("main")->critical("Analysis runs more than one event loop!");
  }

  // as a first testcase, we work on selecting good muons
  Logger::get("main")->info("Overall runtime (real time: {}, CPU time: {})",
                            timer.RealTime(), timer.CpuTime());
}
