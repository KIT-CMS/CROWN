#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TStopwatch.h"
#include "lorentzvectors.hxx"
#include "metfilter.hxx"
#include "pairselection.hxx"
#include "physicsobjects.hxx"
#include "quantities.hxx"
#include "utility/Logger.hxx"
#include <ROOT/RLogger.hxx>

static std::vector<std::string> varSet = {"run", "luminosityBlock", "event"};

int main() {
    TStopwatch timer;
    timer.Start();

    // ROOT logging
    auto verbosity = ROOT::Experimental::RLogScopedVerbosity(
        ROOT::Detail::RDF::RDFLogChannel(),
        ROOT::Experimental::ELogLevel::kInfo);

    // file logging
    ROOT::RDataFrame df(
        "Events",
        "/work/sbrommer/ntuple_prototype/files/DYJetsToLL_SUmmer19_UL18.root");
    // for testing, we limit to 1000 events only
    // 1st stage: Good object selection
    Logger::enableFileLogging("logs/main.txt");
    Logger::setLevel(Logger::LogLevel::INFO);
    Logger::get("main")->info("Starting Setup of Dataframe");
    auto df2 = df.Range(0, 0); // Run on entire file
    // auto df2 = df.Range(100000); // Run on part of file
    // MET Filters
    auto df3 = metfilter::ApplyMetFilter(df2, "Flag_goodVertices",
                                         "METFilter_goodVertices");

    // Electron

    // Taus
    // Steps for good taus in the mt channel:
    // tau_pt > 30 GeV
    // tau_eta < 2.3
    // Tau_decayMode != 5,6,7,8,9
    // Tau_dz < 0.2
    // Tau_idDeepTau2017v2p1VSjet > 4
    // Tau_idDeepTau2017v2p1VSe > 4
    // Tau_idDeepTau2017v2p1VSmu > 1
    auto df4 = physicsobject::CutPt(df3, "Tau_pt", "tau_maskpt", 30);
    auto df5 = physicsobject::CutEta(df4, "Tau_eta", "tau_masketa", 2.3);
    auto df6 = physicsobject::tau::FilterDecayModes(df5, "tau_maskdecaymodes",
                                                    {0, 1, 2, 3, 4, 10, 11});
    auto df7 = physicsobject::tau::FilterTauID(df6, "tau_maskIDvsJet",
                                               "Tau_idDeepTau2017v2p1VSjet", 4);
    auto df8 = physicsobject::tau::FilterTauID(df7, "tau_maskIDvsEle",
                                               "Tau_idDeepTau2017v2p1VSe", 4);
    auto df9 = physicsobject::tau::FilterTauID(df8, "tau_maskIDvsMu",
                                               "Tau_idDeepTau2017v2p1VSmu", 1);
    auto df10 = physicsobject::CutDz(df9, "Tau_dz", "tau_maskdz", 0.2);

    auto df12 = physicsobject::CombineMasks(
        df10, "good_taus_mask",
        {"tau_maskpt", "tau_masketa", "tau_maskdecaymodes", "tau_maskIDvsJet",
         "tau_maskIDvsEle", "tau_maskIDvsMu", "tau_maskdz"});
    // muons
    auto df20 = physicsobject::CutPt(df12, "Muon_pt", "muon_maskpt", 23);
    auto df21 = physicsobject::CutEta(df20, "Muon_eta", "muon_masketa", 2.5);
    auto df22 =
        physicsobject::muon::FilterID(df21, "muon_maskid", "Muon_mediumId");
    auto df23 = physicsobject::muon::FilterIsolation(
        df22, "muon_maskiso", "Muon_pfRelIso04_all", 0.15);
    auto df24 = physicsobject::CombineMasks(
        df23, "good_muons_mask",
        {"muon_maskpt", "muon_masketa", "muon_maskid", "muon_maskiso"});

    // Build the pair !
    auto df24_1 = physicsobject::FilterObjects(df24, "nTau", 1, "NumberOfTaus");
    auto df24_2 =
        physicsobject::FilterObjects(df24_1, "nMuon", 1, "NumberOfMuons");
    auto mt_df = pairselection::mutau::PairSelection(
        df24_2, {"Tau_pt", "Tau_rawDeepTau2017v2p1VSjet", "Muon_pt",
                          "Muon_pfRelIso04_all", "good_taus_mask", "good_muons_mask"}, "mtpair");
    auto mt_df_1 =
        pairselection::filterGoodPairs(mt_df, "mtpair", "GoodMuTauPairs");
    auto mt_df_2 = lorentzvectors::mutau::build(
        mt_df_1, "mtpair", {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"},
        {"Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass"}, "p4_1", "p4_2");
    auto mt_df_3 = quantities::pt(mt_df_2, varSet, "pt_1", "p4_1");
    auto mt_df_4 = quantities::pt(mt_df_3, varSet, "pt_2", "p4_2");
    auto mt_df_5 = quantities::eta(mt_df_4, varSet, "eta_1", "p4_1");
    auto mt_df_6 = quantities::eta(mt_df_5, varSet, "eta_2", "p4_2");
    auto mt_df_7 = quantities::phi(mt_df_6, varSet, "phi_1", "p4_1");
    auto mt_df_8 = quantities::phi(mt_df_7, varSet, "phi_2", "p4_2");
    auto mt_df_9 = quantities::m_vis(mt_df_8, varSet, "m_vis", {"p4_1", "p4_2"});

    auto df_final = mt_df_9;
    auto cutReport = df_final.Report();

    // Logger::get("main")->debug(df_final.Describe()); // <-- starting from
    // ROOT 6.25

    varSet.push_back("good_muons_mask");
    varSet.push_back("good_taus_mask");
    varSet.push_back("mtpair");
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
        Logger::get("main")->critical("Analysis runs more than one event loop!",
                                      timer.RealTime(), timer.CpuTime());
    }

    // as a first testcase, we work on selecting good muons
    Logger::get("main")->info("Overall runtime (real time: {}, CPU time: {})",
                              timer.RealTime(), timer.CpuTime());
}
