#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

#include "physicsobjects.hxx"
#include "metfilter.hxx"

static std::vector<std::string> varSet = {};

int main(){

    ROOT::RDataFrame df("Events", "/ceph/htautau/nanoaod/testing/DYJetsToLL_SUmmer19_UL18.root");
    // for testing, we limit to 1000 events only
    // 1st stage: Good object selection
    auto df2 = df.Range(1000);
    // MET Filters
    auto df3 = metfilter::ApplyMetFilter(df2, "Flag_goodVertices");
    // Muons

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
    auto df5 = physicsobject::CutEta(df4, "Tau_pt", "tau_masketa", 2.3);
    auto df6 = physicsobject::tau::FilterDecayModes(df5, "tau_maskdecaymodes", {1,2,3,4,10});
    auto df7 = physicsobject::tau::FilterTauID(df6, "tau_maskIDvsJet", "Tau_idDeepTau2017v2p1VSjet", 4);
    auto df8 = physicsobject::tau::FilterTauID(df7, "tau_maskIDvsEle", "Tau_idDeepTau2017v2p1VSe", 4);
    auto df9 = physicsobject::tau::FilterTauID(df8, "tau_maskIDvsMu", "Tau_idDeepTau2017v2p1VSmu", 1);
    auto df11 = physicsobject::CutMax(df9, "Tau_dz", "tau_maskdz", 0.2);
    // auto df6 = physicsobject::tau::CutEta(df5, varSet, "tau_cutfoo", "tau_blub");
    // auto df7 = df6.Define("nTaus_filtered","int(good_taus.size())");
    // varSet.push_back("nTaus_filtered");
    // Jets
    varSet.push_back("tau_maskIDvsMu");
    auto df_final = df11;
    df_final.Report();
    df_final.Snapshot("ntuple", "test.root", varSet);
    // as a first testcase, we work on selecting good muons
    return 0;
}