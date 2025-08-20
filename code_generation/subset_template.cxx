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
#include "include/fatjets.hxx"
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

ROOT::RDF::RNode {subsetname} (ROOT::RDF::RNode df0, OnnxSessionManager &onnxSessionManager, correctionManager::CorrectionManager &correctionManager) {

    //    { commands }
}
