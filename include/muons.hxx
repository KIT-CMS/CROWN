#ifndef GUARD_MUONS_H
#define GUARD_MUONS_H

namespace physicsobject {
namespace muon {
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID);
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold);
ROOT::RDF::RNode AntiCutIsolation(ROOT::RDF::RNode df,
                                  const std::string &maskname,
                                  const std::string &isolationName,
                                  const float &Threshold);
ROOT::RDF::RNode CutIsTrackerOrIsGlobal(ROOT::RDF::RNode df,
                                        const std::string &isTracker,
                                        const std::string &isGlobal,
                                        const std::string &maskname);
ROOT::RDF::RNode GenerateRndmRVec(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &objCollection, int seed);
ROOT::RDF::RNode
applyRoccoRData(ROOT::RDF::RNode df, const std::string &outputname,
                const std::string &filename, const int &position,
                const std::string &objCollection,
                const std::string &chargColumn, const std::string &ptColumn,
                const std::string &etaColumn, const std::string &phiColumn,
                int error_set, int error_member);
ROOT::RDF::RNode
applyRoccoRMC(ROOT::RDF::RNode df, const std::string &outputname,
              const std::string &filename, const int &position,
              const std::string &objCollection, const std::string &chargColumn,
              const std::string &ptColumn, const std::string &etaColumn,
              const std::string &phiColumn, const std::string &genPtColumn,
              const std::string &nTrackerLayersColumn,
              const std::string &rndmColumn, int error_set, int error_member);
} // namespace muon
} // namespace physicsobject
#endif /* GUARD_MUONS_H */