#ifndef GUARD_ELECTRONS_H
#define GUARD_ELECTRONS_H

namespace physicsobject {
namespace electron {

ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap);
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correctionManager,
             const std::string &corrected_pt, const std::string &pt,
             const std::string &eta, const std::string &sf_barrel,
             const std::string &sf_endcap, const std::string &sf_file,
             const std::string &jsonESname);
// deprecated version without CorrectionManager
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df, const std::string &corrected_pt,
             const std::string &pt, const std::string &eta,
             const std::string &sf_barrel, const std::string &sf_endcap,
             const std::string &sf_file, const std::string &jsonESname);
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df,
               correctionManager::CorrectionManager &correctionManager,
               const std::string &corrected_pt, const std::string &pt,
               const std::string &eta, const std::string &gain,
               const std::string &ES_sigma_up, const std::string &ES_sigma_down,
               const std::string &era, const std::string &variation,
               const std::string &ES_file);
// deprecated version without CorrectionManager
ROOT::RDF::RNode
PtCorrectionMC(ROOT::RDF::RNode df, const std::string &corrected_pt,
               const std::string &pt, const std::string &eta,
               const std::string &gain, const std::string &ES_sigma_up,
               const std::string &ES_sigma_down, const std::string &era,
               const std::string &variation, const std::string &ES_file);
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID);
ROOT::RDF::RNode CutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const int &IDvalue);
ROOT::RDF::RNode AntiCutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                             const std::string &nameID, const int &IDvalue);
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold);
ROOT::RDF::RNode CutIP(ROOT::RDF::RNode df, const std::string &eta,
                       const std::string &detasc, const std::string &dxy,
                       const std::string &dz, const std::string &maskname,
                       const float &abseta_eb_ee, const float &max_dxy_eb,
                       const float &max_dz_eb, const float &max_dxy_ee,
                       const float &max_dz_ee);

ROOT::RDF::RNode CutGap(ROOT::RDF::RNode df, const std::string &eta,
                        const std::string &detasc, const std::string &maskname,
                        const float &end_eb, const float &start_ee);

} // namespace electron
} // namespace physicsobject
#endif /* GUARD_ELECTRONS_H */