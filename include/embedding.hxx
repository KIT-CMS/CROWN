#ifndef GUARD_EMBEDDING_H
#define GUARD_EMBEDDING_H

namespace embedding {

namespace scalefactor {

ROOT::RDF::RNode
SelectionTrigger(ROOT::RDF::RNode df,
                  correctionManager::CorrectionManager &correction_manager,
                  const std::string &outputname,
                  const std::string &pt_1, const std::string &eta_1,
                  const std::string &pt_2, const std::string &eta_2,
                  const std::string &sf_file,
                  const std::string &sf_name);
ROOT::RDF::RNode
SelectionId(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correction_manager,
             const std::string &outputname, 
             const std::string &pt, const std::string &eta,
             const std::string &sf_file,
             const std::string &sf_name);
} // end namespace scalefactor

namespace muon {

ROOT::RDF::RNode
Scalefactor(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correction_manager,
            const std::string &output,
            const std::string &pt, const std::string &eta,
            const std::string &sf_file, const std::string &sf_name,
            const std::string correction_type,
            const float &extrapolation_factor = 1.0);
} // end namespace muon

namespace electron {

ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap);
ROOT::RDF::RNode
PtCorrection(ROOT::RDF::RNode df,
             correctionManager::CorrectionManager &correction_manager,
             const std::string &outputname, const std::string &pt,
             const std::string &eta, const std::string &es_file,
             const std::string &correction_name,
             const std::string &variation_barrel,
             const std::string &variation_endcap);
ROOT::RDF::RNode
Scalefactor(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correction_manager,
            const std::string &output,
            const std::string &pt, const std::string &eta,
            const std::string &sf_file, const std::string &sf_name,
            const std::string correction_type,
            const float &extrapolation_factor = 1.0);
} // end namespace electron

namespace tau {

ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &pt, const std::string &decay_mode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11);

namespace scalefactor {

ROOT::RDF::RNode
Id_vsJet_lt(ROOT::RDF::RNode df,
            correctionManager::CorrectionManager &correction_manager,
            const std::string &outputname,
            const std::string &pt, const std::string &decay_mode,
            const std::string &gen_match, 
            const std::string &sf_file,
            const std::string &sf_name,
            const std::vector<int> &selected_dms,
            const std::string &wp, const std::string &vsele_wp,
            const std::string &sf_dependence,
            const std::string &sf_vsjet_tau20to25,
            const std::string &sf_vsjet_tau25to30,
            const std::string &sf_vsjet_tau30to35,
            const std::string &sf_vsjet_tau35to40,
            const std::string &sf_vsjet_tau40toInf);
} // end namespace scalefactor
} // end namespace tau
} // end namespace embedding
#endif /* GUARD_EMBEDDING_H */