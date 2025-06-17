#ifndef GUARD_REWEIGHTING_H
#define GUARD_REWEIGHTING_H

namespace event {
namespace reweighting {

ROOT::RDF::RNode
Pileup(ROOT::RDF::RNode df,
          correctionManager::CorrectionManager &correction_manager,
          const std::string &outputname, const std::string &true_pileup_number,
          const std::string &corr_file, const std::string &corr_name,
          const std::string &variation);
ROOT::RDF::RNode PartonShower(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &ps_weights,
                            const float isr, const float fsr);
ROOT::RDF::RNode LHEscale(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &lhe_scale_weights,
                            const float mu_r, const float mu_f);
ROOT::RDF::RNode LHEpdf(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &lhe_pdf_weights,
                        const std::string &variation);
ROOT::RDF::RNode LHEalphaS(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &lhe_pdf_weights,
                        const std::string &variation);
ROOT::RDF::RNode TopPt(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &gen_pdg_id,
                        const std::string &gen_status,
                        const std::string &gen_pt);
ROOT::RDF::RNode ZPtMass(ROOT::RDF::RNode df,
                            const std::string &outputname,
                            const std::string &gen_boson,
                            const std::string &workspace_file,
                            const std::string &functor_name,
                            const std::string &argset);
} // end namespace reweighting
} // end namespace event
#endif /* GUARD_REWEIGHTING_H */
