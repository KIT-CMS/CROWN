#ifndef GUARD_MET_H
#define GUARD_MET_H

#include "../include/utility/utility.hxx"
#include "../include/event.hxx"

namespace met {

ROOT::RDF::RNode genBosonPt(ROOT::RDF::RNode df, const std::string &outputname,
                            const std::string &inputvector);
ROOT::RDF::RNode genBosonEta(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &inputvector);
ROOT::RDF::RNode genBosonPhi(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &inputvector);
ROOT::RDF::RNode genBosonRapidity(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &inputvector);
ROOT::RDF::RNode DefineRecoilsDilep(ROOT::RDF::RNode df,
                                    const std::string &paralell_outputname,
                                    const std::string &perpendicular_outputname,
                                    const std::string &lepton1,
                                    const std::string &lepton2,
                                    const std::string &met_p4);
ROOT::RDF::RNode
DefineRecoilsSinglelep(ROOT::RDF::RNode df,
                       const std::string &paralell_outputname,
                       const std::string &perpendicular_outputname,
                       const std::string &lepton1, const std::string &met_p4);
ROOT::RDF::RNode
propagateLeptonsToMet(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &met,
                      const std::string &p4_1_uncorrected,
                      const std::string &p4_2_uncorrected,
                      const std::string &p4_1, const std::string &p4_2,
                      bool apply_propagation);
ROOT::RDF::RNode propagateLeptonsToMet(
    ROOT::RDF::RNode df, const std::string &outputname, const std::string &met,
    const std::string &p4_1_uncorrected, const std::string &p4_2_uncorrected,
    const std::string &p4_3_uncorrected, const std::string &p4_1,
    const std::string &p4_2, const std::string &p4_3,
    bool apply_propagation);
ROOT::RDF::RNode propagateLeptonsToMet(ROOT::RDF::RNode df,
                                       const std::string &outputname,
                                       const std::string &met,
                                       const std::string &p4_1_uncorrected,
                                       const std::string &p4_1,
                                       bool apply_propagation);
ROOT::RDF::RNode propagateJetsToMet(
    ROOT::RDF::RNode df, const std::string &outputname, const std::string &met,
    const std::string &jet_pt_corrected, const std::string &jet_eta_corrected,
    const std::string &jet_phi_corrected, const std::string &jet_mass_corrected,
    const std::string &jet_pt, const std::string &jet_eta,
    const std::string &jet_phi, const std::string &jet_mass,
    bool apply_propagation, float min_jet_pt);
ROOT::RDF::RNode applyRecoilCorrections( //Run 2
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &met, const std::string &genmet,
    const std::string &jet_pt, 
    const std::string &recoilfile, const std::string &systematicsfile,
    bool applyRecoilCorrections, bool resolution, bool response, bool shiftUp,
    bool shiftDown, bool isWjets);
ROOT::RDF::RNode applyRecoilCorrectionsRun3(
    ROOT::RDF::RNode df, correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname,
    const std::string &p4_met, const std::string &p4_gen_boson,
    const std::string &p4_vis_gen_boson, const std::string &n_jets,
    const std::string &corr_file, const std::string &corr_name,
    const std::string &method, const std::string &order,
    const std::string &variation, bool apply_correction);
ROOT::RDF::RNode applyMetXYCorrections(ROOT::RDF::RNode df,
                                       const std::string &input_p4,
                                       const std::string &npv,
                                       const std::string &run,
                                       const std::string &output_p4,
                                       const std::string &corr_file, bool isMC);
} // end namespace met

namespace lorentzvector {

/**
 * @brief This function propagates object corrections to the MET based on
 * Lorentz vectors. The objects can be e.g. leptons like muons or taus. If the
 * energy of an object is corrected/changed (e.g. via some scale factor) or due
 * to a shift, this change in energy has to be propagated to the MET vector, and
 * the MET vector has to be adapted accordingly. The MET is recalculated via
 *
 * \f[
 *  E_{T,miss,x}^{\text{corrected}} = E_{T,miss,x} + p_{x,\text{object}}
 *        - p_{x,\text{object}}^{\text{corrected}} \\
 *  E_{T,miss,y}^{\text{corrected}} = E_{T,miss,y} + p_{y,\text{object}}
 *        - p_{y,\text{object}}^{\text{corrected}} \\
 *  E_{T,miss}^{\text{corrected}} = \sqrt{E_{T,miss,x}^{\text{corrected}} * E_{T,miss,x}^{\text{corrected}}
 *        + E_{T,miss,y}^{\text{corrected}} * E_{T,miss,y}^{\text{corrected}}}
 * \f]
 *
 * @param df input dataframe
 * @param outputname name of the new column containing the corrected MET Lorentz
 * vector
 * @param p4_met initial/uncorrected MET Lorentz vector
 * @param args the possible arguments need to have a specific structure:
 * - first half: uncorrected lorentz vectors of the objects (e.g. leptons)
 * - second half: corrected lorentz vectors (same number as uncorrected) of
 *   the objects (e.g. leptons)
 * - last argument: boolean indicating whether the propagation should be applied
 *   or just the original MET vector should be returned
 *
 * @return a dataframe with the new column containing the corrected MET Lorentz
 * vector
 */
template <typename... Args>
ROOT::RDF::RNode PropagateToMET(ROOT::RDF::RNode df,
                                       const std::string &outputname,
                                       const std::string &p4_met,
                                       Args... args) {
    using namespace ROOT::VecOps;
    auto argTuple = std::make_tuple(args...);
    auto apply_propagation = utility::extractLastArgument(argTuple);
    std::vector<std::string> LV_list = utility::popLastArgument(argTuple);
    LV_list.push_back(p4_met);

    const auto n_LVs = sizeof...(Args);
    const auto n_Objects = (n_LVs - 1) / 2;
    static_assert((n_Objects > 0) && (n_Objects % 2) == 0,
            "Provide an even, non-zero number of object columns: first half uncorrected, second half corrected.");

    if (apply_propagation) {
        // correct for the objects, store the MET in an intermediate column
        Logger::get("lorentzvector::PropagateToMet")
            ->debug("Setting up correction for {} objects", n_Objects);
        return df.Define(
            outputname,
            utility::PassAsVec<n_LVs, ROOT::Math::PtEtaPhiMVector>(
                [n_Objects](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &LVs) {
                    const ROOT::Math::PtEtaPhiMVector &met =
                        LVs.at(n_Objects * 2);
                    ROOT::Math::PtEtaPhiMVector corrected_met = met;

                    // We propagate the object corrections to the MET by scaling the x
                    // and y component of the MET according to the correction of the objects
                    for (auto i = 0; i < n_Objects; ++i) {
                        float corr_x = LVs.at(i).Px() - LVs.at(n_Objects + i).Px();
                        float corr_y = LVs.at(i).Py() - LVs.at(n_Objects + i).Py();
                        float MetX = corrected_met.Px() + corr_x;
                        float MetY = corrected_met.Py() + corr_y;
                        Logger::get("lorentzvector::PropagateToMET")->debug("corr_x {}", corr_x);
                        Logger::get("lorentzvector::PropagateToMET")->debug("corr_y {}", corr_y);
                        Logger::get("lorentzvector::PropagateToMET")->debug("MetX {}", MetX);
                        Logger::get("lorentzvector::PropagateToMET")->debug("MetY {}", MetY);

                        corrected_met.SetPxPyPzE(MetX, MetY, 0,
                                        std::sqrt(MetX * MetX + MetY * MetY));
                        Logger::get("lorentzvector::PropagateToMET")
                            ->debug("corrected_object pt - {}", LVs.at(n_Objects + i).Pt());
                        Logger::get("lorentzvector::PropagateToMET")
                            ->debug("uncorrected_object pt - {}", LVs.at(i).Pt());
                        Logger::get("lorentzvector::PropagateToMET")
                            ->debug("initial MET {}", met.Pt());
                        Logger::get("lorentzvector::PropagateToMET")
                            ->debug("corrected MET {}", corrected_met.Pt());
                    }
                    return corrected_met;
                }),
            LV_list);
    } else {
        // if we do not apply the propagation, just rename the met column to
        // the new outputname and dont change anything else
        return event::quantity::Rename<ROOT::Math::PtEtaPhiMVector>(df, outputname, p4_met);
    }
}
} // end namespace lorentzvector

namespace physicsobject {

ROOT::RDF::RNode PropagateToMET(
    ROOT::RDF::RNode df, const std::string &outputname, const std::string &p4_met,
    const std::string &pt_corrected, const std::string &eta_corrected,
    const std::string &phi_corrected, const std::string &mass_corrected,
    const std::string &pt, const std::string &eta,
    const std::string &phi, const std::string &mass,
    bool apply_propagation, float min_pt);
} // end namespace physicsobject
#endif /* GUARD_MET_H */
