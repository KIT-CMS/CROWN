#ifndef GUARD_MET_H
#define GUARD_MET_H

#include "../include/utility/utility.hxx"
#include "../include/event.hxx"

namespace met {

/**
 * @brief Calculates the hadronic recoil and its components parallel and perpendicular
 * to a visible system's momentum. The recoil is a crucial observable in analyses with
 * invisible particles (like neutrinos).
 *
 * The hadronic recoil (\f$\vec{u}_T\f$) is defined as the vector sum of the transverse
 * momenta of all hadronic particles in the event. By momentum conservation in the
 * transverse plane, it must balance the momentum of the hard-scatter system (visible
 * leptons/photons and invisible MET).
 *
 * It is calculated as:
 * \f[
 * \vec{u}_T = -(\vec{p}_T^{\text{visible}} + \vec{p}_T^{\text{miss}})
 * \f]
 * where \f$\vec{p}_T^{\text{visible}}\f$ is the vector sum of the transverse momenta
 * of the specified input Lorentz vectors (e.g., \f$\vec{p}_{T, \ell\ell}\f$ for a
 * Z boson event). The function then projects this recoil vector onto the axis defined
 * by the direction of the visible system's transverse momentum,
 * \f$\hat{q}_T = \frac{\vec{p}_T^{\text{visible}}}{|\vec{p}_T^{\text{visible}}|}\f$.
 *
 * - **Parallel Component (\f$u_\parallel\f$)**: The projection of the recoil onto the
 *   visible system's axis. It is a measure of the recoil's response.
 *   \f[ u_\parallel = \vec{u}_T \cdot \hat{q}_T \f]
 * - **Perpendicular Component (\f$u_\perp\f$)**: The component of the recoil orthogonal
 *   to the visible system's axis. It is a measure of the recoil's resolution.
 *   \f[ u_\perp = |\vec{u}_T \times \hat{q}_T| \f] (calculated with a sign).
 *
 * The function returns a vector of two doubles: `{u_parallel, u_perp}`.
 *
 * @tparam Lorentzvectors variadic template parameter pack representing Lorentz vectors
 * @param df input dataframe
 * @param outputname name of the new column containing a `std::vector<double>` with the
 * parallel and perpendicular recoil components
 * @param p4_met name of the column containing the MET Lorentz vector
 * @param lorentzvectors parameter pack of column names containing Lorentz vectors of
 * visible particles (e.g., "Lepton_p4_1", "Lepton_p4_2").
 *
 * @return a dataframe with the new column containing parallel and perpendicular
 * recoil components
 *
 * @warning Implemented calculation was not cross checked.
 */
template <typename... Lorentzvectors>
ROOT::RDF::RNode GetHadronicRecoil(ROOT::RDF::RNode df,
                           const std::string &outputname,
                           const std::string &p4_met,
                           Lorentzvectors... lorentzvectors) {
    auto argTuple = std::make_tuple(lorentzvectors...);
    std::vector<std::string> LorentzvectorList{lorentzvectors...};
    const auto nLVs = sizeof...(Lorentzvectors);

    using namespace ROOT::VecOps;
    auto calc_recoil = [](const ROOT::Math::PtEtaPhiMVector &p4_met,
                          const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &LVs) {
            // Calculate the multilepton vector
            ROOT::Math::PtEtaPhiMVector vis_recoil_vector;
            for (const auto &vector : LVs) {
                vis_recoil_vector += vector;
            }
            double met = p4_met.Et();
            double met_phi = p4_met.Phi();

            double pUX = -met * cos(met_phi) -
                         vis_recoil_vector.Pt() * cos(vis_recoil_vector.Phi());
            double pUY = -met * sin(met_phi) -
                         vis_recoil_vector.Pt() * sin(vis_recoil_vector.Phi());
            double pU = sqrt(pUX * pUX + pUY * pUY);
            double pCos = (pUX * cos(vis_recoil_vector.Phi()) +
                           pUY * sin(vis_recoil_vector.Phi())) /
                          pU;
            double pSin = -(pUX * sin(vis_recoil_vector.Phi()) -
                           pUY * cos(vis_recoil_vector.Phi())) /
                          pU;
            // Return the vector of parallel and perpendicular components
            return std::vector<double>{pU * pCos, pU * pSin};
        };
    return df.Define(outputname,
                     utility::PassAsVec<nLVs, ROOT::Math::PtEtaPhiMVector>(calc_recoil),
                     {p4_met, LorentzvectorList});
};
ROOT::RDF::RNode calculateGenBosonPt(
    ROOT::RDF::RNode df,  const std::string &outputname,
    const std::string &genparticle_pt,
    const std::string &genparticle_eta, const std::string &genparticle_phi,
    const std::string &genparticle_mass, const std::string &genparticle_id,
    const std::string &genparticle_status,
    const std::string &genparticle_statusflag, 
    bool is_data);
ROOT::RDF::RNode RecoilCorrection(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &p4_met, const std::string &p4_gen_boson,
    const std::string &p4_vis_gen_boson, const std::string &jet_pt,
    const std::string &corr_file, const std::string &syst_file,
    const bool apply_correction, const bool resolution, const bool response,
    const bool shift_up, const bool shift_down, const bool is_Wjets);
ROOT::RDF::RNode METPhiCorrection(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &p4_met,
                        const std::string &n_pv,
                        const std::string &run,
                        const std::string &corr_file,
                        const std::string &corr_name);
ROOT::RDF::RNode METPhiCorrection(ROOT::RDF::RNode df,
                        const std::string &outputname,
                        const std::string &p4_met,
                        const std::string &n_pv,
                        const std::string &corr_file,
                        const std::string &corr_name,
                        const std::string &met_type,
                        const std::string &era,
                        const bool is_mc,
                        const std::string &stat_variation,
                        const std::string &pileup_variation);
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
