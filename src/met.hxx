#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "basefunctions.hxx"
#include "utility/Logger.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>

namespace met {
/**
 * @brief Function to get pT from the met Lorentz vector
 *
 * @param df the input dataframe
 * @param metvector the name of the column containing the met lorentz vector
 * @param outputname name of the new column containing the met pT value
 * @return a new dataframe containing the met pT value
 */
auto metPt(auto df, const std::string &metvector,
           const std::string &outputname) {
    auto GetMetPt = [](const ROOT::Math::PtEtaPhiEVector &met) {
        return met.Pt();
    };
    return df.Define(outputname, GetMetPt, {metvector});
}
/**
 * @brief Function to get Phi from the met Lorentz vector
 *
 * @param df the input dataframe
 * @param metvector the name of the column containing the met lorentz vector
 * @param outputname name of the new column containing the met phi value
 * @return a new dataframe containing the met phi value
 */
auto metPhi(auto df, const std::string &metvector,
            const std::string &outputname) {
    auto GetMetPhi = [](const ROOT::Math::PtEtaPhiEVector met) {
        return met.Phi();
    };
    return df.Define(outputname, GetMetPhi, {metvector});
}
/**
 * @brief Function used to propagate lepton corrections to the met. If the
 energy of a lepton is corrected (via some scale factor) or due to a shift, this
 change in energy has to be propagated to the met vector, and the met vector has
 to be adapted accordingly. The met is recalculated via
 @code
  Recalculate MET with corrected lepton energies :
  MetX_corrected = MetX + Px - Px_corrected
  MetY_corrected = MetY + Py - Py_corrected
  MET_corrected = sqrt(MetX_corrected * MetX_corrected + MetY_corrected *
 MetY_corrected)
 @endcode
 * @param df the input dataframe
 * @param met the uncorrected met lorentz vector
 * @param p4_1_uncorrected the uncorrected lorentz vector of the first lepton
 * @param p4_2_uncorrected the uncorrected lorentz vector of the second lepton
 * @param p4_1 the corrected lorentz vector of the first lepton
 * @param p4_2 the corrected lorentz vector of the second lepton
 * @param outputname name of the column containing the corrected met lorentz
 vector
 * @return a new df containing the corrected met lorentz vector
 */
auto propagateLeptons(auto df, const std::string &met,
                      const std::string &p4_1_uncorrected,
                      const std::string &p4_2_uncorrected,
                      const std::string &p4_1, const std::string &p4_2,
                      const std::string &outputname) {
    auto scaleMet = [](const ROOT::Math::PtEtaPhiEVector &met,
                       const ROOT::Math::PtEtaPhiMVector &uncorrected_object,
                       const ROOT::Math::PtEtaPhiMVector &corrected_object) {
        // We propagate the lepton corrections to the MET by scaling the x and y
        // component of the MET according to the correction of the lepton
        // Recalculate MET with corrected lepton energies :
        // MetX_corrected = MetX + Px - Px_corrected
        // MetY_corrected = MetY + Py - Py_corrected
        // MET_corrected = sqrt(MetX_corrected * MetX_corrected + MetY_corrected
        // * MetY_corrected)
        float corr_x = uncorrected_object.Px() - corrected_object.Px();
        float corr_y = uncorrected_object.Py() - corrected_object.Py();
        float MetX = met.Px() - corr_x;
        float MetY = met.Py() - corr_y;
        Logger::get("propagateLeptons")->debug("corr_x {}", corr_x);
        Logger::get("propagateLeptons")->debug("corr_y {}", corr_y);
        Logger::get("propagateLeptons")->debug("MetX {}", MetX);
        Logger::get("propagateLeptons")->debug("MetY {}", MetY);
        ROOT::Math::PtEtaPhiEVector corrected_met;
        corrected_met.SetPxPyPzE(MetX, MetY, 0,
                                 std::sqrt(MetX * MetX + MetY * MetY));
        Logger::get("propagateLeptons")
            ->debug("corrected_object pt - {}", corrected_object.Pt());
        Logger::get("propagateLeptons")
            ->debug("uncorrected_object pt - {}", uncorrected_object.Pt());
        Logger::get("propagateLeptons")->debug("old met {}", met.Pt());
        Logger::get("propagateLeptons")
            ->debug("corrected met {}", corrected_met.Pt());
        return corrected_met;
    };
    // first correct for the first lepton, store the met in an
    // intermediate column
    Logger::get("propagateLeptons")
        ->debug("Setting up correction for first lepton {}", p4_1);
    auto df1 = df.Define(outputname + "_intermidiate", scaleMet,
                         {met, p4_1_uncorrected, p4_1});
    // after the second lepton correction, the correct output column is used
    Logger::get("propagateLeptons")
        ->debug("Setting up correction for second lepton {}", p4_2);
    auto df2 =
        df1.Define(outputname, scaleMet,
                   {outputname + "_intermidiate", p4_2_uncorrected, p4_2});
    return df2;
}
} // end namespace met
