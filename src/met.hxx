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
    auto GetMetPt = [](const ROOT::Math::PtEtaPhiEVector met) {
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
} // end namespace met
