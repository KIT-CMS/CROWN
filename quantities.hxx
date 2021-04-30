#include "utility/Logger.hxx"
#include <Math/Vector4D.h>
/// The namespace that is used to hold the functions for basic quantities that
/// are needed for every event
namespace quantities {

/// Function to calculate the pt from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param varSet - vector of variables with are snapsotted in the end
/// \param outputname name of the new column containing the pt value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

auto pt(auto df, std::vector<std::string> varSet, const std::string &outputname,
        const std::string &inputvector) {
    varSet.push_back(outputname);
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return p4.pt(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param varSet - vector of variables with are snapsotted in the end
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

auto eta(auto df, std::vector<std::string> varSet,
         const std::string &outputname, const std::string &inputvector) {
    varSet.push_back(outputname);
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return p4.eta(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param varSet - vector of variables with are snapsotted in the end
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

auto phi(auto df, std::vector<std::string> varSet,
         const std::string &outputname, const std::string &inputvector) {
    varSet.push_back(outputname);
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return p4.phi(); },
        {inputvector});
}
/// Function to calculate the mass from a pair of lorentz vectors and add it to
/// the dataframe
///
/// \param df the dataframe to add the quantity to
/// \param varSet - vector of variables with are snapsotted in the end
/// \param outputname name of the new column containing the pt value
/// \param inputvector1 name of the column containing the first lorentz vector
/// \param inputvector2 name of the column containing the second lorentz vector
///
/// \returns a dataframe with the new column

auto m_vis(auto df, std::vector<std::string> varSet,
           const std::string &outputname,
           const std::vector<std::string> &inputvectors) {
    varSet.push_back(outputname);
    // build visible mass from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            auto const dileptonsystem = p4_1 + p4_2;
            return dileptonsystem.mass();
        },
        inputvectors);
}
} // end namespace quantities