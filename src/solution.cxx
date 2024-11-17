#ifndef GUARD_SOLUTION_H
#define GUARD_SOLUTION_H

// #include "../include/defaults.hxx"
#include "../include/basefunctions.hxx"
// #include "ROOT/RDFHelpers.hxx"
// #include "ROOT/RDataFrame.hxx"
// #include "ROOT/RVec.hxx"

#include "../include/defaults.hxx"
#include "../include/utility/CorrectionManager.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/Logger.hxx"
#include "utility/RooFunctorThreadsafe.hxx"
#include "utility/utility.hxx"
#include <nlohmann/json.hpp>

namespace solution {

/// Function to select objects with number of Muons
///
/// \param[in] df the input dataframe
/// \param[in] nMuons name of the number of Muons column in the NanoAOD dataframe
/// \param[in] nMuons_req number of required Muons
///
/// \return a dataframe filtered by the cut
ROOT::RDF::RNode CutNMuon(ROOT::RDF::RNode df, const std::string &nMuons, const unsigned int &nMuons_req) {
    auto checkNMuon = [nMuons_req](const unsigned int &nM) {
        bool mask = nM == nMuons_req;
        return mask;
    };
    auto df1 = df.Filter(checkNMuon, {nMuons}, "CutNMuon");
    return df1;
}

/// Function to select objects with number of Muons
///
/// \param[in] df the input dataframe
/// \param[in] CMuons name of the charge of Muons column in the NanoAOD dataframe
/// \param[in] CMuons_req required sign of Muons after multiply
///
/// \return a dataframe filtered by the cut
ROOT::RDF::RNode CutCMuons(ROOT::RDF::RNode df, const std::string &CMuons, const int &CMuons_req) {
    auto checkCMuonMul = [CMuons_req](const ROOT::RVec<int> &charges) {
        bool mask = (charges[0] * charges[1]) == CMuons_req;
        return mask;
    };
    auto df1 = df.Filter(checkCMuonMul, {CMuons}, "CutCMuons");
    return df1;
}

/// Function to calculate invariant mass for objects, using
/// ROOT::VecOps::InvariantMass
///
/// \param[in] df the input dataframe
/// \param[in] new_column name of the new column in the ntuple
/// \param[in] particle_pts name of the pt column in the NanoAOD dataframe
/// \param[in] particle_etas name of the eta column in the NanoAOD dataframe
/// \param[in] particle_phis name of the phi column in the NanoAOD dataframe
/// \param[in] particle_masses name of the mass column in the NanoAOD dataframe
///
/// \return a dataframe containing the new column
ROOT::RDF::RNode GetInvariantMass(
    ROOT::RDF::RNode df, 
    const std::string &new_column, 
    const std::string &particle_pts, 
    const std::string &particle_etas, 
    const std::string &particle_phis, 
    const std::string &particle_masses
) {
    auto InvariantMassLambda = [](
        const ROOT::RVec<float>& pts, 
        const ROOT::RVec<float>& etas, 
        const ROOT::RVec<float>& phis, 
        const ROOT::RVec<float>& masses
    ) {
        return ROOT::VecOps::InvariantMass(pts, etas, phis, masses);
    };
    auto df1 = df.Define(
        new_column, 
        InvariantMassLambda, 
        {particle_pts, particle_etas, particle_phis, particle_masses}
    );
    return df1;
}

/// Function to select objects with number of Muons
///
/// \param[in] df the input dataframe
/// \param[in] new_column name of the new column in the ntuple
/// \param[in] CMuons name of the charge of Muons column in the NanoAOD dataframe
///
/// \return a dataframe filtered by the cut
ROOT::RDF::RNode MuonCSum(ROOT::RDF::RNode df, const std::string &new_column, const std::string &CMuons) {
    auto CSum = [](const ROOT::RVec<int> &charges) {
        return charges[0] + charges[1];
    };
    auto df1 = df.Define(new_column, CSum, {CMuons});
    return df1;
}

} // namespace solution
#endif /* GUARD_SOLUTION_H */
