#ifndef GUARD_REWEIGHTING_H
#define GUARD_REWEIGHTING_H

#include "../include/basefunctions.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TH1.h"
#include <Math/Vector4D.h>

/// namespace used for reweighting related functions
namespace reweighting {
/**
 * @brief Function used to read out pileup weights
 *
 * @param df The input dataframe
 * @param weightname name of the derived weight
 * @param truePUMean name of the column containing the true PU mean of simulated
 * events
 * @param filename path to the rootfile
 * @param histogramname name of the histogram stored in the rootfile
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode puweights(auto &df, const std::string &weightname,
                           const std::string &truePUMean,
                           const std::string &filename,
                           const std::string &histogramname) {

    float bin_density = 1.0;
    std::vector<float> puweights;
    Logger::get("puweights")
        ->debug("Loading pile-up weights from {}", filename);
    {
        TFile inputfile(filename.c_str(), "READ");
        TH1D *puhist = (TH1D *)inputfile.Get(histogramname.c_str());
        for (int i = 0; i <= puhist->GetNbinsX(); ++i) {
            puweights.push_back(puhist->GetBinContent(i));
        }
        bin_density = 1.0 / puhist->GetBinWidth(1);
        delete puhist;
        inputfile.Close();
    }

    auto puweightlambda = [bin_density, puweights](const float pu) {
        size_t puBin = static_cast<size_t>(pu * bin_density);
        return puBin < puweights.size() ? puweights[puBin] : 1.0;
    };
    auto df1 = df.Define(weightname, puweightlambda, {truePUMean});
    return df1;
}

/**
 * @brief Function used to calculate top pt reweighting
 *
 * @param df The input dataframe
 * @param weightname name of the derived weight
 * @param gen_pdgids name of the column containing the PDG-IDs of the generator
 * particles
 * @param gen_status name of the column containing the status flags of the
 * generator particles, where bit 13 contains the isLastCopy flag
 * @param gen_pt name of the column containing the pt of the generator particles
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode topptreweighting(auto &df, const std::string &weightname,
                                  const std::string &gen_pdgids,
                                  const std::string &gen_status,
                                  const std::string &gen_pt) {

    auto ttbarreweightlambda = [](const ROOT::RVec<int> pdgid,
                                  const ROOT::RVec<int> status,
                                  const ROOT::RVec<float> pt) {
        std::vector<float> top_pts;
        for (size_t i = 0; i < pdgid.size(); i++) {
            if (std::abs(pdgid[i]) == 6 && ((status[i] >> 13) & 1) == 1)
                top_pts.push_back(pt[i]);
        }
        if (top_pts.size() != 2) {
            std::cout << top_pts.size();
            Logger::get("topptreweighting")
                ->error("TTbar reweighting applied to event with not exactly "
                        "two top quarks. Probably due to wrong sample type.");
            throw std::runtime_error("Bad number of top quarks.");
        }
        if (top_pts[0] > 472.0)
            top_pts[0] = 472.0;
        if (top_pts[1] > 472.0)
            top_pts[1] = 472.0;
        const float parameter_a = 0.088;
        const float parameter_b = -0.00087;
        const float parameter_c = 0.00000092;
        return sqrt(exp(parameter_a + parameter_b * top_pts[0] +
                        parameter_c * top_pts[0] * top_pts[0]) *
                    exp(parameter_a + parameter_b * top_pts[1] +
                        parameter_c * top_pts[1] * top_pts[1]));
    };
    auto df1 = df.Define(weightname, ttbarreweightlambda,
                         {gen_pdgids, gen_status, gen_pt});
    return df1;
}

/**
 * @brief Function used to evaluate Z pt mass weights
 *
 * @param df The input dataframe
 * @param weightname name of the generated weight
 * @param gen_boson name of the column that contains a pair of Lorentzvectors,
 * where the first one is the one of the genboson
 * @param workspace_file path to the file which contains the workspace to be
 * read
 * @param functor_name name of the function from the workspace
 * @param argset arguments of the function
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode zPtMassReweighting(auto &df, const std::string &weightname,
                                    const std::string &gen_boson,
                                    const std::string &workspace_file,
                                    const std::string &functor_name,
                                    const std::string &argset) {

    // retrieve pt and mass of gen boson reconstructed with the method used by
    // recoil corrections; resulting quantities are only for the purpose of this
    // method
    auto df1 = df.Define(gen_boson + "_pt",
                         [](const std::pair<ROOT::Math::PtEtaPhiMVector,
                                            ROOT::Math::PtEtaPhiMVector> &p4) {
                             return (float)p4.first.pt();
                         },
                         {gen_boson});
    auto df2 = df1.Define(gen_boson + "_mass",
                          [](const std::pair<ROOT::Math::PtEtaPhiMVector,
                                             ROOT::Math::PtEtaPhiMVector> &p4) {
                              return (float)p4.first.mass();
                          },
                          {gen_boson});

    // set up workspace
    Logger::get("zPtMassReweighting")
        ->debug("Setting up functions for zPtMassReweighting");
    Logger::get("zPtMassReweighting")
        ->debug("zPtMassReweighting - Function {} // argset {}", functor_name,
                argset);

    const std::shared_ptr<RooFunctorThreadsafe> weight_function =
        loadFunctor(workspace_file, functor_name, argset);
    auto df3 = basefunctions::evaluateWorkspaceFunction(
        df2, weightname, weight_function, gen_boson + "_mass",
        gen_boson + "_pt");
    return df3;
}
} // namespace reweighting
#endif /* GUARD_REWEIGHTING_H */