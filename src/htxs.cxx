#ifndef GUARDHTXS_H
#define GUARDHTXS_H

#include "../include/basefunctions.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/ggF_qcd_uncertainty_2017.hxx"
#include "../include/utility/qq2Hqq_uncert_scheme.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TGraphErrors.h"
/// namespace used for HTXS related functions
namespace htxs {
/**
 * @brief Function to derive the ggH NNLO weights
 *
 * @param df the input dataframe
 * @param weight_name Name of the derived weight in the dataframe.
 * @param rootfilename Path to the rootfile containing the weight graphs.
 * Corresponding cutoffs are hardcoded in this function.
 * @param generator Generator that was used to simulate the ggH sample, either
 * powheg or amcatnlo.
 * @param htxs_pth Name of the column with pt(H) from the htxs module.
 * @param htxs_njets Name of the column with the number of jets from the htxs
 * module.
 * @returns a dataframe with the weight column included
 */
ROOT::RDF::RNode
ggHNNLOWeights(ROOT::RDF::RNode df, const std::string &weight_name,
               const std::string &rootfilename, const std::string &generator,
               const std::string &htxs_pth, const std::string &htxs_njets) {
    TGraphErrors *WeightsGraphs[4];
    TFile rootFile(rootfilename.c_str(), "READ");
    if (generator == "powheg") {
        WeightsGraphs[0] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_powheg_0jet");
        WeightsGraphs[1] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_powheg_1jet");
        WeightsGraphs[2] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_powheg_2jet");
        WeightsGraphs[3] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_powheg_3jet");
    } else if (generator == "amcatnlo") {
        WeightsGraphs[0] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_mcatnlo_0jet");
        WeightsGraphs[1] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_mcatnlo_1jet");
        WeightsGraphs[2] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_mcatnlo_2jet");
        WeightsGraphs[3] =
            (TGraphErrors *)rootFile.Get("gr_NNLOPSratio_pt_mcatnlo_3jet");
    } else
        Logger::get("ggHNNLOWeights")
            ->critical("WARNING: Invalid ggH generator configured. "
                       "ggHNNLOWeights cannot be determined!");
    rootFile.Close();
    const Float_t cutoff[4] = {125.0, 625.0, 800.0, 925.0};
    auto readout_lambda = [WeightsGraphs, cutoff](const Float_t &htxs_pth,
                                                  const UChar_t &htxs_njets) {
        int njets = std::min(3, int(htxs_njets));
        return WeightsGraphs[njets]->Eval(std::min(htxs_pth, cutoff[njets]));
    };
    auto df1 = df.Define(weight_name, readout_lambda, {htxs_pth, htxs_njets});
    return df1;
}

/**
 * @brief Function to derive the WG1 ggH uncertainty weights.
 *
 * @param df the input dataframe
 * @param weight_names Names of the derived weight in the dataframe in the order
 * given by the WG1 macro.
 * @param htxs_flag Name of the column with the htxs stage1 (NOT 1.1, 1.2 or
 * later!) flag.
 * @param htxs_pth Name of the column with pt(H) from the htxs module.
 * @param htxs_njets Name of the column with the number of jets from the htxs
 * module.
 * @returns a dataframe with the weight column included.
 */
ROOT::RDF::RNode
ggH_WG1_uncertainties(ROOT::RDF::RNode df,
                      const std::vector<std::string> &weight_names,
                      const std::string &htxs_flag, const std::string &htxs_pth,
                      const std::string &htxs_njets) {
    auto df1 = df.Define(
        "ggH_WG1_uncertainties",
        [](const Int_t &flag, const Float_t &pth, const UChar_t &njets) {
            return qcd_ggF_uncertSF_2017(njets, pth, flag);
        },
        {htxs_flag, htxs_pth, htxs_njets});
    auto df2 = basefunctions::UnrollVectorQuantity<double>(
        df1, "ggH_WG1_uncertainties", weight_names);
    return df2;
}

/**
 * @brief Function to derive the WG1 qqH uncertainty weights. The application is
 * explicitly restricted to qqH events according to the STXS flag such that e.g.
 * VH samples can be run with this but VHlep events obtain a weight of 1.0.
 *
 * @param df the input dataframe
 * @param weight_names Names of the derived weight in the dataframe in the order
 * given by the WG1 macro.
 * @param htxs_flag Name of the column with the fine htxs stage1.1 flag.
 * module.
 * @param idx initial uncertainty index. It is 0 by default and should not be
 * set by the user. The parameter is needed for the recursive operation of this
 * function.
 * @returns a dataframe with the weight columns included.
 */
ROOT::RDF::RNode
qqH_WG1_uncertainties(ROOT::RDF::RNode df,
                      const std::vector<std::string> &weight_names,
                      const std::string &htxs_flag, const size_t &idx = 0) {
    if (idx >= weight_names.size()) {
        return df;
    }
    auto dflamda = [idx](const int &stxs1flag) {
        if (stxs1flag >= 200 && stxs1flag < 300) {
            return vbf_uncert_stage_1_1(idx, stxs1flag, 1.0);
        } else {
            return 1.0;
        }
    };
    auto df1 = df.Define(weight_names.at(idx), dflamda, {htxs_flag});
    return qqH_WG1_uncertainties(df1, weight_names, htxs_flag, idx + 1);
}
} // namespace htxs
#endif /* GUARDHTXS_H */