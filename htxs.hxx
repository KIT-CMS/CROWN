#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TGraphErrors.h"
#include "utility/Logger.hxx"
/// namespace used for HTXS related functions
namespace htxs {
/**
 * @brief Function to derive the ggH NNLO weights
 *
 * @param weight_name Name of the derived weight in the dataframe.
 * @param rootfilename Path to the rootfile containing the weight graphs.
 * @param generator Generator that was used to simulate the ggH sample, either powheg or amcatnlo.
 * @param htxs_pth Name of the column with pt(H) from the htxs module.
 * @param htxs_njets Name of the column with the number of jets from the htxs module.
 * @returns a dataframe with the weight column included
 */
auto ggHNLLOWeights(auto &df, const std::string &weight_name,
                    const std::string &rootfilename,
                    const std::string &generator, const std::string &htxs_pth,
                    const std::string &htxs_njets) {
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
        Logger::get("ggHNLLOWeights")
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
} // namespace htxs