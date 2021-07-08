#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TH1.h"
#include "utility/Logger.hxx"

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
auto puweights(auto &df, const std::string &weightname,
               const std::string &truePUMean, const std::string &filename,
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
} // namespace reweighting