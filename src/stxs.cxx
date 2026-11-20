#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TGraphErrors.h>
#include <array>
#include <event.hxx>
#include <memory>
#include <nlohmann/json.hpp>
#include <set>
#include <unordered_map>
#include <utility/CorrectionManager.hxx>
#include <utility/Logger.hxx>
#include <utility/ggF_qcd_uncertainty_2017.hxx>
#include <utility/qq2Hqq_uncert_scheme.hxx>
#include <vector>

/// namespace used for STXS related functions
namespace stxs {
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
 *
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
 *
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
    auto df2 = event::quantity::Unroll<double>(df1, weight_names,
                                               "ggH_WG1_uncertainties");
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
 *
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

namespace {

// `stxs` field of the sample information `{family, bins: {htxs_bin -> [norm,
// R_up, R_down]}}` (written to the sample database from
// https://github.com/KIT-CMS/TAUER/blob/master/analyses/smhtt/common/stxs_tools.py)
struct STXSEntry {
    int family = 0;
    std::unordered_map<int, std::array<double, 3>> bins;

    // [norm, R_up, R_down] of a bin, nullptr if the bin is not in the entry
    const std::array<double, 3> *lookup(int bin) const {
        auto it = bins.find(bin);
        return it == bins.end() ? nullptr : &it->second;
    }
};

std::shared_ptr<const STXSEntry>
LoadSTXSEntry(correctionManager::CorrectionManager &correction_manager) {
    const nlohmann::json &sample_info =
        event::quantity::LoadSampleInfo(correction_manager);
    auto entry = std::make_shared<STXSEntry>();
    if (!sample_info.contains("stxs") || sample_info.at("stxs").is_null()) {
        Logger::get("stxs::LoadSTXSEntry")
            ->warn("No STXS entry in the sample information of {}, family 0 "
                   "and norm 1",
                   sample_info.value("nick", std::string("?")));
        return entry;
    }
    const nlohmann::json &stxs = sample_info.at("stxs");
    entry->family = stxs.at("family").get<int>();
    for (const auto &[bin, factors] : stxs.at("bins").items()) {
        entry->bins[std::stoi(bin)] = {factors.at(0).get<double>(),
                                       factors.at(1).get<double>(),
                                       factors.at(2).get<double>()};
    }
    return entry;
}

} // namespace

/**
 * @brief Defines the STXS normalization of each event from the `stxs` field
 * of the sample information (see `event::quantity::LoadSampleInfo`):
 * `{family, bins: {htxs_bin -> [norm, R_up, R_down]}}`. Samples without the
 * field get family 0 and norm 1, bins missing from it norm 1.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * sample information
 * @param family_output name of the column holding the sample family code
 * @param norm_output name of the column with the per-bin normalization factor
 * @param htxs_bin name of the HTXS stage-1.2 category column
 *
 * @return a dataframe with the two new columns
 */
ROOT::RDF::RNode
Normalization(ROOT::RDF::RNode df,
              correctionManager::CorrectionManager &correction_manager,
              const std::string &family_output, const std::string &norm_output,
              const std::string &htxs_bin) {
    auto entry = LoadSTXSEntry(correction_manager);
    const int family = entry->family;
    auto df1 = df.Define(family_output, [family]() { return family; });
    return df1.Define(norm_output,
                      [entry](int bin) {
                          auto factors = entry->lookup(bin);
                          return factors ? (*factors)[0] : 1.0;
                      },
                      {htxs_bin});
}

/**
 * @brief This function creates the column with the LHEScaleWeight factor of
 * one nuisance of the STXS scheme, selected by `variation`:
 * - `nominal`: 1 for all events
 * - `<nuisance>_up` / `<nuisance>_down`: for the events in the (family, HTXS
 * bin) pairs the nuisance applies to, the LHE scale weight relative to the
 * central one divided by the per-bin ratio `R_up` / `R_down`, and 1 for all
 * others.
 *
 * Dividing by `R_up` / `R_down` removes the effect of the scale variation on
 * the cross section of the HTXS bin, keeping only its effect on the
 * acceptance. The family and the ratios are taken from the `stxs` field of the
 * sample information `{family, bins: {htxs_bin -> [norm, R_up, R_down]}}`,
 * the (family, HTXS bin) pairs of the nuisance from the nuisance table
 * `nuisance -> [[family, HTXS bin], ...]` (both written by TAUER's
 * analyses/smhtt/common/stxs_tools.py). Samples without the `stxs` field get
 * 1, bins missing from it keep the unnormalized LHE scale weight.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * nuisance table and the sample information
 * @param outputname name of the new column
 * @param nuisance_table_path path to the JSON nuisance table
 * @param variation `nominal` or `<nuisance>_up` / `<nuisance>_down`
 * @param htxs_bin name of the HTXS stage-1.2 category column
 * @param lhe_scale_weights of the `LHEScaleWeight` column, its
 * mu_R = mu_F = 0.5 (`down`) is index 0, mu_R = mu_F = 2 (`up`) index 8 with
 * 9 entries and 7 with 8 entries (see event::reweighting::LHEscale)
 *
 * @return a dataframe with the new column of type `double`
 */
ROOT::RDF::RNode
LHEScale(ROOT::RDF::RNode df,
         correctionManager::CorrectionManager &correction_manager,
         const std::string &outputname, const std::string &nuisance_table_path,
         const std::string &variation, const std::string &htxs_bin,
         const std::string &lhe_scale_weights) {
    if (variation == "nominal") {
        return df.Define(outputname, []() { return 1.0; });
    }
    const auto split = variation.rfind('_');
    const std::string nuisance = variation.substr(0, split);
    const bool use_up = variation.substr(split + 1) == "up";
    auto entry = LoadSTXSEntry(correction_manager);
    // HTXS bins of this sample's family the nuisance applies to
    auto bins = std::make_shared<std::set<int>>();
    for (const auto &pair :
         correction_manager.loadjson(nuisance_table_path)->at(nuisance)) {
        if (pair.at(0).get<int>() == entry->family) {
            bins->insert(pair.at(1).get<int>());
        }
    }
    return df.Define(
        outputname,
        [entry, bins, use_up](int bin, const ROOT::RVec<float> &scale_weights) {
            if (!bins->count(bin)) {
                return 1.0;
            }
            if (scale_weights.size() != 8 && scale_weights.size() != 9) {
                Logger::get("stxs::LHEScale")
                    ->error("Invalid number of LHE scale weights: {}",
                            scale_weights.size());
                throw std::runtime_error("Invalid number of LHE scale weights");
            }
            // mu_R = mu_F = 2: index 8 with the nominal entry (9 weights),
            // 7 without it (8 weights); mu_R = mu_F = 0.5: index 0 in both
            const std::size_t index =
                use_up ? (scale_weights.size() == 9 ? 8 : 7) : 0;
            // relative to the central weight as the ratios, not 1 in W+-H
            // MiNNLO
            const double central =
                scale_weights.size() == 9 && scale_weights.at(4) != 0
                    ? scale_weights.at(4)
                    : 1.0;
            const double scale = scale_weights.at(index) / central;
            auto factors = entry->lookup(bin);
            return factors ? scale / (*factors)[use_up ? 1 : 2] : scale;
        },
        {htxs_bin, lhe_scale_weights});
}
} // namespace stxs
