#include <event.hxx>
#include <utility/CorrectionManager.hxx>
#include <utility/Logger.hxx>
#include <utility/ggF_qcd_uncertainty_2017.hxx>
#include <utility/qq2Hqq_uncert_scheme.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TGraphErrors.h>
#include <array>
#include <memory>
#include <nlohmann/json.hpp>
#include <set>
#include <unordered_map>
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

// `samples` section of the STXS table
// `nick -> {family, bins: {htxs_bin -> [norm, R_up, R_down]}}`
// (written by
// https://github.com/KIT-CMS/TAUER/blob/master/analyses/smhtt/common/stxs_tools.py)
struct STXSTable {
    struct Entry {
        int family;
        std::unordered_map<int, std::array<double, 3>> bins;
    };
    std::vector<Entry> entries;
    std::unordered_map<std::string, int> index;

    // index of the nick in `entries`, -1 if the nick is not in the table
    int find(const std::string &nick) const {
        auto it = index.find(nick);
        return it == index.end() ? -1 : it->second;
    }
    // [norm, R_up, R_down] of a bin, nullptr if the nick or the bin is not in
    // the table
    const std::array<double, 3> *lookup(int idx, int bin) const {
        if (idx < 0) {
            return nullptr;
        }
        auto it = entries.at(idx).bins.find(bin);
        return it == entries.at(idx).bins.end() ? nullptr : &it->second;
    }
};

std::shared_ptr<const STXSTable>
LoadSTXSTable(correctionManager::CorrectionManager &correction_manager,
              const std::string &stxs_table_path) {
    const nlohmann::json &json = *correction_manager.loadjson(stxs_table_path);
    auto table = std::make_shared<STXSTable>();
    for (const auto &[nick, entry] : json.at("samples").items()) {
        STXSTable::Entry stxs_entry{entry.at("family").get<int>(), {}};
        for (const auto &[bin, factors] : entry.at("bins").items()) {
            stxs_entry.bins[std::stoi(bin)] = {factors.at(0).get<double>(),
                                               factors.at(1).get<double>(),
                                               factors.at(2).get<double>()};
        }
        table->index[nick] = table->entries.size();
        table->entries.push_back(stxs_entry);
    }
    return table;
}

} // namespace

/**
 * @brief Looks up the STXS normalization of each event in the `samples`
 * section of the STXS table (written by
 * https://github.com/KIT-CMS/TAUER/blob/master/analyses/smhtt/common/stxs_tools.py):
 * `nick -> {family, bins: {htxs_bin -> [norm, R_up, R_down]}}`. The sample
 * nick is read from the `nick` column of the ntuple. Nicks missing from the
 * table (and bins missing from it) get family 0 and norm 1.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * table
 * @param family_output name of the column holding the sample family code
 * @param norm_output name of the column with the per-bin normalization factor
 * @param stxs_table_path path to the STXS JSON table
 * @param nick name of the column containing the sample nick
 * @param htxs_bin name of the HTXS stage-1.2 category column
 *
 * @return a dataframe with the two new columns
 */
ROOT::RDF::RNode
Normalization(ROOT::RDF::RNode df,
              correctionManager::CorrectionManager &correction_manager,
              const std::string &family_output, const std::string &norm_output,
              const std::string &stxs_table_path, const std::string &nick,
              const std::string &htxs_bin) {
    auto table = LoadSTXSTable(correction_manager, stxs_table_path);
    const std::string index_column = family_output + "_nick_index";
    auto df1 = df.Define(
        index_column,
        [table](const std::string &nick) { return table->find(nick); }, {nick});
    auto df2 =
        df1.Define(family_output,
                   [table](int idx) {
                       return idx < 0 ? 0 : table->entries.at(idx).family;
                   },
                   {index_column});
    return df2.Define(norm_output,
                      [table](int idx, int bin) {
                          auto factors = table->lookup(idx, bin);
                          return factors ? (*factors)[0] : 1.0;
                      },
                      {index_column, htxs_bin});
}

/**
 * @brief This function creates the column with the LHEScaleWeight factor of
 * one nuisance of the STXS scheme, selected by `variation`:
 * - `nominal`: 1 for all events
 * - `<nuisance>_up` / `<nuisance>_down`: for the events in the (family, HTXS
 * bin) pairs the nuisance applies to, the LHE scale weight divided by the
 * per-bin ratio `R_up` / `R_down` of the STXS table, and 1 for all others.
 *
 * Dividing by `R_up` / `R_down` removes the effect of the scale variation on
 * the cross section of the HTXS bin, keeping only its effect on the
 * acceptance. Everything is taken from the STXS table (written by TAUER's
 * analyses/smhtt/common/stxs_tools.py): the family and the ratios from its
 * `samples` section `nick -> {family, bins: {htxs_bin -> [norm, R_up,
 * R_down]}}`, using the `nick` column of the ntuple, and the (family, HTXS
 * bin) pairs of the nuisance from its `lhe_scale_nuisances` section
 * `nuisance -> [[family, HTXS bin], ...]`. Events of nicks missing from the
 * table get 1, bins missing from it keep the unnormalized LHE scale weight.
 *
 * @param df input dataframe
 * @param correction_manager correction manager responsible for loading the
 * table
 * @param outputname name of the new column
 * @param stxs_table_path path to the STXS JSON table
 * @param variation `nominal` or `<nuisance>_up` / `<nuisance>_down`
 * @param nick name of the column containing the sample nick
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
         const std::string &outputname, const std::string &stxs_table_path,
         const std::string &variation, const std::string &nick,
         const std::string &htxs_bin, const std::string &lhe_scale_weights) {
    if (variation == "nominal") {
        return df.Define(outputname, []() { return 1.0; });
    }
    const auto split = variation.rfind('_');
    const std::string nuisance = variation.substr(0, split);
    const bool use_up = variation.substr(split + 1) == "up";
    const nlohmann::json &stxs_table =
        *correction_manager.loadjson(stxs_table_path);
    auto members = std::make_shared<std::set<std::pair<int, int>>>();
    for (const auto &pair : stxs_table.at("lhe_scale_nuisances").at(nuisance)) {
        members->insert({pair.at(0).get<int>(), pair.at(1).get<int>()});
    }
    auto table = LoadSTXSTable(correction_manager, stxs_table_path);
    const std::string index_column = outputname + "_nick_index";
    auto df1 = df.Define(
        index_column,
        [table](const std::string &nick) { return table->find(nick); }, {nick});
    return df1.Define(
        outputname,
        [table, members, use_up](int idx, int bin,
                                 const ROOT::RVec<float> &scale_weights) {
            if (idx < 0 ||
                !members->count({table->entries.at(idx).family, bin})) {
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
            const double scale = scale_weights.at(index);
            auto factors = table->lookup(idx, bin);
            return factors ? scale / (*factors)[use_up ? 1 : 2] : scale;
        },
        {index_column, htxs_bin, lhe_scale_weights});
}
} // namespace stxs
