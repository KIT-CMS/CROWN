#pragma once

#include <array>
#include <optional>
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <utility/CorrectionManager.hxx>

namespace trigger {

bool matchParticle(
    const ROOT::Math::PtEtaPhiMVector &particle,
    ROOT::RVec<float> &triggerobject_pts, ROOT::RVec<float> &triggerobject_etas,
    ROOT::RVec<float> &triggerobject_phis, ROOT::RVec<int> &triggerobject_ids,
    ROOT::RVec<int> &triggerobject_filterbits, const float &pt_threshold,
    const float &eta_threshold, const int &trigger_particle_id_value,
    const std::vector<int> &trigger_bit_value, const float &deltaR_threshold);

namespace experimental {
ROOT::RDF::RNode
ObjectFlag(ROOT::RDF::RNode df, const std::string &outputname,
           const std::vector<std::string> &inputs, const std::string &hlt_path,
           const std::vector<float> &pt_thresholds,
           const std::vector<float> &eta_thresholds,
           const std::vector<int> &trigger_particle_id_values,
           const std::vector<std::vector<int>> &trigger_bit_values,
           const float &deltaR_threshold);
} // end namespace experimental

ROOT::RDF::RNode SingleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit, const std::string &hlt_path,
    const float &pt_threshold, const float &eta_threshold,
    const int &trigger_particle_id_value,
    const std::vector<int> &trigger_bit_value, const float &deltaR_threshold);
ROOT::RDF::RNode SingleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit, const float &pt_threshold,
    const float &eta_threshold, const int &trigger_particle_id_value,
    const std::vector<int> &trigger_bit_value, const float &deltaR_threshold);
ROOT::RDF::RNode DoubleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle_1, const std::string &particle_2,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit, const std::string &hlt_path,
    const float &pt_threshold_1, const float &pt_threshold_2,
    const float &eta_threshold_1, const float &eta_threshold_2,
    const int &trigger_particle_id_value_1,
    const int &trigger_particle_id_value_2,
    const std::vector<int> &trigger_bit_value_1,
    const std::vector<int> &trigger_bit_value_2, const float &deltaR_threshold);
ROOT::RDF::RNode DoubleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle_1, const std::string &particle_2,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit, const float &pt_threshold_1,
    const float &pt_threshold_2, const float &eta_threshold_1,
    const float &eta_threshold_2, const int &trigger_particle_id_value_1,
    const int &trigger_particle_id_value_2,
    const std::vector<int> &trigger_bit_value_1,
    const std::vector<int> &trigger_bit_value_2, const float &deltaR_threshold);

ROOT::RDF::RNode
GetPrescaleValues(ROOT::RDF::RNode df,
                  correctionManager::CorrectionManager &correction_manager,
                  const std::string &outputname, const std::string &hlt_path,
                  const std::string &run, const std::string &lumiblock,
                  const std::string &prescale_file);
ROOT::RDF::RNode SelectScaleFactor(ROOT::RDF::RNode df,
                                   const std::string &outputname,
                                   const std::string &pass_flag,
                                   const std::vector<std::string> &sf_inputs,
                                   const int &n_pass_sfs);
} // end namespace trigger
