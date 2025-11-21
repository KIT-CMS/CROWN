#ifndef GUARD_TRIGGERS_H
#define GUARD_TRIGGERS_H

typedef std::bitset<30> IntBits;

namespace trigger {

bool matchParticle(const ROOT::Math::PtEtaPhiMVector &particle,
                   ROOT::RVec<float> &triggerobject_pts,
                   ROOT::RVec<float> &triggerobject_etas,
                   ROOT::RVec<float> &triggerobject_phis,
                   ROOT::RVec<int> &triggerobject_ids,
                   ROOT::RVec<int> &triggerobject_filterbits,
                   const float &pt_threshold, const float &eta_threshold,
                   const int &trigger_particle_id_value, const int &trigger_bit_value,
                   const float &deltaR_threshold);

ROOT::RDF::RNode SingleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &triggerobject_id, const std::string &triggerobject_filterbit,
    const std::string &hlt_path, const float &pt_threshold, const float &eta_threshold,
    const int &trigger_particle_id_value, const int &trigger_bit_value,
    const float &deltaR_threshold);
ROOT::RDF::RNode SingleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &triggerobject_id, const std::string &triggerobject_filterbit,
    const float &pt_threshold, const float &eta_threshold,
    const int &trigger_particle_id_value, const int &trigger_bit_value,
    const float &deltaR_threshold);
ROOT::RDF::RNode DoubleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle_1, const std::string &particle_2,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit, const std::string &hlt_path,
    const float &pt_threshold_1, const float &pt_threshold_2,
    const float &eta_threshold_1, const float &eta_threshold_2,
    const int &trigger_particle_id_value_1, const int &trigger_particle_id_value_2,
    const int &trigger_bit_value_1, const int &trigger_bit_value_2,
    const float &deltaR_threshold);
ROOT::RDF::RNode DoubleObjectFlag(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &particle_1, const std::string &particle_2,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &triggerobject_id,
    const std::string &triggerobject_filterbit,
    const float &pt_threshold_1, const float &pt_threshold_2,
    const float &eta_threshold_1, const float &eta_threshold_2,
    const int &trigger_particle_id_value_1, const int &trigger_particle_id_value_2,
    const int &trigger_bit_value_1, const int &trigger_bit_value_2,
    const float &deltaR_threshold);

ROOT::RDF::RNode GetPrescaleValues(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correction_manager,
    const std::string &outputname, const std::string &hlt_path,
    const std::string &run, const std::string &lumiblock,
    const std::string &prescale_file);
} // end namespace trigger
#endif /* GUARD_TRIGGERS_H */
