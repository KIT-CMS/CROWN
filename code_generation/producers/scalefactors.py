import code_generation.quantities.output as q
from code_generation.producer import Producer, ProducerGroup
from code_generation.producer import ExtendedVectorProducer


############################
# Muon ID, ISO SF
# The readout is done via correctionlib
############################

Muon_1_ID_SF_RooWorkspace = Producer(
    name="MuonID_SF_RooWorkspace",
    call='scalefactor::muon::id_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_1, q.eta_1],
    output=[q.idWeight_old_1],
    scopes=["mt", "mm"],
)
Muon_1_Iso_SF_RooWorkspace = Producer(
    name="MuonIso_SF_RooWorkspace",
    call='scalefactor::muon::iso_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
    input=[q.pt_1, q.eta_1, q.iso_1],
    output=[q.isoWeight_old_1],
    scopes=["mt", "mm"],
)
Muon_1_ID_SF = Producer(
    name="MuonID_SF",
    call='scalefactor::muon::id({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_id_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.idWeight_1],
    scopes=["mt", "mm"],
)
Muon_1_Iso_SF = Producer(
    name="MuonIso_SF",
    call='scalefactor::muon::iso({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_iso_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.isoWeight_1],
    scopes=["mt", "mm"],
)
Muon_2_ID_SF_RooWorkspace = Producer(
    name="MuonID_SF_RooWorkspace",
    call='scalefactor::muon::id_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_2, q.eta_2],
    output=[q.idWeight_old_2],
    scopes=["em", "mm"],
)
Muon_2_Iso_SF_RooWorkspace = Producer(
    name="MuonIso_SF_RooWorkspace",
    call='scalefactor::muon::iso_rooworkspace({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
    input=[q.pt_2, q.eta_2, q.iso_2],
    output=[q.isoWeight_old_2],
    scopes=["em", "mm"],
)
Muon_2_ID_SF = Producer(
    name="MuonID_SF",
    call='scalefactor::muon::id({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_id_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.idWeight_2],
    scopes=["em", "mm"],
)
Muon_2_Iso_SF = Producer(
    name="MuonIso_SF",
    call='scalefactor::muon::iso({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_iso_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.isoWeight_2],
    scopes=["em", "mm"],
)
MuonIDIso_SF = ProducerGroup(
    name="MuonIDIso_SF",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "em", "mm"],
    subproducers={
        "mt": [
            Muon_1_ID_SF,
            Muon_1_Iso_SF,
        ],
        "em": [
            Muon_2_ID_SF,
            Muon_2_Iso_SF,
        ],
        "mm": [
            Muon_1_ID_SF,
            Muon_1_Iso_SF,
            Muon_2_ID_SF,
            Muon_2_Iso_SF,
        ],
    },
)
MuonIDIso_SF_RooWorkspace = ProducerGroup(
    name="MuonIDIso_SF_RooWorkspace",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "em", "mm"],
    subproducers={
        "mt": [
            Muon_1_ID_SF_RooWorkspace,
            Muon_1_Iso_SF_RooWorkspace,
        ],
        "em": [
            Muon_2_ID_SF_RooWorkspace,
            Muon_2_Iso_SF_RooWorkspace,
        ],
        "mm": [
            Muon_1_ID_SF_RooWorkspace,
            Muon_1_Iso_SF_RooWorkspace,
            Muon_2_ID_SF_RooWorkspace,
            Muon_2_Iso_SF_RooWorkspace,
        ],
    },
)

############################
# Tau ID/ISO SF
# The readout is done via correctionlib
############################
Tau_1_VsJetTauID_SF = ExtendedVectorProducer(
    name="VsJetTauID_SF",
    call='scalefactor::tau::id_vsJet({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsjet_tau_id_WP}", "{tau_sf_variation}", "{tau_vsjet_sf_dependence}", {output}, "{tau_sf_file}", "{tau_id_discriminator}")',
    input=[q.pt_1, q.decaymode_1, q.gen_match_1],
    output="tau_1_vsjet_sf_outputname",
    scope=["tt"],
    vec_config="vsjet_tau_id",
)
Tau_1_VsEleTauID_SF = ExtendedVectorProducer(
    name="VsEleTauID_SF",
    call='scalefactor::tau::id_vsEleMu({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsele_tau_id_WP}", "{tau_sf_variation}", {output}, "{tau_sf_file}", "{tau_id_discriminator}")',
    input=[q.eta_1, q.decaymode_1, q.gen_match_1],
    output="tau_1_vsele_sf_outputname",
    scope=["tt"],
    vec_config="vsele_tau_id",
)
Tau_1_VsMuTauID_SF = ExtendedVectorProducer(
    name="VsMuTauID_SF",
    call='scalefactor::tau::id_vsEleMu({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsmu_tau_id_WP}", "{tau_sf_variation}", {output}, "{tau_sf_file}", "{tau_id_discriminator}")',
    input=[q.eta_1, q.decaymode_1, q.gen_match_1],
    output="tau_1_vsmu_sf_outputname",
    scope=["tt"],
    vec_config="vsmu_tau_id",
)
Tau_2_VsJetTauID_SF = ExtendedVectorProducer(
    name="VsJetTauID_SF",
    call='scalefactor::tau::id_vsJet({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsjet_tau_id_WP}", "{tau_sf_variation}", "{tau_vsjet_sf_dependence}", {output}, "{tau_sf_file}", "{tau_id_discriminator}")',
    input=[q.pt_2, q.decaymode_2, q.gen_match_2],
    output="tau_2_vsjet_sf_outputname",
    scope=["mt"],
    vec_config="vsjet_tau_id",
)
Tau_2_VsEleTauID_SF = ExtendedVectorProducer(
    name="VsEleTauID_SF",
    call='scalefactor::tau::id_vsEleMu({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsele_tau_id_WP}", "{tau_sf_variation}", {output}, "{tau_sf_file}", "{tau_id_discriminator}")',
    input=[q.eta_2, q.decaymode_2, q.gen_match_2],
    output="tau_2_vsele_sf_outputname",
    scope=["mt"],
    vec_config="vsele_tau_id",
)
Tau_2_VsMuTauID_SF = ExtendedVectorProducer(
    name="VsMuTauID_SF",
    call='scalefactor::tau::id_vsEleMu({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsmu_tau_id_WP}", "{tau_sf_variation}", {output}, "{tau_sf_file}", "{tau_id_discriminator}")',
    input=[q.eta_2, q.decaymode_2, q.gen_match_2],
    output="tau_2_vsmu_sf_outputname",
    scope=["mt"],
    vec_config="vsmu_tau_id",
)
TauID_SF = ProducerGroup(
    name="TauID_SF",
    call=None,
    input=None,
    output=None,
    scopes=["tt", "mt", "et"],
    subproducers={
        "tt": [
            Tau_1_VsJetTauID_SF,
            Tau_1_VsEleTauID_SF,
            Tau_1_VsMuTauID_SF,
            Tau_2_VsJetTauID_SF,
            Tau_2_VsEleTauID_SF,
            Tau_2_VsMuTauID_SF,
        ],
        "mt": [
            Tau_2_VsJetTauID_SF,
            Tau_2_VsEleTauID_SF,
            Tau_2_VsMuTauID_SF,
        ],
        "et": [
            Tau_2_VsJetTauID_SF,
            Tau_2_VsEleTauID_SF,
            Tau_2_VsMuTauID_SF,
        ],
    },
)