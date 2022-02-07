import code_generation.quantities.output as q
from code_generation.producer import Producer, ProducerGroup
from code_generation.producer import ExtendedVectorProducer


############################
# Muon ID, ISO SF
# The readout is done via RooWorkspaces or correctionlib (with *_UL)
############################

Muon_1_ID_SF_old = Producer(
    name="MuonID_SF_old",
    call='scalefactor::muon::id_old({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_1, q.eta_1],
    output=[q.idWeight_old_1],
    scopes=["mt", "mm"],
)
Muon_1_Iso_SF_old = Producer(
    name="MuonIso_SF_old",
    call='scalefactor::muon::iso_old({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
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
Muon_2_ID_SF_old = Producer(
    name="MuonID_SF_old",
    call='scalefactor::muon::id_old({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_2, q.eta_2],
    output=[q.idWeight_old_2],
    scopes=["em", "mm"],
)
Muon_2_Iso_SF_old = Producer(
    name="MuonIso_SF_old",
    call='scalefactor::muon::iso_old({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
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
MuonIDIso_SF_old = ProducerGroup(
    name="MuonIDIso_SF_old",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "em", "mm"],
    subproducers={
        "mt": [
            Muon_1_ID_SF_old,
            Muon_1_Iso_SF_old,
        ],
        "em": [
            Muon_2_ID_SF_old,
            Muon_2_Iso_SF_old,
        ],
        "mm": [
            Muon_1_ID_SF_old,
            Muon_1_Iso_SF_old,
            Muon_2_ID_SF_old,
            Muon_2_Iso_SF_old,
        ],
    },
)

############################
# Tau ID/ISO SF
# The readout is done via correctionlib
############################
Tau_1_VsJetTauID_SF = ExtendedVectorProducer(
    name="VsJetTauID_SF",
    call='scalefactor::tau::id({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsjet_tau_id_WP}", "{tau_sf_variation}", "{tau_sf_dependence}", {output}, "{tau_sf_file}", "DeepTau2017v2p1VSjet")',
    input=[q.pt_1, q.eta_1, q.decaymode_1, q.gen_match_1],
    output="tau_1_vsjet_sf_outputname",
    scope=["mt"],
    vec_config="vsjet_tau_id",
)
Tau_1_VsEleTauID_SF = ExtendedVectorProducer(
    name="VsEleTauID_SF",
    call='scalefactor::tau::id({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsele_tau_id_WP}", "{tau_sf_variation}", "{tau_sf_dependence}", {output}, "{tau_sf_file}", "DeepTau2017v2p1VSe")',
    input=[q.pt_1, q.eta_1, q.decaymode_1, q.gen_match_1],
    output="tau_1_vsele_sf_outputname",
    scope=["mt"],
    vec_config="vsele_tau_id",
)
Tau_1_VsMuTauID_SF = ExtendedVectorProducer(
    name="VsMuTauID_SF",
    call='scalefactor::tau::id({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsmu_tau_id_WP}", "{tau_sf_variation}", "{tau_sf_dependence}", {output}, "{tau_sf_file}", "DeepTau2017v2p1VSmu")',
    input=[q.pt_1, q.eta_1, q.decaymode_1, q.gen_match_1],
    output="tau_1_vsmu_sf_outputname",
    scope=["mt"],
    vec_config="vsmu_tau_id",
)
Tau_2_VsJetTauID_SF = ExtendedVectorProducer(
    name="VsJetTauID_SF",
    call='scalefactor::tau::id({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsjet_tau_id_WP}", "{tau_sf_variation}", "{tau_sf_dependence}", {output}, "{tau_sf_file}", "DeepTau2017v2p1VSjet")',
    input=[q.pt_2, q.eta_2, q.decaymode_2, q.gen_match_2],
    output="tau_2_vsjet_sf_outputname",
    scope=["mt"],
    vec_config="vsjet_tau_id",
)
Tau_2_VsEleTauID_SF = ExtendedVectorProducer(
    name="VsEleTauID_SF",
    call='scalefactor::tau::id({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsele_tau_id_WP}", "{tau_sf_variation}", "{tau_sf_dependence}", {output}, "{tau_sf_file}", "DeepTau2017v2p1VSe")',
    input=[q.pt_2, q.eta_2, q.decaymode_2, q.gen_match_2],
    output="tau_2_vsele_sf_outputname",
    scope=["mt"],
    vec_config="vsele_tau_id",
)
Tau_2_VsMuTauID_SF = ExtendedVectorProducer(
    name="VsMuTauID_SF",
    call='scalefactor::tau::id({df}, {input}, {vec_open}{tau_dms}{vec_close}, "{vsmu_tau_id_WP}", "{tau_sf_variation}", "{tau_sf_dependence}", {output}, "{tau_sf_file}", "DeepTau2017v2p1VSmu")',
    input=[q.pt_2, q.eta_2, q.decaymode_2, q.gen_match_2],
    output="tau_2_vsmu_sf_outputname",
    scope=["mt"],
    vec_config="vsmu_tau_id",
)
TauID_SF = ProducerGroup(
    name="TauID_SF",
    call=None,
    input=None,
    output=None,
    scopes=["tt", "mt"],
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
    },
)