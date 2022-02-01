import code_generation.quantities.output as q
from code_generation.producer import Producer, ProducerGroup

############################
# Muon ID, ISO SF
# The readout is done via RooWorkspaces or correctionlib (with *_UL)
############################

Muon_1_ID_SF = Producer(
    name="MuonID_SF",
    call='scalefactor::muon::id({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_1, q.eta_1],
    output=[q.idWeight_1],
    scopes=["mt", "mm"],
)
Muon_1_Iso_SF = Producer(
    name="MuonIso_SF",
    call='scalefactor::muon::iso({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
    input=[q.pt_1, q.eta_1, q.iso_1],
    output=[q.isoWeight_1],
    scopes=["mt", "mm"],
)
Muon_1_ID_SF_UL = Producer(
    name="MuonID_SF_UL",
    call='scalefactor::muon::id_ul({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_id_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.idWeightUL_1],
    scopes=["mt", "mm"],
)
Muon_1_Iso_SF_UL = Producer(
    name="MuonIso_SF_UL",
    call='scalefactor::muon::iso_ul({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_iso_sf_name}")',
    input=[q.pt_1, q.eta_1],
    output=[q.isoWeightUL_1],
    scopes=["mt", "mm"],
)
Muon_2_ID_SF = Producer(
    name="MuonID_SF",
    call='scalefactor::muon::id({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_2, q.eta_2],
    output=[q.idWeight_2],
    scopes=["em", "mm"],
)
Muon_2_Iso_SF = Producer(
    name="MuonIso_SF",
    call='scalefactor::muon::iso({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
    input=[q.pt_2, q.eta_2, q.iso_2],
    output=[q.isoWeight_2],
    scopes=["em", "mm"],
)
Muon_2_ID_SF_UL = Producer(
    name="MuonID_SF_UL",
    call='scalefactor::muon::id_ul({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_id_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.idWeightUL_2],
    scopes=["em", "mm"],
)
Muon_2_Iso_SF_UL = Producer(
    name="MuonIso_SF_UL",
    call='scalefactor::muon::iso_ul({df}, {input}, "{muon_sf_year_id}", "{muon_sf_varation}", {output}, "{muon_sf_file}", "{muon_iso_sf_name}")',
    input=[q.pt_2, q.eta_2],
    output=[q.isoWeightUL_2],
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
            #Muon_1_ID_SF,
            #Muon_1_Iso_SF,
            Muon_1_ID_SF_UL,
            Muon_1_Iso_SF_UL,
        ],
        "em": [
            Muon_2_ID_SF,
            Muon_2_Iso_SF,
            Muon_2_ID_SF_UL,
            Muon_2_Iso_SF_UL,
        ],
        "mm": [
            Muon_1_ID_SF,
            Muon_1_Iso_SF,
            Muon_2_ID_SF,
            Muon_2_Iso_SF,
            Muon_1_ID_SF_UL,
            Muon_1_Iso_SF_UL,
            Muon_2_ID_SF_UL,
            Muon_2_Iso_SF_UL,
        ],
    },
)

############################
# Tau ID/ISO SF
# The readout is done via correctionlib
############################

Tau_2_VsJetTauID_SF = Producer(
    name="VsJetTauID_SF",
    call='scalefactor::tau::id({df}, {input}, "{vsJet_WP}", "{tau_sf_variation}", "{tau_sf_dependence}", {output}, "{tau_sf_file}", "DeepTau2017v2p1VSjet", {vec_open}{tau_dms}{vec_close})',
    input=[q.pt_2, q.decaymode_2, q.gen_match_2],
    output=[q.vsJetTauIDWeight_2],
    scopes=["mt"],
)