import code_generation.quantities.output as q
from code_generation.producer import Producer, ProducerGroup

############################
# Muon ID, ISO SF
# The readout is done via RooWorkspaces
############################

MuonID_SF = Producer(
    name="MuonID_SF",
    call='scalefactor::muon::id({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_id_name}", "{muon_sf_id_args}")',
    input=[q.pt_1, q.eta_1],
    output=[q.idWeight_1],
    scopes=["mt"],
)
MuonIso_SF = Producer(
    name="MuonIso_SF",
    call='scalefactor::muon::iso({df}, {input}, {output}, "{muon_sf_workspace}", "{muon_sf_iso_name}", "{muon_sf_iso_args}")',
    input=[q.pt_1, q.eta_1, q.iso_1],
    output=[q.isoWeight_1],
    scopes=["mt"],
)
MuonIDIso_SF = ProducerGroup(
    name="MuonIDIso_SF",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[
        MuonID_SF,
        MuonIso_SF,
    ],
)
