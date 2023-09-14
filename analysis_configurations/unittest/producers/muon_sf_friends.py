from ..quantities import output as q
from code_generation.producer import Producer


MuonIDSF_friends_1 = Producer(
    name="MuonIDSF_friends_1",
    call='scalefactor::embedding::muon_sf({df}, {input}, {output}, "{muon_sf_file}", "emb", "{muon_id_sf}")',
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_mu_friend_1],
    scopes=["mt", "mm"],
)

MuonIsoSF_friends_1 = Producer(
    name="MuonIsoSF_friends_1",
    call='scalefactor::embedding::muon_sf({df}, {input}, {output}, "{muon_sf_file}", "emb", "{muon_iso_sf}")',
    input=[q.pt_1, q.eta_1],
    output=[q.iso_wgt_mu_friend_1],
    scopes=["mt", "mm"],
)

Rename_IDSF = Producer(
    name="Rename_IDSF",
    call="basefunctions::rename<double>({df}, {input}, {output})",
    input=[q.id_wgt_mu_friend_1],
    output=[q.id_wgt_mu_friend_1_renamed],
    scopes=["mt", "mm"],
)
