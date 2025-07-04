from ..quantities import output as q
from code_generation.producer import Producer


MuonIDSF_friends_1 = Producer(
    name="MuonIDSF_friends_1",
    call="""embedding::muon::Scalefactor(
        {df}, 
        correctionManager, 
        {output},
        {input}, 
        "{muon_sf_file}", 
        "{muon_id_sf}",
        "emb")
        """,
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_mu_friend_1],
    scopes=["mt", "mm"],
)

MuonIsoSF_friends_1 = Producer(
    name="MuonIsoSF_friends_1",
    call="""embedding::muon::Scalefactor(
        {df}, 
        correctionManager, 
        {output},
        {input}, 
        "{muon_sf_file}", 
        "{muon_id_sf}",
        "emb")
        """,
    input=[q.pt_1, q.eta_1],
    output=[q.iso_wgt_mu_friend_1],
    scopes=["mt", "mm"],
)

Rename_IDSF = Producer(
    name="Rename_IDSF",
    call="event::quantity::Rename<double>({df}, {output}, {input})",
    input=[q.id_wgt_mu_friend_1],
    output=[q.id_wgt_mu_friend_1_renamed],
    scopes=["mt", "mm"],
)
