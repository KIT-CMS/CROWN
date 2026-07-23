from ..quantities import output as q
from code_generation.producer import Producer
from code_generation.helpers import defaults

with defaults(scopes=["mt", "mm"]):
    with defaults(input=[q.pt_1, q.eta_1]):
        MuonIDSF_friends_1 = Producer(
            call="""embedding::muon::Scalefactor(
                {df}, 
                correctionManager, 
                {output},
                {input}, 
                "{muon_sf_file}", 
                "{muon_id_sf}",
                "emb")
                """,
            output=[q.id_wgt_mu_friend_1],
        )
        MuonIsoSF_friends_1 = Producer(
            call="""embedding::muon::Scalefactor(
                {df}, 
                correctionManager, 
                {output},
                {input}, 
                "{muon_sf_file}", 
                "{muon_id_sf}",
                "emb")
                """,
            output=[q.iso_wgt_mu_friend_1],
        )
    Rename_IDSF = Producer(
        call="event::quantity::Rename<double>({df}, {output}, {input})",
        input=[q.id_wgt_mu_friend_1],
        output=[q.id_wgt_mu_friend_1_renamed],
    )
