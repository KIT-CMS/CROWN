from ..quantities import output as q
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

############################
# Muon RECO, ID, ISO SF
# The readout is done via correctionlib
############################

with defaults(scopes=["mt", "mm"], input=[q.pt_1, q.eta_1]):
    Muon_1_Reco_SF = Producer(
        call="""physicsobject::muon::scalefactor::Reco(
            {df},
            correctionManager,
            {output},
            {input},
            "{muon_sf_file}",
            "{muon_reco_sf_name}",
            "{muon_sf_variation}")
            """,
        output=[q.reco_wgt_mu_1],
    )
    Muon_1_ID_SF = Producer(
        call="""physicsobject::muon::scalefactor::Id(
            {df},
            correctionManager,
            {output},
            {input},
            "{muon_sf_file}",
            "{muon_id_sf_name}",
            "{muon_sf_variation}")
            """,
        output=[q.id_wgt_mu_1],
    )
    Muon_1_Iso_SF = Producer(
        call="""physicsobject::muon::scalefactor::Iso(
            {df},
            correctionManager,
            {output},
            {input},
            "{muon_sf_file}",
            "{muon_iso_sf_name}",
            "{muon_sf_variation}")
            """,
        output=[q.iso_wgt_mu_1],
    )

with defaults(scopes=["mm"], input=[q.pt_2, q.eta_2]):
    Muon_2_Reco_SF = Producer(
        call="""physicsobject::muon::scalefactor::Reco(
            {df},
            correctionManager,
            {output},
            {input},
            "{muon_sf_file}",
            "{muon_reco_sf_name}",
            "{muon_sf_variation}")
            """,
        output=[q.reco_wgt_mu_2],
    )
    Muon_2_ID_SF = Producer(
        call="""physicsobject::muon::scalefactor::Id(
            {df},
            correctionManager,
            {output},
            {input},
            "{muon_sf_file}",
            "{muon_id_sf_name}",
            "{muon_sf_variation}")
            """,
        output=[q.id_wgt_mu_2],
    )
    Muon_2_Iso_SF = Producer(
        call="""physicsobject::muon::scalefactor::Iso(
            {df},
            correctionManager,
            {output},
            {input},
            "{muon_sf_file}",
            "{muon_iso_sf_name}",
            "{muon_sf_variation}")
            """,
        output=[q.iso_wgt_mu_2],
    )

with defaults(call=None, input=None, output=None):
    Muon_SFs = ProducerGroup(
        scopes=["mt", "mm"],
        subproducers={
            "mt": [
                Muon_1_Reco_SF,
                Muon_1_ID_SF,
                Muon_1_Iso_SF,
            ],
            "mm": [
                Muon_1_Reco_SF,
                Muon_1_ID_SF,
                Muon_1_Iso_SF,
                Muon_2_Reco_SF,
                Muon_2_ID_SF,
                Muon_2_Iso_SF,
            ],
        },
    )
