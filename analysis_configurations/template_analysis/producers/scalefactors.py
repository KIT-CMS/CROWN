from ..quantities import output as q
from code_generation.producer import Producer, ProducerGroup

############################
# Muon RECO, ID, ISO SF
# The readout is done via correctionlib
############################

Muon_1_Reco_SF = Producer(
    name="Muon_1_Reco_SF",
    call="""scalefactor::muon::reco(
        {df}, 
        correctionManager, 
        {input}, 
        "{muon_sf_varation}", 
        {output}, 
        "{muon_sf_file}", 
        "{muon_reco_sf_name}")
        """,
    input=[q.pt_1, q.eta_1],
    output=[q.reco_wgt_mu_1],
    scopes=["mt", "mm"],
)
Muon_1_ID_SF = Producer(
    name="Muon_1_ID_SF",
    call="""scalefactor::muon::id(
        {df}, 
        correctionManager, 
        {input}, 
        "{muon_sf_varation}", 
        {output}, 
        "{muon_sf_file}", 
        "{muon_id_sf_name}")
        """,
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_mu_1],
    scopes=["mt", "mm"],
)
Muon_1_Iso_SF = Producer(
    name="Muon_1_Iso_SF",
    call="""scalefactor::muon::iso(
        {df}, 
        correctionManager, 
        {input}, 
        "{muon_sf_varation}", 
        {output}, 
        "{muon_sf_file}", 
        "{muon_iso_sf_name}")
        """,
    input=[q.pt_1, q.eta_1],
    output=[q.iso_wgt_mu_1],
    scopes=["mt", "mm"],
)

Muon_2_Reco_SF = Producer(
    name="Muon_2_Reco_SF",
    call="""scalefactor::muon::reco(
        {df}, 
        correctionManager, 
        {input}, 
        "{muon_sf_varation}", 
        {output}, 
        "{muon_sf_file}", 
        "{muon_reco_sf_name}")
        """,
    input=[q.pt_2, q.eta_2],
    output=[q.reco_wgt_mu_2],
    scopes=["mm"],
)
Muon_2_ID_SF = Producer(
    name="Muon_2_ID_SF",
    call="""scalefactor::muon::id(
        {df}, 
        correctionManager, 
        {input}, 
        "{muon_sf_varation}", 
        {output}, 
        "{muon_sf_file}", 
        "{muon_id_sf_name}")
        """,
    input=[q.pt_2, q.eta_2],
    output=[q.id_wgt_mu_2],
    scopes=["mm"],
)
Muon_2_Iso_SF = Producer(
    name="Muon_2_Iso_SF",
    call="""scalefactor::muon::iso(
        {df}, 
        correctionManager, 
        {input}, 
        "{muon_sf_varation}", 
        {output}, 
        "{muon_sf_file}", 
        "{muon_iso_sf_name}")
        """,
    input=[q.pt_2, q.eta_2],
    output=[q.iso_wgt_mu_2],
    scopes=["mm"],
)

Muon_SFs = ProducerGroup(
    name="Muon_SFs",
    call=None,
    input=None,
    output=None,
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
