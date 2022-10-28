from code_generation.producer import Producer, ProducerGroup, ExtendedVectorProducer
from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD


EmbeddingGenWeight = Producer(
    name="EmbeddingGenWeight",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.genWeight],
    output=[q.emb_genweight],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)

TauEmbeddingInitialMETEt = Producer(
    name="TauEmbeddingInitialMETEt",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_initialMETEt],
    output=[q.emb_initialMETEt],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingInitialMETphi = Producer(
    name="TauEmbeddingInitialMETphi",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_initialMETphi],
    output=[q.emb_initialMETphi],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingInitialPuppiMETEt = Producer(
    name="TauEmbeddingInitialPuppiMETEt",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_initialPuppiMETEt],
    output=[q.emb_initialPuppiMETEt],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingInitialPuppiMETphi = Producer(
    name="TauEmbeddingInitialPuppiMETphi",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_initialPuppiMETphi],
    output=[q.emb_initialPuppiMETphi],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingIsMediumLeadingMuon = Producer(
    name="TauEmbeddingIsMediumLeadingMuon",
    call="basefunctions::rename<Bool_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_isMediumLeadingMuon],
    output=[q.emb_isMediumLeadingMuon],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingIsMediumTrailingMuon = Producer(
    name="TauEmbeddingIsMediumTrailingMuon",
    call="basefunctions::rename<Bool_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_isMediumTrailingMuon],
    output=[q.emb_isMediumTrailingMuon],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingIsTightLeadingMuon = Producer(
    name="TauEmbeddingIsTightLeadingMuon",
    call="basefunctions::rename<Bool_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_isTightLeadingMuon],
    output=[q.emb_isTightLeadingMuon],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingIsTightTrailingMuon = Producer(
    name="TauEmbeddingIsTightTrailingMuon",
    call="basefunctions::rename<Bool_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_isTightTrailingMuon],
    output=[q.emb_isTightTrailingMuon],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingnInitialPairCandidates = Producer(
    name="TauEmbeddingInitialPairCandidates",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_InitialPairCandidates],
    output=[q.emb_InitialPairCandidates],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingSelectionOldMass = Producer(
    name="TauEmbeddingSelectionOldMass",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_SelectionOldMass],
    output=[q.emb_SelectionOldMass],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingSelectionNewMass = Producer(
    name="TauEmbeddingSelectionNewMass",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.TauEmbedding_SelectionNewMass],
    output=[q.emb_SelectionNewMass],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)

EmbeddingQuantities = ProducerGroup(
    name="EmbeddingQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
    subproducers=[
        EmbeddingGenWeight,
        TauEmbeddingInitialMETEt,
        TauEmbeddingInitialMETphi,
        TauEmbeddingInitialPuppiMETEt,
        TauEmbeddingInitialPuppiMETphi,
        TauEmbeddingIsMediumLeadingMuon,
        TauEmbeddingIsMediumTrailingMuon,
        TauEmbeddingIsTightLeadingMuon,
        TauEmbeddingIsTightTrailingMuon,
        TauEmbeddingnInitialPairCandidates,
        TauEmbeddingSelectionOldMass,
        TauEmbeddingSelectionNewMass,
    ],
)

# Selection scalefactor

TauEmbeddingTriggerSelectionSF = Producer(
    name="TauEmbeddingTriggerSelectionSF",
    call='scalefactor::embedding::selection_trigger({df}, {input}, {output}, "{embedding_selection_sf_file}", "{embedding_selection_trigger_sf}")',
    input=[q.gen_pt_1, q.gen_eta_1, q.gen_pt_2, q.gen_eta_2],
    output=[q.emb_triggersel_wgt],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
TauEmbeddingIDSelectionSF_1 = Producer(
    name="TauEmbeddingIDSelectionSF_1",
    call='scalefactor::embedding::selection_id({df}, {input}, {output}, "{embedding_selection_sf_file}", "{embedding_selection_id_sf}")',
    input=[q.gen_pt_1, q.gen_eta_1],
    output=[q.emb_idsel_wgt_1],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)

TauEmbeddingIDSelectionSF_2 = Producer(
    name="TauEmbeddingIDSelectionSF_2",
    call='scalefactor::embedding::selection_id({df}, {input}, {output}, "{embedding_selection_sf_file}", "{embedding_selection_id_sf}")',
    input=[q.gen_pt_2, q.gen_eta_2],
    output=[q.emb_idsel_wgt_2],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)

TauEmbeddingSelectionSF = ProducerGroup(
    name="TauEmbeddingSelectionSF",
    call=None,
    input=None,
    output=None,
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
    subproducers=[
        TauEmbeddingTriggerSelectionSF,
        TauEmbeddingIDSelectionSF_1,
        TauEmbeddingIDSelectionSF_2,
    ],
)

# Muon ID/Iso/Trigger SFS

TauEmbeddingMuonIDSF_1 = Producer(
    name="TauEmbeddingMuonIDSF_1",
    call='scalefactor::embedding::muon_sf({df}, {input}, {output}, "{embedding_muon_sf_file}", "emb", "{embedding_muon_id_sf}")',
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_mu_1],
    scopes=["mt", "mm"],
)

TauEmbeddingMuonIDSF_2 = Producer(
    name="TauEmbeddingMuonIDSF_2",
    call='scalefactor::embedding::muon_sf({df}, {input}, {output}, "{embedding_muon_sf_file}", "emb", "{embedding_muon_id_sf}")',
    input=[q.pt_2, q.eta_2],
    output=[q.id_wgt_mu_2],
    scopes=["mm", "em"],
)

TauEmbeddingMuonIsoSF_1 = Producer(
    name="TauEmbeddingMuonIsoSF_1",
    call='scalefactor::embedding::muon_sf({df}, {input}, {output}, "{embedding_muon_sf_file}", "emb", "{embedding_muon_iso_sf}")',
    input=[q.pt_1, q.eta_1],
    output=[q.iso_wgt_mu_1],
    scopes=["mt", "mm"],
)

TauEmbeddingMuonIsoSF_2 = Producer(
    name="TauEmbeddingMuonIsoSF_2",
    call='scalefactor::embedding::muon_sf({df}, {input}, {output}, "{embedding_muon_sf_file}", "emb", "{embedding_muon_iso_sf}")',
    input=[q.pt_2, q.eta_2],
    output=[q.iso_wgt_mu_2],
    scopes=["mm", "em"],
)
MTGenerateSingleMuonTriggerSF = ExtendedVectorProducer(
    name="MTGenerateSingleMuonTriggerSF",
    call='scalefactor::embedding::muon_sf({df}, {input}, {output}, "{embedding_muon_sf_file}", "emb", "{embedding_trigger_sf}")',
    input=[q.pt_1, q.eta_1],
    output="flagname",
    scope=["mt"],
    vec_config="singlemuon_trigger_sf",
)

# Electron ID/Iso/Trigger SFS

TauEmbeddingElectronIDSF_1 = Producer(
    name="TauEmbeddingElectronIDSF_1",
    call='scalefactor::embedding::electron_sf({df}, {input}, {output}, "{embedding_electron_sf_file}", "emb", "{embedding_electron_id_sf}")',
    input=[q.pt_1, q.eta_1],
    output=[q.id_wgt_mu_1],
    scopes=["et", "ee", "em"],
)

TauEmbeddingElectronIDSF_2 = Producer(
    name="TauEmbeddingElectronIDSF_2",
    call='scalefactor::embedding::electron_sf({df}, {input}, {output}, "{embedding_electron_sf_file}", "emb", "{embedding_electron_id_sf}")',
    input=[q.pt_2, q.eta_2],
    output=[q.id_wgt_mu_2],
    scopes=["ee"],
)

TauEmbeddingElectronIsoSF_1 = Producer(
    name="TauEmbeddingElectronIsoSF_1",
    call='scalefactor::embedding::electron_sf({df}, {input}, {output}, "{embedding_electron_sf_file}", "emb", "{embedding_electron_iso_sf}")',
    input=[q.pt_1, q.eta_1],
    output=[q.iso_wgt_mu_1],
    scopes=["et", "ee", "em"],
)

TauEmbeddingElectronIsoSF_2 = Producer(
    name="TauEmbeddingElectronIsoSF_2",
    call='scalefactor::embedding::electron_sf({df}, {input}, {output}, "{embedding_electron_sf_file}", "emb", "{embedding_electron_iso_sf}")',
    input=[q.pt_2, q.eta_2],
    output=[q.iso_wgt_mu_2],
    scopes=["ee"],
)
ETGenerateSingleElectronTriggerSF = ExtendedVectorProducer(
    name="ETGenerateSingleElectronTriggerSF",
    call='scalefactor::embedding::electron_sf({df}, {input}, {output}, "{embedding_electron_sf_file}", "emb", "{embedding_trigger_sf}")',
    input=[q.pt_1, q.eta_1],
    output="flagname",
    scope=["et"],
    vec_config="singlelectron_trigger_sf",
)