from code_generation.producer import Producer, ProducerGroup
from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD


EmbeddingGenWeight = Producer(
    name="EmbeddingGenWeight",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.genWeight],
    output=[q.emb_genweight],
    scopes=["global"],
)

TauEmbeddingInitialMETEt = Producer(
    name="TauEmbeddingInitialMETEt",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_initialMETEt],
    output=[q.emb_initialMETEt],
    scopes=["global"],
)
TauEmbeddingInitialMETphi = Producer(
    name="TauEmbeddingInitialMETphi",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_initialMETphi],
    output=[q.emb_initialMETphi],
    scopes=["global"],
)
TauEmbeddingInitialPuppiMETEt = Producer(
    name="TauEmbeddingInitialPuppiMETEt",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_initialPuppiMETEt],
    output=[q.emb_initialPuppiMETEt],
    scopes=["global"],
)
TauEmbeddingInitialPuppiMETphi = Producer(
    name="TauEmbeddingInitialPuppiMETphi",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_initialPuppiMETphi],
    output=[q.emb_initialPuppiMETphi],
    scopes=["global"],
)
TauEmbeddingIsMediumLeadingMuon = Producer(
    name="TauEmbeddingIsMediumLeadingMuon",
    call="event::quantity::Rename<Bool_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_isMediumLeadingMuon],
    output=[q.emb_isMediumLeadingMuon],
    scopes=["global"],
)
TauEmbeddingIsMediumTrailingMuon = Producer(
    name="TauEmbeddingIsMediumTrailingMuon",
    call="event::quantity::Rename<Bool_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_isMediumTrailingMuon],
    output=[q.emb_isMediumTrailingMuon],
    scopes=["global"],
)
TauEmbeddingIsTightLeadingMuon = Producer(
    name="TauEmbeddingIsTightLeadingMuon",
    call="event::quantity::Rename<Bool_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_isTightLeadingMuon],
    output=[q.emb_isTightLeadingMuon],
    scopes=["global"],
)
TauEmbeddingIsTightTrailingMuon = Producer(
    name="TauEmbeddingIsTightTrailingMuon",
    call="event::quantity::Rename<Bool_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_isTightTrailingMuon],
    output=[q.emb_isTightTrailingMuon],
    scopes=["global"],
)
TauEmbeddingnInitialPairCandidates = Producer(
    name="TauEmbeddingInitialPairCandidates",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_InitialPairCandidates],
    output=[q.emb_InitialPairCandidates],
    scopes=["global"],
)
TauEmbeddingSelectionOldMass = Producer(
    name="TauEmbeddingSelectionOldMass",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_SelectionOldMass],
    output=[q.emb_SelectionOldMass],
    scopes=["global"],
)
TauEmbeddingSelectionNewMass = Producer(
    name="TauEmbeddingSelectionNewMass",
    call="event::quantity::Rename<Float_t>({df}, {output}, {input})",
    input=[nanoAOD.TauEmbedding_SelectionNewMass],
    output=[q.emb_SelectionNewMass],
    scopes=["global"],
)

EmbeddingQuantities = ProducerGroup(
    name="EmbeddingQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
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
