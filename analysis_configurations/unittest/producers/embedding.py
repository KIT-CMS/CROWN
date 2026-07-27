from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults
from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD

with defaults(scopes=["global"]):
    with defaults(call="event::quantity::Rename<Float_t>({df}, {output}, {input})"):
        EmbeddingGenWeight = Producer(
            input=[nanoAOD.genWeight],
            output=[q.emb_genweight],
        )

        TauEmbeddingInitialMETEt = Producer(
            input=[nanoAOD.TauEmbedding_initialMETEt],
            output=[q.emb_initialMETEt],
        )
        TauEmbeddingInitialMETphi = Producer(
            input=[nanoAOD.TauEmbedding_initialMETphi],
            output=[q.emb_initialMETphi],
        )
        TauEmbeddingInitialPuppiMETEt = Producer(
            input=[nanoAOD.TauEmbedding_initialPuppiMETEt],
            output=[q.emb_initialPuppiMETEt],
        )
        TauEmbeddingInitialPuppiMETphi = Producer(
            input=[nanoAOD.TauEmbedding_initialPuppiMETphi],
            output=[q.emb_initialPuppiMETphi],
        )
        TauEmbeddingnInitialPairCandidates = Producer(
            input=[nanoAOD.TauEmbedding_nInitialPairCandidates],
            output=[q.emb_InitialPairCandidates],
        )
        TauEmbeddingSelectionOldMass = Producer(
            input=[nanoAOD.TauEmbedding_SelectionOldMass],
            output=[q.emb_SelectionOldMass],
        )
        TauEmbeddingSelectionNewMass = Producer(
            input=[nanoAOD.TauEmbedding_SelectionNewMass],
            output=[q.emb_SelectionNewMass],
        )
    with defaults(call="event::quantity::Rename<Bool_t>({df}, {output}, {input})"):
        TauEmbeddingIsMediumLeadingMuon = Producer(
            input=[nanoAOD.TauEmbedding_isMediumLeadingMuon],
            output=[q.emb_isMediumLeadingMuon],
        )
        TauEmbeddingIsMediumTrailingMuon = Producer(
            input=[nanoAOD.TauEmbedding_isMediumTrailingMuon],
            output=[q.emb_isMediumTrailingMuon],
        )
        TauEmbeddingIsTightLeadingMuon = Producer(
            input=[nanoAOD.TauEmbedding_isTightLeadingMuon],
            output=[q.emb_isTightLeadingMuon],
        )
        TauEmbeddingIsTightTrailingMuon = Producer(
            input=[nanoAOD.TauEmbedding_isTightTrailingMuon],
            output=[q.emb_isTightTrailingMuon],
        )

    with defaults(call=None, input=None, output=None):
        EmbeddingQuantities = ProducerGroup(
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
