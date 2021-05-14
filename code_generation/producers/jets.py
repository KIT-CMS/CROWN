import code_generation.quantities as q
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for selection possible good jets
####################
JetPtCut = Producer(
    name="JetPtCut",
    call='physicsobject::CutPt({df}, "Jet_pt", {output}, {min_jet_pt})',
    input=[],
    output=[],
    scopes=["global"],
)
JetEtaCut = Producer(
    name="JetEtaCut",
    call='physicsobject::CutEta({df}, "Jet_eta", {output}, {max_jet_eta})',
    input=[],
    output=[],
    scopes=["global"],
)
JetIDFilter = Producer(
    name="JetIDFilter",
    call='physicsobject::jet::FilterID({df}, {output}, "Jet_jetId", {jet_id})',
    input=[],
    output=[],
    scopes=["global"],
)
GoodJets = ProducerGroup(
    name="GoodJets",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_jets_mask],
    scopes=["global"],
    subproducers=[JetPtCut, JetEtaCut, JetIDFilter],
)

####################
# Set of producers to apply a veto of jets overlapping with ditaupair candidates and ordering jets by their pt
# 1. check all jets vs the two lepton candidates, if they are not within deltaR = 0.5, keep them --> mask
# 2. Combine mask with good_jets_mask
# 3. Generate JetCollection, an RVec containing all indices of good Jets in pt order
# 4. generate jet quantity outputs
####################
VetoOverlappingJets = Producer(
    name="VetoOverlappingJets",
    call="jetveto::VetoOverlappingJets({df}, {output}, {input}",
    input=[q.Jet_eta, q.Jet_phi, q.eta_1, q.eta_2, q.phi_1, q.phi_2],
    output=[],
    scopes=["mt"],
)

GoodJetsWithVeto = ProducerGroup(
    name="GoodJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_jets_with_leptonveto_mask],
    scopes=["mt"],
    subproducers=[GoodJets, VetoOverlappingJets],
)

OrderJetsByPt = Producer(
    name="OrderJetsByPt",
    call="jetveto::OrderJetsByPt({df}, {output}, {input})",
    input=[],
    output=[q.good_jet_collection],
    scopes=["mt"],
)


VetoJets = ProducerGroup(
    name="VetoJets",
    call="",
    input=[],
    output=[],
    scopes=["mt"],
    subproducers=[GoodJetsWithVeto, OrderJetsByPt],
)
