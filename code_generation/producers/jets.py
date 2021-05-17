import code_generation.quantities as q
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for selection possible good jets
####################
JetPtCut = Producer(
    name="JetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_jet_pt})",
    input=[q.Jet_pt],
    output=[],
    scopes=["global"],
)
BJetPtCut = Producer(
    name="BJetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_bjet_pt})",
    input=[q.Jet_pt],
    output=[],
    scopes=["global"],
)
JetEtaCut = Producer(
    name="JetEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_jet_eta})",
    input=[q.Jet_eta],
    output=[],
    scopes=["global"],
)
BJetEtaCut = Producer(
    name="BJetEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_bjet_eta})",
    input=[q.Jet_eta],
    output=[],
    scopes=["global"],
)
JetIDFilter = Producer(
    name="JetIDFilter",
    call='physicsobject::jet::FilterID({df}, {output}, "Jet_jetId", {jet_id})',
    input=[],
    output=[q.jet_id_mask],
    scopes=["global"],
)
BTag = Producer(
    name="BTag",
    call='physicsobject::CutPt({df}, "Jet_btagDeepFlavB", {output}, {btag_cut})',
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
    call="jet::VetoOverlappingJets({df}, {output}, {input}, {deltaR_jet_veto})",
    input=[q.Jet_eta, q.Jet_phi, q.p4_1, q.p4_2],
    output=[q.jet_overlap_veto_mask],
    scopes=["mt"],
)

GoodJetsWithVeto = ProducerGroup(
    name="GoodJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[],
    scopes=["mt"],
    subproducers=[JetPtCut, JetEtaCut, JetIDFilter, VetoOverlappingJets],
)

GoodBJetsWithVeto = ProducerGroup(
    name="GoodBJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.jet_id_mask, q.jet_overlap_veto_mask],
    output=[],
    scopes=["mt"],
    subproducers=[BJetPtCut, BJetEtaCut, BTag],
)

JetCollection = ProducerGroup(
    name="JetCollection",
    call="jet::OrderJetsByPt({df}, {output}, {input})",
    input=[q.Jet_pt],
    output=[q.good_jet_collection],
    scopes=["mt"],
    subproducers=[GoodJetsWithVeto],
)

BJetCollection = ProducerGroup(
    name="BJetCollection",
    call="jet::OrderJetsByPt({df}, {output}, {input})",
    input=[q.Jet_pt],
    output=[q.good_bjet_collection],
    scopes=["mt"],
    subproducers=[GoodBJetsWithVeto],
)

##########################
# Basic Jet Quantities
# njets, pt, eta, phi
##########################

LVJet1 = Producer(
    name="LVJet1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[q.good_jet_collection, q.Jet_pt, q.Jet_eta, q.Jet_phi, q.Jet_mass],
    output=[q.jet_p4_1],
    scopes=["mt"],
)
LVJet2 = Producer(
    name="LVJet2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[q.good_jet_collection, q.Jet_pt, q.Jet_eta, q.Jet_phi, q.Jet_mass],
    output=[q.jet_p4_2],
    scopes=["mt"],
)
NumberOfJets = Producer(
    name="NumberOfJets",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_jet_collection],
    output=[q.njets],
    scopes=["mt"],
)

jpt_1 = Producer(
    name="jpt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jpt_1],
    scopes=["mt"],
)
jpt_2 = Producer(
    name="jpt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jpt_2],
    scopes=["mt"],
)
jeta_1 = Producer(
    name="jeta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jeta_1],
    scopes=["mt"],
)
jeta_2 = Producer(
    name="jeta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jeta_2],
    scopes=["mt"],
)
jphi_1 = Producer(
    name="jphi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jphi_1],
    scopes=["mt"],
)
jphi_2 = Producer(
    name="jphi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jphi_2],
    scopes=["mt"],
)
mjj = Producer(
    name="jphi_2",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.jet_p4_1, q.jet_p4_2],
    output=[q.mjj],
    scopes=["mt"],
)
BasicJetQuantities = ProducerGroup(
    name="BasicJetQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[
        LVJet1,
        LVJet2,
        NumberOfJets,
        jpt_1,
        jeta_1,
        jphi_1,
        jpt_2,
        jeta_2,
        jphi_2,
        mjj,
    ],
)
