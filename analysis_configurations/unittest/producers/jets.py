from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for selection possible good jets
####################
JetPtCorrection = Producer(
    name="JetPtCorrection",
    call="physicsobject::jet::JetPtCorrection({df}, {output}, {input}, {jet_reapplyJES}, {jet_jes_sources}, {jet_jes_shift}, {jet_jer_shift}, {jet_jec_file}, {jet_jer_tag}, {jet_jes_tag}, {jet_jec_algo})",
    input=[
        nanoAOD.Jet_pt,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        nanoAOD.Jet_area,
        nanoAOD.Jet_rawFactor,
        nanoAOD.Jet_ID,
        nanoAOD.GenJet_pt,
        nanoAOD.GenJet_eta,
        nanoAOD.GenJet_phi,
        nanoAOD.rho,
    ],
    output=[q.Jet_pt_corrected],
    scopes=["global"],
)
JetMassCorrection = Producer(
    name="JetMassCorrection",
    call="physicsobject::ObjectMassCorrectionWithPt({df}, {output}, {input})",
    input=[
        nanoAOD.Jet_mass,
        nanoAOD.Jet_pt,
        q.Jet_pt_corrected,
    ],
    output=[q.Jet_mass_corrected],
    scopes=["global"],
)
# in data and embdedded sample, we simply rename the nanoAOD jets to the jet_pt_corrected column
RenameJetPt = Producer(
    name="RenameJetPt",
    call="basefunctions::rename<ROOT::RVec<float>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_pt],
    output=[q.Jet_pt_corrected],
    scopes=["global"],
)
RenameJetMass = Producer(
    name="RenameJetMass",
    call="basefunctions::rename<ROOT::RVec<float>>({df}, {input}, {output})",
    input=[nanoAOD.Jet_mass],
    output=[q.Jet_mass_corrected],
    scopes=["global"],
)
RenameJetsData = ProducerGroup(
    name="RenameJetsData",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[RenameJetPt, RenameJetMass],
)
JetEnergyCorrection = ProducerGroup(
    name="JetEnergyCorrection",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[JetPtCorrection, JetMassCorrection],
)
JetPtCut = Producer(
    name="JetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_jet_pt})",
    input=[q.Jet_pt_corrected],
    output=[],
    scopes=["global"],
)
BJetPtCut = Producer(
    name="BJetPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_bjet_pt})",
    input=[q.Jet_pt_corrected],
    output=[],
    scopes=["global"],
)
JetEtaCut = Producer(
    name="JetEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_jet_eta})",
    input=[nanoAOD.Jet_eta],
    output=[],
    scopes=["global"],
)
BJetEtaCut = Producer(
    name="BJetEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_bjet_eta})",
    input=[nanoAOD.Jet_eta],
    output=[],
    scopes=["global"],
)
JetIDCut = Producer(
    name="JetIDCut",
    call="physicsobject::jet::CutID({df}, {output}, {input}, {jet_id})",
    input=[nanoAOD.Jet_ID],
    output=[q.jet_id_mask],
    scopes=["global"],
)
JetPUIDCut = Producer(
    name="JetPUIDCut",
    call="physicsobject::jet::CutPUID({df}, {output}, {input}, {jet_puid}, {jet_puid_max_pt})",
    input=[nanoAOD.Jet_PUID, q.Jet_pt_corrected],
    output=[q.jet_puid_mask],
    scopes=["global"],
)
BTagCut = Producer(
    name="BTagCut",
    call="physicsobject::jet::CutRawID({df}, {input}, {output}, {btag_cut})",
    input=[nanoAOD.BJet_discriminator],
    output=[],
    scopes=["global"],
)
GoodJets = ProducerGroup(
    name="GoodJets",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_jets_mask],
    scopes=["global"],
    subproducers=[JetPtCut, JetEtaCut, JetIDCut, JetPUIDCut],
)
GoodBJets = ProducerGroup(
    name="GoodBJets",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.jet_id_mask, q.jet_puid_mask],
    output=[q.good_bjets_mask],
    scopes=["global"],
    subproducers=[BJetPtCut, BJetEtaCut, BTagCut],
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
    input=[nanoAOD.Jet_eta, nanoAOD.Jet_phi, q.p4_1, q.p4_2],
    output=[q.jet_overlap_veto_mask],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)

GoodJetsWithVeto = ProducerGroup(
    name="GoodJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.good_jets_mask],
    output=[],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[VetoOverlappingJets],
)

GoodBJetsWithVeto = Producer(
    name="GoodBJetsWithVeto",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[q.good_bjets_mask, q.jet_overlap_veto_mask],
    output=[],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)

JetCollection = ProducerGroup(
    name="JetCollection",
    call="jet::OrderJetsByPt({df}, {output}, {input})",
    input=[q.Jet_pt_corrected],
    output=[q.good_jet_collection],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[GoodJetsWithVeto],
)

BJetCollection = ProducerGroup(
    name="BJetCollection",
    call="jet::OrderJetsByPt({df}, {output}, {input})",
    input=[q.Jet_pt_corrected],
    output=[q.good_bjet_collection],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[GoodBJetsWithVeto],
)

##########################
# Basic Jet Quantities
# njets, pt, eta, phi, b-tag value
##########################

LVJet1 = Producer(
    name="LVJet1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.good_jet_collection,
        q.Jet_pt_corrected,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        q.Jet_mass_corrected,
    ],
    output=[q.jet_p4_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
LVJet2 = Producer(
    name="LVJet2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.good_jet_collection,
        q.Jet_pt_corrected,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        q.Jet_mass_corrected,
    ],
    output=[q.jet_p4_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
NumberOfJets = Producer(
    name="NumberOfJets",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_jet_collection],
    output=[q.njets],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jpt_1 = Producer(
    name="jpt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jpt_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jpt_2 = Producer(
    name="jpt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jpt_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jeta_1 = Producer(
    name="jeta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jeta_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jeta_2 = Producer(
    name="jeta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jeta_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jphi_1 = Producer(
    name="jphi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_1],
    output=[q.jphi_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jphi_2 = Producer(
    name="jphi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.jet_p4_2],
    output=[q.jphi_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jtag_value_1 = Producer(
    name="jtag_value_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.BJet_discriminator, q.good_jet_collection],
    output=[q.jtag_value_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
jtag_value_2 = Producer(
    name="jtag_value_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.BJet_discriminator, q.good_jet_collection],
    output=[q.jtag_value_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
mjj = Producer(
    name="jphi_2",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.jet_p4_1, q.jet_p4_2],
    output=[q.mjj],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
BasicJetQuantities = ProducerGroup(
    name="BasicJetQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[
        LVJet1,
        LVJet2,
        NumberOfJets,
        jpt_1,
        jeta_1,
        jphi_1,
        jtag_value_1,
        jpt_2,
        jeta_2,
        jphi_2,
        jtag_value_2,
        mjj,
    ],
)

##########################
# Basic b-Jet Quantities
# nbtag, pt, eta, phi, b-tag value
##########################

LVBJet1 = Producer(
    name="LVBJet1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.good_bjet_collection,
        q.Jet_pt_corrected,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        q.Jet_mass_corrected,
    ],
    output=[q.bjet_p4_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
LVBJet2 = Producer(
    name="LVBJet2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.good_bjet_collection,
        q.Jet_pt_corrected,
        nanoAOD.Jet_eta,
        nanoAOD.Jet_phi,
        q.Jet_mass_corrected,
    ],
    output=[q.bjet_p4_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
NumberOfBJets = Producer(
    name="NumberOfBJets",
    call="quantities::jet::NumberOfJets({df}, {output}, {input})",
    input=[q.good_bjet_collection],
    output=[q.nbtag],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bpt_1 = Producer(
    name="bpt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.bjet_p4_1],
    output=[q.bpt_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bpt_2 = Producer(
    name="bpt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.bjet_p4_2],
    output=[q.bpt_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
beta_1 = Producer(
    name="beta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.bjet_p4_1],
    output=[q.beta_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
beta_2 = Producer(
    name="beta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.bjet_p4_2],
    output=[q.beta_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bphi_1 = Producer(
    name="bphi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.bjet_p4_1],
    output=[q.bphi_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
bphi_2 = Producer(
    name="bphi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.bjet_p4_2],
    output=[q.bphi_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
btag_value_1 = Producer(
    name="btag_value_1",
    call="quantities::jet::btagValue({df}, {output}, {input}, 0)",
    input=[nanoAOD.BJet_discriminator, q.good_bjet_collection],
    output=[q.btag_value_1],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
btag_value_2 = Producer(
    name="btag_value_2",
    call="quantities::jet::btagValue({df}, {output}, {input}, 1)",
    input=[nanoAOD.BJet_discriminator, q.good_bjet_collection],
    output=[q.btag_value_2],
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
)
BasicBJetQuantities = ProducerGroup(
    name="BasicBJetQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt", "em", "mm", "ee"],
    subproducers=[
        LVBJet1,
        LVBJet2,
        NumberOfBJets,
        bpt_1,
        beta_1,
        bphi_1,
        btag_value_1,
        bpt_2,
        beta_2,
        bphi_2,
        btag_value_2,
    ],
)
