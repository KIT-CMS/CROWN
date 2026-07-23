from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

####################
# Set of producers used for selection possible good jets
####################
with defaults(scopes=["global"]):
    JetPtSmearingSeed = Producer(
        call="event::quantity::GenerateSeed({df}, {output}, {input}, {jet_jer_master_seed})",
        input=[
            nanoAOD.luminosityBlock,
            nanoAOD.run,
            nanoAOD.event,
        ],
        output=[],
    )
    JetPtCorrection = ProducerGroup(
        call="""physicsobject::jet::PtCorrectionMC(
            {df}, 
            correctionManager, 
            {output}, 
            {input}, 
            {jet_jec_file}, 
            {jet_jec_algo}, 
            {jet_jes_tag}, 
            "{jet_jes_source}", 
            {jet_jer_tag}, 
            {jet_reapplyJES}, 
            {jet_jes_shift}, 
            "{jet_jer_shift}",
            "{era}")
            """,
        input=[
            nanoAOD.Jet_pt,
            nanoAOD.Jet_eta,
            nanoAOD.Jet_phi,
            nanoAOD.Jet_area,
            nanoAOD.Jet_rawFactor,
            nanoAOD.Jet_jetId,
            nanoAOD.GenJet_pt,
            nanoAOD.GenJet_eta,
            nanoAOD.GenJet_phi,
            nanoAOD.fixedGridRhoFastjetAll,
        ],
        output=[q.Jet_pt_corrected],
        subproducers=[JetPtSmearingSeed],
    )
    JetMassCorrection = Producer(
        call="physicsobject::MassCorrectionWithPt({df}, {output}, {input})",
        input=[
            nanoAOD.Jet_mass,
            nanoAOD.Jet_pt,
            q.Jet_pt_corrected,
        ],
        output=[q.Jet_mass_corrected],
    )
    # in data and embdedded sample, we simply rename the nanoAOD jets to the jet_pt_corrected column
    with defaults(
        call="event::quantity::Rename<ROOT::RVec<float>>({df}, {output}, {input})"
    ):
        RenameJetPt = Producer(
            input=[nanoAOD.Jet_pt],
            output=[q.Jet_pt_corrected],
        )
        RenameJetMass = Producer(
            input=[nanoAOD.Jet_mass],
            output=[q.Jet_mass_corrected],
        )
    with defaults(call=None, input=None, output=None):
        RenameJetsData = ProducerGroup(
            subproducers=[RenameJetPt, RenameJetMass],
        )
        JetEnergyCorrection = ProducerGroup(
            subproducers=[JetPtCorrection, JetMassCorrection],
        )
    with defaults(output=[]):
        JetPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_jet_pt})",
            input=[q.Jet_pt_corrected],
        )
        BJetPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_bjet_pt})",
            input=[q.Jet_pt_corrected],
        )
        JetEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_jet_eta})",
            input=[nanoAOD.Jet_eta],
        )
        BJetEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_bjet_eta})",
            input=[nanoAOD.Jet_eta],
        )
        BTagCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {btag_cut})",
            input=[nanoAOD.Jet_btagDeepFlavB],
        )
    JetIDCut = Producer(
        call="physicsobject::CutMin<int>({df}, {output}, {input}, {jet_id})",
        input=[nanoAOD.Jet_jetId],
        output=[q.jet_id_mask],
    )
    JetPUIDCut = Producer(
        call="physicsobject::jet::CutPileupID({df}, {output}, {input}, {jet_puid}, {jet_puid_max_pt})",
        input=[nanoAOD.Jet_puId, q.Jet_pt_corrected],
        output=[q.jet_puid_mask],
    )
    with defaults(
        call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")'
    ):
        GoodJets = ProducerGroup(
            input=[],
            output=[q.good_jets_mask],
            subproducers=[JetPtCut, JetEtaCut, JetIDCut, JetPUIDCut],
        )
        GoodBJets = ProducerGroup(
            input=[q.jet_id_mask, q.jet_puid_mask],
            output=[q.good_bjets_mask],
            subproducers=[BJetPtCut, BJetEtaCut, BTagCut],
        )

####################
# Set of producers to apply a veto of jets overlapping with ditaupair candidates and ordering
# jets by their pt
# 1. check all jets vs the two lepton candidates, if they are not within deltaR = 0.5,
#    keep them --> mask
# 2. Combine mask with good_jets_mask
# 3. Generate JetCollection, an RVec containing all indices of good Jets in pt order
# 4. generate jet quantity outputs
####################
with defaults(scopes=["mt", "et", "tt", "em", "mm", "ee"]):
    VetoOverlappingJets = Producer(
        call="physicsobject::jet::VetoOverlappingJets({df}, {output}, {input}, {deltaR_jet_veto})",
        input=[nanoAOD.Jet_eta, nanoAOD.Jet_phi, q.p4_1, q.p4_2],
        output=[q.jet_overlap_veto_mask],
    )
    with defaults(
        call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")', output=[]
    ):
        GoodJetsWithVeto = ProducerGroup(
            input=[q.good_jets_mask],
            subproducers=[VetoOverlappingJets],
        )
        GoodBJetsWithVeto = Producer(
            input=[q.good_bjets_mask, q.jet_overlap_veto_mask],
        )
    with defaults(
        call="physicsobject::OrderByPt({df}, {output}, {input})",
        input=[q.Jet_pt_corrected],
    ):
        JetCollection = ProducerGroup(
            output=[q.good_jet_collection],
            subproducers=[GoodJetsWithVeto],
        )
        BJetCollection = ProducerGroup(
            output=[q.good_bjet_collection],
            subproducers=[GoodBJetsWithVeto],
        )

    ##########################
    # Basic Jet Quantities
    # njets, pt, eta, phi, b-tag value
    ##########################
    with defaults(
        input=[
            q.Jet_pt_corrected,
            nanoAOD.Jet_eta,
            nanoAOD.Jet_phi,
            q.Jet_mass_corrected,
            q.good_jet_collection,
        ]
    ):
        LVJet1 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 0)",
            output=[q.jet_p4_1],
        )
        LVJet2 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 1)",
            output=[q.jet_p4_2],
        )
    NumberOfJets = Producer(
        call="physicsobject::Count({df}, {output}, {input})",
        input=[q.good_jet_collection],
        output=[q.njets],
    )
    with defaults(input=[q.jet_p4_1]):
        jpt_1 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.jpt_1],
        )
        jeta_1 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.jeta_1],
        )
        jphi_1 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.jphi_1],
        )
    with defaults(input=[q.jet_p4_2]):
        jpt_2 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.jpt_2],
        )
        jeta_2 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.jeta_2],
        )
        jphi_2 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.jphi_2],
        )
    with defaults(input=[nanoAOD.Jet_btagDeepFlavB, q.good_jet_collection]):
        jtag_value_1 = Producer(
            call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
            output=[q.jtag_value_1],
        )
        jtag_value_2 = Producer(
            call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
            output=[q.jtag_value_2],
        )
    mjj = Producer(
        call="lorentzvector::GetMass({df}, {output}, {input})",
        input=[q.jet_p4_1, q.jet_p4_2],
        output=[q.mjj],
    )
    with defaults(call=None, input=None, output=None):
        BasicJetQuantities = ProducerGroup(
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
    with defaults(
        input=[
            q.Jet_pt_corrected,
            nanoAOD.Jet_eta,
            nanoAOD.Jet_phi,
            q.Jet_mass_corrected,
            q.good_bjet_collection,
        ]
    ):
        LVBJet1 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 0)",
            output=[q.bjet_p4_1],
        )
        LVBJet2 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 1)",
            output=[q.bjet_p4_2],
        )
    NumberOfBJets = Producer(
        call="physicsobject::Count({df}, {output}, {input})",
        input=[q.good_bjet_collection],
        output=[q.nbtag],
    )
    with defaults(input=[q.bjet_p4_1]):
        bpt_1 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.bpt_1],
        )
        beta_1 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.beta_1],
        )
        bphi_1 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.bphi_1],
        )
    with defaults(input=[q.bjet_p4_2]):
        bpt_2 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.bpt_2],
        )
        beta_2 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.beta_2],
        )
        bphi_2 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.bphi_2],
        )
    with defaults(input=[nanoAOD.Jet_btagDeepFlavB, q.good_bjet_collection]):
        btag_value_1 = Producer(
            call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
            output=[q.btag_value_1],
        )
        btag_value_2 = Producer(
            call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
            output=[q.btag_value_2],
        )
    with defaults(call=None, input=None, output=None):
        BasicBJetQuantities = ProducerGroup(
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
