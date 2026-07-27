from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup, ExtendedVectorProducer
from code_generation.helpers import defaults

####################
# Set of general producers for DiTauPair Quantities
####################

with defaults(scopes=["mt", "et", "tt", "em", "ee", "mm"]):
    with defaults(input=[q.p4_1]):
        pt_1 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.pt_1],
        )
        eta_1 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.eta_1],
        )
        phi_1 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.phi_1],
        )
        mass_1 = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.mass_1],
        )
    with defaults(input=[q.p4_2]):
        pt_2 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.pt_2],
        )
        eta_2 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.eta_2],
        )
        phi_2 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.phi_2],
        )
        mass_2 = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.mass_2],
        )
    with defaults(input=[q.p4_1, q.p4_2]):
        m_vis = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.m_vis],
        )
        pt_vis = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.pt_vis],
        )

####################
# Set of channel specific producers
####################
with defaults(scopes=["mt", "mm"]):
    with defaults(call="event::quantity::Get<float>({df}, {output}, {input}, 0)"):
        muon_dxy_1 = Producer(
            input=[nanoAOD.Muon_dxy, q.dileptonpair],
            output=[q.dxy_1],
        )
        muon_dz_1 = Producer(
            input=[nanoAOD.Muon_dz, q.dileptonpair],
            output=[q.dz_1],
        )
        muon_iso_1 = Producer(
            input=[nanoAOD.Muon_pfRelIso04_all, q.dileptonpair],
            output=[q.iso_1],
        )
    muon_q_1 = Producer(
        call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
        input=[nanoAOD.Muon_charge, q.dileptonpair],
        output=[q.q_1],
    )
with defaults(scopes=["em", "mm"]):
    with defaults(call="event::quantity::Get<float>({df}, {output}, {input}, 1)"):
        muon_dxy_2 = Producer(
            input=[nanoAOD.Muon_dxy, q.dileptonpair],
            output=[q.dxy_2],
        )
        muon_dz_2 = Producer(
            input=[nanoAOD.Muon_dz, q.dileptonpair],
            output=[q.dz_2],
        )
        muon_iso_2 = Producer(
            input=[nanoAOD.Muon_pfRelIso04_all, q.dileptonpair],
            output=[q.iso_2],
        )
    muon_q_2 = Producer(
        call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
        input=[nanoAOD.Muon_charge, q.dileptonpair],
        output=[q.q_2],
    )
with defaults(scopes=["et", "ee", "em"]):
    with defaults(call="event::quantity::Get<float>({df}, {output}, {input}, 0)"):
        electron_dxy_1 = Producer(
            input=[nanoAOD.Electron_dxy, q.dileptonpair],
            output=[q.dxy_1],
        )
        electron_dz_1 = Producer(
            input=[nanoAOD.Electron_dz, q.dileptonpair],
            output=[q.dz_1],
        )
        electron_iso_1 = Producer(
            input=[nanoAOD.Electron_pfRelIso03_all, q.dileptonpair],
            output=[q.iso_1],
        )
    electron_q_1 = Producer(
        call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
        input=[nanoAOD.Electron_charge, q.dileptonpair],
        output=[q.q_1],
    )
with defaults(scopes=["ee"]):
    with defaults(call="event::quantity::Get<float>({df}, {output}, {input}, 1)"):
        electron_dxy_2 = Producer(
            input=[nanoAOD.Electron_dxy, q.dileptonpair],
            output=[q.dxy_2],
        )
        electron_dz_2 = Producer(
            input=[nanoAOD.Electron_dz, q.dileptonpair],
            output=[q.dz_2],
        )
        electron_iso_2 = Producer(
            input=[nanoAOD.Electron_pfRelIso03_all, q.dileptonpair],
            output=[q.iso_2],
        )
    electron_q_2 = Producer(
        call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
        input=[nanoAOD.Electron_charge, q.dileptonpair],
        output=[q.q_2],
    )
with defaults(scopes=["tt"]):
    with defaults(call="event::quantity::Get<float>({df}, {output}, {input}, 0)"):
        tau_dxy_1 = Producer(
            input=[nanoAOD.Tau_dxy, q.dileptonpair],
            output=[q.dxy_1],
        )
        tau_dz_1 = Producer(
            input=[nanoAOD.Tau_dz, q.dileptonpair],
            output=[q.dz_1],
        )
        tau_iso_1 = Producer(
            input=[nanoAOD.Tau_rawDeepTau2017v2p1VSjet, q.dileptonpair],
            output=[q.iso_1],
        )
    with defaults(call="event::quantity::Get<int>({df}, {output}, {input}, 0)"):
        tau_q_1 = Producer(
            input=[nanoAOD.Tau_charge, q.dileptonpair],
            output=[q.q_1],
        )
        tau_decaymode_1 = Producer(
            input=[nanoAOD.Tau_decayMode, q.dileptonpair],
            output=[q.tau_decaymode_1],
        )
    taujet_pt_1 = Producer(
        call="event::quantity::GetJetForObject<float>({df}, {output}, {input}, 0)",
        input=[nanoAOD.Jet_pt, nanoAOD.Tau_jetIdx, q.dileptonpair],
        output=[q.taujet_pt_1],
    )
    VsJetTauIDFlag_1 = ExtendedVectorProducer(
        call="physicsobject::tau::quantity::IDFlag_v9({df}, {output}, {input}, 0, {vsjet_tau_id_WPbit})",
        input=[nanoAOD.Tau_idDeepTau2017v2p1VSjet, q.dileptonpair],
        output="tau_1_vsjet_id_outputname",
        vec_config="vsjet_tau_id",
    )
    VsEleTauIDFlag_1 = ExtendedVectorProducer(
        call="physicsobject::tau::quantity::IDFlag_v9({df}, {output}, {input}, 0, {vsele_tau_id_WPbit})",
        input=[nanoAOD.Tau_idDeepTau2017v2p1VSe, q.dileptonpair],
        output="tau_1_vsele_id_outputname",
        vec_config="vsele_tau_id",
    )
    VsMuTauIDFlag_1 = ExtendedVectorProducer(
        call="physicsobject::tau::quantity::IDFlag_v9({df}, {output}, {input}, 0, {vsmu_tau_id_WPbit})",
        input=[nanoAOD.Tau_idDeepTau2017v2p1VSmu, q.dileptonpair],
        output="tau_1_vsmu_id_outputname",
        vec_config="vsmu_tau_id",
    )
with defaults(scopes=["mt", "et", "tt"]):
    with defaults(call="event::quantity::Get<float>({df}, {output}, {input}, 1)"):
        tau_dxy_2 = Producer(
            input=[nanoAOD.Tau_dxy, q.dileptonpair],
            output=[q.dxy_2],
        )
        tau_dz_2 = Producer(
            input=[nanoAOD.Tau_dz, q.dileptonpair],
            output=[q.dz_2],
        )
        tau_iso_2 = Producer(
            input=[nanoAOD.Tau_rawDeepTau2017v2p1VSjet, q.dileptonpair],
            output=[q.iso_2],
        )
    with defaults(call="event::quantity::Get<int>({df}, {output}, {input}, 1)"):
        tau_q_2 = Producer(
            input=[nanoAOD.Tau_charge, q.dileptonpair],
            output=[q.q_2],
        )
        tau_decaymode_2 = Producer(
            input=[nanoAOD.Tau_decayMode, q.dileptonpair],
            output=[q.tau_decaymode_2],
        )
    taujet_pt_2 = Producer(
        call="event::quantity::GetJetForObject<float>({df}, {output}, {input}, 1)",
        input=[nanoAOD.Jet_pt, nanoAOD.Tau_jetIdx, q.dileptonpair],
        output=[q.taujet_pt_2],
    )
    VsJetTauIDFlag_2 = ExtendedVectorProducer(
        call="physicsobject::tau::quantity::IDFlag_v9({df}, {output}, {input}, 1, {vsjet_tau_id_WPbit})",
        input=[nanoAOD.Tau_idDeepTau2017v2p1VSjet, q.dileptonpair],
        output="tau_2_vsjet_id_outputname",
        vec_config="vsjet_tau_id",
    )
    VsEleTauIDFlag_2 = ExtendedVectorProducer(
        call="physicsobject::tau::quantity::IDFlag_v9({df}, {output}, {input}, 1, {vsele_tau_id_WPbit})",
        input=[nanoAOD.Tau_idDeepTau2017v2p1VSe, q.dileptonpair],
        output="tau_2_vsele_id_outputname",
        vec_config="vsele_tau_id",
    )
    VsMuTauIDFlag_2 = ExtendedVectorProducer(
        call="physicsobject::tau::quantity::IDFlag_v9({df}, {output}, {input}, 1, {vsmu_tau_id_WPbit})",
        input=[nanoAOD.Tau_idDeepTau2017v2p1VSmu, q.dileptonpair],
        output="tau_2_vsmu_id_outputname",
        vec_config="vsmu_tau_id",
    )
with defaults(call="event::quantity::Define<int>({df}, {output}, -1)", input=[]):
    tau_decaymode_1_notau = Producer(
        output=[q.tau_decaymode_1],
        scopes=["et", "mt", "em", "ee", "mm"],
    )
    tau_decaymode_2_notau = Producer(
        output=[q.tau_decaymode_2],
        scopes=["em", "ee", "mm"],
    )
with defaults(call=None, input=None, output=None):
    UnrollMuLV1 = ProducerGroup(
        scopes=["mt", "mm"],
        subproducers=[
            pt_1,
            eta_1,
            phi_1,
            mass_1,
            muon_dxy_1,
            muon_dz_1,
            muon_q_1,
            muon_iso_1,
        ],
    )
    UnrollMuLV2 = ProducerGroup(
        scopes=["mm", "em"],
        subproducers=[
            pt_2,
            eta_2,
            phi_2,
            mass_2,
            muon_dxy_2,
            muon_dz_2,
            muon_q_2,
            muon_iso_2,
        ],
    )
    UnrollElLV1 = ProducerGroup(
        scopes=["et", "ee", "em"],
        subproducers=[
            pt_1,
            eta_1,
            phi_1,
            mass_1,
            electron_dxy_1,
            electron_dz_1,
            electron_q_1,
            electron_iso_1,
        ],
    )
    UnrollElLV2 = ProducerGroup(
        scopes=["ee"],
        subproducers=[
            pt_2,
            eta_2,
            phi_2,
            mass_2,
            electron_dxy_2,
            electron_dz_2,
            electron_q_2,
            electron_iso_2,
        ],
    )
    UnrollTauLV1 = ProducerGroup(
        scopes=["tt"],
        subproducers=[
            pt_1,
            eta_1,
            phi_1,
            mass_1,
            tau_dxy_1,
            tau_dz_1,
            tau_q_1,
            tau_iso_1,
            tau_decaymode_1,
            taujet_pt_1,
            VsJetTauIDFlag_1,
            VsEleTauIDFlag_1,
            VsMuTauIDFlag_1,
        ],
    )
    UnrollTauLV2 = ProducerGroup(
        scopes=["et", "mt", "tt"],
        subproducers=[
            pt_2,
            eta_2,
            phi_2,
            mass_2,
            tau_dxy_2,
            tau_dz_2,
            tau_q_2,
            tau_iso_2,
            tau_decaymode_2,
            taujet_pt_2,
            VsJetTauIDFlag_2,
            VsEleTauIDFlag_2,
            VsMuTauIDFlag_2,
        ],
    )

#####################
# Producer Groups
#####################

with defaults(call=None, input=None, output=None):
    MTDiTauPairQuantities = ProducerGroup(
        scopes=["mt"],
        subproducers=[UnrollMuLV1, UnrollTauLV2, tau_decaymode_1_notau, m_vis, pt_vis],
    )
    MMDiTauPairQuantities = ProducerGroup(
        scopes=["mm"],
        subproducers=[
            UnrollMuLV1,
            UnrollMuLV2,
            tau_decaymode_1_notau,
            tau_decaymode_2_notau,
            m_vis,
            pt_vis,
        ],
    )
    ETDiTauPairQuantities = ProducerGroup(
        scopes=["et"],
        subproducers=[UnrollElLV1, UnrollTauLV2, tau_decaymode_1_notau, m_vis, pt_vis],
    )
    TTDiTauPairQuantities = ProducerGroup(
        scopes=["tt"],
        subproducers=[UnrollTauLV1, UnrollTauLV2, m_vis, pt_vis],
    )
    EMDiTauPairQuantities = ProducerGroup(
        scopes=["em"],
        subproducers=[UnrollElLV1, UnrollMuLV2, m_vis, pt_vis],
    )

## advanced event quantities (can be caluculated when ditau pair and met and all jets are determined)
## leptons: q.p4_1, q.p4_2
## met: met_p4_recoilcorrected
## jets: good_jet_collection (if only the leading two are needed: q.jet_p4_1, q.jet_p4_2
## bjets: gen_bjet_collection

with defaults(scopes=["mt", "et", "tt", "em", "ee", "mm"]):
    LV_dilepton_pair = Producer(
        call="lorentzvector::Sum({df}, {output}, {input})",
        input=[q.p4_1, q.p4_2],
        output=[q.p4_dilepton],
    )
    Pzetamissvis = Producer(
        call="quantities::PzetaMissVis({df}, {output}, {input})",
        input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
        output=[q.pzetamissvis],
    )
    with defaults(call="quantities::TransverseMass({df}, {output}, {input})"):
        mTdileptonMET = Producer(
            input=[q.p4_dilepton, q.met_p4_recoilcorrected],
            output=[q.mTdileptonMET],
        )
        mt_1 = Producer(
            input=[q.p4_1, q.met_p4_recoilcorrected],
            output=[q.mt_1],
        )
        mt_2 = Producer(
            input=[q.p4_2, q.met_p4_recoilcorrected],
            output=[q.mt_2],
        )
        mt_tot = Producer(
            input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
            output=[q.mt_tot],
        )
    with defaults(call="lorentzvector::GetPt({df}, {output}, {input})"):
        pt_tt = Producer(
            input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
            output=[q.pt_tt],
        )
        pt_ttjj = Producer(
            input=[q.p4_1, q.p4_2, q.jet_p4_1, q.jet_p4_2, q.met_p4_recoilcorrected],
            output=[q.pt_ttjj],
        )
    with defaults(call=None, input=None, output=None):
        DiTauPairMETQuantities = ProducerGroup(
            subproducers=[
                LV_dilepton_pair,
                Pzetamissvis,
                mTdileptonMET,
                mt_1,
                mt_2,
                pt_tt,
                pt_ttjj,
                mt_tot,
            ],
        )
with defaults(
    input=[
        q.pt_1,
        q.pt_2,
        q.eta_1,
        q.eta_2,
        q.phi_1,
        q.phi_2,
        q.mass_1,
        q.mass_2,
        q.met,
        q.metphi,
        q.metcov00,
        q.metcov01,
        q.metcov11,
        q.tau_decaymode_1,
        q.tau_decaymode_2,
    ],
    output=[q.p4_fastmtt],
):
    p4_fastmtt_mt = Producer(
        call='quantities::FastMtt({df}, {output}, {input}, "mt")',
        scopes=["mt"],
    )
    p4_fastmtt_et = Producer(
        call='quantities::FastMtt({df}, {output}, {input}, "et")',
        scopes=["et"],
    )
    p4_fastmtt_tt = Producer(
        call='quantities::FastMtt({df}, {output}, {input}, "tt")',
        scopes=["tt"],
    )
    p4_fastmtt_em = Producer(
        call='quantities::FastMtt({df}, {output}, {input}, "em")',
        scopes=["em"],
    )

with defaults(scopes=["mt", "et", "tt", "em"]):
    with defaults(input=[q.p4_fastmtt]):
        pt_fastmtt = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.pt_fastmtt],
        )
        eta_fastmtt = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.eta_fastmtt],
        )
        phi_fastmtt = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.phi_fastmtt],
        )
        m_fastmtt = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.m_fastmtt],
        )
    with defaults(call=None, input=None, output=None):
        FastMTTQuantities = ProducerGroup(
            subproducers={
                "mt": [p4_fastmtt_mt, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
                "et": [p4_fastmtt_et, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
                "tt": [p4_fastmtt_tt, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
                "em": [p4_fastmtt_em, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
            },
        )
