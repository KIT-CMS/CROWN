from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import BaseFilter, Producer, ProducerGroup, VectorProducer
from code_generation.helpers import defaults
from .electrons import DiElectronVeto
from .muons import DiMuonVeto

####################
# Set of general producers for event quantities
####################

with defaults(scopes=["global"]):
    JSONFilter = BaseFilter(
        call="""event::filter::GoldenJSON(
            {df}, 
            correctionManager, 
            "GoldenJSONFilter", 
            {input}, 
            "{golden_json_file}")
            """,
        input=[nanoAOD.run, nanoAOD.luminosityBlock],
    )
    PrefireWeight = Producer(
        call="event::quantity::Rename<float>({df}, {output}, {input})",
        input=[nanoAOD.L1PreFiringWeight_Nom],
        output=[q.prefireweight],
    )
    with defaults(input=[]):
        with defaults(output=None):
            RunLumiEventFilter = VectorProducer(
                call="""event::filter::Quantity<{RunLumiEventFilter_Quantity_Types}>(
                    {df}, 
                    "RunLumiEventFilter", 
                    "{RunLumiEventFilter_Quantities}", 
                    {vec_open}{RunLumiEventFilter_Selections}{vec_close})
                    """,
                vec_configs=[
                    "RunLumiEventFilter_Quantities",
                    "RunLumiEventFilter_Quantity_Types",
                    "RunLumiEventFilter_Selections",
                ],
            )
            MetFilter = VectorProducer(
                call='event::filter::Flag({df}, "{met_filters}", "{met_filters}")',
                vec_configs=["met_filters"],
            )
        is_data = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_data})",
            output=[q.is_data],
        )

        is_embedding = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_embedding})",
            output=[q.is_embedding],
        )
        is_ttbar = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_ttbar})",
            output=[q.is_ttbar],
        )
        is_dyjets = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_dyjets})",
            output=[q.is_dyjets],
        )
        is_wjets = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_wjets})",
            output=[q.is_wjets],
        )
        is_ggh_htautau = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_ggh_htautau})",
            output=[q.is_ggh_htautau],
        )
        is_vbf_htautau = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_vbf_htautau})",
            output=[q.is_vbf_htautau],
        )
        is_diboson = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_diboson})",
            output=[q.is_diboson],
        )
        DiLeptonVeto = ProducerGroup(
            call='event::CombineFlags({df}, {output}, {input}, "any_of")',
            output=[q.dilepton_veto],
            subproducers=[DiElectronVeto, DiMuonVeto],
        )
    with defaults(call=None, input=None, output=None):
        SampleFlags = ProducerGroup(
            subproducers=[
                is_data,
                is_embedding,
                is_ttbar,
                is_dyjets,
                is_wjets,
                is_ggh_htautau,
                is_vbf_htautau,
                is_diboson,
            ],
        )
    Lumi = Producer(
        call="event::quantity::Rename<UInt_t>({df}, {output}, {input})",
        input=[nanoAOD.luminosityBlock],
        output=[q.lumi],
    )
    npartons = Producer(
        call="event::quantity::Rename<UChar_t>({df}, {output}, {input})",
        input=[nanoAOD.LHE_Njets],
        output=[q.npartons],
    )
    PUweights = Producer(
        call="""event::reweighting::Pileup(
            {df}, 
            correctionManager, 
            {output}, 
            {input}, 
            "{PU_reweighting_file}", 
            "{PU_reweighting_era}", 
            "{PU_reweighting_variation}")
            """,
        input=[nanoAOD.Pileup_nTrueInt],
        output=[q.puweight],
    )

with defaults(scopes=["global", "em", "et", "mt", "tt", "mm"]):
    ZPtMassReweighting = Producer(
        call="""event::reweighting::ZPtMass(
            {df}, 
            {output}, 
            {input}, 
            "{zptmass_file}", 
            "{zptmass_functor}", 
            "{zptmass_arguments}")
            """,
        input=[
            q.recoil_genboson_p4,
        ],
        output=[q.ZPtMassReweightWeight],
    )
    TopPtReweighting = Producer(
        call="event::reweighting::TopPt({df}, {output}, {input})",
        input=[
            nanoAOD.GenPart_pdgId,
            nanoAOD.GenPart_statusFlags,
            nanoAOD.GenPart_pt,
        ],
        output=[q.topPtReweightWeight],
    )
    GGH_NNLO_Reweighting = Producer(
        call="""htxs::ggHNNLOWeights(
            {df}, 
            {output}, 
            "{ggHNNLOweightsRootfile}", 
            "{ggH_generator}", 
            {input})
            """,
        input=[nanoAOD.HTXS_Higgs_pt, nanoAOD.HTXS_njets30],
        output=[q.ggh_NNLO_weight],
    )
    GGH_WG1_Uncertainties = Producer(
        call="htxs::ggH_WG1_uncertainties({df}, {output_vec}, {input})",
        input=[
            nanoAOD.HTXS_stage_1_pTjet30,
            nanoAOD.HTXS_Higgs_pt,
            nanoAOD.HTXS_njets30,
        ],  # using non-updated stage1 flag required by the used macro
        output=[
            q.THU_ggH_Mu,
            q.THU_ggH_Res,
            q.THU_ggH_Mig01,
            q.THU_ggH_Mig12,
            q.THU_ggH_VBF2j,
            q.THU_ggH_VBF3j,
            q.THU_ggH_PT60,
            q.THU_ggH_PT120,
            q.THU_ggH_qmtop,
        ],
    )
    QQH_WG1_Uncertainties = Producer(
        call="htxs::qqH_WG1_uncertainties({df}, {output_vec}, {input})",
        input=[
            nanoAOD.HTXS_stage1_1_fine_cat_pTjet30GeV
        ],  # using fine stage1.1 flag required by the used macro
        output=[
            q.THU_qqH_TOT,
            q.THU_qqH_PTH200,
            q.THU_qqH_Mjj60,
            q.THU_qqH_Mjj120,
            q.THU_qqH_Mjj350,
            q.THU_qqH_Mjj700,
            q.THU_qqH_Mjj1000,
            q.THU_qqH_Mjj1500,
            q.THU_qqH_25,
            q.THU_qqH_JET01,
        ],
    )
