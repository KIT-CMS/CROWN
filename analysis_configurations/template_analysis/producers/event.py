from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import BaseFilter, Producer, ProducerGroup, VectorProducer

####################
# Set of general producers for event quantities
####################

RunLumiEventFilter = VectorProducer(
    name="RunLumiEventFilter",
    call='basefunctions::FilterIntSelection<{RunLumiEventFilter_Quantity_Types}>({df}, "{RunLumiEventFilter_Quantities}", {vec_open}{RunLumiEventFilter_Selections}{vec_close}, "RunLumiEventFilter")',
    input=[],
    output=None,
    scopes=["global"],
    vec_configs=[
        "RunLumiEventFilter_Quantities",
        "RunLumiEventFilter_Quantity_Types",
        "RunLumiEventFilter_Selections",
    ],
)

JSONFilter = BaseFilter(
    name="JSONFilter",
    call='basefunctions::JSONFilter({df}, "{golden_json_file}", {input}, "GoldenJSONFilter")',
    input=[nanoAOD.run, nanoAOD.luminosityBlock],
    scopes=["global"],
)

PrefireWeight = Producer(
    name="PrefireWeight",
    call="basefunctions::rename<Float_t>({df}, {input}, {output})",
    input=[nanoAOD.prefireWeight],
    output=[q.prefireweight],
    scopes=["global"],
)

is_data = Producer(
    name="isData",
    input=[],
    call="basefunctions::DefineQuantity({df}, {output}, {is_data})",
    output=[q.is_data],
    scopes=["global"],
)

is_embedding = Producer(
    name="is_embedding",
    call="basefunctions::DefineQuantity({df}, {output}, {is_embedding})",
    input=[],
    output=[q.is_embedding],
    scopes=["global"],
)
is_ttbar = Producer(
    name="is_ttbar",
    call="basefunctions::DefineQuantity({df}, {output}, {is_ttbar})",
    input=[],
    output=[q.is_ttbar],
    scopes=["global"],
)
is_dyjets = Producer(
    name="is_dyjets",
    call="basefunctions::DefineQuantity({df}, {output}, {is_dyjets})",
    input=[],
    output=[q.is_dyjets],
    scopes=["global"],
)
is_wjets = Producer(
    name="is_wjets",
    call="basefunctions::DefineQuantity({df}, {output}, {is_wjets})",
    input=[],
    output=[q.is_wjets],
    scopes=["global"],
)
is_diboson = Producer(
    name="is_diboson",
    call="basefunctions::DefineQuantity({df}, {output}, {is_diboson})",
    input=[],
    output=[q.is_diboson],
    scopes=["global"],
)

SampleFlags = ProducerGroup(
    name="SampleFlags",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[
        is_data,
        is_embedding,
        is_ttbar,
        is_dyjets,
        is_wjets,
        is_diboson,
    ],
)

MetFilter = VectorProducer(
    name="MetFilter",
    call='metfilter::ApplyMetFilter({df}, "{met_filters}", "{met_filters}")',
    input=[],
    output=None,
    scopes=["global"],
    vec_configs=["met_filters"],
)

Lumi = Producer(
    name="Lumi",
    call="basefunctions::rename<UInt_t>({df}, {input}, {output})",
    input=[nanoAOD.luminosityBlock],
    output=[q.lumi],
    scopes=["global"],
)

PUweights = Producer(
    name="PUweights",
    call='reweighting::puweights({df}, {output}, {input}, "{PU_reweighting_file}", "{PU_reweighting_hist}")',
    input=[nanoAOD.Pileup_nTrueInt],
    output=[q.puweight],
    scopes=["global"],
)
