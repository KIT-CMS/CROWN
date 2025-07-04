from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import BaseFilter, Producer, ProducerGroup, VectorProducer

####################
# Set of general producers for event quantities
####################

JSONFilter = BaseFilter(
    name="JSONFilter",
    call='event::filter::GoldenJSON({df}, correctionManager, "GoldenJSONFilter", {input}, "{golden_json_file}")',
    input=[nanoAOD.run, nanoAOD.luminosityBlock],
    scopes=["global"],
)

is_data = Producer(
    name="is_data",
    input=[],
    call="event::quantity::Define({df}, {output}, {is_data})",
    output=[q.is_data],
    scopes=["global"],
)
is_embedding = Producer(
    name="is_embedding",
    call="event::quantity::Define({df}, {output}, {is_embedding})",
    input=[],
    output=[q.is_embedding],
    scopes=["global"],
)
is_ttbar = Producer(
    name="is_ttbar",
    call="event::quantity::Define({df}, {output}, {is_ttbar})",
    input=[],
    output=[q.is_ttbar],
    scopes=["global"],
)
is_dyjets = Producer(
    name="is_dyjets",
    call="event::quantity::Define({df}, {output}, {is_dyjets})",
    input=[],
    output=[q.is_dyjets],
    scopes=["global"],
)
is_wjets = Producer(
    name="is_wjets",
    call="event::quantity::Define({df}, {output}, {is_wjets})",
    input=[],
    output=[q.is_wjets],
    scopes=["global"],
)
is_diboson = Producer(
    name="is_diboson",
    call="event::quantity::Define({df}, {output}, {is_diboson})",
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
    call='event::filter::Flag({df}, "{met_filters}", "{met_filters}")',
    input=[],
    output=None,
    scopes=["global"],
    vec_configs=["met_filters"],
)

Lumi = Producer(
    name="Lumi",
    call="event::quantity::Rename<UInt_t>({df}, {output}, {input})",
    input=[nanoAOD.luminosityBlock],
    output=[q.lumi],
    scopes=["global"],
)
