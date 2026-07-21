from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import BaseFilter, Producer, ProducerGroup, VectorProducer
from code_generation.helpers import defaults

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

    with defaults(input=[]):
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
        is_diboson = Producer(
            call="event::quantity::Define<bool>({df}, {output}, {is_diboson})",
            output=[q.is_diboson],
        )

    SampleFlags = ProducerGroup(
        call=None,
        input=None,
        output=None,
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
        call='event::filter::Flag({df}, "{met_filters}", "{met_filters}")',
        input=[],
        output=None,
        vec_configs=["met_filters"],
    )

    Lumi = Producer(
        call="event::quantity::Rename<UInt_t>({df}, {output}, {input})",
        input=[nanoAOD.luminosityBlock],
        output=[q.lumi],
    )
