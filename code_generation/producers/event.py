import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, VectorProducer

####################
# Set of general producers for event quantities
####################

RunLumiEventFilter = VectorProducer(
    name="RunLumiEventFilter",
    call='basefunctions::FilterIntSelection<{RunLumiEventFilter_Quantity_Types}>({df}, "{RunLumiEventFilter_Quantities}", {vec_open}{RunLumiEventFilter_Selections}{vec_close}, "RunLumiEventFilter")',
    input=[],
    output=None,
    scopes=["general"],
    vec_configs=[
        "RunLumiEventFilter_Quantities",
        "RunLumiEventFilter_Quantity_Types",
        "RunLumiEventFilter_Selections",
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
    call="quantities::rename<UInt_t>({df}, {input}, {output})",
    input=[nanoAOD.luminosityBlock],
    output=[q.lumi],
    scopes=["mt", "et", "tt", "em"],
)
