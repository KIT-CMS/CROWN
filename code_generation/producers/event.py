import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, VectorProducer

####################
# Set of general producers for event quantities
####################

RunLumiEventFilter = VectorProducer(
    name="RunLumiEventFilter",
    call='basefunctions::FilterIntSelection({df}, "{RunLumiEventFilter_Quantities}", std::vector<int>({RunLumiEventFilter_Selections}), "RunLumiEventFilter")',
    input=[],
    output=None,
    scopes=["general"],
    vec_configs=["RunLumiEventFilter_Quantities", "RunLumiEventFilter_Selections"],
)

Lumi = Producer(
    name="Lumi",
    call="quantities::rename<UInt_t>({df}, {input}, {output})",
    input=[nanoAOD.luminosityBlock],
    output=[q.lumi],
    scopes=["mt", "et", "tt", "em"],
)
