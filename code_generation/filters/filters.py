import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import VectorProducer

MetFilter = VectorProducer(
    name="MetFilter",
    call='metfilter::ApplyMetFilter({df}, "{met_filters}", "{met_filters}")',
    input=[],
    output=None,
    scopes=["global"],
    vec_configs=["met_filters"],
)

RequireObjects = VectorProducer(
    name="RequireObjects",
    call='physicsobject::FilterObjects({df}, "{require_candidate}", {require_candidate_number}, "{require_candidate}")',
    input=[],
    output=[],
    scopes=["global"],
    vec_configs=["require_candidate", "require_candidate_number"],
)
