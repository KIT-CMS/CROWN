import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, VectorProducer, ProducerGroup

####################
# Set of producers used for contruction of met related quantities
####################

BuildMetVector = Producer(
    name="BuildMetVector",
    call="lorentzvectors::buildMET({df}, {input}, {output})",
    input=[
        nanoAOD.MET_pt,
        nanoAOD.MET_phi,
    ],
    output=[q.met_p4],
    scopes=["global"],
)

MetPt = Producer(
    name="MetPt",
    call="met::metPt({df}, {input}, {output})",
    input=[q.met_p4],
    output=[q.met],
    scopes=["global"],
)

MetPhi = Producer(
    name="MetPhi",
    call="met::metPhi({df}, {input}, {output})",
    input=[q.met_p4],
    output=[q.metphi],
    scopes=["global"],
)
MetCov00 = Producer(
    name="MetCov00",
    call="quantities::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covXX,
    ],
    output=[q.metcov00],
    scopes=["global"],
)
MetCov01 = Producer(
    name="MetCov01",
    call="quantities::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covXY,
    ],
    output=[q.metcov01],
    scopes=["global"],
)
MetCov10 = Producer(
    name="MetCov10",
    call="quantities::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covXY,
    ],
    output=[q.metcov10],
    scopes=["global"],
)
MetCov11 = Producer(
    name="MetCov11",
    call="quantities::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covYY,
    ],
    output=[q.metcov11],
    scopes=["global"],
)
MetSumEt = Producer(
    name="MetSumEt",
    call="quantities::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_sumEt,
    ],
    output=[q.metSumEt],
    scopes=["global"],
)


Met = ProducerGroup(
    name="Met",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[
        BuildMetVector,
        MetPt,
        MetPhi,
        MetCov00,
        MetCov01,
        MetCov10,
        MetCov11,
        MetSumEt,
    ],
)
