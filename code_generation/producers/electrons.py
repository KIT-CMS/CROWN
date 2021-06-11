import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for selection of veto electrons
####################

VetoElectronPtCut = Producer(
    name="ElectronPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_VetoElectron_pt})",
    input=[nanoAOD.Electron_pt],
    output=[],
    scopes=["global"],
)
VetoElectronEtaCut = Producer(
    name="ElectronEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_VetoElectron_eta})",
    input=[nanoAOD.Electron_eta],
    output=[],
    scopes=["global"],
)
VetoElectronIDCut = Producer(
    name="ElectronIDCut",
    call='physicsobject::electron::CutID({df}, {output}, "{VetoElectron_id}")',
    input=[],
    output=[],
    scopes=["global"],
)
VetoElectronIsoCut = Producer(
    name="ElectronIsoCut",
    call="physicsobject::electron::CutIsolation({df}, {output}, {input}, {max_VetoElectron_iso})",
    input=[nanoAOD.Electron_iso],
    output=[],
    scopes=["global"],
)
VetoElectrons = ProducerGroup(
    name="VetoElectrons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[],
    scopes=["global"],
    subproducers=[
        VetoElectronPtCut,
        VetoElectronEtaCut,
        VetoElectronIDCut,
        VetoElectronIsoCut,
    ],
)
GoodElectronsVeto = ProducerGroup(
    name="GoodElectronsVeto",
    call="physicsobject::LeptonVetoFlag({df}, {output}, {input})",
    input=[],
    output=[q.electron_veto_flag],
    scopes=["global"],
    subproducers=[VetoElectrons],
)
