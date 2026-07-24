from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

####################
# Set of producers used for loosest selection of electrons
####################

with defaults(scopes=["global"]):
    with defaults(output=[]):
        ElectronPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_ele_pt})",
            input=[nanoAOD.Electron_pt],
        )
        ElectronEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_ele_eta})",
            input=[nanoAOD.Electron_eta],
        )
        ElectronDxyCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_ele_dxy})",
            input=[nanoAOD.Electron_dxy],
        )
        ElectronDzCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_ele_dz})",
            input=[nanoAOD.Electron_dz],
        )
        ElectronIDCut = Producer(
            call="physicsobject::CutEqual<bool>({df}, {output}, {input}, true)",
            input=[nanoAOD.Electron_mvaFall17V2noIso_WP90],
        )
        ElectronIsoCut = Producer(
            call="physicsobject::CutMax<float>({df}, {output}, {input}, {max_ele_iso})",
            input=[nanoAOD.Electron_pfRelIso03_all],
        )
    BaseElectrons = ProducerGroup(
        call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
        input=[],
        output=[q.base_electrons_mask],
        subproducers=[
            ElectronPtCut,
            ElectronEtaCut,
            ElectronDxyCut,
            ElectronDzCut,
            ElectronIDCut,
            ElectronIsoCut,
        ],
    )

####################
# Set of producers used for more specific selection of electrons in channels
####################

with defaults(scopes=["em", "et", "ee"]):
    with defaults(output=[]):
        GoodElectronPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_electron_pt})",
            input=[nanoAOD.Electron_pt],
        )
        GoodElectronEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_electron_eta})",
            input=[nanoAOD.Electron_eta],
        )
        GoodElectronIsoCut = Producer(
            call="physicsobject::CutMax<float>({df}, {output}, {input}, {electron_iso_cut})",
            input=[nanoAOD.Electron_pfRelIso03_all],
        )
    GoodElectrons = ProducerGroup(
        call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
        input=[q.base_electrons_mask],
        output=[q.good_electrons_mask],
        subproducers=[
            GoodElectronPtCut,
            GoodElectronEtaCut,
            GoodElectronIsoCut,
        ],
    )
    VetoElectrons = Producer(
        call="physicsobject::VetoSingleObject({df}, {output}, {input}, {electron_index_in_pair})",
        input=[q.base_electrons_mask, q.dileptonpair],
        output=[q.veto_electrons_mask],
    )
VetoSecondElectron = Producer(
    call="physicsobject::VetoSingleObject({df}, {output}, {input}, {second_electron_index_in_pair})",
    input=[q.veto_electrons_mask, q.dileptonpair],
    output=[q.veto_electrons_mask_2],
    scopes=["ee"],
)
ExtraElectronsVeto = Producer(
    call="physicsobject::Veto({df}, {output}, {input})",
    input={
        "em": [q.veto_electrons_mask],
        "et": [q.veto_electrons_mask],
        "mt": [q.base_electrons_mask],
        "tt": [q.base_electrons_mask],
        "mm": [q.base_electrons_mask],
        "ee": [q.veto_electrons_mask_2],
    },
    output=[q.electron_veto_flag],
    scopes=["em", "et", "mt", "tt", "mm", "ee"],
)
NumberOfGoodElectrons = Producer(
    call="physicsobject::Count({df}, {output}, {input})",
    input=[q.good_electrons_mask],
    output=[q.nelectrons],
    scopes=["et", "em", "ee"],
)

####################
# Set of producers used for di-electron veto
####################

with defaults(scopes=["global"]):
    with defaults(output=[]):
        DiElectronVetoPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_dielectronveto_pt})",
            input=[nanoAOD.Electron_pt],
        )
        DiElectronVetoIDCut = Producer(
            call="physicsobject::CutMin<int>({df}, {output}, {input}, {dielectronveto_id_wp})",
            input=[nanoAOD.Electron_cutBased],
        )
        DiElectronVetoElectrons = ProducerGroup(
            call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
            input=ElectronEtaCut.output
            + ElectronDxyCut.output
            + ElectronDzCut.output
            + ElectronIsoCut.output,
            subproducers=[
                DiElectronVetoPtCut,
                DiElectronVetoIDCut,
            ],
        )
    DiElectronVeto = ProducerGroup(
        call="physicsobject::LeptonPairVeto({df}, {output}, {input}, {dileptonveto_dR})",
        input=[
            nanoAOD.Electron_pt,
            nanoAOD.Electron_eta,
            nanoAOD.Electron_phi,
            nanoAOD.Electron_mass,
            nanoAOD.Electron_charge,
        ],
        output=[q.dielectron_veto],
        subproducers=[DiElectronVetoElectrons],
    )
