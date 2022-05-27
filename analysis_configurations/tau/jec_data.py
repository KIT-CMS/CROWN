from code_generation.systematics import SystematicShift
from .producers import jets as jets


def add_jetCorrectionData(configuration, era):
    #########################
    # Jet energy corrections for data
    #########################
    if era == "2018":
        configuration.add_shift(
            SystematicShift(
                name="jec2018A",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL18_RunA_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2018B",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL18_RunB_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2018C",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL18_RunC_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2018D",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL18_RunD_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
    
    if era == "2017":
        configuration.add_shift(
            SystematicShift(
                name="jec2017B",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL17_RunB_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2017C",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL17_RunC_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2017D",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL17_RunD_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2017E",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL17_RunE_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2017F",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL17_RunF_V5_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
    if era == "2016postVFP":
        configuration.add_shift(
            SystematicShift(
                name="jec2016FGHpostVFP",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL16_RunFGH_V7_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
    if era == "2016preVFP":
        configuration.add_shift(
            SystematicShift(
                name="jec2016BCDpreVFP",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL16APV_RunBCD_V7_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
        configuration.add_shift(
            SystematicShift(
                name="jec2016EFpreVFP",
                shift_config={
                    "global": {
                        "jet_jes_tag_data": '"Summer19UL16APV_RunEF_V7_DATA"',
                    },
                },
                producers={"global": jets.JetEnergyCorrection_data},
            ),
            samples=["data"],
        )
