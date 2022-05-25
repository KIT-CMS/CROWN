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
