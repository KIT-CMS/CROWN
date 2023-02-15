import correctionlib.convert
import json
import gzip


syst = ""
# syst = '_combined_syst'

corrs = {}

corrs["2016preVFP_UL"] = {}
corrs["2016preVFP_UL"][
    "trigger_histname"
] = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoVeryTight_abseta_pt"
corrs["2016preVFP_UL"][
    "iso_histname"
] = "NUM_VeryTightRelIso_DEN_TightIDandIPCut_abseta_pt"
corrs["2016preVFP_UL"]["trigger_filename"] = (
    "2016preVFP_UL/" + corrs["2016preVFP_UL"]["trigger_histname"] + ".root"
)
corrs["2016preVFP_UL"]["iso_filename"] = (
    "2016preVFP_UL/" + corrs["2016preVFP_UL"]["iso_histname"] + ".root"
)

corrs["2016postVFP_UL"] = {}
corrs["2016postVFP_UL"][
    "trigger_histname"
] = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoVeryTight_abseta_pt"
corrs["2016postVFP_UL"][
    "iso_histname"
] = "NUM_VeryTightRelIso_DEN_TightIDandIPCut_abseta_pt"
corrs["2016postVFP_UL"]["trigger_filename"] = (
    "2016postVFP_UL/" + corrs["2016postVFP_UL"]["trigger_histname"] + ".root"
)
corrs["2016postVFP_UL"]["iso_filename"] = (
    "2016postVFP_UL/" + corrs["2016postVFP_UL"]["iso_histname"] + ".root"
)

corrs["2017_UL"] = {}
corrs["2017_UL"][
    "trigger_histname"
] = "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoVeryTight_abseta_pt"
corrs["2017_UL"]["iso_histname"] = "NUM_VeryTightRelIso_DEN_TightIDandIPCut_abseta_pt"
corrs["2017_UL"]["trigger_filename"] = (
    "2017_UL/" + corrs["2017_UL"]["trigger_histname"] + ".root"
)
corrs["2017_UL"]["iso_filename"] = (
    "2017_UL/" + corrs["2017_UL"]["iso_histname"] + ".root"
)

corrs["2018_UL"] = {}
corrs["2018_UL"][
    "trigger_histname"
] = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoVeryTight_abseta_pt"
corrs["2018_UL"]["iso_histname"] = "NUM_VeryTightRelIso_DEN_TightIDandIPCut_abseta_pt"
corrs["2018_UL"]["trigger_filename"] = (
    "2018_UL/" + corrs["2018_UL"]["trigger_histname"] + ".root"
)
corrs["2018_UL"]["iso_filename"] = (
    "2018_UL/" + corrs["2018_UL"]["iso_histname"] + ".root"
)


print(corrs.keys())
print(corrs)

for y in corrs.keys():
    corrs[y]["trigger_corr"] = correctionlib.convert.from_uproot_THx(
        path=corrs[y]["trigger_filename"] + ":" + corrs[y]["trigger_histname"] + syst,
        axis_names=["eta", "pt"],
    )
    corrs[y]["trigger_corr"].description = "Muon trigger SFs for " + y
    corrs[y]["trigger_corr"].data.flow = "clamp"

    corrs[y]["iso_corr"] = correctionlib.convert.from_uproot_THx(
        path=corrs[y]["iso_filename"] + ":" + corrs[y]["iso_histname"] + syst,
        axis_names=["eta", "pt"],
    )
    corrs[y]["iso_corr"].description = "Muon iso SFs for " + y
    corrs[y]["iso_corr"].data.flow = "clamp"


file_description = "custom corrections for Muon HLT and VeryTightIso in UL"

cset = {}
for y in corrs.keys():
    cset[y] = correctionlib.schemav2.CorrectionSet(
        schema_version=2,
        description="{} ({})".format(file_description, y),
        corrections=[
            corrs[y]["trigger_corr"],
            corrs[y]["iso_corr"],
        ],
    )

    outfile_name = "{}/trigger_iso_{}{}.json".format(y, y, syst)
    with open(outfile_name, "w") as fout:
        fout.write(cset[y].json(exclude_unset=True))

    with gzip.open(outfile_name + ".gz", "wt") as fout:
        fout.write(cset[y].json(exclude_unset=True))


test = json.dumps(cset["2018_UL"].json(exclude_unset=True), indent=4)
print(test)
