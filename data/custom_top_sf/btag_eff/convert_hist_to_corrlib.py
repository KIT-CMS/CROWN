import correctionlib.convert
import json
import gzip


syst = ""

corrs = {}

corrs["2016preVFP"] = {}
# corrs['2016preVFP']['btag_eff_filename'] = '2016preVFP_UL/btag_eff_2016preVFP.root'.replace('.root', syst + '.root')

corrs["2016postVFP"] = {}
# corrs['2016postVFP']['btag_eff_filename'] = '2016postVFP_UL/btag_eff_2016postVFP.root'.replace('.root', syst + '.root')

corrs["2017"] = {}
# corrs['2017']['btag_eff_filename'] = '2017_UL/btag_eff_2017.root'.replace('.root', syst + '.root')

corrs["2018"] = {}
# corrs['2018']['btag_eff_filename'] = '2018_UL/btag_eff_2018.root'.replace('.root', syst + '.root')

print(corrs.keys())
print(corrs)

for y in corrs.keys():
    for f in ["b", "c", "l"]:
        corrs[y]["btag_eff_" + f + "_top"] = correctionlib.convert.from_uproot_THx(
            path=y + "_UL/btageff__top_" + y + ".root:h_eff_" + f,
            axis_names=["eta", "pt"],
        )
        corrs[y]["btag_eff_" + f + "_top"].description = (
            "btag eff (" + f + ") for top in " + y
        )
        corrs[y]["btag_eff_" + f + "_top"].data.flow = "clamp"
        corrs[y]["btag_eff_" + f + "_top"].name = "top_" + f

        corrs[y]["btag_eff_" + f + "_ewk"] = correctionlib.convert.from_uproot_THx(
            path=y + "_UL/btageff__ewk_" + y + ".root:h_eff_" + f,
            axis_names=["eta", "pt"],
        )
        corrs[y]["btag_eff_" + f + "_ewk"].description = (
            "btag eff (" + f + ") for ewk in " + y
        )
        corrs[y]["btag_eff_" + f + "_ewk"].data.flow = "clamp"
        corrs[y]["btag_eff_" + f + "_ewk"].name = "ewk_" + f


file_description = "b tag eff maps for UL"

cset = {}


for y in corrs.keys():
    cset[y] = correctionlib.schemav2.CorrectionSet(
        schema_version=2,
        description="{} ({})".format(file_description, y),
        corrections=[
            corrs[y]["btag_eff_b_top"],
            corrs[y]["btag_eff_c_top"],
            corrs[y]["btag_eff_l_top"],
            corrs[y]["btag_eff_b_ewk"],
            corrs[y]["btag_eff_c_ewk"],
            corrs[y]["btag_eff_l_ewk"],
        ],
    )

    outfile_name = "{}_UL/btag_eff_{}{}.json".format(y, y, syst)
    with open(outfile_name, "w") as fout:
        fout.write(cset[y].json(exclude_unset=True))

    with gzip.open(outfile_name + ".gz", "wt") as fout:
        fout.write(cset[y].json(exclude_unset=True))


test = json.dumps(cset["2018"].json(exclude_unset=True), indent=4)
print(test)
