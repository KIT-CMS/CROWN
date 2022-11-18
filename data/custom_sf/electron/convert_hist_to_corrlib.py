import correctionlib.convert
import json
import gzip


syst = ''
#syst = '_combined_syst'

corrs = {}

corrs['2016preVFP'] = {}
corrs['2016preVFP']['trigger_histname'] = 'EGamma_SF2D'
corrs['2016preVFP']['trigger_filename'] = '2016preVFP_UL/trig_2016preVFP.root'

corrs['2016postVFP'] = {}
corrs['2016postVFP']['trigger_histname'] = 'EGamma_SF2D'
corrs['2016postVFP']['trigger_filename'] = '2016postVFP_UL/trig_2016postVFP.root'

corrs['2017'] = {}
corrs['2017']['trigger_histname'] = 'EGamma_SF2D'
corrs['2017']['trigger_filename'] = '2017_UL/trig_2017.root'

corrs['2018'] = {}
corrs['2018']['trigger_histname'] = 'EGamma_SF2D'
corrs['2018']['trigger_filename'] = '2018_UL/trig_2018.root'


print(corrs.keys())
print(corrs)

for y in corrs.keys():
    corrs[y]['trigger_corr'] = correctionlib.convert.from_uproot_THx(path = corrs[y]['trigger_filename'] + ':' + corrs[y]['trigger_histname'] + syst, axis_names = ['eta','pt'])
    corrs[y]['trigger_corr'].description = 'Electron trigger SFs for ' + y
    corrs[y]['trigger_corr'].data.flow = "clamp"



file_description = 'custom corrections for Electron HLT in UL'

cset = {}
for y in corrs.keys():
    cset[y] = correctionlib.schemav2.CorrectionSet(
        schema_version=2,
        description='{} ({})'.format(file_description, y),
        corrections=[
            corrs[y]['trigger_corr'],
        ],
    )

    outfile_name = '{}_UL/trigger_{}{}.json'.format(y,y,syst)
    with open(outfile_name, "w") as fout:
        fout.write(cset[y].json(exclude_unset=True))

    with gzip.open(outfile_name + ".gz" , "wt") as fout:
        fout.write(cset[y].json(exclude_unset=True))


test = json.dumps(cset['2018'].json(exclude_unset=True),
                          indent = 4)
print(test)
