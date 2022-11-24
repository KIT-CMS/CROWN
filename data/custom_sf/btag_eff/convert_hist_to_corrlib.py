import correctionlib.convert
import json
import gzip




syst = ''

corrs = {}

corrs['2016preVFP'] = {}
# corrs['2016preVFP']['btag_eff_filename'] = '2016preVFP_UL/btag_eff_2016preVFP.root'.replace('.root', syst + '.root')

corrs['2016postVFP'] = {}
# corrs['2016postVFP']['btag_eff_filename'] = '2016postVFP_UL/btag_eff_2016postVFP.root'.replace('.root', syst + '.root')

corrs['2017'] = {}
# corrs['2017']['btag_eff_filename'] = '2017_UL/btag_eff_2017.root'.replace('.root', syst + '.root')

corrs['2018'] = {}
# corrs['2018']['btag_eff_filename'] = '2018_UL/btag_eff_2018.root'.replace('.root', syst + '.root')











print(corrs.keys())
print(corrs)

for y in corrs.keys():
        for f in ['b', 'c', 'udsg']:
                corrs[y]['btag_eff_' + f + '_ttbar'] = correctionlib.convert.from_uproot_THx(path = y + '_UL/btag_ttbar.root:h_eff_' + f + 'jets_m_deepjet', axis_names = ['eta','pt'])
                corrs[y]['btag_eff_' + f + '_ttbar'].description = 'btag eff (' + f +') for ttbar in ' + y
                corrs[y]['btag_eff_' + f + '_ttbar'].data.flow = "clamp"
                corrs[y]['btag_eff_' + f + '_ttbar'].name = 'ttbar_' + f

                corrs[y]['btag_eff_' + f + '_wjets'] = correctionlib.convert.from_uproot_THx(path = y + '_UL/btag_wjets.root:h_eff_' + f + 'jets_m_deepjet', axis_names = ['eta','pt'])
                corrs[y]['btag_eff_' + f + '_wjets'].description = 'btag eff (' + f +') for wjets in ' + y
                corrs[y]['btag_eff_' + f + '_wjets'].data.flow = "clamp"
                corrs[y]['btag_eff_' + f + '_wjets'].name = 'wjets_' + f



file_description = 'b tag eff maps for UL'

cset = {}


for y in corrs.keys():

    cset[y] = correctionlib.schemav2.CorrectionSet(
        schema_version=2,
        description='{} ({})'.format(file_description, y),
        corrections=[
            corrs[y]['btag_eff_b_ttbar'],
            corrs[y]['btag_eff_c_ttbar'],
            corrs[y]['btag_eff_udsg_ttbar'],
            corrs[y]['btag_eff_b_wjets'],
            corrs[y]['btag_eff_c_wjets'],
            corrs[y]['btag_eff_udsg_wjets'],
        ],
    )

    outfile_name = '{}_UL/btag_eff_{}{}.json'.format(y,y,syst)
    with open(outfile_name, "w") as fout:
        fout.write(cset[y].json(exclude_unset=True))

    with gzip.open(outfile_name + ".gz" , "wt") as fout:
        fout.write(cset[y].json(exclude_unset=True))


test = json.dumps(cset['2018'].json(exclude_unset=True),
                          indent = 4)
print(test)
