import json
import argparse
import copy

parser = argparse.ArgumentParser(description="Add NMSSM entries to datasets.json")
parser.add_argument("--mhigh", type=int, required=True, help="Heavy higgs mass")
parser.add_argument("--mlow", type=int, required=True, help="Light higgs mass")
parser.add_argument("--era", type=int, required=True, help="era")


args = parser.parse_args()
mhigh = args.mhigh
mlow = args.mlow
era = args.era
new_nick  = ""

if era == 2018:
    new_nick = "NMSSMM{}h1M125tautauh2M{}_RunIIAutumn18MiniAOD_102X_13TeV_MINIAOD_madgraph-pythia8_v1".format(mhigh,mlow)
elif era == 2017:
    new_nick = "NMSSMM{}h1M125tautauh2M{}_RunIIFall17MiniAODv2_PU2017_13TeV_MINIAOD_madgraph-pythia8_v1".format(mhigh,mlow)
elif era == 2016:
    new_nick = "NMSSMM{}h1M125tautauh2M{}_RunIISummer16MiniAODv3_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_v1".format(mhigh,mlow)
else:
    print "Select era from 2016, 2017 or 2018"
    exit(1)


with open("datasets.json", "r") as ds_file:
    datasets = json.load(ds_file)

if era == 2018:
    example_dataset = "NMSSMM320h1M125tautauh2M60_RunIIAutumn18MiniAOD_102X_13TeV_MINIAOD_madgraph-pythia8_v1"
elif era == 2017:
    example_dataset = "NMSSMM320h1M125tautauh2M60_RunIIFall17MiniAODv2_PU2017_13TeV_MINIAOD_madgraph-pythia8_v1"
elif era == 2016:
    example_dataset = "NMSSMM320h1M125tautauh2M60_RunIISummer16MiniAODv3_PUMoriond17_13TeV_MINIAOD_madgraph-pythia8_v1"


base_entry = copy.deepcopy(datasets[example_dataset])

with open("miniAOD_lists_{}/M{}_h1_M125_tautau_h2_M{}_bb_miniAOD.txt".format(era,mhigh,mlow),"r") as filelist:
    for line in filelist.readlines():
        if line[:9]=="events = ":
            n_events = int(line[9:])
            break
with open("miniAOD_lists_{}/M{}_h1_M125_tautau_h2_M{}_bb_miniAOD.txt".format(era,mhigh,mlow),"r") as filelist:
    n_files = len(filelist.readlines())-5

datasets[new_nick] = {}
for key in base_entry.keys():
    if key=="process":
        datasets[new_nick][key] = "NMSSMM{}h1M125tautauh2M{}".format(mhigh,mlow)
    elif key=="n_events_generated":
        datasets[new_nick][key] = n_events
    elif key=="n_files":
        datasets[new_nick][key] = n_files
    else:
        datasets[new_nick][key] = base_entry[key]

with open("datasets.json", 'w') as out_json:
    out_json.write(json.dumps(datasets, sort_keys=True, indent=2))
    out_json.close()
