#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import logging
import numpy as np
logger = logging.getLogger("sum_training_weights")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

### YAML + ORDERED DICT MAGIC
from collections import OrderedDict
import yaml
from yaml import Loader, Dumper
from yaml.representer import SafeRepresenter

def dict_representer(dumper, data):
   return dumper.represent_dict(data.items())
Dumper.add_representer(OrderedDict, dict_representer)

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))
Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Sum training weights of classes in training dataset.")
    parser.add_argument("--era", required=True, help="Experiment era")
    parser.add_argument("--channel", required=True, help="Analysis channel")
    parser.add_argument("--dataset",nargs='+', type=str,required=True, help="Training datasets.")
    parser.add_argument("--dataset-config-file", type=str,required=True, help="Specifies the config file created by ml/create_training_dataset.sh calling ml/write_dataset_config.py")
    parser.add_argument("--training-template", type=str,required=False, help="Specifies the config file setting the model, used variables...")
    parser.add_argument("--write-weights", type=bool, default=True, help="Overwrite inverse weights to ml/$era_$channel_training.yaml")
    parser.add_argument("--weight-branch", default="training_weight", type=str, help="Branch with weights.")
    return parser.parse_args()


def dictToString(exdict):
    return str(["{} : {}".format(key, value) for key, value in sorted(exdict.items(), key=lambda x: str(x[1]))])


def main(args):
    logger.info("Process training datasets %s.", args.dataset)
    f = [ROOT.TFile(filename) for filename in args.dataset]

    dsConfDict= yaml.load(open(args.dataset_config_file, "r"), Loader=yaml.SafeLoader)
    ### use the classes that have processes mapped to them
    classes = set([dsConfDict["processes"][key]["class"] for key in list(dsConfDict["processes"].keys())])

    if args.training_template == None:
        args.training_template= "ml/templates/{}_{}_training.yaml".format(args.era, args.channel)
    trainingTemplateDict=yaml.load(open(args.training_template, "r"), Loader=yaml.SafeLoader)

    ### Weight Calculation
    counts = []
    sum_all = 0.0
    for name in classes:
        logger.debug("Process class %s.", name)
        sum_ = 0.0
        for f_ in f:
            tree = f_.Get(name)
            if tree == None:
                logger.fatal("Tree %s does not exist in file.", name)
                raise Exception
            for event in tree:
                evw=getattr(event, args.weight_branch)
                if np.isnan( evw):
                    logger.fatal("Fatal: no event weight in class {} with ID {}".format(name,getattr(event, "event")))
                    raise Exception
                else:
                    sum_ += evw
        sum_all += sum_
        counts.append(sum_)

    ### Weight printing
    for i, name in enumerate(classes):
        logger.info(
            "Class {} (sum, fraction, inverse): {:g}, {:g}, {:g}".format(
                name, counts[i], counts[i] / sum_all, sum_all / counts[i]))

    logger.info( "{}-{}: Class weights before update: {}".format(args.era, args.channel,dictToString(trainingTemplateDict["class_weights"])))

    newWeightsDict={}
    for i, name in enumerate(classes):
        newWeightsDict[name]=sum_all / counts[i]

    ### Warning for big changes
    if set(list(trainingTemplateDict["class_weights"].keys()))==set(list(newWeightsDict.keys())):
        for i, name in enumerate(classes):
            oldweight=trainingTemplateDict["class_weights"][name]
            newweight=newWeightsDict[name]
            if newweight/oldweight > 2 or newweight/oldweight < .5:
                logger.warning( "{}-{}: Class weights for {} changing by more than a factor of 2".format(args.era, args.channel,name))
    else:
        logger.warning("Training classes in {} and {} differ".format(args.dataset_config_file,args.training_template))

    ## Sort the clases, so testing plots/... are easierer to compare
    priolist=["qqh","ggh","emb","ztt","tt","db","misc","zll","w","noniso","ss","ff"]
    odDict=OrderedDict({})
    for key in priolist:
        if key in newWeightsDict:
            odDict[key]=newWeightsDict[key]
            del newWeightsDict[key]
    ## add classes that are not in the priolist at the end
    for key in sorted(newWeightsDict):
        odDict[key]=newWeightsDict[key]
    del newWeightsDict

    ###
    # attach the weights dict with the classes to the dsConfDict
    dsConfDict["classes"]=list(odDict.keys())
    dsConfDict["class_weights"]=odDict


    ############ Logic for merging the configs
    if "classes" in list(trainingTemplateDict.keys()):
        if set(trainingTemplateDict["classes"])!=set(classes):
            logger.warning("Training classes in {} and {} differ".format(args.dataset_config_file,args.training_template))
        #exit 1

    mergeddict=OrderedDict({})
    ## check if there are relevant overwrites:
    ## dsConfDict has priority over trainingTemplateDict
    ## nothing that is provided by the write_dataset_config.py is overwritten
    for d in [dsConfDict, trainingTemplateDict]:
        for key in ["classes","class_weights"]+list(d.keys()):
            #if True:
            # making sure processes appear at the end
            if key not in ["processes"]:
                if key in mergeddict and mergeddict[key]!=d[key]:
                    logger.warning("Key overlap for key {}, {} should overwrite {}".format(key,d[key],mergeddict[key]))
                else: mergeddict[key]=d[key]
    mergeddict["processes"]=dsConfDict["processes"]
    with open(mergeddict["output_path"]+"/dataset_config.yaml","w") as f:
        yaml.dump(mergeddict, f,Dumper=Dumper, default_flow_style=False)

    logger.info( "{}-{}: Dict after merge: without processes".format(args.era, args.channel))
    print(dictToString({key: value for key, value in mergeddict.items() if key != "processes"}))

    logger.info( "{}-{}: Class weights after update: {}".format(args.era, args.channel,dictToString(trainingTemplateDict["class_weights"])))


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
