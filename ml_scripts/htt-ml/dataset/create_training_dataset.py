#!/usr/bin/env python

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True  # disable ROOT internal argument parser
#  ROOT.ROOT.EnableImplicitMT(12); # Tell ROOT you want to go parallel
import argparse
import yaml
import os
import subprocess
from array import array
from multiprocessing import Pool

import logging

logger = logging.getLogger("create_training_dataset")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)


def parse_arguments():
    logger.debug("Parse arguments.")
    parser = argparse.ArgumentParser(description="Create training dataset")
    parser.add_argument("config", help="Datasets config file")
    return parser.parse_args()


def parse_config(filename):
    logger.debug("Load YAML config: {}".format(filename))
    return yaml.load(open(filename, "r"), Loader=yaml.SafeLoader)


def generateRootFiles(jobconfig):
    process, num_fold, config = jobconfig
    logger.info(process)
    logger.debug("Collect events of process {} for fold {}.".format(process, num_fold))

    # Create output file
    output_filename = os.path.join(
        config["output_path"], "merge_fold{}_{}.root".format(num_fold, process)
    )

    # Collect all files for this process in a chain. Create also chains for friend files
    chain = ROOT.TChain(config["tree_path"])  ## "mt_nominal/ntuple"
    friendchains = {}
    for friendPath in config["friend_paths"]:  ####/ceph/htautau/2017/nnscore_friends/
        friendTreeName = os.path.basename(os.path.normpath(friendPath))
        friendchains[friendTreeName] = ROOT.TChain(config["tree_path"])

    # for each file, add ntuple TTree to the chain and do the same for the the friendTrees
    for filename in config["processes"][process]["files"]:
        path = os.path.join(config["base_path"], filename)
        if not os.path.isfile(path):
            logger.fatal("File does not exist: {}".format(path))
            raise Exception
        chain.AddFile(path)
        # Make sure, that friend files are put in the same order together
        for friendPath in config["friend_paths"]:
            friendFileName = os.path.join(friendPath, filename)
            if not os.path.isfile(friendFileName):
                logger.fatal("File does not exist: {}".format(friendFileName))
                raise Exception
            friendTreeName = os.path.basename(os.path.normpath(friendPath))
            logger.debug(
                "Attaching friendtree for {}, filename{}".format(
                    friendTreeName, friendFileName
                )
            )
            friendchains[friendTreeName].AddFile(friendFileName)

    logger.debug("Joining TChains")
    for friendTreeName in friendchains.keys():
        logger.debug("Adding to mainchain: {}".format(friendTreeName))
        chain.AddFriend(friendchains[friendTreeName], friendTreeName)

    logger.debug("Calculationg number of events")
    # Disable branch "nickname" of type "string" to prevent futile searching
    chain.SetBranchStatus("nickname",0)
    rdf = ROOT.RDataFrame(chain)
    chain_numentries = rdf.Count().GetValue()
    if chain_numentries == 0:
        logger.fatal("Chain (before skimming) does not contain any events.")
        raise Exception
    logger.info("Found {} events for process {}.".format(chain_numentries, process))

    # Skim the events with the cut string
    cut_string = "({EVENT_BRANCH}%2=={NUM_FOLD})&&({CUT_STRING})".format(
        EVENT_BRANCH=config["event_branch"],
        NUM_FOLD=num_fold,
        CUT_STRING=config["processes"][process]["cut_string"],
    )
    logger.debug("Skim events with cut string: {}".format(cut_string))

    rdf = rdf.Filter(cut_string)

    chain_skimmed_numentries = rdf.Count().GetValue()
    if not chain_skimmed_numentries > 0:
        logger.fatal("Chain (after skimming) does not contain any events.")
        raise Exception
    logger.debug(
        "Found {} events for process {} after skimming.".format(
            chain_skimmed_numentries, process
        )
    )

    # # Write training weight to new branch
    logger.debug(
        "Add training weights with weight string: {}".format(
            config["processes"][process]["weight_string"]
        )
    )
    rdf = rdf.Define(
        config["training_weight_branch"],
        "(float)(" + config["processes"][process]["weight_string"] + ")",
    )

    opt = ROOT.ROOT.RDF.RSnapshotOptions()
    opt.fMode = "RECREATE"

    rdf.Snapshot(config["processes"][process]["class"], output_filename, "^((?!nickname).)*$", opt)
    logger.info("snapshot created for process {}!".format(process))
    return (num_fold, output_filename)


def main(args, config):
    jobL = [
        (process, fold, config) for process in config["processes"] for fold in range(2)
    ]
    pool = Pool(12)
    created_files = pool.map(generateRootFiles, jobL)
    for num_fold in range(2):
        logger.info("Merge input files for fold {}.".format(num_fold))
        foldfile = [fn for fold, fn in created_files if fold == num_fold]
        # Combine all skimmed files using `hadd`
        logger.debug(
            "Call `hadd` to combine files of processes for fold {}.".format(num_fold)
        )
        output_file = os.path.join(
            config["output_path"],
            "fold{}_{}".format(num_fold, config["output_filename"]),
        )
        subprocess.call(["hadd", "-f", output_file] + foldfile)
        logger.info("Created output file: {}".format(output_file))


if __name__ == "__main__":
    args = parse_arguments()
    config = parse_config(args.config)
    main(args, config)
