#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import logging
logger = logging.getLogger("sum_training_weights")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

import yaml
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description="Sum training weights of classes in training dataset.")
    parser.add_argument("--tag", required=True, help="The tag of the folder, usually mc, emb, ff etc.")
    parser.add_argument("--tag_name", required=False, default="", help="The tag of the specific analyis folder.")
    parser.add_argument("--input_base_path", required=False, help="Base path of inputs.")
    parser.add_argument("--channel", required=True, help="Analysis channel")
    parser.add_argument("--output_dir", type=str,required=True, help="Output directory")
    return parser.parse_args()

def main(args):
    eras = ["2016", "2017", "2018"]
    configs = []
    for era in eras:
        config_path = "{}{}/{}_{}_{}/dataset_config.yaml".format(args.input_base_path, args.tag_name, era, args.channel,args.tag)
        logger.info("Try to open {}".format(config_path))
        config = yaml.load(open(config_path, 'r'), Loader =yaml.SafeLoader)
        configs.append(config)

    all_era_template = yaml.load(open('ml/templates/all_eras_training_{}.yaml'.format(args.channel)), Loader=yaml.SafeLoader)

    for i_era, era in enumerate(eras):
        all_era_template["datasets_{}".format(era)] = configs[i_era]["datasets"]
        all_era_template["class_weights_{}".format(era)] = configs[i_era]["class_weights"]

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    all_era_template["classes"] = configs[0]["classes"]
    all_era_template["variables"] = configs[0]["variables"]

    all_era_template["output_path"] = args.output_dir

    output_file = args.output_dir + "/dataset_config.yaml"

    logger.info("Writing new dataset config for all eras to {}".format(output_file))
    yaml.dump(all_era_template, open(output_file, 'w'), default_flow_style=False)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
