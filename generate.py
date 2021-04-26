# Dummy Python generation script, just to showcase the cmake integration
# This script just copies the template around and does not modify it

import argparse
from shutil import copyfile
from os import path

parser = argparse.ArgumentParser(description='Generate the C++ code for a given config')
parser.add_argument('template', type=str, help='Path to the template')
parser.add_argument('output', type=str, help='Path to the output directory')
parser.add_argument('config', type=str, help='Name of the config')
args = parser.parse_args()

files = [f'analysis_{i}.cxx' for i in range(20)]

with open(path.join(args.output, 'files.txt'), 'w') as f:
    for filename in files:
        copyfile(args.template, path.join(args.output, filename))
        f.write(filename + '\n')
