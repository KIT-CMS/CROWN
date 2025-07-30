"""
This script produces a Python file containing NanoAOD branch definitions.
Each branch is written as a variable assignment to a NanoAODQuantity object,
accompanied by a documentation string that contains the branch data type and 
description. Branches are grouped by collections.

This script can be used if a new nanoAOD version or run specific version is 
needed for an analysis. The automatically generated file can be used to
import the branches of the nanoAOD version in the analysis code.
"""

import ROOT
import argparse
from collections import defaultdict


def get_all_unique_branch_info(input_files):
    """
    Collect all unique branch names and their info from multiple ROOT files.

    Args:
        input_files (list): List of input ROOT file paths.

    Returns:
        dict: Mapping branch name -> (btype, btitle)
    """
    branch_info = {}
    for input_file in input_files:
        f = ROOT.TFile.Open(input_file)
        tree = f.Get("Events")
        for branch in tree.GetListOfBranches():
            name = branch.GetName()
            if name not in branch_info:
                leaf = branch.GetLeaf(name)
                btype = branch.GetClassName() or (leaf.GetTypeName() if leaf else "unknown")
                btitle = branch.GetTitle()
                branch_info[name] = (btype, btitle)
        f.Close()
    return branch_info


def get_collection(branch_names):
    """
    Group branch names into collections based on their names.

    Args:
        branch_names (iterable): Iterable of branch names.

    Returns:
        tuple: (collections, nbranches)
    """
    collections = defaultdict(list)
    nbranches = {}
    for name in branch_names:
        if name.startswith("n") and len(name) > 1:
            collection = name[1:]
            nbranches[collection] = name
        elif "_" in name:
            collection = name.split("_")[0]
            collections[collection].append(name)
        else:
            collections["__other__"].append(name)
    return collections, nbranches

def add_branch_to_lines(name, branch_info, used, lines):
    """
    Add a branch definition to the lines list if it hasn't been used yet.

    Args:
        name (str): Branch name.
        branch_info (dict): Mapping branch name -> (btype, btitle).
        used (set): Set of already used branch names.
        lines (list): List to append the branch definition to.
    """
    if name in used:
        return
    btype, btitle = branch_info.get(name, ("unknown", "unknown"))
    left = f'{name} = NanoAODQuantity("{name}")'
    comment = f'\n"""dtype: {btype}; description: {btitle}"""'
    lines.append((left, comment))
    used.add(name)

def dump_nanoaod_collections(branch_info, output_file):
    """
    Write a Python file containing the definitions of unique branches.

    Args:
        branch_info (dict): Mapping branch name -> (btype, btitle)
        output_file (str): Path to the output Python file.
    """
    collections, nbranches = get_collection(branch_info.keys())

    lines = []
    used = set()
    for collection in sorted(collections):
        nbranch = nbranches.get(collection)
        if nbranch:
            add_branch_to_lines(nbranch, branch_info, used, lines)
        for name in sorted(collections[collection]):
            add_branch_to_lines(name, branch_info, used, lines)
        lines.append(("", ""))  # blank line between collections

    # Add any branches not already used (e.g. uncategorized)
    for name in branch_info:
        add_branch_to_lines(name, branch_info, used, lines)

    maxlen = max(len(left) for left, _ in lines if left)

    with open(output_file, "w") as out:
        out.write("from code_generation.quantity import NanoAODQuantity\n\n")
        for left, comment in lines:
            if left == "" and comment == "":
                out.write("\n")
            else:
                out.write(f"{left.ljust(maxlen)}  {comment}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Produce a Python file with all unique NanoAOD branches from multiple ROOT files."
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="Input NanoAOD ROOT files (space separated list)."
    )
    parser.add_argument(
        "output_file",
        help="Output Python file."
    )
    args = parser.parse_args()

    branch_info = get_all_unique_branch_info(args.input_files)
    dump_nanoaod_collections(branch_info, args.output_file)