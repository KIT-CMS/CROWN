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
import sys
from collections import defaultdict

def get_collection(branches):
    """
    Group branches into collections based on their names.

    A branch whose name starts with "n" is assumed to be a counter for a collection.
    Other branches that contain an underscore ("_") are grouped by the part of the name
    before the underscore. Branches that don't match these rules are grouped in '__other__'.

    Args:
        branches (list): List of ROOT.TBranch objects.

    Returns:
        tuple: A tuple (collections, nbranches) where:
            - collections is a defaultdict(list) grouping branches by collection name.
            - nbranches is a dict mapping collection names to their corresponding nCOLLECTIONNAME branch.
    """
    collections = defaultdict(list)
    nbranches = {}
    for branch in branches:
        name = branch.GetName()
        if name.startswith("n") and len(name) > 1:
            # For collection counter branches, e.g. nElectron
            collection = name[1:]
            nbranches[collection] = branch
        elif "_" in name:
            # For branches named like Electron_pt, group by 'Electron'
            collection = name.split("_")[0]
            collections[collection].append(branch)
        else:
            # Uncategorized branches go into __other__
            collections["__other__"].append(branch)
    return collections, nbranches

def dump_nanoaod_collections(input_file, output_file):
    """
    Process a NanoAOD ROOT file and writes a Python file containing the definitions of branches.

    Each branch is output as a variable assignment to NanoAODQuantity("<branch_name>"),
    followed by a triple-quoted documentation string that contains the branch's data type and description.
    Branches are grouped by collection. If a counter branch (nCOLLECTIONNAME) exists, it is output first
    in that collection.

    Args:
        input_file (str): Path to the input NanoAOD ROOT file.
        output_file (str): Path to the output Python file.
    """
    # Open the ROOT file and get the 'Events' tree
    f = ROOT.TFile.Open(input_file)
    tree = f.Get("Events")
    branches = list(tree.GetListOfBranches())
    collections, nbranches = get_collection(branches)

    # Prepare list of lines with assignments and documentation strings.
    lines = []
    used = set()
    for collection in sorted(collections):
        # Add nCOLLECTIONNAME variable at the top if it exists
        nbranch = nbranches.get(collection)
        if nbranch:
            name = nbranch.GetName()
            leaf = nbranch.GetLeaf(name)
            btype = nbranch.GetClassName() or (leaf.GetTypeName() if leaf else "unknown")
            btitle = nbranch.GetTitle()
            left = f'{name} = NanoAODQuantity("{name}")'
            # The documentation string that will appear on hover in the IDE.
            comment = f'\n"""dtype: {btype}; description: {btitle}"""'
            lines.append((left, comment))
            used.add(name)
        # Add collection branches (excluding the nCOLLECTIONNAME branch, if already added).
        for branch in sorted(collections[collection], key=lambda b: b.GetName()):
            name = branch.GetName()
            if name in used:
                continue
            leaf = branch.GetLeaf(name)
            btype = branch.GetClassName() or (leaf.GetTypeName() if leaf else "unknown")
            btitle = branch.GetTitle()
            left = f'{name} = NanoAODQuantity("{name}")'
            comment = f'\n"""dtype: {btype}; description: {btitle}"""'
            lines.append((left, comment))
            used.add(name)
        lines.append(("", ""))  # Insert a blank line between collections

    # Process branches not assigned to any collection.
    for branch in branches:
        name = branch.GetName()
        if name in used:
            continue
        leaf = branch.GetLeaf(name)
        btype = branch.GetClassName() or (leaf.GetTypeName() if leaf else "unknown")
        btitle = branch.GetTitle()
        left = f'{name} = NanoAODQuantity("{name}")'
        comment = f'\n"""dtype: {btype}; description: {btitle}"""'
        lines.append((left, comment))

    # Compute maximum left-hand side length for alignment
    maxlen = max(len(left) for left, _ in lines if left)

    # Write the output file with aligned assignment and documentation strings.
    with open(output_file, "w") as out:
        out.write("from code_generation.quantity import NanoAODQuantity\n\n")
        for left, comment in lines:
            if left == "" and comment == "":
                out.write("\n")
            else:
                out.write(f"{left.ljust(maxlen)}  {comment}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python produce_nanoAOD_branches_file.py input.root output.py")
        sys.exit(1)
    dump_nanoaod_collections(sys.argv[1], sys.argv[2])