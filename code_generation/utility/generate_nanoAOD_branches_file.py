"""
Generates a Python file containing `NanoAODQuantity` definitions for all
unique branches found in a set of NanoAOD ROOT files.

This script scans ROOT files, finds all branches, automatically groups them
into collections (e.g., 'Jet', 'Muon') based on their naming conventions,
and writes them into a structured, importable Python module. This can be used
for creating (run-specific) branch definitions for new nanoAOD versions.
"""

import ROOT
import os
import argparse
from collections import defaultdict
from typing import List, Dict, Tuple, Set


def get_all_unique_branch_info(input_files: List[str]) -> Dict[str, Tuple[str, str]]:
    """
    Scans multiple NanoAOD ROOT files and extracts information about all
    unique branches.

    It iterates through the 'Events' TTree of each input file. For each
    branch, it records its name, C++ data type, and its description.
    Branches with the same name across different files are only recorded once.

    Args:
        input_files: A list of paths to the input ROOT files.

    Returns:
        A dictionary where keys are unique branch names and values are tuples
        containing the branch's data type (str) and title (str).
        Example: {'Jet_pt': ('Float_t', 'pt for jet')}
    """
    branch_info = {}
    print("--- Scanning input files for branches ---")
    for input_file in input_files:
        print(f"  - Processing {input_file}...")
        f = ROOT.TFile.Open(input_file)
        tree = f.Get("Events")
        for branch in tree.GetListOfBranches():
            name = branch.GetName()
            if name not in branch_info:
                leaf = branch.GetLeaf(name)
                btype = leaf.GetTypeName()
                bdescription = branch.GetTitle()
                branch_info[name] = (btype, bdescription)
        f.Close()
    print(f"Found {len(branch_info)} unique branches.\n")
    return branch_info


def get_collections(
    branch_names: List[str],
) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """
    Groups a list of branch names into collections based on NanoAOD naming conventions.

    It identifies collections using two rules:
    1.  'n-branches': A branch like 'nJet' defines a collection named 'Jet'.
    2.  Shared Prefixes: Branches sharing a common prefix up to the last
        underscore (e.g., 'Jet_pt', 'Jet_eta') are grouped into a 'Jet' collection.

    Args:
        branch_names: A flat list of all unique branch names from the ROOT files.

    Returns:
        A tuple containing two dictionaries:
        1.  A map from collection names to a list of their corresponding branch names.
            Example: {'Jet': ['Jet_pt', 'Jet_eta', ...]}
        2.  A map from a collection name to its corresponding 'n' branch name.
            Example: {'Jet': 'nJet'}
    """
    collections = set()
    underscore_branches = []
    nbranches = {}

    # Identify 'n-branch' collections and gather other branches for prefix analysis
    for branch in branch_names:
        if branch.startswith("n"):
            collection_name = branch[1:]
            nbranches[collection_name] = branch
            collections.add(collection_name)
        elif "_" in branch:
            underscore_branches.append(branch)

    # Find collections from shared prefixes
    underscore_branches.sort()
    for i in range(len(underscore_branches) - 1):
        common_prefix = os.path.commonprefix(
            [underscore_branches[i], underscore_branches[i + 1]]
        )
        if "_" in common_prefix:
            collection_candidate = common_prefix[: common_prefix.rfind("_")]
            if collection_candidate:
                collections.add(collection_candidate)

    # Assign all branches to their identified collections
    collection_map = defaultdict(list)
    sorted_collections = sorted(list(collections), key=len, reverse=True)

    for branch in branch_names:
        # Check for 'n-branch' case first for direct assignment
        if branch.startswith("n") and len(branch) > 1:
            collection_map[branch[1:]].append(branch)
            continue

        # Check for shared prefix, matching longest collections first
        for collection_name in sorted_collections:
            if branch.startswith(collection_name + "_"):
                collection_map[collection_name].append(branch)
                break
            # Handle cases where the branch name is the collection name (e.g., event)
            elif branch == collection_name:
                collection_map[collection_name].append(branch)
                break

    return (dict(collection_map), nbranches)


def add_branch_to_lines(
    name: str,
    branch_info: Dict[str, Tuple[str, str]],
    used: Set[str],
    lines: List[Tuple[str, str]],
) -> None:
    """
    Formats a single branch into a Python definition and adds it to the output list.

    This is a helper function that checks if a branch has already been processed.
    If not, it formats the branch as a `NanoAODQuantity` definition string and a
    corresponding docstring comment with its type and description. It then appends
    this pair to the `lines` list and marks the branch as used.

    Args:
        name: The name of the branch to process.
        branch_info: The main dictionary with all branch data.
        used: A set to track which branches have already been written to avoid
              duplicates. This argument is modified in-place.
        lines: The list where the formatted (code, comment) tuples are appended.
               This argument is modified in-place.
    """
    if name in used:
        return
    btype, bdescription = branch_info.get(name)
    quantity = f'{name} = NanoAODQuantity("{name}")'
    comment = f'\n"""dtype: {btype}; description: {bdescription} """'
    lines.append((quantity, comment))
    used.add(name)


def dump_nanoaod_collections(
    branch_info: Dict[str, Tuple[str, str]], output_file: str
) -> None:
    """
    Generates and writes the final Python file with all NanoAOD branch definitions.

    This function orchestrates the file generation. It first calls `get_collections`
    to group branches. It then iterates through the sorted collections, writing out
    the `NanoAODQuantity` definition for each branch. Uncategorized branches are
    added at the end. The final output is written to the specified file.

    Args:
        branch_info: The dictionary of all unique branches and their metadata.
        output_file: The path to the Python file to be created.
    """
    collections, nbranches = get_collections(branch_info.keys())

    lines = []
    used = set()
    print("--- Grouping branches into collections ---")
    for collection in sorted(collections):
        # The 'n' branch should always come first in its collection.
        nbranch = nbranches.get(collection)
        if nbranch:
            add_branch_to_lines(nbranch, branch_info, used, lines)

        # Sort the rest of the branches in the collection alphabetically.
        for name in sorted(collections[collection]):
            add_branch_to_lines(name, branch_info, used, lines)
        lines.append(("", ""))  # blank line between collections

    # Add any remaining branches that were not categorized into a collection
    for name in branch_info:
        add_branch_to_lines(name, branch_info, used, lines)

    print(f"--- Writing output to analysis_configurations/quantities/{output_file} ---")
    with open("analysis_configurations/quantities/" + output_file, "w") as out:
        out.write("from code_generation.quantity import NanoAODQuantity\n\n")
        for quantity, comment in lines:
            if quantity == "" and comment == "":
                out.write("\n")
            else:
                out.write(f"{quantity}")
                out.write(f"{comment}\n")
    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Produce a Python file with all unique NanoAOD branches from multiple ROOT files."
    )
    parser.add_argument(
        "-i",
        "--input-files",
        nargs="+",
        required=True,
        help="One or more input NanoAOD ROOT files (space-separated list).",
        metavar="IN_FILES",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        required=True,
        help="Name of the output Python file. The file is saved in 'analysis_configurations/quantities/'",
        metavar="OUT_FILE",
    )

    args = parser.parse_args()

    branch_info = get_all_unique_branch_info(args.input_files)
    dump_nanoaod_collections(branch_info, args.output_file)
