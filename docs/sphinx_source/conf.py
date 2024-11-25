# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../../"))

## for read the docs
import subprocess, os


def configureDoxyfile(input_dir, output_dir):
    with open("../Doxyfile.in", "r") as file:
        filedata = file.read()

    filedata = filedata.replace("@DOXYGEN_INPUT_DIR@", input_dir)
    filedata = filedata.replace("@DOXYGEN_OUTPUT_DIR@", output_dir)

    with open("Doxyfile", "w") as file:
        file.write(filedata)


# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"

breathe_projects = {}

if read_the_docs_build:
    input_dir = "../../src/"
    output_dir = "build"
    configureDoxyfile(input_dir, output_dir)
    subprocess.call("doxygen", shell=True)
    breathe_projects["CROWN"] = output_dir + "/xml"

# ...

# -- Project information -----------------------------------------------------

project = "CROWN"
copyright = "2023, Sebastian Brommer, Maximilian Burkart, Artur Gottmann, Sebastian Wozniewski, Stefan Wunsch"
author = "Sebastian Brommer, Maximilian Burkart, Artur Gottmann, Sebastian Wozniewski, Stefan Wunsch"

# The full version, including alpha/beta/rc tags
release = "0.3"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = []

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_logo = "../logos/crown_logo_bw.svg"
# html_static_path = ["../_static"]
html_theme_options = {
    "logo_only": True,
    "display_version": False,
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

# ...
extensions = [
    "breathe",
    "sphinx.ext.viewcode",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
    "sphinx.ext.coverage",
]
# ...

autodoc_typehints = "description"
autoclass_content = "both"

# Breathe Configuration
breathe_default_project = "CROWN"

exhale_args = {
    # Suggested optional arguments
    "createTreeView": True,
    "treeViewIsBootstrap": True,
    "pageLevelConfigMeta": ":github_url: https://github.com/KIT-CMS/CROWN",
}
breathe_domain_by_extension = {"h": "cpp", "hxx": "cpp"}
pygments_style = "sphinx"
