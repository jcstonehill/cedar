# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys, os

# Determine if we're on Read the Docs server
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# On Read the Docs, we need to mock a few third-party modules so we don't get
# ImportErrors when building documentation
from unittest.mock import MagicMock

MOCK_MODULES = [
    'openmoc', 'openmc.data.reconstruct',
]
sys.modules.update((mod_name, MagicMock()) for mod_name in MOCK_MODULES)


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Cedar'
copyright = '2025, Jacob Stonehill'
author = 'Jacob Stonehill'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
]

autosummary_generate = True

templates_path = ['_templates']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'venv', '.venv', '_templates']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"

html_static_path = ['_static']
html_css_files = [
    "custom.css"
]