# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'cardioemulator'
copyright = '2021, Francesco Regazzoni'
author = 'Francesco Regazzoni'

import cardioemulator
version = cardioemulator.__version__
release = cardioemulator.__version__

# -- General configuration ---------------------------------------------------

extensions = [
	"sphinx.ext.autodoc",
	"sphinx.ext.napoleon",
	"sphinx.ext.autosummary",
	"sphinx.ext.viewcode",
    'nbsphinx',
]

autosummary_generate = True  # Turn on sphinx.ext.autosummary
autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries
html_show_sourcelink = False  # Remove 'view source code' from top of page (for html, not python)
autodoc_inherit_docstrings = True  # If no class summary, inherit base class summary

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Copy example notebooks into docs/_examples ------------------------------

import shutil
print("Copy example notebooks into docs/_examples")
shutil.rmtree('_examples', ignore_errors=True)
os.mkdir('_examples')
shutil.copyfile('../example/example.ipynb', '_examples/example.ipynb')

# -- Options for HTML output -------------------------------------------------

html_theme = 'pydata_sphinx_theme'
"_static/pandas.svg"

html_theme_options = {
    "github_url": "https://github.com/FrancescoRegazzoni/cardioemulator",
    "show_toc_level": 1,
	"navbar_start": ["navbar-logo", "version"],
	"icon_links": [
        {
            "name": "Paper",
            "url": "https://www.mate.polimi.it/biblioteca/add/qmox/35-2021.pdf",
            "icon": "fa fa-book",
        },
    ],
}

html_static_path = ['_static']