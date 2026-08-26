"""Sphinx configuration for the LiionTR documentation."""

from __future__ import annotations

import sys
from pathlib import Path


# Make the LiionTR source package importable by Sphinx.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = PROJECT_ROOT / "src"

sys.path.insert(0, str(SOURCE_ROOT))


# -- Project information -----------------------------------------------------

project = "LiionTR"
copyright = "2026, François Dion"
author = "François Dion"
release = "0.1.0"


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
]

templates_path = ["_templates"]

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]


# -- Autodoc configuration ---------------------------------------------------

autodoc_member_order = "bysource"

autodoc_typehints = "description"

autodoc_class_signature = "separated"

autosummary_generate = True


# -- Napoleon configuration --------------------------------------------------

napoleon_google_docstring = False
napoleon_numpy_docstring = True

napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True

napoleon_use_param = True
napoleon_use_rtype = True


# -- Intersphinx configuration -----------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}


# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_static_path = ["_static"]

html_title = "LiionTR Documentation"

html_show_sourcelink = True

# -- Options for Bibtex -------------------------------------------------
bibtex_bibfiles = ["../references.bib"]
bibtex_reference_style = "author_year"