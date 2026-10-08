Installation
============

Python Version
--------------

The package metadata currently declares Python 3.12 or newer.

Development Installation
------------------------

Clone or unpack the LiionTR repository, create a virtual environment, and
install the package in editable mode.

Linux or macOS:

.. code-block:: bash

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install -e .

Windows PowerShell:

.. code-block:: powershell

   python -m venv .venv
   .\.venv\Scripts\Activate.ps1
   python -m pip install -e .

Development Tools
-----------------

The project defines a ``dev`` extra for testing and code-quality tools:

.. code-block:: bash

   python -m pip install -e ".[dev]"

The current ``dev`` extra includes pytest, coverage support, Ruff, Black, mypy,
and pre-commit.

Documentation Tools
-------------------

Sphinx dependencies are currently installed separately:

.. code-block:: bash

   python -m pip install sphinx sphinx-rtd-theme sphinxcontrib-bibtex

Build the documentation from the project root with:

.. code-block:: bash

   sphinx-build -b html -W docs docs/_build/html

Optional Dependencies
---------------------

Cantera
~~~~~~~

Cantera is an optional dependency used by the implemented equilibrium
thermochemistry backend. It is required for ``CanteraEquilibriumBackend`` and
for tests exercising elemental equilibrium, coupled gas energy, and elemental
venting.

.. code-block:: bash

   python -m pip install cantera

Without Cantera, the empirical species-yield gas path and the rest of the core
package remain usable; Cantera-specific tests are expected to skip.

OpenFOAM
~~~~~~~~

OpenFOAM is not a Python dependency of LiionTR. It is used as an external
spatial transport/CFD backend. The current coupling work exchanges tabulated
volumetric heat-generation histories rather than importing OpenFOAM into the
Python package.
