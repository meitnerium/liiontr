Installation
============

Development Installation
------------------------

Clone the LiionTR repository and create a Python virtual environment.

From the project root, install the package in editable mode:

.. code-block:: powershell

   python -m pip install -e .

Install the development dependencies required for testing and
documentation:

.. code-block:: powershell

   python -m pip install pytest ruff interrogate sphinx sphinx-rtd-theme


Optional Dependencies
---------------------

Some LiionTR capabilities may rely on optional scientific packages.

Cantera
~~~~~~~

Cantera is intended to provide thermochemical and chemical-equilibrium
capabilities for gas-phase thermal runaway products.

.. code-block:: powershell

   python -m pip install cantera

Cantera integration is currently under development.