LiionTR Documentation
=====================

LiionTR is a scientific Python framework for reduced-order modeling of
lithium-ion battery thermal runaway.

The framework couples reaction kinetics, heat release, gas generation,
thermochemistry, pressure evolution, and venting through explicit and testable
interfaces. Optional external tools complement LiionTR: Cantera provides an
implemented equilibrium-thermochemistry backend, while OpenFOAM is being used
as a spatial thermal/CFD backend through a verified tabulated heat-source
interface.

.. note::

   LiionTR is under active development. Numerical verification is relatively
   mature for the implemented reduced-order subsystems, but experimental
   validation and calibration remain ongoing work.

Project Status
--------------

.. toctree::
   :maxdepth: 2

   project_status

Getting Started
---------------

.. toctree::
   :maxdepth: 2

   installation
   quickstart

Theory
------

.. toctree::
   :maxdepth: 2

   theory/index

Architecture
------------

.. toctree::
   :maxdepth: 2

   architecture/index

Verification
------------

.. toctree::
   :maxdepth: 2

   verification/index

API Reference
-------------

.. toctree::
   :maxdepth: 2

   api/index

References
----------

.. toctree::
   :maxdepth: 1

   references

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
