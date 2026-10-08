# LiionTR

LiionTR (Lithium-Ion Thermal Runaway) is a modular scientific Python framework
for reduced-order simulation of lithium-ion battery thermal runaway.

The project is designed around explicit, testable couplings between reaction
kinetics, heat release, gas generation, thermochemistry, pressure evolution,
and venting. LiionTR owns the thermal-runaway physics and orchestration layer;
specialized packages such as Cantera and OpenFOAM are used as complementary
backends when higher-fidelity thermochemistry or spatial transport is needed.

## Current capabilities

The current codebase includes:

- configurable cell geometry and effective material properties;
- Arrhenius, threshold, power-law, autocatalytic, and inhibited reaction models;
- context-dependent and multi-channel reaction networks;
- a lumped thermal model with convective heat loss;
- SciPy BDF integration with temperature, pressure, and vent-opening events;
- empirical species-yield gas generation;
- species-resolved gas inventories and ideal-gas pressure reconstruction;
- choked and unchoked compressible venting with irreversible vent opening;
- elemental inventories and reaction-derived element generation;
- energy-scaled elemental generation calibrated from vent-gas data;
- first-order delayed release of generated elements to the gas volume;
- an optional Cantera equilibrium backend supporting constant-`T,V` and
  constant-`U,V` equilibrium states;
- coupled cell/gas thermal-energy accounting and energy carried by vented gas;
- literature models/data based on Hu et al. (2020) and Howard et al. (2025).

A companion OpenFOAM verification effort has also demonstrated a one-way
thermal-source interface in which LiionTR supplies volumetric heat generation
`q'''(t)` in W/m³ and OpenFOAM resolves spatial heat transport. The verified
V0–V4 cases are not yet packaged in this repository snapshot.

## Model hierarchy

LiionTR currently supports two gas-modeling paths:

1. **Empirical species-yield path**

   Reaction progress → prescribed species yields → gas inventory → ideal-gas
   pressure → compressible venting.

2. **Elemental thermochemistry path**

   Reaction/energy release → elemental inventory → optional release dynamics →
   Cantera equilibrium → thermodynamic state/pressure → element-consistent
   venting and vent-energy removal.

The empirical path remains useful for reduced-order calibration and regression
verification. The elemental path provides the foundation for composition-aware
thermochemistry.

## Quick start

Create an environment and install LiionTR in editable mode:

```bash
python -m venv .venv
source .venv/bin/activate   # Linux/macOS
python -m pip install -e .
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e .
```

Then run a bundled example:

```bash
python examples/first_simulation.py
python examples/hu2020_simulation.py
python examples/vented_thermal_runaway.py
```

Cantera is optional:

```bash
python -m pip install cantera
```

## Documentation

The Sphinx documentation is under `docs/` and covers:

- installation and usage;
- physical and mathematical models;
- software architecture;
- thermochemical and gas modeling;
- verification status;
- API reference.

Build it with the documentation dependencies installed:

```bash
python -m pip install sphinx sphinx-rtd-theme sphinxcontrib-bibtex
sphinx-build -b html -W docs docs/_build/html
```

## Verification and validation status

LiionTR follows a verification-first development strategy. The repository
contains unit, coupled, and regression tests for the reduced-order model.
Companion OpenFOAM manufactured-solution cases have additionally confirmed
first-order temporal convergence, second-order spatial convergence, and the
external tabulated heat-source interface.

These results establish **numerical verification**, not complete physical
validation. Experimental calibration and validation against ARC/DSC, pressure,
venting, and gas-composition data remain active development priorities.

## Near-term roadmap

The next priorities are:

1. package a real LiionTR thermal-source exporter for OpenFOAM;
2. version the V0–V4 OpenFOAM verification cases with the repository;
3. extend/calibrate elemental thermochemistry beyond the current mechanisms and
   carrier-species assumptions;
4. strengthen validation against experimental thermal, pressure, venting, and
   composition datasets;
5. add sensitivity/uncertainty workflows;
6. investigate PyBaMM coupling for electrochemical initial state and abuse
   initiation;
7. progress from lumped cells toward spatial cell and multi-cell propagation
   models.

## Project status

LiionTR is a research codebase under active development. Interfaces, parameter
sets, and physical assumptions may evolve as verification and experimental
validation proceed.
