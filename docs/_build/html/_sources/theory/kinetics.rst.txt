Reaction Kinetics
=================

Overview
--------

LiionTR separates reaction kinetics from reaction progress models.

The kinetic model determines the temperature dependence of a reaction
rate coefficient, while the progress model determines how the current
reaction state modifies the effective reaction rate.

For a generic thermal runaway reaction, LiionTR uses the conceptual form

.. math::

   \frac{d\alpha}{dt}
   =
   k(T) f(\alpha, \mathrm{context}),

where :math:`\alpha` is the reaction conversion, :math:`k(T)` is the
temperature-dependent kinetic rate, and :math:`f` is a progress
function that may depend on conversion or on the state of other
reactions.


Reaction Conversion
-------------------

Reaction progress is represented using a dimensionless conversion
variable

.. math::

   0 \leq \alpha \leq 1.

The interpretation is:

.. math::

   \alpha = 0

for an unreacted state and

.. math::

   \alpha = 1

for a fully reacted state.

The remaining reactant fraction is therefore

.. math::

   1 - \alpha.


Arrhenius Kinetics
------------------

The standard kinetic model in LiionTR is the Arrhenius law

.. math::

   k(T)
   =
   A
   \exp\left(
       -\frac{E_a}{R T}
   \right),

where

.. list-table::
   :header-rows: 1
   :widths: 20 45 20

   * - Symbol
     - Description
     - SI unit

   * - :math:`k(T)`
     - Temperature-dependent kinetic rate
     - depends on kinetic model

   * - :math:`A`
     - Pre-exponential factor
     - typically s\ :sup:`-1`

   * - :math:`E_a`
     - Activation energy
     - J mol\ :sup:`-1`

   * - :math:`R`
     - Universal gas constant
     - J mol\ :sup:`-1` K\ :sup:`-1`

   * - :math:`T`
     - Absolute temperature
     - K

The Arrhenius law introduces the strong temperature dependence that is
characteristic of thermal runaway reactions.

As temperature increases, the exponential term increases rapidly and
can produce strong positive feedback between heat generation and
reaction rate.


Temperature Threshold Kinetics
------------------------------

Some thermal runaway reactions are activated only above a prescribed
temperature threshold.

LiionTR provides a threshold kinetic model that returns zero below the
activation temperature and evaluates an underlying kinetic model above
the threshold.

Conceptually,

.. math::

   k_{\mathrm{threshold}}(T)
   =
   \begin{cases}
   0,
   &
   T < T_{\mathrm{threshold}},
   \\[6pt]
   k(T),
   &
   T \geq T_{\mathrm{threshold}}.
   \end{cases}

This representation is useful for empirical literature models where a
reaction is assumed to begin only after a specified onset temperature.


Progress Models
===============

Power-Law Progress
------------------

A basic progress model is based on the remaining reactant fraction,

.. math::

   f(\alpha)
   =
   (1-\alpha)^n,

where :math:`n` is the reaction order.

The resulting progress rate is

.. math::

   \frac{d\alpha}{dt}
   =
   k(T)(1-\alpha)^n.

For :math:`n=1`, this corresponds to a first-order reaction with respect
to the unreacted fraction.


Autocatalytic Progress
----------------------

LiionTR also supports autocatalytic progress laws.

A generic autocatalytic form may be written as

.. math::

   f(\alpha)
   =
   \alpha^m
   (1-\alpha)^n,

where :math:`m` controls the dependence on reacted material and
:math:`n` controls depletion of remaining reactant.

This type of expression can reproduce reaction behavior in which the
reaction rate initially increases as products or active sites are
formed, before decreasing as the reactant is depleted.


Threshold Progress
------------------

A progress model can also be conditioned on an external state variable.

LiionTR uses this mechanism to represent interactions between thermal
runaway reactions.

The generic form is

.. math::

   f_{\mathrm{threshold}}
   =
   \begin{cases}
   0,
   &
   x > x_{\mathrm{threshold}},
   \\[6pt]
   f(\alpha),
   &
   x \leq x_{\mathrm{threshold}},
   \end{cases}

where :math:`x` may represent the conversion or remaining fraction of
another reaction.

This allows one reaction to activate only after another reaction has
progressed sufficiently.


Exponential Inhibition
----------------------

LiionTR provides an exponential inhibition model of the form

.. math::

   f_{\mathrm{inhibition}}
   =
   f(\alpha)
   \exp(-x),

where :math:`x` is a state-dependent inhibition variable.

The inhibition factor approaches unity as :math:`x` approaches zero and
suppresses reaction progress for larger values of :math:`x`.

This model is particularly useful for thermal runaway mechanisms in
which a protective layer or remaining material inhibits another
reaction.


Reaction Context
----------------

Progress models may depend on quantities other than their own reaction
conversion.

LiionTR therefore provides a reaction context that exposes the current
state of the complete reaction network.

Context-dependent variables may include:

* the conversion of another reaction;
* the remaining fraction of another reaction;
* a normalized or transformed reaction state;
* quantities derived from multiple reaction states.

This architecture avoids embedding mutable state directly inside
reaction objects and keeps the numerical solver responsible for the
actual ODE state vector.


Linear Conversion Variables
---------------------------

A reaction-context variable may scale another reaction conversion
relative to a reference conversion.

Conceptually,

.. math::

   x
   =
   \frac{\alpha}{\alpha_{\mathrm{ref}}},

or an equivalent normalized relation depending on the configured
variable definition.

These variables are useful when a literature model defines reaction
coupling in terms of a normalized conversion quantity.


Remaining-Fraction Ratio Variables
----------------------------------

LiionTR also supports context variables defined from the remaining
fraction of another reaction.

The remaining fraction is

.. math::

   r
   =
   1-\alpha.

A normalized remaining-fraction quantity may then be constructed using
a reference value,

.. math::

   x
   =
   \frac{r}{r_{\mathrm{ref}}}.

This representation is used by models in which a protective layer or
reactant inventory controls another reaction rate.


Separation of Kinetics and Progress
-----------------------------------

The separation between kinetic and progress models is intentional.

The kinetic model answers:

   How strongly does temperature affect the intrinsic reaction rate?

The progress model answers:

   How does the current physical or chemical state modify that rate?

The complete reaction rate is obtained only after combining both
contributions.

This design makes it possible to reuse the same Arrhenius kinetics with
different progress laws and to reproduce complex literature models
without creating a separate monolithic class for every reaction.


Units and Conventions
---------------------

LiionTR uses SI units internally.

Typical kinetic quantities are:

.. list-table::
   :header-rows: 1
   :widths: 35 25

   * - Quantity
     - Unit

   * - Temperature
     - K

   * - Activation energy
     - J mol\ :sup:`-1`

   * - Pre-exponential factor
     - typically s\ :sup:`-1`

   * - Reaction conversion
     - dimensionless

   * - Reaction progress rate
     - s\ :sup:`-1`

Care is required when importing parameters from literature sources,
because activation energies, reaction enthalpies, and reactant contents
are often reported using non-SI units.


Current Limitations
-------------------

The present kinetic framework is primarily intended for reduced-order
thermal runaway models.

Current limitations include:

* no uncertainty propagation for kinetic parameters;
* no automatic parameter estimation from ARC or DSC experiments;
* no explicit pressure dependence in the standard kinetic models;
* no transport-limited reaction formulation;
* no automatic selection of competing kinetic mechanisms;
* no direct electrochemical state dependence such as local lithium
  concentration or electrode potential.

Future extensions may include parameter calibration, more general
kinetic laws, uncertainty analysis, and coupling to electrochemical
models such as PyBaMM.


Implementation
--------------

The kinetic model interface is defined by

:class:`liiontr.kinetics.model.KineticModel`.

Arrhenius kinetics are implemented by

:class:`liiontr.kinetics.arrhenius.Arrhenius`.

Temperature-threshold kinetics are implemented by

:class:`liiontr.kinetics.threshold.TemperatureThresholdKinetics`.

Reaction progress models are implemented in

:mod:`liiontr.kinetics.progress`.

Context-dependent reaction state variables are provided by

:mod:`liiontr.reactions.context_variable`.