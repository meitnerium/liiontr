Reaction Models
===============

Overview
--------

LiionTR represents thermal runaway chemistry as a network of reaction
models coupled through a common state vector.

Each reaction model contributes two quantities to the governing
equations:

* a reaction progress rate;
* a mass-specific heat-generation rate.

For a reaction :math:`i`, the generic progress equation is

.. math::

   \frac{d\alpha_i}{dt}
   =
   k_i(T)
   f_i(\alpha_i, \mathrm{context}),

where :math:`\alpha_i` is the reaction conversion, :math:`k_i(T)` is
the kinetic rate, and :math:`f_i` is the reaction progress function.

The reaction framework is intentionally independent from the numerical
solver. Reaction objects evaluate rates from a supplied thermodynamic
state but do not internally advance their own state.


Reaction State
--------------

Each reaction is represented by a dimensionless conversion variable

.. math::

   0 \leq \alpha_i \leq 1.

The corresponding remaining reactant fraction is

.. math::

   r_i = 1 - \alpha_i.

The numerical solver stores all reaction conversions in its state
vector. Reaction objects therefore remain stateless with respect to
integration history.


Reaction Rate
-------------

For a standard reaction model, LiionTR combines a kinetic model and a
progress model.

The progress rate is

.. math::

   \dot{\alpha}
   =
   \frac{d\alpha}{dt}
   =
   k(T) f(\alpha, \mathrm{context}).

The kinetic model determines the temperature dependence, while the
progress model determines the dependence on reaction state and
potentially on the state of other reactions.


Reaction Enthalpy
-----------------

A reaction is assigned a specific reaction enthalpy

.. math::

   \Delta H,

expressed in J/kg of reacting material.

For an exothermic thermal runaway model, LiionTR uses the reaction
enthalpy as the amount of heat released as the reaction progresses.

If a reaction occupies a mass fraction :math:`w` of the total cell,
the mass-specific heat-generation rate relative to total cell mass can
be written as

.. math::

   \dot{q}
   =
   w \Delta H
   \frac{d\alpha}{dt},

with units of W/kg of cell.

For a cell of mass :math:`m_{\mathrm{cell}}`, the corresponding total
heat-generation rate is

.. math::

   \dot{Q}
   =
   m_{\mathrm{cell}}
   \dot{q}.

Thus,

.. math::

   \dot{Q}
   =
   m_{\mathrm{cell}}
   w \Delta H
   \frac{d\alpha}{dt}.


Mass Fraction
-------------

The reaction mass fraction

.. math::

   w

represents the fraction of total cell mass associated with the material
participating in a given reaction model.

This quantity is dimensionless.

For literature models that report reactant content per unit cell
volume, LiionTR provides helper classes that convert volumetric content
to an effective cell mass fraction.

If a reacting material has volumetric content

.. math::

   W_r

in kg/m\ :sup:`3`, and the effective cell density is :math:`\rho`,
then the corresponding mass fraction is

.. math::

   w
   =
   \frac{W_r}{\rho}.


Reaction Context
----------------

Some thermal runaway reactions depend on the state of other reactions.

LiionTR represents these dependencies through a
:class:`liiontr.reactions.context.ReactionContext`.

The reaction context contains the current conversion vector and exposes
quantities such as:

* conversion of another reaction;
* remaining fraction of another reaction;
* normalized context variables derived from reaction state.

This architecture allows reaction coupling without storing mutable
cross-reaction state inside individual reaction objects.


Context Variables
-----------------

Context variables transform reaction-network state into quantities used
by progress models.

Examples include normalized conversion variables and remaining-fraction
ratios.

A generic normalized remaining-fraction variable may be written as

.. math::

   x
   =
   \frac{1-\alpha_j}
        {r_{\mathrm{ref}}},

where :math:`\alpha_j` is the conversion of another reaction and
:math:`r_{\mathrm{ref}}` is a reference remaining fraction.

Such variables can represent effects including protective-layer
thickness, depletion of an inhibitor, or activation of a subsequent
reaction.


Reaction Object
---------------

The standard :class:`liiontr.reactions.reaction.Reaction` combines:

* a reaction name;
* a kinetic model;
* a reaction enthalpy;
* a reacting-material mass fraction;
* a progress model.

Its two principal operations are:

.. math::

   \dot{\alpha}
   =
   k(T) f(\alpha, \mathrm{context})

and

.. math::

   \dot{q}
   =
   w \Delta H \dot{\alpha}.

This makes the reaction object the basic physical unit used by a
reaction network.


Multi-Channel Reactions
-----------------------

Some physical reactions are represented by multiple parallel kinetic
channels that share a single reaction conversion.

LiionTR represents these models using
:class:`liiontr.reactions.multichannel.MultiChannelReaction`.

For :math:`N` kinetic channels sharing conversion :math:`\alpha`, each
channel has its own rate coefficient and reaction enthalpy.

A conceptual heat-generation expression is

.. math::

   \dot{q}
   =
   w
   f(\alpha, \mathrm{context})
   \sum_{j=1}^{N}
   \Delta H_j k_j(T).

The channels therefore contribute separately to heat generation while
sharing the same reaction state.

This structure is useful for literature models in which one physical
material undergoes multiple simultaneous decomposition pathways.


Shared Conversion
-----------------

For a multi-channel reaction, only one conversion variable is integrated.

The progress rate must therefore represent the combined progression of
the physical reaction rather than introducing one independent state
variable for every kinetic channel.

This distinction is important because parallel heat-generation channels
do not necessarily imply independent reactant inventories.


Reaction Network
----------------

A :class:`liiontr.reactions.network.ReactionNetwork` groups multiple
reaction models into a single thermal runaway mechanism.

For a network containing :math:`N` reactions, the conversion state is

.. math::

   \boldsymbol{\alpha}
   =
   \left[
       \alpha_1,
       \alpha_2,
       \ldots,
       \alpha_N
   \right].

The network evaluates the progress-rate vector

.. math::

   \frac{d\boldsymbol{\alpha}}{dt}
   =
   \left[
       \dot{\alpha}_1,
       \dot{\alpha}_2,
       \ldots,
       \dot{\alpha}_N
   \right].

Each reaction receives the current temperature, its own conversion, and
a reaction context representing the complete network state.


Network Heat Generation
-----------------------

The total mass-specific heat-generation rate is the sum of the
individual reaction contributions,

.. math::

   \dot{q}_{\mathrm{total}}
   =
   \sum_{i=1}^{N}
   \dot{q}_i.

For a cell of mass :math:`m_{\mathrm{cell}}`, the total heat-generation
rate supplied to the thermal model is

.. math::

   \dot{Q}_{\mathrm{gen}}
   =
   m_{\mathrm{cell}}
   \dot{q}_{\mathrm{total}}.

This value is then used in the lumped thermal energy balance.


Reaction Network Backend
------------------------

LiionTR provides a chemistry backend that connects a reaction network to
the thermal problem.

The backend evaluates:

* reaction progress rates;
* total reaction heat generation.

This layer separates the generic thermal problem from the specific
reaction-network implementation.

The architecture also leaves room for alternative chemistry backends in
the future.


State Coupling
--------------

The complete thermal-reaction subsystem may be represented as

.. math::

   \frac{dT}{dt}
   =
   \frac{
       \dot{Q}_{\mathrm{gen}}
       -
       \dot{Q}_{\mathrm{loss}}
   }{
       C_{\mathrm{th}}
   },

together with

.. math::

   \frac{d\alpha_i}{dt}
   =
   k_i(T)
   f_i(\boldsymbol{\alpha}),
   \qquad
   i = 1,\ldots,N.

The coupling is bidirectional:

1. temperature controls reaction rates through the kinetic models;
2. reactions release heat;
3. heat increases temperature;
4. increased temperature accelerates subsequent reactions.

This positive feedback is the fundamental mechanism of thermal runaway
in the reduced-order model.


Initial Conversions
-------------------

A reaction network requires one initial conversion for each reaction.

The initial state is

.. math::

   \boldsymbol{\alpha}_0
   =
   \left[
       \alpha_{1,0},
       \ldots,
       \alpha_{N,0}
   \right].

Literature models sometimes specify an initial remaining reactant
fraction rather than an initial conversion.

LiionTR uses

.. math::

   \alpha_0
   =
   1-r_0

to convert between these representations.


Validation
----------

The reaction network validates the consistency of reaction definitions
and conversion vectors before evaluating rates.

Examples of required consistency include:

* unique reaction names;
* one conversion value per reaction;
* physically meaningful conversion values;
* valid reaction mass fractions and kinetic parameters.

These checks help detect model-definition errors before numerical
integration begins.


Current Assumptions
-------------------

The current reaction framework assumes:

* each integrated reaction is represented by a scalar conversion;
* reaction state is spatially uniform throughout the cell;
* heat-release parameters are effective lumped quantities;
* reaction enthalpies do not vary explicitly with temperature;
* reaction rates do not directly depend on gas pressure unless encoded
  through a custom model;
* reaction products are not automatically determined by chemical
  equilibrium;
* transport limitations are not explicitly resolved.

These assumptions are appropriate for reduced-order thermal runaway
models but may not capture all mechanisms occurring inside a real
battery cell.


Relationship to Gas Generation
------------------------------

Reaction progress can also be used to drive gas-generation models.

In the current empirical framework, gas generation may be associated
with reaction progress through prescribed species yields.

Conceptually,

.. math::

   \dot{n}_s
   =
   Y_s
   m_{\mathrm{cell}}
   w
   \frac{d\alpha}{dt},

where :math:`Y_s` is the yield of gas species :math:`s` in mol/kg of
reacted material.

This direct-yield approach is intentionally simple.

A future thermochemical architecture may instead convert reaction
progress into elemental inventories and determine equilibrium gas
composition using a backend such as Cantera.


Implementation
--------------

The standard reaction model is

:class:`liiontr.reactions.reaction.Reaction`.

Multi-channel reactions are implemented by

:class:`liiontr.reactions.multichannel.MultiChannelReaction`.

The reaction-network container is

:class:`liiontr.reactions.network.ReactionNetwork`.

Context information is provided by

:class:`liiontr.reactions.context.ReactionContext`.

Context variables are implemented in

:mod:`liiontr.reactions.context_variable`.

The reaction-model interface is defined by

:class:`liiontr.reactions.model.ReactionModel`.

The reaction-network chemistry backend is

:class:`liiontr.chemistry.reaction_backend.ReactionNetworkBackend`.