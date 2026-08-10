# ADR-0001: Problem-oriented simulation architecture

## Status

Accepted

## Date

2026-08-04


# Context

LiionTR aims to become a research framework for the simulation of
lithium-ion battery thermal runaway.

The framework must support:

- cylindrical cells (18650, 21700, 4680);
- pouch cells;
- prismatic cells;
- ARC experiments;
- DSC experiments;
- gas generation;
- venting;
- thermal propagation between cells;
- future extension to 1D, 2D and 3D models.


A simple battery-centric architecture is not sufficient because the
same physical model may apply to different battery geometries.


# Decision

LiionTR adopts a problem-oriented architecture.


The main abstraction becomes:

