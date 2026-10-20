# Solver documentation

Documentation is being consolidated here incrementally. Coverage is preliminary
and limited to the guides listed below.

## Models

- [Reduced 1D transport](models/transport-1d.md): configuration, Bohm /
  gyro-Bohm diffusion, pinch, scientific conventions and output.

## Inputs

- [Magnetic equilibrium input](inputs/equilibrium.md): fixed/moving equilibria,
  HDF5 formats, flux reconstruction, COCOS conventions and legacy conversion.
- [External heating input](inputs/heating.md): ion/electron power-density
  profiles, HDF5 formats, time sampling, normalization and output.
- [Source schedules and density feedback](inputs/controls.md): prescribed
  sources, measurement chords, controller selection, timing and restart history.

## Simulation

- [Simulation configuration](simulation/configuration.md): parameter templates,
  geometry and input preparation, initialization and fixed-diffusion scans.
- [Physical time stepping](simulation/time-stepping.md): mode selection,
  units, failure recovery, growth and optional Newton prediction.
- [Restart and solution output](simulation/restart-and-output.md): saved-solution
  initialization, physical continuation, output naming and checkpoint metadata.

## Numerical utilities

- [One-dimensional splines](numerics/splines.md): interpolation conventions,
  Fortran usage and current limitations of the retained utility.

## Other documentation

Installation, building and basic runs remain in the [root README](../README.md).
The [regression guide](../regression_tests/README.md), its operational tutorials,
and the [unit-test guide](../unit_tests/README.md) remain beside their harnesses.
