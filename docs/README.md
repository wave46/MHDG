# Solver documentation

Documentation is being consolidated here incrementally. Coverage is preliminary
and limited to the guides listed below.

## Models

- [Physics cutoff conventions](models/physics-cutoffs.md): floors, caps,
  invalid-state extensions and the adopted Ohmic resistivity rule.
- [Parallel heat-flux limiter](models/heat-flux-limiter.md): free-streaming
  caps, floors and consistent flux linearization.
- [AMJUEL rates and source linearization](models/amjuel-rates.md): physical
  coefficients, continuation, derivatives and their source-assembly roles.
- [Wall recycling: physical law and HDG implementation](models/wall-recycling.md):
  element/trace states, recycled diffusion, stabilization and balance interpretation.
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

- [Solver flow](simulation/solver-flow.md): initialization, Newton assembly and
  solve, acceptance, recovery and outer continuation.
- [Simulation configuration](simulation/configuration.md): parameter templates,
  geometry and input preparation, initialization and fixed-diffusion scans.
- [Physical time stepping](simulation/time-stepping.md): mode selection,
  units, failure recovery, growth and optional Newton prediction.
- [Restart and solution output](simulation/restart-and-output.md): saved-solution
  initialization, physical continuation, output naming and checkpoint metadata.

## Numerical utilities

- [Mesh adaptivity](numerics/adaptivity.md): optional element-count limit and
  its interaction with Newton recovery and accepted timesteps.

- [Newton formulation](numerics/newton-formulation.md): equivalent absolute and
  incremental solves, local elimination and boundary handling.
- [HDG residuals and Newton convergence](numerics/hdg-residuals.md): residual
  evaluation, normalization, reporting and experimental stopping criteria.
- [Diffusion in HDG](numerics/diffusion.md): element/trace coefficient states,
  gradients and lagged stabilization.
- [One-dimensional splines](numerics/splines.md): interpolation conventions,
  Fortran usage and current limitations of the retained utility.

## Development

- [Intentional exact real comparisons](development/real-comparisons.md): reviewed
  compiler warnings, their numerical and input-contract rationale, and review policy.

## Other documentation

Installation, building and basic runs remain in the [root README](../README.md).
The [regression guide](../regression_tests/README.md), its operational tutorials,
and the [unit-test guide](../unit_tests/README.md) remain beside their harnesses.
