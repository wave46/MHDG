# Focused unit tests

This directory contains small deterministic checks for numerical machinery that
can be exercised without launching a full MHDG case.  Production-equilibrium
and end-to-end coverage remains in `regression_tests`.

Run the magnetic-topology checks from `lib/` after loading the build
environment:

```bash
source Make.inc/init_vars_libs.sh
make check-magnetic-topology
make check-magnetic-geometry
make check-transport-region-policy
make check-transport-taper
make check-transport-pinch
```

The topology target also checks that reversing the psi convention preserves the
outward normalized-flux direction.  The second target checks all four supported
poloidal-field conventions, noncanonical and percent-level field perturbations,
the same-alpha Grad--Shafranov toroidal current and its fitted `+1`/`-1` sign
convention, exact
nodal/quadrature topology caches, and cache generation refresh.  The
region-policy target checks
parsing, region inclusion, signed pinch orientation, null-normal suppression,
and the explicit legacy magnetic-field fallback.  Test executables are generated
in `lib/` and removed by `make clean`.  The taper target checks the disabled
slope, the default `0.7` factor beyond the LCFS, and lower-bound clipping.
The pinch target checks production namelist selection, particle-only defaults,
mixed/all/disabled selections, duplicate and zero entries, invalid IDs and
overlong lists, signed velocity normalization, untouched neutral rows, radial
and topology suppression, configuration preservation across profile
reinitialization, reset defaults, and HDF5 selection metadata. It links the
current NGammaTiTe-family solver objects, needs Python 3 with h5py and NumPy,
and does not launch a solver simulation or require external case data.

For all neutral unit checks, use the serial 2D
`NGammaTiTeNeutral` or `NGammaTiTeNeutralGamma` build configuration:

```bash
cd lib
source Make.inc/init_vars_libs.sh
make check-neutral-flux-limiter
make check-neutral-wall-flux
make check-neutral-wall-input
make check-neutral-wall-diagnostics
```

Keep the build cache for the same configuration; run `make clean` before
switching model or build configuration. The input and diagnostics drivers use
serial solver objects and do not initialize MPI. They need Python 3;
the diagnostics check also needs h5py and NumPy.

- `check-neutral-flux-limiter`: diffusion/pressure flux, perpendicular projection,
  supplied-speed cap, epsilon regularization and diffusion/cap floors.
- `check-neutral-wall-flux`: production Maxwellian normalization, limited Ti,
  finite-difference derivatives (`2e-7` scaled tolerance), Euler homogeneity and
  common scaling (`5e-13`), five/six-component states, small positive plasma
  density, zero neutral density and inactive `Re_n=1`. It also checks Ti/fixed
  limiter speed selection and wall independence from that selection.
- `check-neutral-wall-input`: production input/adimensionalization calls in
  temporary directories; 15 cases cover defaults, inheritance, pump overrides
  and range/nonfinite rejection.
- `check-neutral-wall-diagnostics`: production accumulator/HDF5 writer; nine
  cases cover zero/active loss, modes, merges, signs, particle totals and
  separation from plasma recycling and relocated puff/pump.

All targets print `PASS` on success and fail with a nonzero exit status. They
need no external case data or solver simulation. Warm reference comparisons
belong to the regression harness: existing neutral workflows cover inactive
compatibility, and `neutral_wall` covers `Re_n=0.99` with inherited pump recycling
and with `Re_n_pump=0.95`. See [regression instructions](../regression_tests/README.md).
