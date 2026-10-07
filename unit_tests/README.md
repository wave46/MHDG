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

For a 2D `NGammaTiTeNeutral` build, `make check-me-restart` reads small synthetic
HDF5 checkpoints through the production input, initialization, and solution
readers. It checks static bootstrap initialization, real ME steps 0/1/50/99/100,
feedback-state restoration, rejection before out-of-range history access, and
inactive ME compatibility. It needs Python 3 with h5py and NumPy, uses no external
equilibrium or case data, and does not launch a solver simulation.

In ME runs, `nts` is the final global step number: with `nts=100`, restart steps
50, 99, and 100 leave 50, 1, and 0 advances. The magnetic filename remains
`equilibrium_(step+1)` when another advance remains. With nts=100, restarting
at 99 loads frame 0100, and restarting at 100 reads no magnetic files. The
final advance skips the next field/control update. A restart beyond
the final step is invalid. The non-ME step-count behavior is unchanged.

For the executable loop and MPI exit paths, an optional focused check consumes
a prepared external ME case (the compact synthetic case has the required layout):

```bash
python3 unit_tests/check_me_step_limit.py --case-directory /path/to/prepared/case \
  --solver /path/to/MHDG-NGammaTiTeNeutral-parall-2D \
  --run-root /tmp/new-me-step-check --mpi-ranks 2
```

Load the solver build environment first. The script leaves the source case
unchanged, freezes its first magnetic pair and control values into new fixtures,
and checks fresh/resumed advances, final-step zero-work MPI exit, beyond-limit
rejection, and uninterrupted versus resumed numbered checkpoints with prescribed
puff and density feedback. It compares solution and feedback arrays, and saves
logs/reports in the new run directory. The step-50 fixture advances to step 52
to exercise the same bound with only two solves; step 99 advances to 100.
This verifies step machinery, not the moving-equilibrium physics trajectory.
