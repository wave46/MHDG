# Focused unit tests

This directory contains small deterministic checks for numerical machinery that
can be exercised without launching a full MHDG case.  Production-equilibrium
and end-to-end coverage remains in `regression_tests`.

`make check-linear-algebra` checks the production BLAS/LAPACK wrappers using
rectangular products, pivoted single/multiple-RHS solves, inversion, LU factors
and real right-eigenvector residuals. It also checks cross-product orientation
and antisymmetry. It needs only `LinearAlgebra.o` and the configured BLAS/LAPACK
libraries; no solver run, mesh or MPI launch is needed.

`make check-bdf-coefficients` from `lib/` checks constant-step BDF1-6 polynomial
differentiation, startup order ramping, unused history and unchanged BDF1 weights.
It also checks the production history initialization/update routines at fresh
and mature restart indices, constant-state preservation and single-state reseeding.
It extracts the internal production coefficient routine verbatim and compiles it
with the real global time state; no solver run or external data is needed. This
target uses the current build configuration without a clean rebuild.

`make check-strong-dirichlet` uses the current 2D NGammaTiTe-family build.
It checks interleaved face ordering and matching frozen stabilization for
constant, variable diagonal and matrix tau. The stabilization and prescribed
contributions are extracted from production assembly. With PaStiX, it also
assembles and solves a small trace system with prescribed faces. The MPI
executable can run with `mpirun.openmpi -np 2 ./strong_dirichlet_check` to check
ghost metadata, remote owner indices, prescribed ghost faces and trace exchange.
It uses generated algebraic fixtures and needs no external case data.

`python3 unit_tests/check_bohm_momentum.py` from the repository root checks
the 2D temperature-model Bohm trace-density target, independent element/trace
derivatives, diagonal and full stabilization, and signed outgoing-flow selection.
It checks sound-speed sensitivity in both stabilization branches, grazing
incidence and suppression. It compiles extracted production code in a temporary directory
with gfortran; no solver build, mesh, MPI launch or external data is needed.

`make check-amjuel-extrapolation` uses the current 2D NGammaTiTe-family build
with `AMJUELSPLINES`. It compares the production rate Jacobian with finite
differences in density and electron energy below, inside, at and above the
temperature interval, with in-range and clamped densities. It also checks the
existing zero-rate cutoff. It needs no MPI launch or external case data.

`make check-energy-conductivity` checks the production ion/electron conductivity
values and derivatives at negative, zero, floored, ordinary and capped specific
energies. It covers active-branch derivatives at the floor/cap, the convergence
test's cap bypass, density/momentum couplings and degree-zero homogeneity.
It uses the current NGammaTiTe-family configuration without mesh or case data.

`make check-heat-flux-limiter` uses the current NGammaTiTe-family configuration.
It compares ion/electron limited-flux derivatives with finite differences and
checks reconstruction from the Jacobian and affine term in ordinary,
temperature-floor and cap-floor states, with either flux sign and zero flux.
It uses the production value/derivative helpers without MPI or external data.

`make check-neutral-gamma-transport` requires a NeutralGamma configuration.
It checks velocity, convection and pressure-curvature derivatives and flux
reconstruction at small positive and ordinary neutral densities, including both
convection modes, signed/zero neutral momentum, and physical/conservative
velocity conversion. No mesh or MPI launch is needed.

`make check-neutral-energy-source` uses a temperature model with
`NEUTRALGAMMA`. It checks that the production source matrix and affine term
recover the physical neutral kinetic-energy transfer at small positive densities,
including zero neutral momentum and state-dependent rate coefficients.
It extracts the internal source assembly verbatim and uses the production
energy-factor helpers; no MPI launch or external case data is needed.

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
  temporary directories; cases cover defaults, inheritance, pump overrides,
  range/nonfinite rejection, retired ITER-core BC 6 and testcases 80–89,
  supported source/circular presets, ignored legacy `diffmin` and core sources,
  `bohmtypebc=0` with rejection of the removed alternate Bohm selectors, and
  legacy `saveTau=.false.` with rejection of removed `saveTau=.true.` output.
  It also checks the k-build testcase restriction to 50–69 and its absence
  in non-k builds; the compiled test reports its build kind to the checker.
- `check-neutral-wall-diagnostics`: production accumulator/HDF5 writer; nine
  cases cover zero/active loss, modes, merges, signs, particle totals and
  separation from plasma recycling and relocated puff/pump.

All targets print `PASS` on success and fail with a nonzero exit status. They
need no external case data or solver simulation. Warm reference comparisons
belong to the regression harness: existing neutral workflows cover inactive
compatibility, and `neutral_wall` covers `Re_n=0.99` with inherited pump recycling
and with `Re_n_pump=0.95`. See [regression instructions](../regression_tests/README.md).

`make check-equilibrium-series MODE=parall` checks the isolated equilibrium
provider with small generated HDF5 inputs. It covers non-square orientation,
nonuniform physical timestamps, exact/intermediate sampling, jumps and reverse
sampling across more frames than fit in the cache, legacy and single-frame stationary equivalence,
optional references, coverage rejection, provider reinitialization and invalid
input shapes/nonfinite values. It needs h5py and NumPy, initializes no MPI and
does not run the plasma solver. The trajectory layout is documented in
[equilibrium input guide](../docs/inputs/equilibrium.md).
The same fixtures check Python conversion of stationary, numbered and paired
current inputs: exact values/orientation and attributes, untouched sources,
non-overwrite, late grid/time rejection without a partial output, and sampling
converted files with the Fortran reader.

The same target checks the [heating provider](../docs/inputs/heating.md#implementation-notes)
with its own non-square grid and timestamps: both power channels, forward and
reverse sampling, stationary/one-frame equivalence, endpoint coverage, cache
ownership, dimension inference and malformed/nonfinite inputs. These are raw
provider checks; solver integration is covered separately by `check-me-restart`.

For a Python-only converter change, reuse the already-built check executable
without running make:

```bash
python unit_tests/check_equilibrium_series.py lib/equilibrium_series_check
```

For a 2D `NGammaTiTeNeutral` build, `make check-me-restart` reads small synthetic
HDF5 checkpoints through the production input, initialization, and solution
readers. It checks static bootstrap initialization, real ME steps 0/1/50/99/100,
feedback-state restoration, rejection before out-of-range history access,
endpoint interpolation/heating selection, and inactive ME compatibility.
The same driver exercises production time-input parsing and seconds-to-solver
conversion, normalized-unit equivalence, the `tfi` alias, independent adaptation
options, and invalid/nonfinite configuration rejection. Physical initialization
uses `t_start` and allocates initial history capacity independently of `nts`.
Endpoint checks use the production clock helper for constant/nonuniform steps,
shortened final intervals, continuation between supplied frames and the last
timestamp roundoff case. History growth preserves accepted entries and does not
change their count. The solver contracts checked here are described in
[physical time stepping](../docs/simulation/time-stepping.md) and
[restart and output](../docs/simulation/restart-and-output.md).
Provider integration checks infer fixed/moving magnetic input without `ME`, sample
between frames independently of numerical counters, and use heating's own time
grid. They also rebuild current and heating values on changed mesh nodes from
the retained providers. Numbered heating is checked with its own timestamps,
including fixed magnetics, missing-time rejection and a lone fixed frame.
Numbered magnetic input remains covered alongside single-file trajectories.
Fixed-field source checks retain constants with empty paths, combine scheduled
puff/concentration with fixed Zeff, and sample the final control timestamp.
Existing density/recycling/XPR controllers are exercised with two physical
intervals and repeated evaluation from restored controller values. These checks
reuse the same synthetic input and executable to exercise control updates directly.
The controller check also exercises linear ramps through repeated timestep
halving and 1.3-times growth for puff/recycling/XPR, replays candidate updates
after restoring state, and compares uninterrupted versus saved/restored feedback
updates. Synthetic ramps cover the timestep-change spike without storing case traces.
The same fixtures check explicit restart initialization from steady/pseudo-time
states, fixed-field physical step-1 continuation, completed checkpoints, and
moving-field continuation between frames with a numerical index larger than the
frame count and initial history capacity. Clock conversion, feedback restoration,
BDF1 accepted history, next-endpoint selection, and nonfinite-time rejection are
checked through the production reader/helpers. Startup sampling also checks that
prescribed values are initialized without changing restored feedback state.
It also checks HDF5 grid dimensions/orientation on a non-square rectangular grid,
stored-current interpolation at initial/endpoint frames, full filenames and
legacy basenames, optional namelist dimensions, dataset-shape rejection and
missing required density-chord fields.
It also checks that a caller requesting optional reference Jtor is told when
the input contains none, without attempting interpolation of missing data.
The same cases load a small Gmsh mesh and checkpoint mesh/solution using stems
and full .msh/.h5 names, and verify that full input names preserve output stems.
The same driver checks physical output indices independently of save frequency,
initial and final names, and indices above 9999 while preserving legacy output
numbering and Newton diagnostic overwriting.
Control readers infer series lengths and preserve puff, Zeff and impurity
interpolation at initial, intermediate and final times, including repeated
loads; mismatched value/time lengths and invalid vector shapes are rejected.
The same target checks global nodal input on one rank and, for parallel builds,
four ranks with four OpenMP threads: reordered/shared nodes, stationary and ME
endpoint data, global flux extrema and axis coordinates, and mesh-size mismatch.
Parallel builds also gather signed nodal samples on one and four ranks, checking
shared-node values, exclusion of ghost-only copies, and all-rank/root-only output.
It needs Python 3 with h5py and NumPy, uses no external
equilibrium or case data, and does not launch a solver simulation. The nodal
checks initialize MPI; parallel builds need a working local MPI launcher.

For completed-restart executable checks with prepared external inputs, use:

```bash
python3 unit_tests/check_me_completed_restart.py --solver /path/to/solver \
  --parameters /path/to/param.txt --checkpoint /path/to/checkpoint.h5 \
  --ranks 4 --omp-threads 4
```

It checks legacy ME and explicit fixed/prescribed completion with deliberately
absent equilibrium files and no physical advances or new output files.

Legacy moving-equilibrium step-limit checks cover fresh and resumed advances
through the final global `nts` bound, zero-work completion and out-of-range
rejection. The [legacy frame/step contract](../docs/inputs/equilibrium.md#legacy-time-stepping)
is separate from explicit physical modes, which stop at `t_final`.

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

The magnetic-geometry checks cover all 16 COCOS identifiers, independent field/current
coefficient fits, and selection with missing references or conflicting metadata.
The existing series checks also cover COCOS preservation during conversion and sampling.


`make check-divergence-refinement` also checks reusable mesh transfer: P1–P3
polynomials on affine/curved meshes, solution and gradient layout, repeated
application, rebuilding for another destination mesh, and BDF-history transfer.
MPI builds additionally check owned source values and destination ghost points.
It uses generated meshes and the existing reference-node file, with no plasma
solve or external case data. Pass the same `MODE`, model and build options as
the current solver cache; no clean is needed when those options are unchanged.


`make check-linear-solve-status` checks collective failure propagation, PETSc KSP
reason classification without requiring PETSc, and (with PaStiX enabled) real
diagonal solves followed by NaN-matrix/Inf-RHS rejection and backend rebuilding.
It verifies that failure leaves the saved trace untouched and releases borrowed
matrix storage. Use the existing build configuration; the target runs one MPI
rank. The resulting `lib/linear_solve_status_check` can also run under MPI.
PETSc runtime/library-API coverage requires a PETSc build on an equipped machine.
