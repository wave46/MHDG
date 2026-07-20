# MHDG regression tests

This directory contains the tracked, non-sensitive regression-test contract.
Meshes, equilibria, parameter files, restarts, reference solutions, logs, and
machine-specific paths remain in external bundles and private run directories.

The harness can build the solver, execute isolated warm or staged cold runs,
compare their HDF5 outputs, collect accepted references, and check a new build
against a golden bundle.

## Quick start with a golden bundle

Install the Python dependencies once:

```bash
python -m pip install -r regression_tests/requirements.txt
```

Copy `settings.example.env` to the ignored `golden.local.env` and fill in its
absolute data, run, executable, launcher, and environment-script paths. Then
run from the repository root:

```bash
# Lightweight canonical warm restart.
regression_tests/regression.sh golden-check

# Build optimized serial and parallel executables first.
regression_tests/regression.sh golden-check --build

# Fixed and adaptive cold workflows in the canonical 4x4 layout.
regression_tests/regression.sh golden-check cold --build

# Periodic full workflow-by-layout matrix.
regression_tests/regression.sh golden-check cold_matrix --build
```

Set `MHDG_REGRESSION_GOLDEN_SETTINGS=/path/to/settings.env` or pass
`--settings FILE` to use another private settings file. `golden.local.env` is
ignored by Git.

## Cases, workflows, and suites

The tracked case is `legacy_case`, pinned to `develop` commit
`29f442db67bac169b2616f2cbe399289d137993c`. Its external-data identifier
remains `legacy_fixed` so bundles created before the rename still work.

The archived `historical_feature` case records evidence from
`feautre/neutrals_pressure`; it is not a routine run or a correctness oracle.

| Workflow | Purpose |
| --- | --- |
| `warm` | Reconverge an unchanged steady restart on a fixed mesh. |
| `cold_fixed` | Start analytically on the refined mesh and complete seven continuation stages without adaptivity. |
| `cold_adaptive` | Start analytically on the coarse mesh, adapt during the first two stages, then complete the same continuations. |

Both cold workflows run `time_init`, `diffusion_reduction`, and five numbered
continuations. Each stage waits for its predecessor and uses its selected HDF5
output as the next restart.

| Suite | Workflows and layouts | Intended use |
| --- | --- | --- |
| `warm` | `warm`, `mpi4_omp4` | Fast routine check. |
| `warm_parallelism` | `warm`, all five layouts | MPI/OpenMP characterization. |
| `cold` | Both cold workflows, `mpi4_omp4` | Routine feature-integration check. |
| `cold_matrix` | Both cold workflows, all five layouts | Periodic or overnight evidence. |

The characterized warm run takes seconds. The canonical cold suite is a
longer routine check, while the initial ten-cell cold matrix took about six
hours on the reference workstation. Runtime is recorded but is not currently
a pass/fail metric.

Tracked layouts are `serial_omp1`, `mpi2_omp1`, `mpi2_omp4`, `mpi4_omp1`, and
`mpi4_omp4`. MPI layouts use Open MPI core binding and assign exclusive cores
to each rank's OpenMP threads.

## External bundle contract

A bundle is a self-contained directory with a validated `manifest.json`. The
manifest records each artifact once using a generic role, bundle-relative
path, size, media type, and SHA-256 checksum. Checksums are generated and
checked automatically; users do not write them by hand.

The basic prepared directory uses these generic names:

```text
legacy_case/
├── param.txt
├── mesh.msh
├── geometry.geo
├── transport_model.nml
├── equilibrium.h5
├── current_density.h5
├── restart.h5
└── reference_mpi4_omp4.h5
```

Cold workflows additionally use:

```text
mesh_adaptive_initial.msh
param_cold_fixed_time_init.txt
param_cold_fixed_diffusion_reduction.txt
param_cold_fixed_continuation_01.txt ... continuation_05.txt
transport_cold_fixed_initial.nml
transport_cold_fixed_continuation_01.nml ... continuation_05.nml
```

The prepared entries may be files or symlinks to a private archive. Bundle
creation copies their contents, producing a portable bundle without external
symlinks. Physical names and paths must never appear in tracked files.

Create and validate a candidate bundle with:

```bash
regression_tests/regression.sh bundle create \
  --case legacy_case \
  --source /private/path/legacy_case \
  --output /private/path/candidate_bundle

regression_tests/regression.sh \
  --settings /private/path/candidate-settings.env \
  check-data
```

`bundle create` refuses to replace an existing output. A `candidate` and a
`golden` bundle use the same data contract; the class records whether its
references have been accepted for regression checks.

## Local settings and builds

Start from `settings.example.env`. The important settings are:

```text
MHDG_REGRESSION_DATA_ROOT=/absolute/path/to/bundle
MHDG_REGRESSION_RUN_ROOT=/absolute/path/to/scratch
MHDG_SERIAL_EXECUTABLE=/absolute/path/to/serial/solver
MHDG_PARALLEL_EXECUTABLE=/absolute/path/to/parallel/solver
MHDG_MPI_LAUNCHER=mpirun.openmpi
MHDG_ENVIRONMENT_SCRIPT=/absolute/path/to/lib/Make.inc/init_vars_libs.sh
```

Build clean optimized serial and parallel executables with:

```bash
regression_tests/regression.sh \
  --settings /private/path/settings.env \
  build --jobs 8
```

The two builds are sequential because they share `.o` and `.mod` files.
Their commands, logs, Git revision and dirty state, toolchain versions,
environment-script checksum, executable checksums, and generated settings are
stored under the private build root. `--build` on a suite performs this step
and uses the generated executables automatically.

Build and run provenance is stored in JSON sidecars. Writing solver commit or
build metadata directly into HDF5 remains deferred.

## Preparing and executing runs

Preparation creates an isolated run directory and renders private paths into a
copy of `param.txt`; it never changes the bundle:

```bash
regression_tests/regression.sh \
  --settings /private/path/settings.env \
  prepare legacy_case warm --layout mpi4_omp4
```

Use `run` to prepare and execute:

```bash
regression_tests/regression.sh \
  --settings /private/path/settings.env \
  run legacy_case warm --layout mpi4_omp4

regression_tests/regression.sh \
  --settings /private/path/settings.env \
  run legacy_case cold_fixed cold_adaptive \
  --layout mpi4_omp4 --run-id cold-01
```

Run directories contain read-only input links, writable `outputs/` and `res/`
directories, rendered parameters, `stdout.log`, `stderr.log`, `run_plan.json`,
and `run_metadata.json`. Staged workflows contain one isolated directory per
stage. A failed stage stops that workflow, while a suite continues with its
next workflow/layout cell.

The generic `positionFeketeNodesTri2D.h5` file must be beside the selected
executable. It is solver runtime data, not physical case data.

## Comparison contract

Fixed-mesh comparison reads HDF5 directly and supports both current grouped
files and older flat `develop` files. It checks:

- final Newton error and finite values;
- exact mesh connectivity and tolerance-based coordinates;
- per-equation `u`, `q`, and `u_tilde` errors;
- reference transport-1D coefficients and profiles when present.

Adaptive comparison uses `HDG_postprocess` to interpolate both solutions at
deterministic interior points. The selected Python environment must therefore
be able to import `hdg_postprocess`; use `PYTHON` and, when necessary,
`PYTHONPATH` to select it.

A golden matrix stores a reference for every accepted
`(workflow, layout, stage)` tuple. New staged runs compare each stage with the
matching tuple, stop comparison at the first divergence, and write a compact
`matrix_comparison.json`. Warm runs and bundles without a stage matrix retain
the final-state comparison path.

The initial numerical limits in `tolerances.json` are:

| Comparison | Relative L2 | Normalized Linf |
| --- | ---: | ---: |
| Warm, same layout | `1e-10` | `1e-9` |
| Warm, cross layout | `5e-8` | `1e-6` |
| Matching fixed cold stage | `1e-8` | `1e-7` |
| Fixed cold fallback against warm reference | `1e-5` | `1e-5` |
| Adaptive solution | `0.05` | `0.1` |
| Adaptive gradient | `0.25` | `0.2` |

All profiles require final Newton error at most `2e-4`. Fixed-mesh coordinates
use absolute tolerance `1e-12`. Normalized Linf is the largest absolute
pointwise difference divided by the largest absolute reference value.

These are regression tolerances for the characterized legacy workflow, not
general physical-accuracy targets.

Compare any completed run manually with:

```bash
regression_tests/regression.sh compare /path/to/completed/run
```

The command reads the run plan and selects fixed-mesh, adaptive-mesh, or staged
reference-matrix comparison automatically.

The command prints a short summary and writes `comparison.json`, or
`matrix_comparison.json` for a staged reference matrix.

## Suites, saved verification, and promotion

Run a suite against a candidate bundle with explicit settings:

```bash
regression_tests/regression.sh \
  --settings /private/path/settings.env \
  suite warm
```

For the expensive matrix, save all results before deciding on references:

```bash
regression_tests/regression.sh \
  --settings /private/path/settings.env \
  suite cold_matrix --run-only --run-id overnight-01
```

The suite summary is updated after every cell. Reuse the same command with
`--resume` after an interruption. Re-run comparison logic without launching
the solver using:

```bash
regression_tests/regression.sh suite-verify /path/to/suite_summary.json
```

Promotion is always explicit and never overwrites an existing bundle:

```bash
regression_tests/regression.sh \
  --settings /private/path/candidate-settings.env \
  bundle promote /path/to/suite_summary.json \
  --output /private/path/golden_bundle \
  --bundle-version 1.0.0-golden.1
```

Promoting a canonical warm suite replaces the warm reference and records its
run evidence. Promoting a deferred cold matrix copies every stage output into
a reference matrix and records its suite/run provenance. Promotion validates
bundle integrity and recorded run completion, but physical acceptance remains
a human decision. Ordinary tests never trigger promotion.

The short `golden-check` interface and the explicit `--settings FILE suite`
interface share the same runner and comparator.

## Adding cases or parameter variants

The current extension points are deliberately generic:

1. Put new private parameter or model files in the external prepared data.
2. Give them generic artifact roles in a tracked case definition.
3. Define the workflow stages, layouts, tolerance profile, and suite selection.
4. Run and review a candidate suite before promoting new references.

Do not encode machine names, experiment identifiers, physical values, or
private paths in tracked identifiers. The follow-up harness audit will focus
on making common parameter variants declarative, reducing the amount of schema
and Python knowledge required from a new user.

## Known limitations and follow-up

- Adaptive meshes can diverge after OpenMP-dependent reduction ordering crosses
  refinement thresholds. Mesh-independent comparison is required; identical
  adapted meshes are not promised.
- The legacy `a_minor` estimate is mesh-sensitive and belongs to the
  Bohm/gyro-Bohm follow-up, not this harness.
- Selected conservation and balance diagnostics belong to the diagnostics PR.
- Historical feature outputs remain reference evidence rather than golden
  truth.
- Golden promotion currently treats one accepted summary at a time.
- A follow-up audit will minimize routine suites, simplify the harness tests and
  fixtures, split multi-purpose functions, reduce unstructured dictionaries,
  and make new parameter variants easier to add.

Run `regression_tests/regression.sh help` for the compact command reference.
