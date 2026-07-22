# MHDG regression tests

This directory contains the tracked test definitions and tools. Physical
inputs and results stay outside Git: meshes, equilibria, parameter files,
restarts, accepted solutions, logs, executables, and local paths belong in
external bundles and private run directories.

The harness builds MHDG, prepares isolated runs, executes and compares solver
workflows, resumes interrupted suites, and publishes reviewed golden data.

## Concepts

| Term | Meaning |
| --- | --- |
| **Case** | Input contract: available workflows and external file roles. |
| **Workflow** | One simulation procedure, such as a warm restart or staged cold start. |
| **Layout** | Executable type, MPI ranks, and OpenMP threads. |
| **Suite** | Selected workflows and layouts executed as one resumable job. |
| **Bundle** | Validated external files plus a generated manifest. A **golden** bundle contains accepted references. |

```text
case + workflow + layout -> run
suite                    -> selected runs and comparisons
accepted suite results   -> golden bundle
```

The routine case is `legacy_case`. `historical_feature` only archives earlier
feature evidence and is not a correctness oracle.

## First-time setup

Run commands from the repository root. Install the Python dependencies:

```bash
python -m pip install -r regression_tests/requirements.txt
```

Adaptive comparisons also need `hdg_postprocess`. Set
`PYTHON=/path/to/python` if the required environment is not `python3`.

Create the ignored local settings file:

```bash
cp regression_tests/settings.example.env regression_tests/golden.local.env
```

Fill in absolute paths for the external bundle, private run root, solver
executables, Open MPI launcher, and environment script. Executable paths are
only needed for prebuilt runs; `--build` generates them automatically.

Validate the configured bundle and run the smallest accepted-reference check:

```bash
regression_tests/regression.sh bundle validate \
  --settings regression_tests/golden.local.env

regression_tests/regression.sh suite check
```

`suite check` defaults to the `warm` suite and
`regression_tests/golden.local.env`. Override the file with `--settings FILE`
or `MHDG_REGRESSION_GOLDEN_SETTINGS`. It must point to a golden-class bundle.

## Test matrix

Warm workflows restart an existing solution; cold workflows start
analytically. Fixed and adaptive indicate whether the mesh can change.

| Workflow | Start and mesh | Work |
| --- | --- | --- |
| `warm` | Existing steady restart; fixed mesh | Reconverge the same state. |
| `cold_fixed` | Analytical start; refined fixed mesh | `time_init`, `diffusion_reduction`, then five continuations. |
| `cold_adaptive` | Analytical start; coarse mesh | Same seven stages; adapt in the first two. |
| `cold_step_fixed` | Analytical start; coarse fixed mesh | One time step and two Newton iterations. |
| `cold_step_adaptive` | Analytical start; coarse adaptive mesh | One time step, two Newton iterations, and one adaptation pass. |

Full cold stages run sequentially and restart from their predecessor. The
one-step workflows probe mesh construction and races, not convergence.

Tracked layouts are:

| Layout | Execution |
| --- | --- |
| `serial_omp1` | Serial, one OpenMP thread. |
| `serial_omp16` | Serial, sixteen OpenMP threads. |
| `mpi4_omp1` | Four MPI ranks, one thread each. |
| `mpi4_omp4` | Four MPI ranks, four threads each; canonical layout. |

Identifiers matching `serial_ompN` or `mpiM_ompN` determine the executable,
ranks, and threads. MPI runs bind each rank to exclusive cores.

| Suite | Coverage | Use |
| --- | --- | --- |
| `warm` | `warm`, `mpi4_omp4` | Fast routine golden check. |
| `race` | Both one-step workflows, `serial_omp1` vs `serial_omp16` | Routine OpenMP race check. |
| `cold` | Both full cold workflows, `mpi4_omp4` | Canonical integration check. |
| `warm_parallelism` | `warm`, all layouts | Periodic layout characterization. |
| `race_matrix` | Both one-step workflows, serial and MPI pairs | Periodic race check. |
| `cold_matrix` | Both full cold workflows, all layouts | Overnight evidence. |

Warm and race suites are short. Full cold workflows are longer, and
`cold_matrix` can take hours. Runtime is recorded but is not a pass criterion.

## Common tasks

### Check accepted references

```bash
# Existing executables; warm suite by default.
regression_tests/regression.sh suite check

# Build clean serial and MPI executables first.
regression_tests/regression.sh suite check --build --build-jobs 8

# Canonical fixed and adaptive cold workflows.
regression_tests/regression.sh suite check cold --build --build-jobs 8
```

### Run candidate or race evidence

`suite run` accepts candidate or golden bundles and requires explicit settings:

```bash
# Direct serial_omp1 versus serial_omp16 race comparison; no golden output.
regression_tests/regression.sh suite run race \
  --settings /private/path/settings.env

# Save an expensive matrix before comparison or acceptance.
regression_tests/regression.sh suite run cold_matrix \
  --settings /private/path/settings.env \
  --run-only --run-id overnight-01
```

Resume an interrupted matrix with the same suite, settings, and run ID:

```bash
regression_tests/regression.sh suite run cold_matrix \
  --settings /private/path/settings.env \
  --run-only --run-id overnight-01 --resume
```

Recorded cells are skipped. An incomplete run is preserved and retried as
`RUN_ID-resume-N`. Resume rejects changed settings, bundle data, catalogs,
executables, or launcher. `--build` cannot be combined with `--resume`; after
an initial `--build`, use the generated `settings.env` printed by that build.

### Build reusable executables

```bash
regression_tests/regression.sh build \
  --settings /private/path/settings.env --jobs 8
```

Serial and MPI builds run sequentially because they share objects. The command
prints the generated settings path and records commands, logs, Git state,
toolchain versions, environment checksum, and executable checksums.

### Inspect or run one workflow

Preparation validates and renders an isolated run but does not launch MHDG:

```bash
regression_tests/regression.sh prepare legacy_case cold_step_adaptive \
  --layout serial_omp16 --settings /private/path/settings.env
```

Execute one or more workflows directly when debugging:

```bash
regression_tests/regression.sh run legacy_case cold_fixed cold_adaptive \
  --layout mpi4_omp4 --run-id investigation-01 \
  --settings /private/path/settings.env
```

Prefer suites for routine work because they preserve one resumable summary.

### Recompare saved results

Neither command below launches MHDG:

```bash
# One completed run; policy comes from run_plan.json.
regression_tests/regression.sh compare /path/to/completed/run

# Every recorded cell or layout pair in a suite.
regression_tests/regression.sh suite compare \
  /path/to/suites/cold_matrix/overnight-01/suite_summary.json
```

### Publish accepted references

After human review, create a new golden bundle:

```bash
regression_tests/regression.sh bundle promote \
  /path/to/suite_summary.json \
  --settings /private/path/candidate-settings.env \
  --output /private/path/new_golden_bundle \
  --bundle-version 1.0.0-golden.1
```

Promotion validates and copies references and provenance. It never overwrites
an output, and tests never promote automatically.

## Outputs and provenance

Runs are stored below `MHDG_REGRESSION_RUN_ROOT`. Builds use
`MHDG_REGRESSION_BUILD_ROOT`, which defaults to `RUN_ROOT/builds`:

```text
runs/
├── builds/.../                    executables, logs, metadata, settings.env
├── suites/SUITE/RUN_ID/           suite_summary.json
└── CASE/WORKFLOW/LAYOUT/RUN_ID/
    ├── run_plan.json              requested inputs and comparison policy
    ├── run_metadata.json          execution result and provenance
    ├── stdout.log / stderr.log
    ├── outputs/
    └── stages/...                 staged-workflow runs
```

Comparisons write `comparison.json`, `matrix_comparison.json`, or
`verification_summary.json`. Suite summaries are updated after every cell.

Each newly built solver writes automatic compile-time identity into its HDF5
solutions:

```text
/provenance/git_commit
/provenance/git_dirty
/provenance/build_id
```

Do not edit these fields. One regression build ID is shared by its serial and
MPI executables. Candidate and golden commit hashes normally differ, so this
provenance is informational rather than a comparison criterion.

## External bundles

The V2 `manifest.json` is generated. Prepare external files using the
filenames declared in `regression_tests/cases/CASE.json`, then create a
candidate bundle:

```bash
regression_tests/regression.sh bundle create \
  --case legacy_case \
  --source /private/path/prepared_legacy_case \
  --output /private/path/candidate_bundle
```

Creation copies files under `inputs/` and generates role mappings, media
types, sizes, and SHA-256 checksums. It refuses to replace an output. Candidate
and golden bundles share this contract; their class records acceptance.

`positionFeketeNodesTri2D.h5` must be beside each executable. It is generic
solver runtime data, not case-specific bundle data.

| Edit directly | Generated; do not edit |
| --- | --- |
| Private settings and prepared physical files | Bundle manifest, identifiers, sizes, checksums |
| `cases/*.json` | Rendered parameter files, input links, run/build metadata |
| `suites.json`, `layouts.json`, `tolerances.json` | Summaries, comparison reports, HDF5 provenance |

## Adding a parameter variant

Common variants require JSON, not Python. In `cases/CASE.json`, inherit the
closest workflow and override only changed parameters:

```json
"cold_step_fixed_two_nr": {
  "extends": "cold_step_fixed",
  "description": "Repeat the fixed coarse-mesh probe with two Newton iterations",
  "parameter_overrides": {
    "nrp": 2
  }
}
```

This inherits inputs, mesh, stages, layout, comparison policy, and the other
parameter overrides. Add the new workflow identifier to an existing or new
suite in `suites.json`.

Overrides accept booleans, finite numbers, and strings. The named assignment
must occur exactly once in the selected MHDG parameter file. Preparation
formats it in a private copy and records the effective values in
`run_plan.json`; the bundle is unchanged.

Workflow overrides apply to every stage. Stage overrides take precedence for
that stage. With `extends`, `parameter_overrides` merge with the parent;
another supplied field replaces the complete parent field.

For a variant requiring another physical file:

1. Add a generic role and filename to the case catalog.
2. Put the file in the external prepared directory and reference its role.
3. Create and validate a candidate bundle; checksums are automatic.
4. Run the smallest relevant suite and review it before promotion.

Parameter values, matching layout identifiers, suite selections, and
tolerance profiles need no Python. Python is reserved for new preparation
behavior, workflow kinds, or comparison algorithms. Keep physical values,
private paths, machine names, and credentials out of tracked identifiers.

## Comparison behavior

Fixed-mesh comparison checks finite values and Newton error, exact
connectivity, tolerance-based coordinates, each equation in `u`, `q`, and
`u_tilde`, and transport-1D data when present.

Adaptive meshes may differ across layouts. Adaptive comparison uses
`HDG_postprocess` to interpolate both solutions at deterministic interior
points instead of requiring equal connectivity.

The adaptive race probes intentionally require identical connectivity. They
currently expose a known OpenMP defect: repeated multi-threaded runs can
build different meshes, while the corresponding one-thread runs reproduce
exactly. Fixed-mesh race probes pass. Keep this diagnostic failure visible
until the adaptation path is made deterministic.

Golden matrices keep a reference for each workflow, layout, and stage. A
staged comparison stops at the first divergent stage. Race suites instead
compare layout pairs produced by the same build directly.

| Comparison | Relative L2 | Normalized Linf |
| --- | ---: | ---: |
| Warm, same layout | `1e-10` | `1e-9` |
| Warm, cross layout | `5e-8` | `1e-6` |
| One-step race probe | `5e-8` | `1e-6` |
| Matching fixed cold stage | `2e-7` | `3e-7` |
| Fixed cold final state against warm reference | `1e-5` | `1e-5` |
| Adaptive solution | `0.05` | `0.1` |
| Adaptive gradient | `0.25` | `0.3` |

Normalized Linf divides the largest pointwise difference by the largest
absolute reference value. Fixed coordinates use absolute tolerance `1e-12`.
Except for transient initialization and race probes, final Newton error must
not exceed `2e-4`. These are
regression limits for `legacy_case`, not physical-accuracy targets.

## Focused synthetic tests

These tests use temporary bundles, small arrays, and fake executables. They do
not launch MHDG or need physical data.

| Changed area | Test group |
| --- | --- |
| Build | `tests.test_build` |
| Bundles and case loading | `tests.test_bundles` |
| Preparation and parameters | `tests.test_preparation` |
| Execution and run metadata | `tests.test_execution` |
| Fixed comparison | `tests.test_fixed_comparison` |
| Adaptive comparison | `tests.test_adaptive_comparison` |
| Reference matrices | `tests.test_matrix_comparison` |
| Suites, layout pairs, resume | `tests.test_suites` |
| Promotion | `tests.test_reference_publication` |

Run only the affected group:

```bash
cd regression_tests
PYTHONPATH=tools python -m unittest tests.test_preparation
```

For a cross-cutting change, run all synthetic tests:

```bash
cd regression_tests
PYTHONPATH=tools python -m unittest discover -s tests -p 'test_*.py'
```

## Limitations

- Adaptive refinement may cross different thresholds after reduction-order
  changes; identical adapted meshes are not promised.
- Runtime is recorded without a timing threshold.
- Promotion verifies technical evidence but physical acceptance remains human.
- One promotion consumes one accepted suite summary.

Run `regression_tests/regression.sh help` for the command synopsis.
