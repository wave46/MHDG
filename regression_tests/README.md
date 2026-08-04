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

The routine case is `legacy_case`.

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
| `cold_step_impurity_off` | Analytical start; coarse fixed mesh | Short disabled-impurity lifecycle check. |
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
| `impurity_scalar_baseline` | Impurity off and N, `mpi4_omp4` | Focused compatibility check. |
| `impurity_mixture` | Impurity off, W, N, and N+W, `mpi4_omp4` | Manual mixture-reference check. |
| `initialization_smoke` | Disabled-impurity analytical start, `mpi4_omp4` | Execution-only initialization evidence. |
| `race` | Both one-step workflows, `serial_omp1` vs `serial_omp16` | Routine OpenMP race check. |
| `cold` | Both full cold workflows, `mpi4_omp4` | Canonical integration check. |
| `warm_parallelism` | `warm`, all layouts | Periodic layout characterization. |
| `race_matrix` | Both one-step workflows, every pair of tracked layouts | Periodic race check. |
| `cold_matrix` | Both full cold workflows, all layouts and all layout pairs | Overnight golden and reproducibility evidence. |
| `stored_field_compatibility` | Limited two-Newton-step scratch start, stored `Br/Bz` | Compatibility smoke for `compute_from_flux=false`. |
| `diverted_warm` | Diverted warm restart, `mpi4_omp4` | Routine diverted topology/transport check. |
| `diverted_warm_parallelism` | Diverted warm restart, all layouts | Diverted layout characterization. |
| `diverted_cold_adaptive` | Full diverted adaptive cold start, `mpi4_omp4` | Canonical diverted reference producer. |
| `diverted_race_matrix` | Diverted two-step fixed/adaptive starts, every layout pair | Lightweight diverted race and cache-refresh evidence. |

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

# Every same-layout golden check and declared layout pair in a suite.
regression_tests/regression.sh suite compare \
  /path/to/suites/cold_matrix/overnight-01/suite_summary.json
```

### Publish accepted references

For a complete refresh, one command starts or continues the persisted golden
campaign and runs every producer and verification stage. It stops once, after
the complete candidate is ready for review:

```bash
regression_tests/regression.sh golden update legacy_case \
  --settings /private/path/source.env --run-id refresh-01 \
  --output /private/path/new_golden_bundle --bundle-version 2.5.0
```

Run the same command again with `--accept campaign` after reviewing the final
candidate. This second command only publishes; it does not rerun completed
stages. `golden status WORKSPACE` shows progress; without `--workspace`, the
workspace is `MHDG_REGRESSION_RUN_ROOT/golden_campaigns/RUN_ID`. Resumes reject
changed inputs and never overwrite a workspace, candidate, or output.

The default refresh includes cold matrices, warm and mixture references,
initialization/stored-field/race evidence, warm layout checks, and final
warm/mixture verification. Repeat
`--only` to select `cold_matrix`, `warm`, or `impurity_mixture`. A cold-matrix
refresh also updates the canonical warm restart, while `warm` updates the warm
reference. Updating only warm or only mixture references records a consistency
warning; select both together to avoid it.

The verified candidate is published atomically as a golden bundle. Campaign
state, declaration, build metadata, suite summaries, run plans/metadata, and
comparison reports are registered below `provenance/golden_campaign/`.

### PR03 limited and diverted refresh

PR03 keeps the exhaustive limited campaign and adds an independent diverted
campaign. Run them sequentially: each campaign performs clean serial and MPI
builds in the same worktree.

The limited source remains the last accepted `legacy_case` golden bundle. Its
full campaign runs the fixed/adaptive cold layout matrix, impurity references,
the stored-field compatibility smoke, warm layouts, the two-step race matrix,
and final verification. The tracked workflows render `compute_from_flux=true`;
the compatibility smoke alone renders it false.

The first diverted source is a candidate bundle, because no diverted golden
exists yet. Prepare the filenames declared by `cases/diverted_case.json`. All
transport namelists used by the cold and warm workflows must contain:

```text
transport_region_policy = 'core_and_main_sol'
c_bohm_n_rho_slope = 0.7
```

Use `puff = 5.0e21` in every staged parameter file and
`diff_n_min_phys = 0.1` in the final continuation and warm transport files.
The case workflows explicitly render `compute_from_flux=true`, even if a
prepared parameter file still contains false. Create the source candidate:

```bash
regression_tests/regression.sh bundle create \
  --case diverted_case \
  --source /private/path/prepared_diverted_case \
  --output /private/path/diverted_case_pr03_source \
  --bundle-version pr03-source.1
```

Point a private settings file at that candidate and start the resumable cold
overnight in `tmux` or another persistent shell:

```bash
regression_tests/regression.sh golden update diverted_case \
  --settings /private/path/diverted-source.env \
  --run-id pr03-diverted-refresh-01 \
  --workspace /private/path/campaigns/pr03-diverted-refresh-01 \
  --output /private/path/golden_bundles/diverted_case_pr03_golden \
  --bundle-version 1.0.0-pr03-golden.1 \
  --bootstrap-candidate \
  --build-jobs 8
```

The first invocation runs the full adaptive `mpi4_omp4` cold producer, maps its
final output to `warm_restart` and `warm_reference`, generates a reconverged
warm reference, runs the warm-layout and race checks, and finishes with final
warm verification. It then stops at `awaiting_acceptance` without publishing.
Inspect status and the reports below the workspace, then repeat the identical
update command with `--accept campaign`:

```bash
regression_tests/regression.sh golden status \
  /private/path/campaigns/pr03-diverted-refresh-01

# Add this to the identical `golden update` command:
--accept campaign
```

If a run is interrupted, repeat the identical command without `--accept`; the
recorded stage resumes. `--bootstrap-candidate` is only for the first diverted
promotion; future refreshes start from the accepted golden and omit it. Never
delete or reuse the workspace or output path.

If a stage records a failure, correct its tracked workflow or source bundle,
then repeat the identical command with `--retry-failed`. The failed attempt is
retained in campaign provenance and the corrected stage receives a distinct
`-retry-N` run ID. Completed producer stages are not rerun.

If the correction affects an earlier campaign declaration, use
`--retry-from STAGE` instead. The campaign records the declaration amendment,
restores the preceding candidate, archives superseded downstream attempts, and
uses distinct run, candidate, report, and settings paths for the replacements.

Refresh the limited golden with the same lifecycle but `legacy_case`, an
accepted golden source settings file, and distinct run/workspace/output names:

```bash
regression_tests/regression.sh golden update legacy_case \
  --settings /private/path/limited-golden-source.env \
  --run-id pr03-limited-refresh-01 \
  --workspace /private/path/campaigns/pr03-limited-refresh-01 \
  --output /private/path/golden_bundles/legacy_case_pr03_golden \
  --bundle-version 2.5.0-pr03-golden.1 \
  --build-jobs 8
```

It runs through the cold matrix, warm reference, impurity references, smoke and
race checks, and final verification in one invocation. Review the composed
candidate once, then publish it with the identical command plus
`--accept campaign`.

The `legacy_case` order was checked against the accepted PR 02 record: clean
builds at `7ce486f`, its full cold-matrix refresh, the passing disabled
initialization probe, final mixture verification, and the validated 118-artifact
`legacy_case_2.4.0-pr02-mixture-golden.1` bundle. Those results defined the
campaign shape; they were not rerun or treated as resumable campaign state
because the historical record does not retain every suite path and fingerprint.

For a smaller manually reviewed promotion, create a golden bundle directly:

```bash
regression_tests/regression.sh bundle promote \
  /path/to/cold_matrix/suite_summary.json \
  /path/to/warm/suite_summary.json \
  --settings /private/path/candidate-settings.env \
  --output /private/path/new_golden_bundle \
  --bundle-version 1.0.0-golden.1
```

Promotion validates and applies one or more accepted summaries in the given
order, including their references and provenance. It never overwrites an
output, and tests never promote automatically.

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

Adaptive golden references may come from a different mesh after an intentional
adaptivity change. Their comparison uses `HDG_postprocess` to interpolate both
solutions at deterministic interior points instead of requiring equal
connectivity.

Race probes require identical connectivity arrays, including node, element,
and boundary ordering. The generated adaptive `temp.msh` files must also be
byte-identical, which rejects tag swaps even when the physical topology is
unchanged. The full race matrix applies that check to every pair of tracked
layouts before comparing `u`, `q`, and `u_tilde` with the race tolerances.
Cold-matrix layout pairs likewise require exact final and retained meshes, then
compare the final HDF5 solution and transport data with the cold cross-layout
tolerances.

Golden matrices keep a reference for each workflow, layout, and stage. A
staged comparison stops at the first divergent stage. Race suites instead
compare layout pairs produced by the same build directly.

| Comparison | Relative L2 | Normalized Linf |
| --- | ---: | ---: |
| Warm, same layout | `1e-10` | `1e-9` |
| Warm, cross layout | `5e-8` | `1e-6` |
| One-step race probe | `5e-8` | `1e-6` |
| Matching fixed cold stage | `2e-7` | `3e-7` |
| Converged cold, cross layout | `3.5e-7` | `3e-7` |
| Fixed cold final state against warm reference | `1e-5` | `1e-5` |
| Adaptive solution | `0.05` | `0.1` |
| Adaptive gradient | `0.25` | `0.3` |

Normalized Linf divides the largest pointwise difference by the largest
absolute reference value. Fixed and race HDF5 coordinates use absolute
tolerance `1e-12`; generated adaptive mesh files are compared byte-for-byte.
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

- Adaptive refinement may cross different thresholds between code revisions;
  adaptive golden-reference comparisons therefore do not promise identical
  meshes. Layout pairs within one race or cold matrix do require exact meshes.
- Runtime is recorded without a timing threshold.
- Promotion verifies technical evidence but physical acceptance remains human.
- Promotion inputs must all describe accepted results from the configured source
  bundle.

Run `regression_tests/regression.sh help` for the command synopsis.
