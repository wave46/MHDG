# MHDG regression tests

`python -m regression_tests` is the entry point for the MHDG scientific regression
harness. It builds selected solver variants, runs real limited and diverted cases in
isolated directories, checks convergence and output identity, and compares results
with accepted solutions or selected parallel layouts. Physical inputs, executables,
runs and goldens stay outside Git.

## Quick start: one warm check

From the repository root, install the Python requirements:

```bash
python -m pip install -r regression_tests/requirements.txt
```

Create the ignored `regression_tests/settings.local.json` with paths to an
existing build manifest and the accepted diverted golden:

```json
{
  "defaults": {
    "build": "/path/to/build/build_metadata.json",
    "bundles": {"diverted_case": "/path/to/diverted/golden_bundle"}
  }
}
```

```bash
python -m regression_tests check warm
```

This runs the five-equation 2D `NGammaTiTeNeutral` model. It restarts the
accepted diverted bootstrap solution on its saved mesh and reconverges it to
steady state with transport 1D off, W radiation on and diffusion at 16 m²/s.
The harness checks output identity and Newton convergence, then
compares the final mesh and solution fields with the accepted diverted warm
reference. For a read-only setup check before launching MHDG, use
`python -m regression_tests doctor warm`. Results are written under
`~/.cache/mhdg-regression` unless you set `run_root` as described below.

## How checks are organized

| Term | Meaning in this harness |
| --- | --- |
| **Case** | A real physical setup, such as limited or diverted, that names its input roles and available workflows. The files stay outside Git. |
| **Workflow** | A solver recipe: starting state, inputs, parameter changes, launches, expected outputs and comparison policy. |
| **Layout** | Executable type, MPI ranks and OpenMP threads for a launch. |
| **Suite** | Workflows for one case, their layouts, reference checks and selected cross-run comparisons. |
| **Profile** | An ordered selection of suites, possibly spanning cases. |
| **Bundle** | External input files mapped to case roles, with a manifest and checksums. |

The current catalog covers 2D only. Both real topologies use the five-equation
`NGammaTiTeNeutral` model for their warm, cold and feature checks. The limited
case also has one focused six-equation `NGammaTiTeNeutralGamma` check: an
analytical cold step with two Newton iterations, OpenMP agreement and detailed
diagnostics. That probe has no NeutralGamma golden reference. Other models and
3D are outside the current harness.

For example, `diverted_case` + `baseline_warm` + `mpi4_omp4` specifies one
workflow run; the `warm` suite checks it against the diverted golden. The terms
used in the coverage table mean:

| Term | Meaning in this harness |
| --- | --- |
| **Warm / cold** | Warm restarts a saved plasma solution; cold initializes plasma without one. A cold *bootstrap* uses successive launches, or stages, to reach a converged state. |
| **Golden** | A published bundle containing accepted reference results for regression comparison. |
| **Fixed / adaptive mesh** | Fixed keeps the selected mesh; adaptive permits refinement. Matching meshes are compared directly; different meshes use sampled comparisons. |
| **Feature checks** | Selected physics changes compared with their own references:<br>• **Neutrals:** source relocation, pressure, perpendicular diffusion and flux limiters.<br>• **Transport 1D:** activation followed by staged diffusion-floor reduction.<br>• **Impurities:** radiation off, N and N+W radiation. |
| **Diagnostics** | Optional balance output checked for presence, consistency and, in selected parallel suites, agreement of totals. |

## Scientific coverage

| Selection | Purpose |
| --- | --- |
| `routine` | Diverted warm golden and short adaptive/source golden; diagnostics off. Target about 30 seconds on the development machine, excluding builds. |
| `routine-extended` | Diverted warm, short adaptive/source with detailed diagnostics, and short neutral/impurity feature goldens. |
| `full` | Both real topologies, complete cold bootstraps and transport activation, distinct diverted features, and focused serial/OpenMP/MPI/hybrid comparisons. |

`full` selects its own evidence; it does not run the other profiles first. Its
cold evidence comprises adaptive three-stage bootstraps for limited and diverted
cases, plus an independent limited fixed-mesh bootstrap on the accepted adaptive
final mesh. Each bootstrap uses `time_init -> diffred -> steady`; transport is off
and the steady diffusion is 16 m²/s. Initialization requires valid finite output.
The later stages must meet the declared Newton limit. The fixed run starts from
analytical plasma, not the current adaptive solution. The limited suite also
compares the two completed bootstrap endpoints.

Transport is a separate feature on **both** topologies. Its five stages lower the
transport diffusion floor through 8, 4, 2, 0.5 and 0.1 m²/s. The diverted case
also checks source relocation, neutral pressure, perpendicular neutral diffusion,
two neutral limiters, radiation off, N and N+W impurity radiation, and stored
magnetic/current fields against its own reference. Limited adds analytical
initialization and the focused NeutralGamma probe described above. There is no
synthetic magnetic geometry.

Parallel suites compare serial OMP1 with serial OMP16 for OpenMP, and serial OMP1
with MPI4×OMP1 for MPI. The diverted short-transport suite additionally compares
MPI4×OMP1 with production MPI4×OMP4. Selected short transport, adaptive/source
and NeutralGamma runs carry detailed diagnostics; the long cold chains do not.
`routine-extended` has no layout-pair comparison: choose a focused parallel suite
or `full` when parallel agreement is the question. Suite relations named `openmp`,
`mpi`, `hybrid` and `all_pairs` generate layout pairs from the available layouts.
Run `list suites` to see current focused selections and their descriptions.

## Settings and builds

`settings.local.json` is the usual place for machine paths and default
selections. The quick start needs only the diverted bundle. For limited checks
or `full`, add the real limited case, `legacy_case`, under `defaults.bundles`.
`settings.example.json` shows both cases and an explicit `run_root`.

`run_root` holds runs and reports and defaults to `~/.cache/mhdg-regression`.
`build_root` holds newly built executables and defaults to `run_root/builds`.
Set these paths once for your machine if the defaults are unsuitable. To use a
different settings file, pass `--settings FILE` for one command or set
`MHDG_REGRESSION_SETTINGS` in the shell. For a single-case command,
`--bundle DIR` or `--build-manifest FILE` can override a saved default. A
multi-case profile uses the case-to-bundle map in `defaults.bundles`; one
`--bundle` path cannot represent both physical cases. The harness finds its
repository catalogs automatically, so ordinary commands need no catalog paths.

If no build is selected, use `check warm --build` or build reusable executables:

```bash
python -m regression_tests build routine --jobs 8
```

The command prints its `build_metadata.json`; select that file as `defaults.build`
for later checks. Builds derive the required serial/MPI/model variants from the
selection, rather than using a fixed list. A new build needs a usable compiler
environment in the shell or the optional `environment_script`. The generic
`positionFeketeNodesTri2D.h5` must be beside each executable and is not a case
bundle input. Adaptive checks also need
[HDG_postprocess](https://github.com/wave46/HDG_postprocess) importable in the
selected Python environment. Plain `doctor` and `check` select the default
`routine` profile, including its adaptive check; `doctor warm` and `check warm`
select only the focused warm suite.

## Commands and results

```bash
python -m regression_tests list cases
python -m regression_tests list workflows diverted_case
python -m regression_tests list suites
python -m regression_tests list layouts

python -m regression_tests check warm --case legacy_case
python -m regression_tests check transport_hybrid
python -m regression_tests check full
```

`run CASE WORKFLOW` executes one workflow for investigation without asserting a
regression comparison. `prepare CASE WORKFLOW` renders its inputs without
launching MHDG. Both use the workflow's default layout unless `--layout` is
supplied. A focused suite uses `check SUITE`; a profile uses `check PROFILE`.
`--diagnostics off|summary|equations|detailed` can override the selected mode
for a focused investigation. `check` requires a golden bundle by default;
`--allow-candidate` is available when deliberately validating an unpublished
candidate that already contains the selected references.

```bash
python -m regression_tests run diverted_case baseline_warm --run-id investigation-01
python -m regression_tests check full --run-id review-01
python -m regression_tests check full --run-id review-01 --resume
```

Suite/profile resume reuses completed runs only while their recorded outputs and
selected inputs still match. It reassesses the scientific checks; failed,
incomplete or changed workflows are preserved and rerun in new directories.
Resume needs the original build and cannot be combined with `--build`.

Runs live under the configured `run_root`:

```text
profiles/PROFILE/RUN_ID/profile_summary.json
suites/SUITE/CASE/RUN_ID/suite_summary.json
CASE/WORKFLOW/LAYOUT/RUN_ID/{run_plan.json,run_metadata.json,stdout.log,outputs/,...}
```

The profile summary links suite summaries; a suite summary lists every workflow,
layout and pair result. Comparison reports contain per-equation norms, mesh and
Newton evidence, failures, and the selected tolerance profile. Build metadata
records the solver revision, selected variants and artifact hashes; run metadata
records what was actually launched and the resulting output hashes. To recheck
saved evidence without launching MHDG:

```bash
python -m regression_tests compare --suite /path/to/suite_summary.json
python -m regression_tests compare /path/to/completed/run
```

A matching adaptive mesh is compared directly, including fields and gradients.
Only valid *different* meshes use deterministic interior-point interpolation via
`hdg_postprocess`, with required point coverage. Parallel mesh checks require
matching numbering and generated-mesh identity. Fixed comparisons check the
mesh, every equation in `u`, `q` and `u_tilde`, and present transport/magnetic
data. Every run checks its selected executable, runtime-file and output identity,
including model, feature switches and build provenance. Intermediate cold stages
must be valid; the final reference cannot conceal a bad earlier stage.

The numerical limits live in `tolerances.json`. The default converged Newton
limit is `2e-4`; initialization and short race probes require a finite recorded
Newton value without that bound. Matching warm states use tighter field limits
than cross-layout probes. Different-mesh adaptive sampling has separate solution
and gradient limits. These are regression limits, not a proof of physical
accuracy, positivity or mesh convergence. Enabled balance diagnostics also check
output presence, finite values, units and terminal/HDF5 agreement; selected
parallel suites compare diagnostic totals. They do not independently establish
that every physical source term is correct.

## External bundles and golden references

Use `bundle readiness` on a prepared source directory to see user-supplied,
workflow-producible, optional, present and missing files. It reports producer
prerequisites but does not run them. A missing required file makes readiness fail
even if a producer is declared. `--workflow` adds that workflow's required files,
including any reference it declares.

```bash
python -m regression_tests bundle readiness diverted_case \
  --source /path/to/prepared-inputs --workflow baseline_warm
python -m regression_tests bundle create --case diverted_case \
  --source /path/to/prepared-inputs --output /path/to/candidate-bundle
python -m regression_tests bundle validate /path/to/candidate-bundle
```

Creation copies available declared source files into a self-contained candidate,
following source symlinks into physical copies. It records roles, hashes and
sizes, validates the result, and refuses to overwrite an output. A new reference
need not exist to create a base candidate; a workflow-specific readiness or
validation request will identify that missing reference. During runs, immutable
bundle inputs are symlinked into temporary directories; rendered parameter and
transport files are private run copies. Source bundles are not modified.

A reference change requires a separate reviewed golden action. `golden refresh`
builds the case's required variants once, runs the ordered producers in
`golden.json`, collects a candidate and runs validation suites. It records old
reference comparisons for review when possible; old agreement is not a condition
for producing a valid new reference. Refresh can be expensive and has no campaign
resume. Keep a failed workspace for inspection and use a new workspace after a
correction.

```bash
python -m regression_tests golden refresh diverted_case \
  --bundle /path/to/source-bundle --workspace /path/to/new-refresh --jobs 8
# Review new-refresh/refresh.json and the linked suite/comparison reports.
python -m regression_tests golden publish /path/to/new-refresh \
  --output /path/to/new-golden --bundle-version 3.0.0 \
  --reason "Accepted physical or numerical change" \
  --provenance "Scientific review record"
```

Publication is a separate explicit command. It revalidates producer success,
convergence, candidate files and validation evidence, records the reason and
provenance, and refuses an existing destination. The published golden contains
physical copies of its inputs and references; it has no dependency on refresh
workspace symlinks. Update the private settings default only after reviewing
and publishing it.

## Extending and maintaining the harness

The repository owns `cases/*.json` (physical roles and case workflow choices),
`workflows.json` (shared procedures and stages), `suites.json` (suites and
profiles), `layouts.json` (available layouts), `tolerances.json` (numerical
limits), and `golden.json` (reference producers and validation suites). Add a
feature to the nearest real case, derive a workflow where possible, then add a
focused suite and a profile selection only if it contributes distinct evidence.
Physical input values and machine paths belong in external sources/settings.
Changing scalar parameters or namelists usually needs JSON overrides, not Python.
The tutorials build from [one warm feature check](tutorials/01-warm-feature.md)
to [layout comparisons](tutorials/02-layout-comparison.md), a
[multi-workflow suite](tutorials/03-scientific-suite.md), a
[cross-case profile](tutorials/04-profile.md), and a
[real case with an input bundle](tutorials/05-real-case-bundle.md).

`python -m regression_tests suggest --base develop` groups advisory checks for
changed files; it reads Git and catalogs but executes nothing. Suggestions do
not replace scientific judgment about cold convergence or diagnostic formats.
The Python harness tests use temporary fixtures and fake executables, not MHDG:

```bash
python -m pytest -q regression_tests/tests
```

Inspect external storage before any removal. `clean` only previews unless given
explicit inventoried directories **and** `--delete`; configured default bundles
and builds, unfinished work, unknown historical data and retained dependencies
are protected. Use a storage root that contains the related runs and summaries,
and run cleanup when no regression work is active:

```bash
python -m regression_tests clean /path/to/regression-data
python -m regression_tests clean /path/to/regression-data runs/CASE/WORKFLOW/LAYOUT/RUN_ID
# Only after reviewing that exact selection and its dependencies:
python -m regression_tests clean /path/to/regression-data runs/CASE/WORKFLOW/LAYOUT/RUN_ID --delete
```

The harness never cleans runs, candidates, historical goldens or builds as a side
effect of `check`, refresh or publication. External storage cleanup is an
explicit user action.
