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
| **Suite** | A focused selection of workflows and layouts. |
| **Profile** | An ordered set of suites, potentially using different cases. |
| **Bundle** | Validated external files plus a generated manifest. A **golden** bundle contains accepted references. |

```text
case + workflow + layout -> run
suite                    -> selected runs and comparisons
golden refresh           -> validated candidate for review
golden publish           -> accepted golden bundle
```

The default `routine` profile uses `diverted_case`. Limited coverage is retained
in `full`. The Python package is the public entry point; catalog paths are
discovered automatically.

```bash
python -m regression_tests --help
python -m regression_tests list cases
python -m regression_tests list workflows diverted_case
python -m regression_tests list suites
python -m regression_tests list layouts
```

Commands return 0 for success, 1 for failed checks/runtime input errors, 2 for
invalid command usage, and 130 for interruption. Solver completion from `run`
is distinct from a passing scientific comparison from `check`.

## First-time setup

Run commands from the repository root. Install the Python dependencies:

```bash
python -m pip install -r regression_tests/requirements.txt
```

Adaptive comparisons also need `hdg_postprocess`. Activate that Python environment
or invoke its interpreter directly instead of `python`.

Machine preferences live in one ignored JSON file. Copy the example and select
an existing build manifest and accepted bundles, or pass those selections on the
command line:

```bash
cp regression_tests/settings.example.json regression_tests/settings.local.json
python -m regression_tests doctor
python -m regression_tests check
python -m regression_tests check routine-extended

# Explicit selections override local defaults.
python -m regression_tests check warm --case diverted_case \
  --bundle /path/to/diverted/golden_bundle \
  --build-manifest /path/to/build/build_metadata.json
```

| Owner | Settings |
| --- | --- |
| Repository | Automatically discovered cases, workflows, suites, layouts and tolerances |
| Machine | `run_root`, optional `build_root`, `build_jobs`, `mpi_launcher`, `environment_script` |
| Invocation/defaults | `--bundle`, `--build-manifest`; optional `defaults.bundles` per case and `defaults.build` |
| Generated build manifest | Executables, solver revision, binary/runtime checksums; never enter these manually |

Use `--settings FILE` or `MHDG_REGRESSION_SETTINGS` to select a different machine
JSON file. Legacy KEY=VALUE `.env` settings are no longer accepted. Relative paths
inside JSON are relative to that file; command-line paths
are relative to the current directory. Without a machine file, the scratch root
is `~/.cache/mhdg-regression`; build output defaults to its `builds/` directory.
No bundle or build is chosen implicitly. `check` prints its effective selections.

Open MPI is discovered from the configured shell environment. Set `mpi_launcher`
only when discovery is insufficient. `environment_script` is optional when the
current shell already provides the runtime libraries/toolchain. A supplied build
needs no compiler setup. `doctor [PROFILE_OR_SUITE]` verifies Python prerequisites, catalogs,
bundle checksums, build artifacts, launcher and scratch accessibility without
creating runs. It does not execute the solver or prove numerical convergence.

`build` and `check --build` produce a manifest; select the printed
`--build-manifest` path for later runs/resume, or save it as `defaults.build`.
A new build ignores the old default build. Build manifests must describe completed
2D NGammaTiTeNeutral or NGammaTiTeNeutralGamma builds; selected binaries and Fekete data must match their
recorded checksums. Runtime environment setup remains machine-owned rather than
being copied from another machine's build provenance.

## Scientific coverage

Warm workflows restart an existing solution; cold workflows start
analytically. Fixed and adaptive indicate whether the mesh can change.

The retained cold recipe has seven stages: `time_init`, `diffusion_reduction`,
then five transport continuations. `cold_fixed` uses the refined fixed mesh;
`cold_adaptive` starts coarse and adapts in the first two stages. Short
`cold_step_fixed` and `cold_step_adaptive` workflows use one time step and two
Newton iterations; the adaptive variant also performs one adaptation pass.
Use `list workflows CASE` for the complete current selection of feature variants.

Full cold stages run sequentially and restart from their predecessor. The
one-step workflows probe mesh construction and races, not convergence.
The limited neutral-feature golden workflows all use the same accepted source-relocated restart.
Each then enables exactly one feature: source relocation alone, neutral
pressure, perpendicular diffusion, fixed-`Tn` limiting, or `Tn=Ti` limiting.
The active limiter declarations always retain
`neutral_wall_sources_in_elements=true`.
The source-only golden keeps the strict fixed-layout tolerance. Repeated
MPI4-by-OMP4 feature solves showed converged endpoint variation up to
`2.90e-7` relative L2, so the other four use the measured
`neutral_feature_same_layout` profile (`5e-7` for relative L2 and normalized
L-infinity). A reference refresh updates the source reference but deliberately
keeps its accepted restart fixed, so every feature continues from the same
state.

The two-Newton-iteration feature race uses a separate `1e-7` relative-L2
profile. Repeating an identical nonlinear `mpi4_omp1` limiter cell varied by
`7.40e-8`, while the source-only control remained near `1e-11`; the larger
threshold avoids treating amplified parallel solver ordering as an assembly
race. It remains five times tighter than the converged-feature profile.

Tracked layouts are:

| Layout | Execution |
| --- | --- |
| `serial_omp1` | Serial, one OpenMP thread. |
| `serial_omp16` | Serial, sixteen OpenMP threads. |
| `mpi4_omp1` | Four MPI ranks, one thread each. |
| `mpi4_omp4` | Four MPI ranks, four threads each; canonical layout. |

Identifiers matching `serial_ompN` or `mpiM_ompN` determine the executable,
ranks, and threads. MPI runs bind each rank to exclusive cores.

`suites.json` supplies common case/layout defaults. A suite can override them,
select `layouts: "all"`, or request generated `relations`: `openmp` compares
serial thread counts against one thread; `mpi` compares serial against MPI with
one thread per rank; `hybrid` compares MPI thread counts at the same rank count;
`all_pairs` compares every selected pair. Relations use the selected layouts
(all registered layouts if omitted), deduplicate pairs and run each participating
layout once. A relation with no matching pair is an error. Reference comparisons
default to enabled for ordinary suites and disabled for relation suites;
`reference_comparisons: true` enables both kinds of evidence together.

| Profile | Coverage |
| --- | --- |
| `routine` (default) | Diverted hybrid warm golden, with transport_1d enabled, and one short adaptive execution smoke; diagnostics off |
| `routine-extended` | Diverted warm golden plus short fixed/adaptive OpenMP, MPI and hybrid comparisons; detailed diagnostics on the parallel runs |
| `full` | Extended checks, full seven-stage adaptive cold runs for both cases, limited fixed cold, and distinct limited neutral/impurity/initialization features |

On the current machine with the selected existing PR08 build, one complete
`routine` command took 24.31 seconds, excluding build time. This is measured
evidence, not a hard performance threshold.

The adaptive smoke checks solver completion and output presence. Quantitative
field/mesh parity is covered by `parallel` in the extended profile. Full retains
source-relocated cold initialization, pressure, perpendicular projection, both
neutral limiter modes, NeutralGamma, stored-field compatibility and off/N/N+W
impurity references. The ordinary warm case supplies the W reference.
Detailed source-placement checks run on the existing short fixed/adaptive source
suite; long cold chains default to diagnostics off. All seven cold stages remain.
Diverted feature goldens are still a later scientific baseline task.

Profiles share focused suite definitions; `full` includes `routine-extended`
once. Use `list suites` to see the selections. `--case` selects another topology
for a focused suite such as `warm`, `parallel`, `cold` or `fixed_cold`; profiles
declare their own cases. `full` gets both bundles from `defaults.bundles`.
`--bundle` is supported for a selection using one case.

`cold_matrix` and `warm_parallelism` remain focused characterization checks.
Normal parallel checks isolate OpenMP, MPI and hybrid effects with three
relations. Golden refresh selects producer workflows directly; reference and
restart files do not require their own suite wrappers.

## Common tasks

```bash
# Fast diverted regression, using the selected existing build.
python -m regression_tests check
python -m regression_tests check routine-extended
python -m regression_tests check full

# Build the required serial/MPI variants once for the whole selection.
python -m regression_tests check routine-extended --build --build-jobs 8

# Focused limited checks or parallel evidence from a candidate bundle.
python -m regression_tests check warm --case legacy_case
python -m regression_tests check parallel --allow-candidate --bundle /path/to/candidate

# Resume the same profile and build after an interruption.
python -m regression_tests check full --run-id overnight-01
python -m regression_tests check full --run-id overnight-01 --resume
```

Completed runs are reused only while their recorded output hashes still match;
comparisons are rerun. Failed, incomplete or changed runs are preserved and the
whole workflow starts in a fresh `RUN_ID-resume-N` directory. A profile stops at
its first failed suite. Resume uses the ordinary suite records for completed and
unfinished work. Changed selected inputs or declarations are rejected.
`--build` cannot be combined with `--resume`; reuse its printed build manifest.

Use `--diagnostics off|summary|equations|detailed` to override the selected mode
for a focused investigation. Four-format and on/off comparisons are focused
checks for diagnostics changes, rather than duplicate routine solver chains.
Suites without golden comparisons validate each stage's mesh, field sizes,
finiteness and declared Newton acceptance. Initialization and short probes require
finite Newton values without a convergence bound; `source_cold` enforces the cold
stage thresholds. Compact validation results are stored in the suite summary.
Golden producers share these checks independently of old-reference agreement.

`--run-only` defers reference, convergence and diagnostic checks; build/output identity
is still checked on every execution. `compare --suite` checks a saved suite summary,
including output validity for suites without reference comparisons.

### Build reusable executables

```bash
python -m regression_tests build full \
  --settings /private/path/settings.json --jobs 8
```

`build [PROFILE_OR_SUITE]` (default `routine`) and `check --build` derive the
required model/execution combinations from the selected workflows and layouts.
Each combination is built once; OpenMP thread counts share an executable. A
workflow defaults to NGammaTiTeNeutral and can declare `model` explicitly;
`cold_step_neutralgamma` requires NGammaTiTeNeutralGamma. Currently `routine`
needs ordinary-neutral MPI, while `full` also needs ordinary-neutral serial and
Gamma serial. These are selection results, not fixed build counts.

Builds run sequentially because they share objects in `lib/`. Existing objects
have an unknown manual configuration, so the harness cleans before building and
between model/execution configurations. An empty tree needs no initial clean.
The command records commands, logs, Git state, toolchain and artifact identities.

Without `--build`, missing required executables fail preflight before a profile
starts. The harness does not silently build midway through a check. Version-3
manifests identify each artifact by model/execution; existing version-2 ordinary
build manifests remain readable so accepted runs can still be investigated.

Before every launch, the selected binary and runtime data must match the build
record. After execution, the output must report the expected model, dimension,
equations and build provenance. Relevant feature switches, neutral parameters
and active impurity species/concentrations are checked against the run inputs.
Zero exit status with an incorrect output contract fails the run and stops a cold
chain. Saved comparisons apply the same identity checks when execution records
are available; standalone debugging comparisons can lack those records.

### Inspect or run one workflow

Preparation validates and renders an isolated run but does not launch MHDG:

```bash
python -m regression_tests prepare legacy_case cold_step_adaptive \
  --layout serial_omp16 --settings /private/path/settings.json
```

Execute one workflow directly when debugging:

```bash
python -m regression_tests run legacy_case cold_fixed \
  --layout mpi4_omp4 --run-id investigation-01 \
  --settings /private/path/settings.json
```

Prefer suites for routine work because they preserve one resumable summary.

### Recompare saved results

Neither command below launches MHDG:

```bash
# One completed run; policy comes from run_plan.json.
python -m regression_tests compare /path/to/completed/run

# Every same-layout golden check and declared layout pair in a suite.
python -m regression_tests compare --suite \
  /path/to/suites/cold_matrix/legacy_case/overnight-01/suite_summary.json
```

Enable diagnostics on an existing suite; its normal checks also validate the
saved diagnostic output:

```bash
python -m regression_tests check warm --diagnostics detailed --run-id warm-detailed
python -m regression_tests check parallel --diagnostics detailed --run-id race-detailed
```

The first command validates output presence, finite values, core units and
terminal/HDF5 content agreement. The second also compares diagnostic scalars
across the existing layout pairs, after the solution and mesh checks pass.
For equation output, rate differences are bounded by `max(1e-12, 1e-9 * scale)`,
where `scale` is the largest absolute rate in that equation across both files;
content is scaled separately. Summary aggregates are scaled individually. The
absolute floor uses each quantity's reported physical units.

`off`, `summary`, `equations`, and `detailed` are supported. Off requires absence
of balance output. The other modes require their core output and readable finite
quantities; detailed source-relocation runs additionally compare the integrated
puff with the configured input and require zero relocated wall puff/pump.
These checks do not establish physical convergence or independently validate all
source and BC terms. There are no exhaustive arithmetic-identity checks.

To check saved diagnostic output alone:

```bash
python -m regression_tests compare --suite --diagnostics /path/to/suite_summary.json
```

This writes a compact `balance_diagnostics_check.json`; full history remains in
`stdout.log`. Ordinary `compare --suite` includes diagnostic output validation and
any declared parallel comparisons in its verification summary.

For a focused real mode check, reuse the same `warm` workflow, build, bundle and
layout with distinct run IDs for the four modes. Compare enabled solutions directly
with the off result using the existing comparator:

```bash
python -m regression_tests compare /path/to/warm-detailed/run \
  --reference /path/to/warm-off/final-solution.h5 --tolerance-profile fixed_same_layout
```

Mode changes do not require dedicated workflow aliases or a cold-chain matrix.

### Refresh and publish goldens

Refresh builds the required executable variants once, runs the case's ordered
producers, collects one candidate, and runs validation suites against it:

```bash
python -m regression_tests golden refresh diverted_case \
  --settings /path/settings.local.json \
  --bundle /path/source_bundle --workspace /path/refresh-workspace --jobs 8
```

The source may be an accepted golden or a candidate containing the required
physical inputs. `--bundle` can be omitted when the case has a configured default.
Refresh does not require old reference files to generate new results. Each
producer must complete, pass its declared Newton convergence checks, and produce
structurally valid, finite fields. Available old-reference comparisons are saved
for review; differences from an old golden do not excuse a failed producer.

Limited refresh retains the full fixed/adaptive cold matrix across four layouts
and all six pairs. Its canonical fixed cold result feeds the warm restart;
subsequent warm, impurity and neutral workflows produce the respective files.
Diverted refresh retains the full adaptive hybrid cold producer followed by warm
reconvergence. Both retain their existing seven-stage recipes. These choices live
in repository-owned `golden.json`; the CLI never takes a catalog path.

Producers pass output files directly to subsequent workflows. There are no
intermediate candidate copies or generated settings files. After collection,
existing warm/parallel/feature suites validate the new candidate. The shared build
serves all producers and validation runs. Run refreshes sequentially in a worktree
because the build uses its solver object directory.

Review `WORKSPACE/refresh.json`, the old-reference comparisons, and the validation
suite reports. The workspace contains `candidate/`, `runs/`, and `build.json`.
A successful refresh reports `ready`; it has not published a golden. Publication
is a separate, explicit command:

```bash
python -m regression_tests golden publish /path/refresh-workspace \
  --output /path/new_golden --bundle-version 2.9.0 \
  --reason "Accepted change to the solver" \
  --provenance "Scientific review record or change identifier"
```

Publication rechecks producer and validation evidence, recorded file identities,
convergence, catalog identities and candidate integrity. A changed output or
failed producer blocks it. It requires a version, reason and review provenance;
solver/build and source-bundle identity are recorded automatically. It refuses an
existing destination and creates a standalone bundle with physical copies of its
inputs, references and evidence. Old campaign provenance is not accumulated into
each new golden; the source bundle identity records its predecessor.

A failed or interrupted refresh retains its workspace and logs for inspection.
Start a new refresh in a new workspace after correcting the cause. Golden refresh
has no campaign resume, stage retries, rewind, catalog amendments or partial
component updates. Ordinary `check --resume` remains available for suites.

## Outputs and provenance

Runs are stored below the configured `run_root`; `build_root` defaults to its
`builds/` directory:

```text
runs/
├── builds/.../                    executables, logs, build_metadata.json
├── profiles/PROFILE/RUN_ID/       profile_summary.json
├── suites/SUITE/CASE/RUN_ID/       suite_summary.json
└── CASE/WORKFLOW/LAYOUT/RUN_ID/
    ├── run_plan.json              requested inputs and comparison policy
    ├── run_metadata.json          execution result and provenance
    ├── stdout.log / stderr.log
    ├── outputs/
    └── stages/...                 staged-workflow runs
```

Comparisons write `comparison.json`, `matrix_comparison.json`, or
`verification_summary.json`. Suite summaries are updated after every cell.

Metadata has four owners:

| Record | Contents |
| --- | --- |
| Build metadata | Revision, dirty state and changed-file list, build configuration, toolchain and produced artifact identities |
| Executed run/stage metadata | Actual command, environment, observed executable/runtime hashes, build-manifest identity, outcome and output hashes |
| Staged workflow metadata | Stage order, directories, outcomes and selected outputs; detailed observations remain in each stage |
| Suite summary | Selected checks, results and file identities needed for resume |

Run observations retain binary hashes because they identify what was actually
launched. Revision and build configuration are read from the build manifest.
Suite summaries do not copy the machine settings dictionary; unused executables
and build-only preferences do not affect resume. Selected runtime files are
checked independently of executable checksums. Executables are selected through
a build manifest.
Published canonical references from staged runs include a local copy of the
final execution metadata.

Each newly built solver writes automatic compile-time identity into its HDF5
solutions:

```text
/provenance/git_commit
/provenance/git_dirty
/provenance/build_id
```

Do not edit these fields. One build ID is shared by the artifacts in that build
set. Each run's output must match its selected build provenance. Candidate and
golden revisions may differ; they are not required to match each other. Per-run
metadata references the build record rather than copying its revision/configuration.

## External bundles

The V2 `manifest.json` is generated. Inspect the files expected for a case and
workflow before creating a candidate bundle:

```bash
python -m regression_tests bundle readiness diverted_case \
  --source /private/path/prepared_diverted_case --workflow warm
```

The report lists filenames, presence, whether each file is required or optional
for the selection, and whether a declared workflow can produce it. Missing
producible files list producer workflows and their prerequisites. Producers are
not launched or chosen automatically. A producer declaration does not mean the
file is present or accepted as a golden.

Readiness is for prepared sources before creation. Use `doctor` to check an
existing bundle as part of the selected regression setup. To validate a bundle
independently of machine settings, use the command below; it checks recorded
sizes, checksums, contained paths and reference-matrix integrity:

```bash
python -m regression_tests bundle validate /private/path/candidate_bundle \
  --workflow warm
```

`--workflow` can be repeated on readiness, create and validate. Without it, these
commands require the base bundle files; with it, they additionally require every
selected workflow's inputs, including its currently declared reference files.
Optional files become required when a selected workflow uses them. Readiness
returns exit code 1 when required files are missing, even if producers exist.

Prepare external files using the reported names, then create a candidate:

```bash
python -m regression_tests bundle create \
  --case legacy_case \
  --source /private/path/prepared_legacy_case \
  --output /private/path/candidate_bundle --workflow warm
```

Creation copies all available declared files under `inputs/`, follows source
symlinks into physical copies, and generates role mappings, media types, sizes
and SHA-256 checksums. It validates before publishing the directory and refuses
to replace an output. Candidate and golden bundles share this contract; their
class records acceptance. `doctor` also checks the selected profile's workflow
requirements, so a valid base bundle alone does not imply readiness for a run.

`positionFeketeNodesTri2D.h5` must be beside each executable. It is generic
solver runtime data, not case-specific bundle data.

| Edit directly | Generated; do not edit |
| --- | --- |
| Private settings and prepared physical files | Bundle manifest, identifiers, sizes, checksums |
| `workflows.json`, `cases/*.json` | Rendered parameter files, input links, run/build metadata |
| `suites.json`, `layouts.json`, `tolerances.json` | Summaries, comparison reports, HDF5 provenance |

## Adding a parameter variant

`workflows.json` defines the shared warm, fixed/adaptive cold and short-step
procedures. Its `cold_bootstrap` and `transport_continuation` sequences compose
the unchanged seven-stage cold recipe. Each case declares its own physical file
roles and selects workflows by name; `{}` uses the shared definition unchanged.
Only workflows selected in the case are exposed, but `extends` can also name an
unselected shared parent. Case overrides are applied before parent resolution,
so derived workflows inherit that case's changes. Descriptions are inherited.

To adjust one stage without repeating the recipe, use `stage_overrides`:

```json
"cold_fixed": {
  "stage_overrides": {
    "continuation_05": {"parameter_overrides": {"tNR": 1e-5}}
  }
}
```

This limited-case override also reaches its adaptive and diagnostic variants.
Stage overrides are keyed by existing stage IDs and can change parameter/transport
roles, parameter values, or the Newton check. Unknown IDs fail during loading.
Sequences contain explicit stages; workflows may combine sequence names and
inline stages. Sequence expansion does not change restart order or execution.

Common variants require JSON, not Python. In `cases/CASE.json`, inherit the
closest workflow and override only changed parameters:

```json
"cold_step_fixed_three_nr": {
  "extends": "cold_step_fixed",
  "description": "Fixed coarse-mesh probe with three Newton iterations",
  "parameter_overrides": {
    "nrp": 3
  }
}
```

This inherits inputs, mesh, stages, layout, comparison policy, and the other
parameter overrides. Add the new workflow identifier to an existing or new
suite in `suites.json`.

Overrides accept booleans, finite numbers, and strings. The named assignment
must occur exactly once, or declare its insertion namelist:

```json
"parameter_overrides": {"new_coefficient": 0.25},
"parameter_namelists": {"new_coefficient": "phys_lst"}
```

`parameter_namelists` can supply shared defaults in `workflows.json`, case
defaults, or workflow/stage overrides. It only inserts missing assignments;
existing assignments retain their location. A missing/ambiguous assignment or
insertion group fails. No Python physics-key whitelist is needed. The renderer
handles one complete scalar assignment per line; keep arrays and more complex
Fortran syntax in the supplied input files. Comments and untouched lines are
preserved. Preparation renders a private copy and records effective values in
`run_plan.json`; immutable bundle inputs remain symlinked and unchanged.

Workflow parameter overrides apply to every stage; stage parameter values take
precedence there. With `extends`, `parameter_overrides` and `parameter_namelists` merge with the parent,
and `stage_overrides` merge by stage ID (including their parameter values).
Other supplied fields, including the `stages` list, replace the parent field.

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

Adaptive comparisons first validate both meshes and inspect coordinates,
connectivity, polynomial order and face numbering. Matching discrete meshes use
direct comparison of `u`, `q`, `u_tilde`, transport and magnetic data. Matching
cold stages use `fixed_stage_reference`; final cold-versus-warm checks use
`cold_fixed_reference`. Shared workflow declarations own these selections through
`comparison.direct_profile` and `comparison.direct_stage_profile`.

Only differing valid meshes use `HDG_postprocess` to interpolate conservative
fields and gradients at deterministic interior points, with `adaptive_reference`
tolerances and full required point coverage. Malformed input is an error; a failed
direct comparison never retries through interpolation. Reports record the selected
method, profile and reason. A debugging `--tolerance-profile` override must fit the
selected method: fixed-field limits for matching meshes, sampled-field limits for
different meshes.

Race probes require identical connectivity arrays, including node, element,
face and boundary ordering, and matching discretization metadata. The generated
adaptive `temp.msh` files must also be
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

## Harness implementation and tests

`catalog.py` owns case/workflow inheritance, layouts and suite/profile resolution;
`config.py` owns machine JSON, build/runtime selection and environment setup. Preparation
and suites accept resolved settings, so they do not parse configuration files.
`bundles.py` owns manifest paths, staged-reference indexes, bundle integrity and
publication-class checks; comparison code consumes its validated references. Shared document/
schema I/O, file/path contracts and errors/record identifiers live in `documents.py`,
`files.py` and `support.py`. All callers use package imports; no `tools/` import
bridge or test path injection is needed.

These tests use temporary bundles, small arrays, and fake executables. They do
not launch MHDG or need physical data.

Run from the repository root:

```bash
# Focused example; substitute the affected test module.
python -m pytest -q regression_tests/tests/test_preparation.py
# Whole Python harness; no solver build or real case inputs needed.
python -m pytest -q regression_tests/tests
```

Tests cover composition and generated layout relations, bundle integrity,
preparation, execution, direct/interpolated/parallel comparisons, diagnostics,
resume, CLI, build and refresh/publication failure handling. Cleanup and suggestions
have focused read-only/protection checks. Passing pytest establishes harness
behavior; use real-data profiles to establish solver regression evidence.

## Limitations and scientific acceptance

- Adaptive golden meshes can differ between revisions; within-build layout pairs
  require exact meshes. Matching adaptive meshes use direct field comparisons.
- Existing limited fixed/adaptive references have different mesh lineage (844 vs
  846 elements). Align the fixed input with the final adaptive mesh during the
  later accepted golden redesign, not by loosening comparison tolerances.
- Diverted has a final golden but lacks historical per-stage references. Its cold
  stages can be compared across layouts; a historical stage comparison needs new
  accepted reference data. Diverted feature goldens remain future work.
- Timing is recorded without a pass/fail threshold. Golden publication verifies
  technical evidence; scientific acceptance requires human review.

Before redesigning the cold recipe or refreshing references, run `full` plus both
cases' `cold_matrix` and `warm_parallelism` checks with unchanged accepted bundles
and the same selected builds. Investigate failures before changing baselines.

### Inspect and clean stored runs/data

Use one encompassing external storage root so retained runs, suites and refresh
workspaces can protect their dependencies. Inventory and selection are read-only
by default:

```bash
python -m regression_tests clean /path/regression-data
python -m regression_tests clean /path/regression-data runs/old-run suites/old-suite
# After reviewing the same explicit paths:
python -m regression_tests clean /path/regression-data runs/old-run suites/old-suite --delete
```

Paths are absolute or relative to the storage root. The report shows directory
kind, file sizes in bytes, selection and protection reasons. Sizes count stored
file lengths without following symlinks; they are not a promise of reclaimed disk
blocks. Select whole inventoried runs, suite/profile records, builds, refresh
workspaces or candidate/golden bundles. There is no automatic age-based selection.

Active/unfinished work, unknown historical data and configured default bundles
or builds are protected. Retained records and symlinks protect referenced data;
select dependent runs and their suite/profile records together to remove them.
`--keep DIRECTORY` additionally protects a directory and inspects its dependencies;
use it for retained work outside the storage root. Use the encompassing root and
list external consumers: the command cannot discover unrelated storage elsewhere.
Published bundle provenance is historical and does not keep old workspaces alive.

Run cleanup while regression jobs and refresh/publish/resume commands are stopped.
`--delete` rescans the tree on each invocation and refuses protected selections
before removing anything. Removal does not follow symlinks. Sources and unrecognized
legacy campaigns stay for manual review. An old golden requires explicit selection
and must no longer be a configured default or a retained run's dependency. Nothing
is removed by `check`, refresh or publication, and no cleanup database is created.

### Suggest checks from code changes

```bash
python -m regression_tests suggest                 # staged, unstaged, untracked
python -m regression_tests suggest --base develop  # also include branch commits
```

With `--base`, committed changes are measured from the merge base of that revision
and HEAD. Local changes are added separately, including both sides of renames;
Git-ignored files are excluded. The command reads Git and the repository catalogs,
prints suggestions with the changed files/reasons, and executes nothing. It needs
no machine settings, solver build or external scientific bundle.

A small explicit mapping in `suggest.py` connects current code owners to checks.
For example, parallel/adaptivity changes suggest routine-extended; shared physics
or initialization changes suggest full, including the cold chains. Shared harness
workflow/comparison/catalog changes suggest pytest plus full real-data regression;
helper-only changes (cleanup, suggestions, reporting, doctor, build-helper tests)
suggest pytest. Harness diagnostics suggest pytest and the existing parallel/source
diagnostic checks. Reasons are printed once per command, with matching files grouped
underneath, and Python tests appear first. Profile membership comes from the suite
catalog so a broader suggested profile absorbs already-covered focused checks.
Unmapped files need manual assessment; Markdown/reStructuredText needs document review.

Harness-only changes use existing executables for their real-data checks. Changes
under solver sources or build files add `--build`. For several such checks, reuse
the first build by replacing subsequent `--build` options with its printed
`--build-manifest PATH`. Suggestions are not a complete verification plan:
file names cannot establish numerical impact. Assess long cold convergence and
focused diagnostic format/on–off solver checks when those behaviors change;
Python diagnostic-checker tests do not establish solver format coverage.
