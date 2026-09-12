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
accepted suite results   -> golden bundle
```

The default `routine` profile uses `diverted_case`. Limited coverage is retained
in `full`. The Python package is the public entry point; catalog paths are
discovered automatically.
Golden campaigns temporarily remain a standalone tool;
the new CLI does not expose the campaign recovery/publication state machine.

```bash
python -m regression_tests --help
python -m regression_tests list cases
python -m regression_tests list workflows diverted_case
python -m regression_tests list suites
python -m regression_tests list layouts
```

Catalogs are discovered inside the package; their paths are not CLI options.
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
file. Relative paths inside JSON are relative to that file; command-line paths
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
NGammaTiTeNeutral 2D builds; selected binaries and Fekete data must match their
recorded checksums. Runtime environment setup remains machine-owned rather than
being copied from another machine's build provenance.

Existing `.env` files remain accepted through explicit `--settings FILE` during
migration. Their old executable-path contract is unchanged. Golden campaigns
still use that format until their replacement; no new `.env` file is needed for
the package CLI. Existing suite summaries made before this settings refactor
cannot be resumed: start a new run ID; their output remains available to compare.

Use `--diagnostics off|summary|equations|detailed` to run the same suite with
an explicit balance-diagnostics mode. The selected override is recorded in
the suite summary and each run plan.

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
| `cold_step_neutralgamma` | Analytical start; coarse fixed mesh | Short NeutralGamma race check with detailed particle diagnostics. |
| `cold_step_adaptive` | Analytical start; coarse adaptive mesh | One time step, two Newton iterations, and one adaptation pass. |
| `warm_neutral_sources_in_elements` | Accepted source-relocated restart; fixed mesh | Reconverge with puff and pump in elements. |
| `warm_neutral_pressure` | Source-relocated warm restart; fixed mesh | Exercise `neutralp_lambda=0.05`. |
| `warm_neutral_perpendicular` | Source-relocated warm restart; fixed mesh | Exercise projected perpendicular neutral diffusion. |
| `warm_neutral_limiter_fixed` | Source-relocated warm restart; fixed mesh | Exercise active limiting with fixed `Tn=2.5 eV`. |
| `warm_neutral_limiter_ti` | Source-relocated warm restart; fixed mesh | Exercise active limiting with `Tn=Ti`. |
| `cold_step_fixed_neutral_sources_in_elements` | Analytical start; coarse fixed mesh | Short relocated-source race and detailed-balance check. |
| `cold_step_adaptive_neutral_sources_in_elements` | Analytical start; coarse adaptive mesh | Short relocated-source adaptation and detailed-balance check. |
| `cold_adaptive_neutral_sources_in_elements` | Analytical start; coarse adaptive mesh | Complete relocated-source cold workflow. |

Full cold stages run sequentially and restart from their predecessor. The
one-step workflows probe mesh construction and races, not convergence.
The PR06 golden workflows all use the same accepted source-relocated restart.
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

The old golden campaign temporarily also uses `cold_matrix`, `warm_parallelism`,
`source_reference`, `source_restart`, `neutral_references` and `impurity_restarts`.
Its full six-pair cold matrix remains available for reference production. Normal
parallel checks isolate OpenMP, MPI and hybrid effects with three relations.

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
Suite summaries now include the case in their directory name; earlier summaries
remain available for comparison but require a new run ID for new checks.

Use `--diagnostics off|summary|equations|detailed` to override the selected mode
for a focused investigation. Four-format and on/off comparisons are focused
checks for diagnostics changes, rather than duplicate routine solver chains.
`--run-only` defers comparisons; `compare --suite` checks a saved suite summary.

### Build reusable executables

```bash
python -m regression_tests build \
  --settings /private/path/settings.json --jobs 8
```

`build` produces both serial and MPI executables. `check --build` builds only
the variants required by the selected suite, once each. Builds run sequentially
because they share objects in `lib/`. Existing objects have an unknown manual
build configuration, so the harness cleans them before building; it also cleans
between serial and MPI variants. An empty build tree needs no initial clean.

The command prints the generated manifest path and records commands, logs, Git
state, toolchain versions, environment checksum, and executable checksums. It
writes no `settings.env`: the manifest selects the new executables directly.
A manifest may contain just one variant; a suite needing another reports the
missing build. Use `build` to produce both for subsequent mixed-layout suites.

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

### Publish accepted references

For a complete refresh, one command starts or continues the persisted golden
campaign and runs every producer and verification stage. It stops once, after
the complete candidate is ready for review:

```bash
python regression_tests/tools/golden_update.py update legacy_case \
  --settings /private/path/source.env --run-id refresh-01 \
  --output /private/path/new_golden_bundle --bundle-version 2.5.0
```

Run the same command again with `--accept campaign` after reviewing the final
candidate. This second command only publishes; it does not rerun completed
stages. `golden status WORKSPACE` shows progress; without `--workspace`, the
workspace is `MHDG_REGRESSION_RUN_ROOT/golden_campaigns/RUN_ID`. Resumes reject
changed inputs and never overwrite a workspace, candidate, or output.

The legacy refresh still produces cold matrices and warm/impurity/neutral
references, then runs focused initialization, parallel and verification checks.
Short source-placement checks use diagnostics; the long source-relocated cold
workflow keeps diagnostics off. Repeat `--only` to select `cold_matrix`, `warm`,
`impurity_mixture`, or `neutral_features`. A
cold-matrix refresh also updates the canonical warm restart, while `warm`
updates the warm reference. Updating only warm or only mixture references
records a consistency warning; select both together to avoid it.
Reference producers may differ from the old fields under review, but golden
promotion stops automatically if any producer fails its declared Newton
convergence check.

The legacy campaign does not carry specialized restart states directly from
the source bundle. Its analytical cold matrix first creates the canonical warm
restart. Dedicated producer stages then derive the impurity-off, nitrogen,
nitrogen-tungsten, and relocated-source restarts from that current warm state.
Reference and race stages consume those regenerated restarts. The previous
golden therefore supplies immutable inputs and comparison references, but not
the solution-state lineage published by the new campaign.

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
python -m regression_tests bundle create \
  --case diverted_case \
  --source /private/path/prepared_diverted_case \
  --output /private/path/diverted_case_pr03_source \
  --bundle-version pr03-source.1
```

Point a private settings file at that candidate and start the resumable cold
overnight in `tmux` or another persistent shell:

```bash
python regression_tests/tools/golden_update.py update diverted_case \
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
python regression_tests/tools/golden_update.py status \
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

Because multiple cases share the campaign catalog, a change confined to another
case does not invalidate an in-progress campaign. Resume records the old and new
catalog identities in `campaign_catalog_refreshes` after confirming that the
selected case declaration is unchanged. A change to the selected case still
requires an explicit `--retry-from STAGE`.

Refresh the limited golden with the same lifecycle but `legacy_case`, an
accepted golden source settings file, and distinct run/workspace/output names:

```bash
python regression_tests/tools/golden_update.py update legacy_case \
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

## Outputs and provenance

Runs are stored below `MHDG_REGRESSION_RUN_ROOT`. Builds use
`MHDG_REGRESSION_BUILD_ROOT`, which defaults to `RUN_ROOT/builds`:

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
checked even when using prebuilt executables without a build manifest.
Summaries from before this metadata change require a new run ID; saved outputs
remain available for comparison. Published canonical references from staged runs
include a local copy of the final execution metadata.

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

## Focused synthetic tests

These tests use temporary bundles, small arrays, and fake executables. They do
not launch MHDG or need physical data.

| Changed area | Test group |
| --- | --- |
| Settings/build selection and doctor | `tests/test_settings.py` (pytest) |
| Public package CLI smoke | `tests/test_cli.py` (pytest) |
| Shared workflows and inheritance | `tests/test_catalog.py` (pytest) |
| Generated layout relations and suite defaults | `tests/test_layouts.py` (pytest) |
| Build | `tests/test_build.py` (pytest) |
| Bundles and case loading | `tests.test_bundles` |
| Preparation and parameters | `tests/test_preparation.py` (pytest) |
| Execution and run metadata | `tests/test_execution.py` (pytest) |
| Fixed comparison | `tests/test_fixed_comparison.py` (pytest) |
| Adaptive comparison | `tests/test_adaptive_comparison.py` (pytest) |
| Reference matrices | `tests/test_matrix_comparison.py` (pytest) |
| Suites, layout pairs, resume | `tests.test_suites` |
| Promotion | `tests.test_reference_publication` |

Run the focused catalog tests with pytest (from the repository root):

```bash
python -m pytest regression_tests/tests/test_catalog.py
```

The remaining test groups are being migrated incrementally. Run only the
affected group:

```bash
cd regression_tests
python -m pytest tests/test_preparation.py tests/test_execution.py
```

For a cross-cutting change, run all synthetic tests:

```bash
cd regression_tests
python -m pytest tests
```

## Limitations

- Adaptive refinement may cross different thresholds between code revisions;
  adaptive golden-reference comparisons therefore do not promise identical
  meshes. Layout pairs within one race or cold matrix do require exact meshes.
- Runtime is recorded without a timing threshold.
- Promotion verifies technical evidence but physical acceptance remains human.
- Promotion inputs must all describe accepted results from the configured source
  bundle.

Run `python -m regression_tests --help` for the command synopsis.
