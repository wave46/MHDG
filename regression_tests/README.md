# MHDG regression harness

This directory contains the tracked, non-sensitive contract for solver
regression tests. Physical case data, complete parameter files, restart
solutions, and golden outputs are distributed separately in an external case
bundle.

Status: external-data contract, candidate and golden bundle tools, isolated
warm and staged cold execution, fixed-mesh HDF5 comparison, and tracked
workflow-by-layout suites. Clean serial/parallel build orchestration is also
available.

## Repository boundary

Tracked here:

- generic case identifiers and workflow definitions;
- the local-settings example;
- external-bundle schemas and examples;
- generic serial, MPI, and OpenMP layouts;
- comparison tolerances and a library-independent fixed-mesh comparator;
- tracked suite definitions and orchestration;
- clean serial/parallel build orchestration and provenance recording;
- staged fixed-mesh and adaptive-mesh workflows.

Not tracked here:

- equilibria, meshes, geometry names, or experiment identifiers;
- solver parameter and model files containing physical case specifications;
- restart solutions, reference solutions, or run logs;
- machine-specific executable, data, or scratch paths.

## Initial cases

`legacy_fixed` is the routine integration case for branch `develop` at commit
`29f442db67bac169b2616f2cbe399289d137993c`. Its first workflow is a same-state
warm restart on a fixed mesh. The characterized run takes about ten seconds
with four MPI ranks and four OpenMP threads.

`historical_feature` identifies archived output from branch
`feautre/neutrals_pressure` at commit
`224fe346acb82c52a3efe5ee1fd6e35de43e3e4e`. It is comparison evidence, not a
routine solver run and not an unquestioned correctness oracle.

As cleaned features are accepted on `develop`, new case definitions and goldens
may be added without changing the external-data contract.

## Workflow vocabulary

- `warm_same_state`: restart from a converged solution without changing the
  physical or numerical settings.
- `staged_fixed_mesh`: start analytically on the refined mesh, then pass each
  stage output to the next stage without adaptivity.
- `staged_adaptive_mesh`: start analytically on the coarse mesh, refine during
  the first two stages, then finish the same continuations on that mesh.
- `archived_compare`: inspect or compare a stored historical result without
  routinely rerunning its source branch.

`legacy_fixed/warm`, `legacy_fixed/cold_fixed`, and
`legacy_fixed/cold_adaptive` are executable. The `cold_matrix` suite records
both cold workflows across every tracked layout; numerical verification is a
separate step so an overnight run is never repeated merely to compare it.

## External bundle contract

The bundle root contains a `manifest.json` conforming to
`schemas/bundle-manifest.schema.json`. Each physical file is recorded once in a
top-level artifact registry, with a bundle-relative path, SHA-256 checksum, and
size. Entries under `case_data` map the generic roles required by a case to
those reusable artifact identifiers.
A SHA-256 checksum is a 64-character fingerprint of a file's exact contents;
it detects a missing, substituted, or changed external artifact. The bundle
creator calculates checksums and sizes automatically; users do not enter them
manually.

The complete `legacy_fixed` case package contains these roles:

- `mesh`;
- `geometry`;
- `equilibrium_magnetic_field`;
- `equilibrium_current_density`;
- `warm_parameters`;
- `transport_configuration`;
- `warm_restart`;
- `warm_reference`.

The same artifact may be referenced by multiple cases or workflows. The
manifest can therefore reuse a mesh, equilibrium, or reference solution
without requiring a special `shared` directory or another physical copy.
The case bundle is created once; serial, MPI, and OpenMP layouts all read these
same artifacts and write only their own run outputs.

The bundle can be unpacked anywhere. A user supplies its root through a local
settings file based on `settings.example.env`.

The optional manifest field `bundle_class` distinguishes a hand-prepared
`candidate` from an accepted `golden` bundle. Both obey the same physical-file
contract and are executed by the same harness commands.

### Preparing case data

Prepare one private, flat directory outside the repository. Use these generic
names for the initial warm case:

```text
legacy_fixed/
├── param.txt
├── mesh.msh
├── geometry.geo
├── transport_model.nml
├── equilibrium.h5
├── current_density.h5
├── restart.h5
└── reference_mpi4_omp4.h5
```

The parameter, mesh, geometry, transport-model, equilibrium, and
current-density files specify the physical case. `restart.h5` is the converged
solution supplied to the solver. `reference_mpi4_omp4.h5` is the accepted
result produced after reconverging that restart with four MPI ranks and four
OpenMP threads. They are kept separate because a warm run can make small but
measurable changes to the starting solution.

Physical filenames before preparation may contain case-specific names, but the
prepared directory and all tracked examples use only these generic names. The
prepared entries may be regular files or symlinks to a private source archive;
the creator copies symlink targets into the final bundle.

To add the fixed-mesh and adaptive cold workflows, prepare these optional files
in the same directory:

```text
mesh_adaptive_initial.msh
param_cold_fixed_time_init.txt
param_cold_fixed_diffusion_reduction.txt
param_cold_fixed_continuation_01.txt ... continuation_05.txt
transport_cold_fixed_initial.nml
transport_cold_fixed_continuation_01.nml ... continuation_05.nml
```

The tracked stage names deliberately omit physical diffusion values. The
external parameter and transport files contain those numerical choices. The
sequence is analytical `time_init`, restart `diffusion_reduction`, then five
restart continuations ending at the warm reference state. Both workflows reuse
the parameter and transport files. `cold_fixed` starts on the refined
`mesh.msh`, while `cold_adaptive` starts analytically on the coarse
`mesh_adaptive_initial.msh` and carries each resulting mesh through the HDF5
restart. The renderer forces `rest_adapt = .false.` throughout `cold_fixed`.
For `cold_adaptive`, it forces `.true.` in `time_init` and
`diffusion_reduction`, then `.false.` in the continuations. All other
adaptivity and physical settings remain external.

These cold files are optional so existing warm-only candidate and golden
bundles remain valid. Either cold workflow requires the complete cold set.

### Creating a bundle

From the repository root, run:

```bash
regression_tests/regression.sh bundle create \
  --case legacy_fixed \
  --source /private/path/legacy_fixed \
  --output /private/path/mhdg_case_bundle
```

The command creates a candidate bundle and will:

- recognize the conventional filenames above;
- require every file listed above;
- refuse to replace an existing output path;
- copy the physical artifacts into a self-contained output bundle;
- calculate every size and SHA-256 checksum;
- generate `manifest.json` with the required role mappings;
- validate the completed bundle before publishing it at the output path;
- leave the prepared source directory unchanged.

The default bundle version is `1.0.0`. Supply another label when needed with
`--bundle-version VERSION`.

The generated physical layout inside the bundle is an implementation detail;
users only prepare the flat source directory. The example manifest documents
the generated data contract and is not intended for manual checksum editing.
`param.txt` is preserved exactly at this stage. The `prepare` command renders
active data paths and `save_folder` into a private run copy; it does not modify
the bundled parameter file.

After creation, copy `settings.example.env` to a private location and set:

```text
MHDG_REGRESSION_SETTINGS_VERSION=1
MHDG_REGRESSION_DATA_ROOT=/private/path/mhdg_case_bundle
```

Then verify the bundle from the repository root:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  check-data
```

The checker reads settings as data rather than sourcing them as shell code. It
verifies manifest fields, safe bundle-relative paths, artifact sizes and
SHA-256 checksums, role mappings, and the roles required by tracked case
definitions. It does not modify the bundle.

## Building regression executables

Build the optimized `NGammaTiTeNeutral` 2D serial and parallel executables from
the current checkout with:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  build
```

The build sources `lib/Make.inc/init_vars_libs.sh` by default. Set
`MHDG_ENVIRONMENT_SCRIPT` only when another absolute path is needed.
`MHDG_REGRESSION_BUILD_JOBS` defaults to 8 and `--jobs N` overrides it once.

Serial and parallel modes share the same `.o` and `.mod` files in `lib/`, so
the harness performs two sequential clean builds. It stores both executables,
the generic Fekete-node file, logs, `build_metadata.json`, and generated
`settings.env` under `MHDG_REGRESSION_BUILD_ROOT`. If that setting is omitted,
the default is `MHDG_REGRESSION_RUN_ROOT/builds`.

The metadata records the Git revision and dirty state, build commands,
environment-script checksum, toolchain versions, and executable checksums.
The generated settings preserve the original private data and run paths while
selecting the new executables.

## Preparing an isolated run

Add the private run root and the executable required by the selected layout.
MPI layouts also require the launcher:

```text
MHDG_REGRESSION_RUN_ROOT=/private/path/regression_runs
MHDG_SERIAL_EXECUTABLE=/absolute/path/to/serial/solver
MHDG_PARALLEL_EXECUTABLE=/absolute/path/to/parallel/solver
MHDG_MPI_LAUNCHER=mpirun.openmpi
```

Then prepare the characterized 4-by-4 layout:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  prepare legacy_fixed warm --layout mpi4_omp4
```

Preparation creates a timestamped run directory under
`MHDG_REGRESSION_RUN_ROOT`. It contains:

```text
legacy_fixed/warm/mpi4_omp4/<run-id>/
├── inputs/             # symlinks to read-only bundle artifacts
├── outputs/            # writable solver output directory
├── res/                # writable mesh/adaptivity workspace
├── positionFeketeNodesTri2D.h5  # symlink to generic solver data
├── param.txt           # rendered private run copy
└── run_plan.json       # command, environment, layout, and provenance
```

The rendered parameter file replaces `transport_model_path`, `field_path`,
`jtor_path`, `geometry_path`, and `save_folder`. The bundled parameter file is
not changed. The solver command is reported and recorded but is not executed.

The generic Fekete-node file is not physical case data and is therefore not
stored in the external bundle. The preparer requires
`positionFeketeNodesTri2D.h5` beside the selected executable and links it into
the run directory, matching the solver's fixed runtime filename.

The tracked layouts are `serial_omp1`, `mpi2_omp1`, `mpi2_omp4`, `mpi4_omp1`,
and `mpi4_omp4`. Every layout reuses the same bundle files.

For the complete fixed-mesh cold workflow, replace `warm` with `cold_fixed`:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  prepare legacy_fixed cold_fixed --layout mpi4_omp4
```

Its run root contains the shared reference and seven isolated stage directories:

```text
legacy_fixed/cold_fixed/mpi4_omp4/<run-id>/
├── inputs/reference.h5
├── run_plan.json
└── stages/
    ├── 01_time_init/
    ├── 02_diffusion_reduction/
    └── 03_continuation_01/ ... 07_continuation_05/
```

Each stage has its own rendered `param.txt`, inputs, outputs, writable `res/`
mesh workspace, runtime file, and `run_plan.json`. Restart links are deliberately
absent during preparation because they are created only after the preceding
stage has completed.

## Executing runs

If the executable depends on library paths exported by the build setup, add
its absolute path to the private settings file:

```text
MHDG_ENVIRONMENT_SCRIPT=/absolute/path/to/lib/Make.inc/init_vars_libs.sh
```

The script is sourced in a child Bash process. Its exported environment is used
for MHDG, after which the selected layout sets `OMP_NUM_THREADS`, `OMP_PLACES`,
and `OMP_PROC_BIND`. If the script setting is omitted, the runner otherwise
inherits the environment from the calling shell.

Use `run` instead of `prepare` to create the run directory and execute its
recorded command:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  run legacy_fixed warm --layout mpi4_omp4
```

The solver runs from the isolated directory with the layout's OpenMP thread
count. MPI commands use Open MPI's `--bind-to core --map-by slot:PE=THREADS` so
each rank receives exclusive cores for its threads. `OMP_PLACES=cores` and
`OMP_PROC_BIND=spread` distribute those threads within the assigned cores.
Standard output and error are written to `stdout.log` and `stderr.log`.
`run_metadata.json` records the command, environment, exit status, runtime,
executable checksum, optional revision/build description, and checksums of all
files produced under `outputs/`.

A run is `completed` only when the solver exits successfully, produces at least
one HDF5 file, and does not report a fatal source/destination file-opening error.
The last condition catches serial Fortran stops that return exit code zero.
Generic Newton `Error:` diagnostics are not treated as failures. This status
describes execution only; comparison records numerical acceptance separately.
Existing run directories are never reused or overwritten.

Execute the fixed-mesh cold sequence with:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  run legacy_fixed cold_fixed --layout mpi4_omp4
```

To launch both complete cold workflows sequentially for an overnight run, list
both workflow IDs in one command:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  run legacy_fixed cold_fixed cold_adaptive \
  --layout mpi4_omp4 --run-id overnight-01
```

The adaptive workflow still starts analytically, refines on the `time_init` to
`diffusion_reduction` restart, and then completes the same five continuations
on the resulting mesh. If one workflow fails, its later stages stop, but the
next workflow in the command is still attempted.

For the complete two-workflow by five-layout overnight matrix, use the tracked
suite instead:

```bash
regression_tests/regression.sh \
  --settings /private/path/regression-settings.env \
  suite cold_matrix --run-only --run-id overnight-01
```

`--run-only` deliberately defers numerical comparison. All ten result paths
are retained in the suite summary for the later verification pass. The summary
is updated after every cell; after an interruption, repeat the same command
with `--resume` to skip all cells already recorded there.

`time_init` runs without a restart. The runner waits for each solver process,
selects its final HDF5 output, links that result as the next stage's restart,
and only then starts the next stage. Any failed stage stops the sequence; its
logs and all earlier stage results remain available, while later stages are
recorded as `not_run`. The workflow-level `stdout.log` and `stderr.log` point
to the last attempted stage, and `run_metadata.json` summarizes the sequence.

The completed workflow root is compatible with the fixed-mesh comparator. Until
cold tolerances are characterized, select a profile explicitly for inspection:

```bash
regression_tests/regression.sh compare \
  /path/to/completed/cold_fixed/run --profile fixed_same_layout
```

The adaptive output is collected with full provenance but is not passed to the
fixed-mesh HDF5 comparator. Its mesh-independent comparison remains a separate
follow-up.

## Initial comparison contract

The named profiles are in `tolerances.json`. For the repeated same-build,
four-MPI-by-four-OpenMP warm case, the initial requirements are:

- final Newton error no larger than `2e-4`;
- finite selected solution and transport datasets;
- exact mesh connectivity;
- mesh-coordinate absolute tolerance `1e-12`;
- per-equation relative L2 tolerance `1e-10`;
- per-equation normalized Linf tolerance `1e-9`.

Normalized Linf is the largest absolute pointwise difference divided by the
largest absolute reference value for that equation.

Cross-layout comparisons initially allow relative L2 `5e-8` and normalized
Linf `1e-6`. These provisional limits cover the characterized serial, two-rank,
and four-rank reduction orderings and should be revisited with more runs and
solver builds.

Runtime is recorded but is not initially a correctness failure. Suite-level
runtime warnings require a characterized reference median and remain planned.

The core comparator uses HDF5 directly and supports both grouped and older flat
solution/mesh layouts. `HDG_postprocess` remains an optional richer layer for
mesh-independent adaptive comparisons.

The provisional cold profiles are based on the initial serial and two-rank
overnight runs. Fixed cold runs allow relative L2 and normalized Linf errors up
to `1e-5`. Adaptive runs require full common-point coverage, solution errors of
at most `0.05` relative L2 and `0.1` normalized Linf, and gradient errors of at
most `0.25` and `0.2`, respectively. Both require a final Newton error no larger
than `2e-4`. These tolerances characterize the legacy workflow, including its
known mesh sensitivity; they are not general accuracy targets.

## Verifying a saved matrix

Verification consumes the run-only suite summary and never launches the solver:

```bash
PYTHON=/path/to/hdg-postprocess-python \
PYTHONPATH=/path/to/HDG_postprocess \
regression_tests/regression.sh suite-verify \
  /path/to/suites/cold_matrix/overnight-01/suite_summary.json
```

Fixed workflows use the library-independent HDF5 comparator. Adaptive
workflows use `HDG_postprocess` to evaluate the reference and candidate at four
interior points per reference triangle. Each run receives its own
`comparison.json`; the matrix result is written beside the source summary as
`verification_summary.json`. A failed comparison does not prevent later cells
from being checked.

## Comparing a completed run

Run the comparator on the directory printed by `run`:

```bash
regression_tests/regression.sh compare /private/path/to/completed/run
```

It selects the final HDF5 save reported in `stdout.log`, falling back to the
unique output without a `_NNNN` time-save suffix. `--candidate` and
`--reference` can resolve an ambiguous or manual comparison. The selected
tolerance profile comes from the tracked workflow; non-default layouts use
`fixed_cross_layout`.

The comparator checks:

- the final `Error:` value from `stdout.log`;
- exact `T`, `Tlin`, and `Tb` connectivity and tolerance-based `X` coordinates;
- finite, per-equation `u`, `q`, and `u_tilde` values;
- transport-1D coefficient and profile datasets present in the reference.

It prints a short pass/fail and worst-error summary, writes `comparison.json` in
the run directory, and returns zero only when all checks pass. Use
`--report /path/report.json` to write it elsewhere.

## Available user interface

```bash
regression_tests/regression.sh help
regression_tests/regression.sh --help
regression_tests/regression.sh bundle create --case legacy_fixed --source /path/to/case --output /path/to/bundle
regression_tests/regression.sh --settings /path/to/settings.env build
regression_tests/regression.sh --settings /path/to/settings.env check-data
regression_tests/regression.sh --settings /path/to/settings.env prepare legacy_fixed warm --layout mpi4_omp4
regression_tests/regression.sh --settings /path/to/settings.env run legacy_fixed warm --layout mpi4_omp4
regression_tests/regression.sh --settings /path/to/settings.env run legacy_fixed cold_fixed --layout mpi4_omp4
regression_tests/regression.sh --settings /path/to/settings.env run legacy_fixed cold_fixed cold_adaptive --layout mpi4_omp4 --run-id overnight-01
regression_tests/regression.sh compare /path/to/completed/run
regression_tests/regression.sh --settings /path/to/settings.env suite warm
regression_tests/regression.sh --settings /path/to/settings.env suite warm --build
regression_tests/regression.sh --settings /path/to/settings.env suite warm_parallelism
regression_tests/regression.sh --settings /path/to/settings.env suite cold_matrix --run-only --run-id overnight-01
regression_tests/regression.sh --settings /path/to/settings.env suite cold_matrix --run-only --run-id overnight-01 --resume
regression_tests/regression.sh suite-verify /path/to/suite_summary.json
regression_tests/regression.sh --settings /path/to/settings.env bundle promote /path/to/suite_summary.json --output /path/to/golden --bundle-version VERSION
regression_tests/regression.sh golden-check
regression_tests/regression.sh golden-check --build
regression_tests/regression.sh golden-check warm_parallelism
```

The regression commands require Python 3 and the packages in `requirements.txt`;
install them into the selected Python environment with:

```bash
python -m pip install -r regression_tests/requirements.txt
```

The existing `hdg-postprocess-py312` environment already provides `jsonschema`,
NumPy, and h5py. The comparator does not depend on `HDG_postprocess` itself. Set
the `PYTHON` environment variable when a non-default interpreter should run the
command:

```bash
PYTHON=/path/to/environment/bin/python \
  regression_tests/regression.sh --settings /path/to/settings.env check-data
```

## Running suites

The canonical warm suite runs and compares the characterized 4-by-4 layout:

```bash
regression_tests/regression.sh --settings /path/to/settings.env suite warm
```

The `warm_parallelism` suite exercises the warm workflow across all tracked
serial, MPI, and OpenMP layouts:

```bash
regression_tests/regression.sh --settings /path/to/settings.env suite warm_parallelism
```

The `cold_matrix` suite executes `cold_fixed` and `cold_adaptive` for each of
the five layouts. Run it without comparisons while producing the golden
overnight evidence:

```bash
regression_tests/regression.sh \
  --settings /path/to/settings.env \
  suite cold_matrix --run-only --run-id overnight-01
```

Each suite validates the bundle once, then gives every workflow/layout cell a
distinct isolated run directory. A failed cell does not prevent later cells
from running. The command returns nonzero if any requested run or comparison
fails and writes
`suite_summary.json` under
`MHDG_REGRESSION_RUN_ROOT/suites/SUITE/RUN_ID/`. Use `--run-id ID` when a
stable label is useful; otherwise a UTC timestamp is generated. The summary is
updated after every completed cell. Reuse that identifier with `--resume` to
continue without repeating recorded cells.

After the run-only matrix finishes, pass its summary to `suite-verify`. This is
a separate command specifically so comparison logic or tolerances can be rerun
without repeating the overnight solver calculations.

Without `--build`, the executor uses the prebuilt executable paths from the
settings file. Add `--build` to a suite command to build the current
checkout first. It never switches Git branches or overwrites reference files.

## Candidate and golden test flows

There are two bundle sources, but only one execution and comparison path.

1. A candidate bundle is created from the private hand-prepared directory with
   `bundle create`. Use it to reproduce manual regressions, characterize
   layouts, and decide whether the inputs and canonical result are suitable.
2. A golden bundle is an immutable copy of that complete candidate package in
   which the passing canonical result has become `warm_reference`. Use it for
   routine comparisons of newly built executables.

After the candidate's full suite passes and has been checked against the
manual regressions, promote it explicitly:

```bash
regression_tests/regression.sh \
  --settings /path/to/candidate-settings.env \
  bundle promote /path/to/suite_summary.json \
  --output /path/to/golden-bundles/legacy_fixed_v1 \
  --bundle-version 1.0.0-golden.1
```

Promotion requires a passing suite containing the workflow's canonical
layout, currently `mpi4_omp4`. It verifies that the canonical solution and its
comparison inputs have not changed, copies the complete source bundle, replaces
only the reference artifact, updates its checksum, and validates the finished
bundle. It never modifies the candidate bundle or replaces an existing output.

The golden bundle therefore contains the warm starting solution together with
the mesh, equilibrium, geometry, parameters, transport configuration, and
accepted reference. Canonical run metadata, comparison, logs, suite summary,
executable checksum, and optional solver revision/build description are stored
under `provenance/golden_reference/` and registered in the manifest.

For a convenient routine check, copy the settings example to the default
untracked location:

```bash
cp regression_tests/settings.example.env regression_tests/golden.local.env
```

Fill its executable, run-root, launcher, and environment-script paths, and set:

```text
MHDG_REGRESSION_DATA_ROOT=/path/to/golden-bundles/legacy_fixed_v1
```

Then the short command validates that the selected bundle is golden and runs
the canonical four-by-four warm suite:

```bash
regression_tests/regression.sh golden-check
```

Build the current checkout first and automatically use its generated settings
and provenance with:

```bash
regression_tests/regression.sh golden-check --build
```

Run all characterized layouts with:

```bash
regression_tests/regression.sh golden-check warm_parallelism
```

Append `--build` to rebuild before the full-layout suite.

`regression_tests/golden.local.env` is ignored by Git. A settings file stored
elsewhere can be selected with `MHDG_REGRESSION_GOLDEN_SETTINGS` or the usual
`--settings FILE` option.

The explicit `--settings FILE suite SUITE` interface remains available for
handmade candidate bundles and other development runs. Both interfaces use the
same preparation, execution, and comparison implementation.

Compiled binaries remain in the private build root; they are not stored in the
portable golden bundle.

## Run outputs and provenance

Execution adds stdout, stderr, produced HDF5 files, and `run_metadata.json` to
the prepared directory. Comparison adds `comparison.json`. Automated builds
supply the solver revision, build description, and build-manifest checksum;
prebuilt executables may still provide revision and description manually in
the private settings file.
