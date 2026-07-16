# MHDG regression harness

This directory contains the tracked, non-sensitive contract for solver
regression tests. Physical case data, complete parameter files, restart
solutions, and golden outputs are distributed separately in an external case
bundle.

Status: contract only. The runner and HDF5 comparator are not implemented yet.

## Repository boundary

Tracked here:

- generic case identifiers and workflow definitions;
- the local-settings example;
- external-bundle schemas and examples;
- comparison tolerances;
- later, the runner, comparator, and their synthetic tests.

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
- `continuation_fixed`: restart on the same mesh while enabling a feature or
  changing one continuation parameter.
- `fixed_mesh_bootstrap`: start from the analytical initialization on the
  already refined mesh, without adaptivity.
- `cold_adaptive`: start from the analytical initialization on the coarse mesh
  and exercise refinement and projection.
- `archived_compare`: inspect or compare a stored historical result without
  routinely rerunning its source branch.

Only `legacy_fixed/warm` is defined as a routine executable workflow in this
first contract.

## External bundle contract

The bundle root contains a `manifest.json` conforming to
`schemas/bundle-manifest.schema.json`. Each physical file is recorded once in a
top-level artifact registry, with a bundle-relative path, SHA-256 checksum, and
size. Entries under `case_data` map the generic roles required by a case to
those reusable artifact identifiers.
A SHA-256 checksum is a 64-character fingerprint of a file's exact contents;
it detects a missing, substituted, or changed external artifact. It can be
generated with `sha256sum <file>`.

The minimum `legacy_fixed/warm` roles are:

- `mesh`;
- `equilibrium_magnetic_field`;
- `equilibrium_current_density`;
- `warm_parameters`;
- `transport_configuration`;
- `warm_restart`;
- `warm_reference`.

The same artifact may be referenced by multiple cases or workflows. Shared
meshes, equilibria, and reference solutions therefore do not need to be copied
or described more than once.

The bundle can be unpacked anywhere. A user supplies its root through a local
settings file based on `settings.example.env`.

## Initial comparison contract

The named profiles are in `tolerances.json`. For the repeated same-build,
four-MPI-by-four-OpenMP warm case, the initial requirements are:

- final Newton error no larger than `2e-4`;
- finite selected solution and transport datasets;
- exact mesh connectivity;
- mesh-coordinate absolute tolerance `1e-12`;
- per-equation relative L2 tolerance `1e-10`;
- per-equation normalized Linf tolerance `1e-9`.

Runtime is reported and may warn when it exceeds twice the reference median; it
is not initially a correctness failure.

The core comparator will use HDF5 directly and support both grouped and older
flat solution layouts. `HDG_postprocess` remains an optional richer layer for
mesh-independent adaptive comparisons.

## Planned user interface

These commands document the intended interface; they are not available until
the runner is implemented:

```bash
regression_tests/regression.sh help
regression_tests/regression.sh --help
regression_tests/regression.sh --settings /path/to/settings.env check-data
regression_tests/regression.sh --settings /path/to/settings.env suite warm
regression_tests/regression.sh --settings /path/to/settings.env suite parallelism
regression_tests/regression.sh --settings /path/to/settings.env run legacy_fixed fixed-bootstrap
regression_tests/regression.sh --settings /path/to/settings.env run legacy_fixed cold-adaptive
```

The `warm` suite is the canonical same-state warm restart. The `parallelism`
suite runs that inexpensive case using the selected serial, MPI, and OpenMP
layouts. The help command will list all suites, workflows, and layouts.

The runner will use isolated run directories and prebuilt executables. It will
not switch Git branches, rebuild the solver, or overwrite reference files.

## Planned run outputs

Each run directory will contain rendered input copies, the executed command,
stdout, the accepted HDF5 result, a machine-readable comparison report, and a
companion `run_metadata.json` file. This metadata file will record the
executable checksum, supplied solver revision/build description, MPI/OpenMP
layout, input/data checksums, and output checksum without embedding
machine-specific paths in tracked files.
