# Tutorial 5: add a real case and input bundle

[Tutorial 4](04-profile.md) selected suites across existing cases. To add a
different physical setup, declare its file roles, provide real inputs outside
Git, then select only the workflows whose physics apply. This example reduces
the existing [diverted case](../cases/diverted_case.json) to a 2D bootstrap and
warm check. It creates a template named `new_case`; supply your own real
geometry, equilibrium and parameter files. It does not introduce 3D or another
solver model.

## Declare the smallest case

Create `regression_tests/cases/new_case.json`:

```json
{
  "schema_version": 2,
  "description": "Real 2D configuration with adaptive bootstrap and warm baseline",
  "files": {
    "required": {
      "geometry": "geometry.geo",
      "equilibrium_magnetic_field": "equilibrium.h5",
      "equilibrium_current_density": "current_density.h5",
      "coarse_mesh": "mesh_adaptive_initial.msh",
      "cold_fixed_time_init_parameters": "param_cold_fixed_time_init.txt",
      "impurity_configuration": "impurity_model_w.nml"
    },
    "optional": {
      "bootstrap_reference": "reference_bootstrap.h5",
      "baseline_reference": "reference_baseline.h5"
    }
  },
  "workflows": {
    "bootstrap_adaptive": {},
    "baseline_warm": {}
  }
}
```

The key on the left of each mapping is a **role** used by the shared
[workflows](../workflows.json); the value is the filename in your prepared
source directory and eventual bundle. Empty workflow objects reuse the shared
definitions without case-specific changes. `required` means every bundle for
this case needs all six supplied inputs, including the coarse mesh needed by
`bootstrap_adaptive`. The bundle copies and records that mesh under the
`coarse_mesh` role; the workflow selects the role. The two `optional`
reference roles are produced later, so they are absent from the initial source
directory.

This example assumes the current 2D `NGammaTiTeNeutral` workflow, W radiation,
transport off and the shared cold recipe are appropriate for the new physical
case. Check those assumptions against your inputs. Use case-specific JSON
parameter overrides or a separate workflow where physics differs; do not copy
diverted input values merely to fit this template.

## Create a self-contained source bundle

Prepare a directory containing your real files under the declared filenames:

```text
prepared-inputs/
  geometry.geo
  equilibrium.h5
  current_density.h5
  mesh_adaptive_initial.msh
  param_cold_fixed_time_init.txt
  impurity_model_w.nml
```

The generic `positionFeketeNodesTri2D.h5` comes from the solver build, not
this case. The prepared files may be links to immutable source inputs;
`bundle create` follows them and copies their contents into a standalone
candidate. Do not place generated reference solutions in this source directory.

```bash
python -m regression_tests bundle readiness new_case --source /path/to/prepared-inputs
python -m regression_tests bundle create --case new_case \
  --source /path/to/prepared-inputs --output /path/to/new-source-bundle
python -m regression_tests bundle validate /path/to/new-source-bundle
```

Base readiness identifies supplied and workflow-producible roles and requires
all six source files. A workflow-specific readiness request also requires that
workflow's reference, so it can report `missing` before golden production.
The candidate manifest records roles and checksums. Bundle creation and later
publication refuse existing output directories.

## Produce and accept references

Add the minimal case campaign to [`golden.json`](../golden.json):

```json
"new_case": {
  "producers": ["bootstrap_adaptive", "baseline_warm"],
  "checks": ["warm"]
}
```

The order matters: bootstrap produces `bootstrap_reference`, which warm uses
as its restart; warm produces `baseline_reference`. The `warm` suite is then
rerun against the new candidate to validate the reference. Refresh also checks
producer outputs and convergence.

```bash
python -m regression_tests golden refresh new_case \
  --bundle /path/to/new-source-bundle --workspace /path/to/new-refresh --jobs 8
# Review new-refresh/refresh.json and its linked run/comparison reports.
python -m regression_tests golden publish /path/to/new-refresh \
  --output /path/to/new-golden --bundle-version 1.0.0 \
  --reason "Reviewed baseline for the new physical case" \
  --provenance "Your scientific review record"
```

Refresh builds the selected solver variant once. Publication is separate,
requires successful producers and validation, and records your actual reason
and provenance. The published golden contains physical copies of its inputs
and references. See the [README](../README.md#external-bundles-and-golden-references)
for the full review contract.

## Select the new case

Map `new_case` to the published bundle in your ignored `settings.local.json`;
reuse a suitable 2D build manifest as in the README. Then:

```bash
python -m regression_tests doctor warm --case new_case
python -m regression_tests check warm --case new_case
```

The focused `warm` suite is reused with `--case`; no duplicate suite definition
is needed. To include this case in a profile, add
`{"suite": "warm", "case": "new_case"}` to that profile's `checks` only if it
adds distinct evidence. Add transport, neutral, impurity or fixed-mesh roles
only when their workflows are actually selected; the case need not copy the
whole diverted catalog.
