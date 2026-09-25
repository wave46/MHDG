# Tutorial 1: one warm feature check

This tutorial follows the existing diverted `stored_field` check. It needs the
build and published diverted bundle from the
[README quick start](../README.md#quick-start-one-warm-check), but requires no
catalog edits. It is the first level: one workflow in a one-workflow suite. Later
levels add layout comparisons, combine workflows, and assemble a profile.

The check starts from the accepted bootstrap solution, reconverges the 2D
`NGammaTiTeNeutral` model, and compares the result with its own accepted warm
reference. Its one physical change from the baseline warm workflow is
`compute_from_flux: false`: MHDG uses the stored magnetic and current fields.

## How the check is declared

In [`workflows.json`](../workflows.json), `stored_field` inherits the restart,
model, inputs and default `mpi4_omp4` layout from `baseline_warm`. It declares
only its changed behavior and its reference:

```json
"stored_field": {
  "extends": "baseline_warm",
  "description": "Converge using the stored magnetic/current fields with an independent reference",
  "reference": "stored_field_reference",
  "outputs": ["stored_field_reference"],
  "comparison": {
    "method": "fixed_hdf5",
    "profile": "feature_converged",
    "cross_layout_profile": "fixed_cross_layout"
  },
  "parameter_overrides": {"compute_from_flux": false}
}
```

| Key | Meaning |
| --- | --- |
| `extends` | Inherit the bootstrap restart, inputs, model, default layout and warm-run settings from `baseline_warm`. |
| `description` | Explain the workflow in catalog listings. |
| `reference` | Input role for the accepted solution that `check` compares with this run. |
| `outputs` | Roles that golden refresh fills from this producer's output. |
| `comparison.method` | Use direct HDF5 comparison on the fixed mesh. |
| `comparison.profile` | Use `feature_converged` limits from `tolerances.json` for the reference comparison. |
| `comparison.cross_layout_profile` | Limits for a run on a nondefault layout against its golden; unused by this default-layout suite. |
| `parameter_overrides.compute_from_flux` | Use stored magnetic/current fields instead of reconstructing them from flux. |

`reference` and `outputs` repeat the same role here because the producer writes
the very file that a later check reads. They describe opposite directions of the
golden flow, not two copies of the solution. Other producers can declare more than
one output role, such as a solution and its mesh.

In [`cases/diverted_case.json`](../cases/diverted_case.json),
`files.optional.stored_field_reference` maps this role to the external
`reference_stored_field.h5`. It may be absent before golden production, but
`check` requires it. `workflows.stored_field: {}` enables the shared workflow
for the diverted case without changing its settings.

In [`suites.json`](../suites.json), the focused suite is just:

```json
"stored_field": {
  "description": "Stored magnetic/current fields: convergence and own golden",
  "workflows": ["stored_field"]
}
```

The suite's `description` appears in `list suites`; its `workflows` array
selects the one workflow to run. It adds no second layout or pair comparison.
In [`golden.json`](../golden.json), `producers` runs `stored_field` and collects
its declared `outputs`; `checks` then runs the `stored_field` suite against the
candidate. Once published, that bundle supplies the accepted reference for an
ordinary `check`.

## Run it

From the repository root, with the quick-start settings in place:

```bash
python -m regression_tests doctor stored_field
python -m regression_tests prepare diverted_case stored_field
python -m regression_tests check stored_field
```

`doctor` checks prerequisites. `prepare` renders the inputs and `run_plan.json`
without launching MHDG; inspect the generated parameter file to confirm the
switch and restart. `check` launches the workflow, verifies output identity and
Newton convergence, and compares its final mesh and fields with its own golden.
It prints the path to `suite_summary.json`, which links its comparison report.
For an exploratory launch without golden comparison, use
`python -m regression_tests run diverted_case stored_field`.

## Reuse the pattern for a new warm feature

1. Choose a saved starting state and the one physical change to exercise. Derive
   a workflow from `baseline_warm` only if its restart and other inherited
   physics are appropriate. Give the result a distinct reference role.
2. Declare that role and enable the workflow in the relevant case. Add a
   one-workflow suite with a comparison method and tolerance supported by the
   expected numerical behavior.
3. Add the workflow to that case's producers and the suite to its validation
   checks in `golden.json`. Use `prepare` and `run` to inspect the new behavior
   before accepting a reference. A new reference needs a reviewed
   [golden refresh and publish](../README.md#external-bundles-and-golden-references);
   refresh runs the case's declared producers, not just the new workflow.
4. Run `check YOUR_SUITE` against the published golden. Keep the suite focused;
   layout comparisons and multi-workflow suites are separate next steps.
