# Tutorial 3: group related workflows into a suite

[Tutorial 1](01-warm-feature.md) built one warm check;
[Tutorial 2](02-layout-comparison.md) ran one workflow on several layouts. This
level groups physically related workflows. The existing diverted `impurities`
suite checks three distinct radiation choices on one `mpi4_omp4` layout.

## Three independent warm runs

All three workflows inherit the accepted W-bootstrap restart from
`baseline_warm`. Each starts from that same saved plasma solution; none starts
from another workflow's result.

| Workflow in [`workflows.json`](../workflows.json) | Physical choice | Its own reference |
| --- | --- | --- |
| `impurity_off` | Set `impurity_radiation: false` and remove the impurity configuration. | `impurity_off_reference` |
| `impurity_n` | Use `impurity_n_configuration` for N radiation. | `impurity_n_reference` |
| `impurity_nw` | Inherit the N workflow, but select `impurity_nw_configuration` for N+W radiation. | `impurity_nw_reference` |

Each workflow lists its reference role in `outputs` for golden refresh and
uses `fixed_hdf5` with the `feature_converged` tolerance for its ordinary check.
[`cases/diverted_case.json`](../cases/diverted_case.json) maps the two
configuration roles and three reference roles to files in the external bundle.
The reference files may be absent before production; the configuration files
are required when their workflows run.

## One scientific suite

In [`suites.json`](../suites.json), the suite is:

```json
"impurities": {
  "description": "Radiation off, N and N+W: convergence and own goldens",
  "workflows": ["impurity_off", "impurity_n", "impurity_nw"]
}
```

| Entry | Meaning |
| --- | --- |
| `impurities` | Suite name passed to `check`. |
| `description` | Text shown by `list suites`. |
| `workflows` | Run all three workflows and assess each result. |
| No `case` | Use the catalog default, `diverted_case`. |
| No `layouts` or `relations` | Use the suite catalog's default `mpi4_omp4` layout, with no layout pairs. |
| No `reference_comparisons` | With no pair relation, golden comparison is enabled for each workflow. |

The three solutions represent **different physics**, so this suite does not
compare them to one another. It checks convergence and each solution against
its own accepted reference. In [`golden.json`](../golden.json), the diverted
campaign lists all three workflows as producers and `impurities` as a validation
suite. The producers run after the bootstrap producer that supplies their
common restart.

## Run and inspect

With the published diverted bundle and MPI build from the
[README quick start](../README.md#quick-start-one-warm-check):

```bash
python -m regression_tests doctor impurities
python -m regression_tests check impurities
```

The printed `suite_summary.json` has three `results`, each with a comparison
report. There are no layout-pair `comparisons`. If one workflow fails, inspect
its own run directory and comparison report; the other two results remain
separate evidence.

## Apply the pattern

For a new feature family, give each physically distinct workflow a justified
starting state and its own reference. Reuse an inherited procedure when its
physics and restart still apply. Put the related workflows in one suite; add
layout relations only when parallel agreement adds distinct evidence. Declare
their input and reference roles in the case, add their workflows under
`golden.json` producers and the suite under checks, and order dependent
producers after their prerequisites. [Tutorial 4](04-profile.md) composes
suites into a profile across cases.
