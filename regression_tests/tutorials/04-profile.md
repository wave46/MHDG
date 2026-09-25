# Tutorial 4: compose a profile across cases

[Tutorial 3](03-scientific-suite.md) grouped related workflows in one suite.
A profile selects suites in order. It adds no solver recipe or numerical
tolerance of its own. The existing `full` profile combines limited and diverted
evidence, so it shows how one selection spans cases.

## Read the `full` profile

In [`suites.json`](../suites.json), the profile is:

```json
"full": {
  "description": "Both topologies, full bootstraps, feature convergence and focused parallel checks",
  "checks": [
    {"suite": "bootstrap"}, {"suite": "bootstrap_limited"},
    {"suite": "warm"}, {"suite": "warm", "case": "legacy_case"},
    {"suite": "transport"}, {"suite": "transport", "case": "legacy_case"},
    {"suite": "transport_hybrid"}, {"suite": "transport_parallel", "case": "legacy_case"},
    {"suite": "adaptive_parallel"}, {"suite": "neutrals"}, {"suite": "neutral_parallel"},
    {"suite": "impurities"}, {"suite": "mixture_parallel"}, {"suite": "stored_field"},
    {"suite": "neutralgamma"}, {"suite": "initialization"}
  ]
}
```

| Entry | Meaning |
| --- | --- |
| `full` | Profile name passed to `doctor` or `check`. |
| `description` | Purpose shown by `list suites`. |
| `checks` | Ordered suite selections. Each suite retains its own workflows, layouts and comparison rules. |
| `suite` | Name of an existing suite; the profile does not redefine it. |
| `case` | Select the physical case for this entry when needed. |

For example, `{"suite": "warm"}` uses the catalog default `diverted_case`;
`{"suite": "warm", "case": "legacy_case"}` runs the same warm suite on limited
data. `bootstrap_limited`, `neutralgamma` and `initialization` already declare
`legacy_case` in their suite definitions, so their profile entries need no
`case` key. A repeated suite on a different case is a separate check.

`full` explicitly selects its evidence. It does not first run `routine` or
`routine-extended`, and it is not a union of every suite in the catalog. The
first entries exercise cold bootstraps and warm baselines on both topologies;
later entries cover transport, selected parallel paths, and distinct features.
The `bootstrap_limited` suite itself owns the limited adaptive/fixed endpoint
comparison.

## Select machine inputs and run

Because this profile spans two cases, `settings.local.json` must map both to
published golden bundles, as shown in
[`settings.example.json`](../settings.example.json). One `--bundle` argument
cannot select both cases. The selected build manifest must contain the solver
variants needed by the selected suites: serial and MPI 2D
`NGammaTiTeNeutral`, plus serial 2D `NGammaTiTeNeutralGamma`.

```bash
python -m regression_tests doctor full
python -m regression_tests check full
```

`doctor` checks catalog, bundle, build and machine prerequisites without
launching MHDG. `check` runs the suites in profile order and stops after a failed
suite. It prints a `profile_summary.json` path; that file records each
completed suite, its case, status and `suite_summary.json` path. Read a suite
summary for its workflow and comparison reports. `full` includes long cold
runs, so use focused suites while developing and run the profile when the
complete evidence is needed.

The profile reads published goldens; it does not produce or update them.
Golden refresh is a separate per-case campaign described in the
[README](../README.md#external-bundles-and-golden-references).

## Apply the pattern

Add an existing suite to a profile only when it contributes distinct evidence.
Use an explicit `case` for a suite reused on another topology; leave it out
when the suite or catalog default is correct. Keep the profile as an ordered
selection of suites rather than copying their workflow or layout definitions.
[Tutorial 5](05-real-case-bundle.md) uses the diverted case and bundle as
a template for authoring a new real-data case.
