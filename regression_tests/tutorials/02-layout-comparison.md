# Tutorial 2: compare one workflow across layouts

[Tutorial 1](01-warm-feature.md) used one warm workflow on its default layout.
This level keeps one workflow and adds a comparison between parallel layouts.
The existing diverted `mixture_parallel` suite is the smallest example: it runs
`impurity_nw_short` with one and sixteen OpenMP threads. Both runs start from the
same saved plasma state and use the same N+W impurity inputs; only the layout
changes. The short workflow performs two Newton iterations.

## Suite declaration

In [`suites.json`](../suites.json), the entire suite is:

```json
"mixture_parallel": {
  "description": "Short N+W golden and OpenMP agreement",
  "workflows": ["impurity_nw_short"],
  "relations": ["openmp"],
  "reference_comparisons": true,
  "tolerance_profile": "race_step"
}
```

| Key | Meaning |
| --- | --- |
| `description` | Text shown by `list suites`. |
| `workflows` | Run the same short warm workflow on each selected layout. |
| `relations` | Generate layout pairs from [`layouts.json`](../layouts.json); `openmp` selects `serial_omp1` → `serial_omp16`. |
| `reference_comparisons` | Compare **each** layout's result with the accepted `impurity_nw_short` golden. |
| `tolerance_profile` | Use `race_step` from [`tolerances.json`](../tolerances.json) for the **pair** comparison. |

The suite names no explicit layout array. With a relation present, it considers
the declared layouts and runs only those in generated pairs. Here that means two
serial-executable runs and one OpenMP pair. The first run has one MPI rank and
one OpenMP thread; the second has one rank and sixteen threads.

There are two distinct comparisons: each run against the **golden**, and the
`serial_omp16` result against `serial_omp1`. The workflow's inherited
`comparison.cross_layout_profile` sets the golden tolerance when a run uses a
layout other than that workflow's default `mpi4_omp4`. The suite's
`tolerance_profile` sets the run-to-run pair tolerance. Neither creates another
golden file.

## Run and inspect

Use the published diverted bundle from the [README quick start](../README.md#quick-start-one-warm-check).
From the repository root:

```bash
python -m regression_tests check mixture_parallel --build
```

`--build` builds only the solver variant this suite needs: 2D
`NGammaTiTeNeutral/serial`. It needs the compiler environment described in the
README. If your selected build manifest already contains that variant, omit
`--build`. The command prints a `suite_summary.json` path. Its `results` contain
the two runs and their golden reports; `comparisons` contains the OpenMP pair
and its report. A failure in either kind of comparison fails the suite.

## Apply this to another feature

Keep the workflow and reference from Tutorial 1. Add a focused suite selecting
that workflow and the relation that isolates the parallel effect you want:

| Relation | Generated comparison in the current catalog | Existing example |
| --- | --- | --- |
| `openmp` | `serial_omp1` → `serial_omp16` | `mixture_parallel` |
| `mpi` | `serial_omp1` → `mpi4_omp1` | `transport_parallel` |
| `hybrid` | `mpi4_omp1` → `mpi4_omp4` | `transport_hybrid` |

Choose a pair tolerance for the behavior being probed, and decide whether each
run also needs a golden check. The three relations can be combined in one suite,
as `transport_hybrid` does, without writing explicit pairs or duplicating the
workflow. [Tutorial 3](03-scientific-suite.md) groups several workflows into
one scientific suite.
