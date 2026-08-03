# Focused unit tests

This directory contains small deterministic checks for numerical machinery that
can be exercised without launching a full MHDG case.  Production-equilibrium
and end-to-end coverage remains in `regression_tests`.

Run the magnetic-topology checks from `lib/` after loading the build
environment:

```bash
source Make.inc/init_vars_libs.sh
make check-magnetic-topology
```

The test executable is generated in `lib/` and removed by `make clean`.
