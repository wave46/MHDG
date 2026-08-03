# Focused unit tests

This directory contains small deterministic checks for numerical machinery that
can be exercised without launching a full MHDG case.  Production-equilibrium
and end-to-end coverage remains in `regression_tests`.

Run the magnetic-topology checks from `lib/` after loading the build
environment:

```bash
source Make.inc/init_vars_libs.sh
make check-magnetic-topology
make check-magnetic-geometry
make check-transport-region-policy
make check-transport-taper
```

The topology target also checks that reversing the psi convention preserves the
outward normalized-flux direction.  The second target checks the signed
poloidal-field fit, percent-level field perturbations, exact nodal/quadrature
topology caches, and cache generation refresh.  The region-policy target checks
parsing, region inclusion, signed pinch orientation, null-normal suppression,
and the explicit legacy magnetic-field fallback.  Test executables are generated
in `lib/` and removed by `make clean`.  The taper target checks the disabled
slope, the default `0.7` factor beyond the LCFS, and lower-bound clipping.
