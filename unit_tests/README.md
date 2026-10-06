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
outward normalized-flux direction.  The second target checks all four supported
poloidal-field conventions, noncanonical and percent-level field perturbations,
the same-alpha Grad--Shafranov toroidal current and its fitted `+1`/`-1` sign
convention, exact
nodal/quadrature topology caches, and cache generation refresh.  The
region-policy target checks
parsing, region inclusion, signed pinch orientation, null-normal suppression,
and the explicit legacy magnetic-field fallback.  Test executables are generated
in `lib/` and removed by `make clean`.  The taper target checks the disabled
slope, the default `0.7` factor beyond the LCFS, and lower-bound clipping.

For neutral thermal speed, wall absorption and flux limiting, use the 2D
`NGammaTiTeNeutral` or `NGammaTiTeNeutralGamma` build configuration:

```bash
cd lib
source Make.inc/init_vars_libs.sh
make check-neutral-flux-limiter
make check-neutral-wall-flux
```

The wall check calls the production temperature and wall-flux routines. It
checks Maxwellian normalization, temperature-floor behavior, finite-difference
derivatives at multiple steps, Euler homogeneity, five/six-component states,
small positive plasma density, zero neutral density, `Re_n=1` inactivity and
independence from the limiter's temperature selector. Function-section linking
keeps this check independent of solver runs and external physical case data.

The limiter check verifies complete diffusion/pressure flux, perpendicular
projection, the supplied-speed cap, epsilon regularization and diffusion/cap
floors. The wall check also verifies the production limiter's Ti/fixed speed
selection. Both targets print `PASS` on success and exit nonzero on failure;
wall assertion failures print the check name, actual value and expected value.
Finite-difference derivatives are checked with a scaled tolerance of `2e-7`;
Euler and common-scaling identities use `5e-13`. Temperature samples include
both sides of the soft limiter but avoid crossing its piecewise joins.

These are local kernel tests, not boundary assembly or reconvergence tests.
The existing `neutral_limiter_fixed` and `neutral_limiter_ti` harness workflows
check warm compatibility of limiter changes against their accepted references.
Wall `Re_n=1` compatibility and an active `Re_n=0.99` solver run are planned
when the boundary term is connected. Keep the existing build cache for the same
configuration; run `make clean` before switching model or build configuration.
