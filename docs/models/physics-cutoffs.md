# Physics cutoff conventions

Conventions for reviewing floors, caps and extensions of constitutive laws.
The adopted Ohmic rule is recorded below; other model-specific limits remain
in their respective guides.

## General conventions

- State each cutoff's physical quantity, units, threshold and affected law.
  A denominator tolerance is not a physical temperature or density floor.
- Retain the original law in its valid range. Define finite behavior outside
  that range, explicitly including zero and negative states. A small negative
  energy must not become an effectively zero temperature inside a singular law.
- Apply a coefficient cutoff to that coefficient's evaluation. It does not
  change the conserved state or establish physical admissibility.
- Preserve continuity at floors and caps. At the exact threshold, use the
  derivative from the frozen side. Treat invalid-state extensions separately.
- Differentiate the regularized law consistently. A frozen coefficient has
  zero derivative in its frozen arguments; other state-dependent factors in
  the full source or flux still contribute to its Jacobian.
- Define temperature cutoffs on temperature, including its density dependence,
  rather than on conserved energy alone. Equivalent normalizations must give
  the same physical threshold.

## Adopted Ohmic convention

For prescribed current density and effective charge,
$Q_\Omega=\eta(T_e)j_\phi^2$, with $\eta\propto Z_{\mathrm{eff}}T_e^{-3/2}$.
The electron-temperature floor for resistivity is **1 eV**.

| Local state | Ohmic heating | State Jacobian |
| --- | --- | --- |
| Density $\le 0$ | Zero | Zero |
| Positive density, $T_e\le1$ eV, including negative electron energy | $\eta(1\,\mathrm{eV})j_\phi^2$ | Zero |
| Positive density, $T_e>1$ eV | Original temperature-dependent law | Derivative of that law |

This freezes resistivity; it neither switches off heating in cold plasma nor
imposes a limit on total Ohmic power. The nonpositive-density rule is a
numerical extension outside the admissible plasma state.

Related conventions: [heat-flux limiter](heat-flux-limiter.md) and
[AMJUEL rates](amjuel-rates.md).
