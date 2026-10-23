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
The source switches off at **0.1 eV** and recovers the original law at **1 eV**.
These are physical electron temperatures, converted once to thresholds on
`U4/U1` using the run's reference temperature and `Mref`.

| Local state | Ohmic heating | State Jacobian |
| --- | --- | --- |
| Density $\le 0$ | Zero | Zero |
| Positive density, $T_e\le0.1$ eV, including negative electron energy | Zero; resistivity is not evaluated | Zero |
| Positive density, $0.1<T_e<1$ eV | $w(T_e)\eta(T_e)j_\phi^2$ | Derivative of both the weight and resistivity |
| Positive density, $T_e\ge1$ eV | Original temperature-dependent law | Derivative of that law |

For $x=(T_e-0.1\,\mathrm{eV})/(0.9\,\mathrm{eV})$, the transition weight is
$w=3x^2-2x^3$. The shared `smooth_transition` utility clips this cubic to zero
below $x=0$ and one above $x=1$, with zero slope at both endpoints. The source
and its state Jacobian are continuous at both temperature thresholds.

There is no resistivity floor or cap: the zero-source branch returns before
evaluating the singular law. The transition can exceed the former 1 eV capped
source within part of the interval. This is a low-temperature numerical
safeguard, not a limit on total Ohmic power or a model of cold current-carrying
plasma. It does not establish positivity of the conserved state. The
nonpositive-density rule remains an extension outside the admissible state.

Related conventions: [heat-flux limiter](heat-flux-limiter.md) and
[AMJUEL rates](amjuel-rates.md).
