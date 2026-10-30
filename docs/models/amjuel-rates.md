# AMJUEL rates and source linearization

## Physical coefficients

The NGammaTiTe family with `AMJUELSPLINES` uses AMJUEL's rate coefficient
$\langle\sigma v\rangle(T_e,n_e)$ and energy-loss coefficient
$\langle\sigma vE\rangle(T_e,n_e)$. Their tabulated units are
cm³/s and eV cm³/s, respectively. The H.4 fit is

$$
\ln\langle\sigma v\rangle
=\sum_{i=0}^{8}\sum_{j=0}^{8}\alpha_{ij}
 (\ln T_e)^i(\ln\widetilde n_e)^j,
\qquad \widetilde n_e=\frac{n_e}{10^{14}\,\mathrm{m}^{-3}},
$$

with $T_e$ inserted in eV. H.10 uses the same fit form for
$\langle\sigma vE\rangle$, with separate coefficients. See
[AMJUEL §§II.3.4 and II.3.10](https://www.eirene.de/Documentation/amjuel.pdf).
MHDG converts cm³ to m³ and shifts the polynomial coefficients once to use
normalized density, temperature and rate units; constant unit conversions
leave the logarithmic slopes unchanged.

The particle source for singly charged hydrogenic plasma is
$\mathcal S_p=n_e n_n\langle\sigma v\rangle_{\mathrm{iz}}
-n_e n_p\langle\sigma v\rangle_{\mathrm{rec}}$.
The energy fits enter the corresponding density-weighted electron losses;
the recombination potential-energy correction is assembled separately.

## Continuation outside the fit interval

MHDG currently uses $0.1\le T_e\le20000$ eV and
$10^{14}\le n_e\le10^{22}$ m⁻³. Density is clamped at its nearest boundary.
Outside the temperature interval, the coefficient is continued as

$$
\langle\sigma v\rangle_{\mathrm{ext}}(T_e,n_e)
=\langle\sigma v\rangle(T_b,n_e)
 \left(\frac{T_e}{T_b}\right)^{s_b(n_e)},
\qquad
s_b=\left.\frac{\partial\ln\langle\sigma v\rangle}
                       {\partial\ln T_e}\right|_{T_b}.
$$

Thus, at fixed temperature and density inside its interval,

$$
\frac{\partial\ln\langle\sigma v\rangle_{\mathrm{ext}}}{\partial\ln n_e}
=\left.\frac{\partial\ln\langle\sigma v\rangle}{\partial\ln n_e}\right|_{T_b}
+\ln\!\left(\frac{T_e}{T_b}\right)
 \left.\frac{\partial^2\ln\langle\sigma v\rangle}
 {\partial\ln T_e\,\partial\ln n_e}\right|_{T_b}.
$$

Both the boundary value and its temperature slope depend on density. The
same expressions apply to the energy-loss coefficient. Outside density
bounds, the direct density derivative is zero; temperature dependence remains.

## Helpers and linearization terms

The [physics helpers](../../src/Models/NGammaTiTe/physics.f90) map as follows;
the polynomial derivatives in the middle rows are evaluated before continuation.

| Helper | Quantity |
| --- | --- |
| `compute_2D_eirene_rate` | Coefficient with continuation and the existing zero-rate cutoff |
| `compute_2D_logeirene_rate` | $\ln\langle\sigma v\rangle$ polynomial |
| `compute_dlogeirene_2D_dlogte_rate` | $\partial\ln\langle\sigma v\rangle/\partial\ln T_e$ |
| `compute_dlogeirene_2D_dlogne_rate` | $\partial\ln\langle\sigma v\rangle/\partial\ln n_e$ |
| `compute_d2logeirene_2D_dlogtedlogne_rate` | Mixed derivative in the continuation formula |
| `compute_2D_eirene_rate_du` | Continued coefficient differentiated with respect to $U$ |

For positive states away from floors, $n_e\propto U_1$ and
$T_e\propto U_4/U_1$. The conservative-variable chain rule is

$$
\frac{\partial\langle\sigma v\rangle}{\partial U_1}
=\frac{\langle\sigma v\rangle}{U_1}
 \left(\frac{\partial\ln\langle\sigma v\rangle}{\partial\ln n_e}
      -\frac{\partial\ln\langle\sigma v\rangle}{\partial\ln T_e}\right),
\qquad
\frac{\partial\langle\sigma v\rangle}{\partial U_4}
=\frac{\langle\sigma v\rangle}{U_4}
 \frac{\partial\ln\langle\sigma v\rangle}{\partial\ln T_e}.
$$

The `compute_sigmaviz` / `compute_sigmavrec` wrappers select the particle fits;
`compute_sigmavEiz` / `compute_sigmavErec` select the energy fits. Their
`compute_d..._dU` partners return the derivatives above, with the appropriate
rate normalization.

For a source $\mathcal S=N(U)\langle\sigma v\rangle(U)$, both factors enter:
$\mathcal J_{\mathcal S}=\langle\sigma v\rangle\,\partial_U N
+N\,\partial_U\langle\sigma v\rangle$.
In [`assemblyNeutral`](../../src/HDG/hdg_ComputeJacobian.f90), these derivatives
enter `Sn` with the left-hand-side sign, $\mathrm{Sn}=-\mathcal J_{\mathcal S}$;
the affine term is $\mathrm{Sn0}=\mathcal J_{\mathcal S}U-\mathcal S$ and
is subtracted from the RHS. Together they represent the linearized source,
retaining its value at the current iterate.

Atomic prefactors use a $10^{-20}$ dimensionless floor on density and ion-energy
arguments as a numerical zero safeguard.
