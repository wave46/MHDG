# Wall recycling

## Physical model

With outward wall normal $\boldsymbol{\nu}$, the physical recycling condition is

$$
-\boldsymbol{\Gamma}_n\cdot\boldsymbol{\nu}
=R\,\boldsymbol{\Gamma}_p\cdot\boldsymbol{\nu}
+S-\kappa n_n-F_{\mathrm{abs}}.
$$

$R$ is the plasma recycling fraction, including TRIM reflection when enabled;
$S$ is inward puff, $\kappa n_n$ is pumping, and $F_{\mathrm{abs}}$ is neutral
wall absorption. Pump boundaries have their own recycling fraction.

The plasma particle flux includes parallel flow, pinch and perpendicular diffusion:

$$
\boldsymbol{\Gamma}_p
=n_pu_{\parallel p}\boldsymbol b+n_p\boldsymbol v_{\mathrm{pinch}}
-D_p\nabla_\perp n_p,
\qquad \nabla_\perp=(I-\boldsymbol b\boldsymbol b^T)\nabla.
$$

The neutral flux includes diffusion, optional pressure transport and, for
NeutralGamma, parallel flow:

$$
\boldsymbol{\Gamma}_n
=n_nu_{\parallel n}\boldsymbol b
-\phi D_n\mathsf P_n\left[
\nabla n_n+\lambda\frac{n_n}{T_i^{\mathrm{lim}}+T_s}\nabla T_i
\right].
$$

Here $u_{\parallel n}=0$ in the diffusive neutral model;
$\mathsf P_n=I$ for isotropic diffusion or $I-\boldsymbol b\boldsymbol b^T$
for perpendicular diffusion. $\phi=1$ when the neutral limiter is off;
$\lambda=0$ when pressure transport is off. $T_i^{\mathrm{lim}}$ and $T_s$
are the limited ion temperature and configured temperature suppression.
Neutral absorption is
$F_{\mathrm{abs}}=(1-r_n)n_nc_n/\sqrt{2\pi}$, where $r_n$ is the neutral
reflection fraction and $c_n$ the thermal speed based on $T_i^{\mathrm{lim}}$.
Plasma recycling and neutral reflection are separate settings.

## HDG realization

The element state is $U$, the independent face trace is $\widehat U$, and
$Q$ supplies element gradients. Define outward physical fluxes
$J_a=\boldsymbol{\Gamma}_a\cdot\boldsymbol{\nu}$ for $a=p,n$. The implemented
neutral boundary equation is

$$
J_n(\widehat U,Q)+\tau_n(n_n-\widehat n_n)
+R J_p(\widehat U,Q)+S-\kappa\widehat n_n-F_{\mathrm{abs}}(\widehat U)=0.
$$

Both the element face equation and BC use $\widehat U$ for wall convection,
state-dependent coefficients, limiter and losses, and $Q$ for gradients.
The limiter and stabilization coefficients remain lagged; stabilization uses
unlimited diffusion. Plasma diffusion is recycled; plasma stabilization
$R\tau_p(n_p-\widehat n_p)$ is excluded.

The general volume/face coefficient convention is described in
[Diffusion in HDG](../numerics/diffusion.md).

These are weak face equations. $U$ and $\widehat U$ need not coincide at
convergence. At unit recycling without puff or losses, the combined numerical
wall flux can retain the plasma stabilization correction; exact numerical
cancellation is not imposed by this physical-flux recycling condition.

If puff/pumping is relocated into elements, its wall terms are removed and
volume operators use $U$. Absorption remains at the wall. The Bohm momentum
target uses trace-derived sound speed times element density on its active
branch; its linearization follows those arguments, with suppression lagged.

## Reading the diagnostics

The **physical balance** uses prescribed wall exchange: the plasma particle
wall dataset contains convection plus pinch, whereas the neutral dataset
includes the full recycling/puff/loss expression above. The **discrete
residual** uses element numerical fluxes, including stabilization; the
**BC residual** measures the boundary equation mismatch. A small discrete
residual alone does not establish the physical particle balance. Saved
balances describe the last pre-solve assembly; saved fields are post-update.

Scope: 2D NGammaTiTe neutral/NeutralGamma/k family, diagonal
stabilization and Bohm/BohmPump/BohmPuff walls; shock diffusion is excluded.
Source: [physical flux coefficients](../../src/Models/NGammaTiTe/physics.f90),
[wall BC](../../src/HDG/hdg_BC.f90),
[element face flux](../../src/HDG/hdg_ComputeJacobian.f90),
[diagnostics](../../src/Utils/Diagnostics/balance_diagnostics_bc.f90).
