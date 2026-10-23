# Wall recycling

## Physical model

With outward wall normal $\boldsymbol{\nu}$, the physical recycling condition is

$$
-\boldsymbol{\Gamma}_n\cdot\boldsymbol{\nu}
=R\,\Gamma_{\mathrm{rec}}
+S-\kappa n_n-F_{\mathrm{abs}}.
$$

$R$ is the plasma recycling fraction, including TRIM reflection when enabled;
$S$ is inward puff, $\kappa n_n$ is pumping, and $F_{\mathrm{abs}}$ is neutral
wall absorption. $R$ uses `Re` at ordinary/puff walls and `Re_pump` at pump
walls, multiplied by the TRIM coefficient when enabled.

The volume/raw reconstructed plasma flux includes parallel flow, pinch and
perpendicular diffusion:

$$
\boldsymbol{\Gamma}_p
=n_pu_{\parallel p}\boldsymbol b+n_p\boldsymbol v_{\mathrm{pinch}}
-D_p\nabla_\perp n_p,
\qquad \nabla_\perp=(I-\boldsymbol b\boldsymbol b^T)\nabla.
$$

For Bohm/BohmPump/BohmPuff walls, define the outward normal convection/pinch
and raw diffusive fluxes

$$
\Gamma_{p,\mathrm{adv}}
=(n_pu_{\parallel p}\boldsymbol b+n_p\boldsymbol v_{\mathrm{pinch}})
\cdot\boldsymbol\nu,\qquad
\Gamma_{p,\mathrm{diff}}=-D_p\nabla_\perp n_p\cdot\boldsymbol\nu.
$$

The ion flux driving recycling is selected by `recycle_plasma_diffusion`
in `PHYS_LST`:

$$
\Gamma_{\mathrm{rec}}=
\begin{cases}
\Gamma_{p,\mathrm{adv}},
&\texttt{.false.}\ \text{(default)},\\
\Gamma_{p,\mathrm{adv}}+\Gamma_{p,\mathrm{diff}},
&\texttt{.true.}\ \text{(legacy)}.
\end{cases}
$$

The default recycles the prescribed Bohm convection/pinch outflow; `.true.`
restores the previous raw-diffusion recycling law. This choice changes neutral
recycling, not plasma or neutral transport.

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
Plasma recycling and neutral reflection are separate settings: $r_n$ uses
`Re_n`, or `Re_n_pump` at pump walls.

## HDG realization

The element state is $U$, the independent face trace is $\widehat U$, and
$Q$ supplies element gradients. Define outward normal fluxes
$\Gamma_{a,\nu}=\boldsymbol\Gamma_a\cdot\boldsymbol\nu$ and stabilization
corrections $\Gamma_{a,\mathrm{stab}}=\tau_a(n_a-\widehat n_a)$.
The numerical ion outflow and density wall condition are

$$
\widehat\Gamma_{p,\nu}
=\Gamma_{p,\mathrm{adv}}+\Gamma_{p,\mathrm{diff}}+\Gamma_{p,\mathrm{stab}},
\qquad
\Gamma_{p,\mathrm{diff}}+\Gamma_{p,\mathrm{stab}}=0
\quad\text{weakly on the wall}.
$$

Thus the prescribed ion outflow is $\Gamma_{p,\mathrm{adv}}$.
The neutral boundary equation is

$$
\Gamma_{n,\nu}(\widehat U,Q)+\Gamma_{n,\mathrm{stab}}
+R\,\Gamma_{\mathrm{rec}}+S-\kappa\widehat n_n-F_{\mathrm{abs}}(\widehat U)=0.
$$

Both the element face equation and BC use $\widehat U$ for wall convection,
state-dependent coefficients, limiter and losses, and $Q$ for gradients.
The limiter and stabilization coefficients remain lagged; stabilization uses
unlimited diffusion. The optional recycled diffusion coupling and its gradient
Jacobian are enabled together. Plasma stabilization is not recycled.

The general volume/face coefficient convention is described in
[Diffusion in HDG](../numerics/diffusion.md). These are weak face equations;
$U$ and $\widehat U$ need not coincide at convergence.

With the default `recycle_plasma_diffusion = .false.`, the net ion loss
after neutral recycling is

$$
\widehat\Gamma_{p,\nu}-R\,\Gamma_{p,\mathrm{adv}}
=(1-R)\Gamma_{p,\mathrm{adv}}
+\Gamma_{p,\mathrm{diff}}+\Gamma_{p,\mathrm{stab}}
\simeq(1-R)\Gamma_{p,\mathrm{adv}}.
$$

With the legacy setting `.true.`, it instead becomes

$$
\widehat\Gamma_{p,\nu}-R(\Gamma_{p,\mathrm{adv}}+\Gamma_{p,\mathrm{diff}})
=(1-R)(\Gamma_{p,\mathrm{adv}}+\Gamma_{p,\mathrm{diff}})
+\Gamma_{p,\mathrm{stab}}
\simeq(1-R)\Gamma_{p,\mathrm{adv}}-R\,\Gamma_{p,\mathrm{diff}}.
$$

Consequently, positive raw ion diffusion can increase neutral return without
increasing numerical ion outflow. A small
$\Gamma_{p,\mathrm{diff}}/\Gamma_{p,\mathrm{adv}}$ can materially change the
net loss when $1-R$ is small. At unit recycling without puff or losses, the
default law cancels the combined wall exchange weakly; the legacy law can
retain $\Gamma_{p,\mathrm{stab}}$.

If puff/pumping is relocated into elements, its wall terms are removed and
volume operators use $U$. Absorption remains at the wall. The Bohm momentum
target uses trace-derived sound speed times trace density on its active
branch; its linearization follows those arguments, with suppression lagged.
Only outward supersonic flow selects the free-outflow branch. Momentum
diffusion remains present on both branches, so the active target does not
impose exactly sonic trace velocity when its diffusive contribution is nonzero.

## Reading the diagnostics

The **physical balance** uses prescribed wall exchange: the plasma particle
wall dataset contains convection plus pinch, whereas the neutral dataset
includes the selected recycling/puff/loss expression above. The **discrete
residual** uses element numerical fluxes, including stabilization; the
**BC residual** measures the boundary equation mismatch. A small discrete
residual alone does not establish the physical particle balance. Saved
balances describe the last pre-solve assembly; saved fields are post-update.
A closed balance checks the selected boundary law; it does not establish that
the legacy return equals $R$ times the numerical ion outflow.

The choice is recorded in
`/simulation_parameters/physics/recycle_plasma_diffusion` in HDF5 output.
On restart, `param.txt` selects it; it is not restored from the checkpoint.
Existing inputs omitting the flag adopt the new default. Set `.true.` explicitly
when reproducing older recycling results or regression baselines.

Scope: 2D NGammaTiTe neutral/NeutralGamma/k family, diagonal
stabilization and Bohm/BohmPump/BohmPuff walls; shock diffusion is excluded.
Source: [physical flux coefficients](../../src/Models/NGammaTiTe/physics.f90),
[wall BC](../../src/HDG/hdg_BC.f90),
[element face flux](../../src/HDG/hdg_ComputeJacobian.f90),
[diagnostics](../../src/Utils/Diagnostics/balance_diagnostics_bc.f90).
