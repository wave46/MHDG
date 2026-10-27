# Smooth Bohm momentum boundary condition

At 2D temperature-model Bohm, BohmPump and BohmPuff faces, angle and energy
reduce the target parallel velocity; signed outward Mach controls release to
free flow. The trace state is $\widehat U$, the element state is $U$, and $Q$
contains element gradients. These are weak face conditions.

## Target and weights

With $s_n=\operatorname{sign}(b\cdot\nu)$, define

$$
M=\frac{s_n\widehat U_2}{\widehat U_1c_s(\widehat U)},\qquad
E=\max(\widehat U_3,10^{-20}),\qquad
\beta=\frac{E}{E+E_{\rm th}}.
$$

`bohm_energy_thresh` is **already dimensionless**: it is read directly in
the same units as $\widehat U_3$, without further input scaling. It is an ion
**total energy-density** scale, including parallel kinetic energy, not an eV
cutoff or a density cutoff separately.

With the current physical normalization, $n_{\rm ref}=10^{19}$ m$^{-3}$,
$u_{\rm ref}\simeq1.384\times10^4$ m/s and
$E_{\rm ref}=m_i u_{\rm ref}^2/e\simeq4.00$ eV,

$$
\widehat U_3=\frac{n_i}{n_{\rm ref}}\left[
\frac{3T_i[\mathrm{eV}]}{2E_{\rm ref}[\mathrm{eV}]}
+\frac12\left(\frac{u_\parallel}{u_{\rm ref}}\right)^2\right].
$$

For example, $E_{\rm th}=2\times10^{-4}$ corresponds at zero flow to
$n_iT_i\simeq5.34\times10^{15}$ eV m$^{-3}$: about **0.053 eV at
$10^{17}$ m$^{-3}$**, or **0.1 eV at $5.3\times10^{16}$ m$^{-3}$**.
At $10^{17}$ m$^{-3}$ the total energy is **0.080 eV per ion** at this
threshold. At zero flow it equals $3T_i/2$, explaining why this energy
estimate and the corresponding ion temperature differ.
At this energy density $\beta=1/2$; parallel flow increases the energy used
by the suppression. These estimates apply to dimensioned plasma cases, not
MMS cases which bypass physical normalization.

The threshold must be finite and nonnegative; `bohmth` must be finite.
The numerical floor $10^{-20}$ only prevents a zero denominator and has
zero derivative at and below its clipping point. $E_{\rm th}=0$ gives $\beta=1$.

For the clipped cubic $S(z)=z^2(3-2z)$ on $(0,1)$, zero below and one above,

$$
A=S\!\left(\frac{|b\cdot\nu|}{\texttt{bohmth}}\right),\qquad
w=1-S\!\left(\frac{M-1/\beta}{0.20}\right),
$$

$$
p_0=w(1-\beta A),\qquad p_s=w\beta A,\qquad p_F=1-w.
$$

The coefficients are dimensionless weights, not independent constants:

| Symbol | Physical or boundary-condition meaning |
|---|---|
| $A$ | Angular reduction of the target speed: zero at tangency, one above the angular threshold |
| $\beta$ | Energy-density reduction of the target speed; approaches zero in depleted, low-energy plasma and one at high energy |
| $w$ | Strength of the imposed condition: one below Mach release, decreasing to zero across the outward supersonic band |
| $p_0$ | Zero-parallel-momentum part of the imposed stabilization target |
| $p_s$ | Outward sonic part of that target, referenced to $s_nc_s\widehat U_1$ |
| $p_F$ | Free-flow part, referenced to element momentum $U_2$ |

$p_0$ and $p_s$ together lower the imposed speed to $s_n\beta A c_s$.
They are fractions of the stabilization condition, not particle populations.
The retained diffusion term belongs to the complete boundary equation below.

These nonnegative weights sum to one. Angle lowers the target without
weakening its enforcement. For `bohmth<=0`, $A=1$ at nonzero incidence;
exact tangency has $A=0$. The Mach transition has fixed width 0.20 and
starts at $1/\beta$, so energy suppression also extends the imposed region.
Inward flow cannot select the free branch. At zero incidence the Mach sign
convention is positive for either floating-point sign of zero. Geometry is
fixed during Newton.

## Boundary equation and cases

For diagonal stabilization, the momentum equation is

$$
B_2=\tau_{22}\left[p_ss_nc_s\widehat U_1+p_FU_2-\widehat U_2\right]
-g_2(\widehat U,Q)=0.
$$

$g_2$ includes direct perpendicular momentum diffusion and the complete
$W_2$ density-gradient contribution. Full momentum diffusion is retained in
all branches. Consequently, diffusion can offset the target velocity;
the imposed branch is a mixed condition, not exact trace Mach equality.

| State | $p_0$ | $p_s$ | $p_F$ | Stabilization target for momentum |
|---|---|---|---|---|
| $M\le1/\beta$, $\beta=A=1$ | 0 | 1 | 0 | $s_nc_s\widehat U_1$ |
| $M\le1/\beta$, arbitrary angle/energy | $1-\beta A$ | $\beta A$ | 0 | $s_n\beta A c_s\widehat U_1$ |
| $M\le1/\beta$, exact tangency | 1 | 0 | 0 | Zero parallel momentum |
| $\beta\to0$ at finite Mach | $\to1$ | $\to0$ | 0 | Tends to zero parallel momentum |
| $1/\beta<M<1/\beta+0.20$ | $w(1-\beta A)$ | $w\beta A$ | $1-w$ | $ws_n\beta A c_s\widehat U_1+(1-w)U_2$ |
| $M\ge1/\beta+0.20$ | 0 | 0 | 1 | Element momentum $U_2$; free branch |

Even in the transition, angle and energy reduce the imposed target; they do
not prescribe a larger Mach.
The same $-g_2$ appears in every row. This is a numerical wall regularization,
not a kinetic sheath model. Energy sheath conditions and recycling retain
their existing flux laws. Isothermal models keep their previous hard selector
and target.

## Linearization and ownership

`bohm_boundary` evaluates weights and their trace derivatives through separate
angle, suppression and Mach helpers. HDG assembly owns matrix insertion.
Let $c=p_ss_nc_s$ and $f=p_F$. With tau and transport held fixed, the added
momentum terms are

$$
\frac{\partial B_2}{\partial U_j}=\tau_{22}f\delta_{j2},\qquad
\frac{\partial B_2}{\partial Q}=-\frac{\partial g_2}{\partial Q},
$$

$$
\frac{\partial B_2}{\partial\widehat U_j}
=\tau_{22}\left[c\delta_{j1}+\widehat U_1\frac{\partial c}{\partial\widehat U_j}
+U_2\frac{\partial f}{\partial\widehat U_j}-\delta_{j2}\right]
-\frac{\partial g_2}{\partial\widehat U_j}.
$$

| Contribution | Destination |
|---|---|
| Free element momentum | `Alu` |
| Trace momentum, target density, sound speed and weight derivatives | `All` |
| Direct and $W_2$ gradient derivatives | `Alq` |
| $W_2$ trace derivatives | `All` |
| $\tau_{22}(\widehat U_1\,\partial c/\partial\widehat U+U_2\,\partial f/\partial\widehat U)\cdot\widehat U$ | Affine `fh` |

The absolute Newton RHS is $Jx-B$. Sound speed has degree zero, so
$(\partial c_s/\partial\widehat U)\cdot\widehat U=0$. Suppression has
$(\partial\beta/\partial\widehat U)\cdot\widehat U=\beta(1-\beta)$ above
the floor: its Jacobian and affine terms remain. The Mach threshold is also
differentiated.

The optional NeutralGamma target uses the same active weight $w$, including
its derivative and affine term; its velocity coefficient remains lagged.
Balance diagnostics use the assembled nonlinear target, without affine terms.
