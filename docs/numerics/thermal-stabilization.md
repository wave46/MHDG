# Thermal stabilization

For the 2D temperature model, automatic diagonal stabilization (`stab=5`)
combines an advective estimate, perpendicular diffusion and parallel thermal
conduction. For each energy equation,

$$
\tau_s=a+\frac{p}{h_e}D_{\perp,s}
+\frac{p}{h_e}|b\cdot n|\,\|b_{1:2}\|
\frac{\kappa_{\parallel,s}(\min(T_s,T_{\max,s}))}{\widehat U_1}.
$$

Here $p$ is polynomial degree, $h_e$ is the stored element-length estimate
(the shortest vertex edge for triangles), and $n$ is the outward unit normal.
The wrapper passes $p/h_e$ as `length_factor`; it stays fixed during each
Newton assembly. Temperatures and transport coefficients
are in solver units. The advective estimate $a$ is the existing maximum of two characteristic-speed
estimates. The prescribed diagonal option (`stab=1`) uses the input values.
Here $\kappa_{\parallel,s}$ is the unlimited parallel conductivity used by
the stabilization, with the existing $T^{5/2}$ law. Physical conductivity
uses `epn`.

The heat-flux limiter is **excluded** from stabilization. It limits the
physical parallel conductive flux separately. Thus saturation of conduction
does not simultaneously weaken control of the element–trace temperature
jump. Perpendicular diffusion and the prescribed sheath heat loss are
unchanged. Neutral stabilization likewise uses its unlimited diffusion
coefficient, supplied before applying the neutral physical limiter.

An unlimited penalty can carry numerical heat flux larger than the physical
limited conduction. It is an HDG stabilization contribution, not an extra
physical conductivity or a guarantee that the total numerical heat flux is
limited. Its effect on wall heat loss still depends on resolution and the
element–trace jump; convergence alone does not establish wall-flux accuracy.

The model-local `physics_stabilization` submodule owns coefficient evaluation;
HDG assembly owns jumps, matrix placement and boundary conditions. Values are
recomputed from the current trace state without a cache or restart state.
The shared value evaluator has no gradient or limiter dependence. The
[physical heat limiter](../models/heat-flux-limiter.md) uses the same frozen
limiter policy in volume and face tangents.

For the five-equation 2D plasma–neutral model, the tangent differentiates all
five diagonal coefficients with respect to $\widehat U$: the common advective
estimate, both thermal terms and the unlimited neutral diffusion coefficient.
Supplied plasma diffusion profiles remain fixed during this assembly.
There is no coefficient derivative with respect to the element gradient $Q$.
Other model families retain their existing frozen-coefficient tangents.

For an ordinary face, stabilization contributes
$\tau_i(\widehat U)(U_i-\widehat U_i)$. Its trace derivative is

$$
\frac{\partial}{\partial\widehat U_j}
  [\tau_i(\widehat U)(U_i-\widehat U_i)]
=-\tau_i\delta_{ij}
 +(U_i-\widehat U_i)\frac{\partial\tau_i}{\partial\widehat U_j}.
$$

The existing assembly already adds $+\tau_i$ in the element-state column
and $-\tau_i$ in the trace-state column. The new helper adds only the second
term above, with shape-function products and the face integration weight.
It also adds the matching absolute-Newton affine RHS:

$$
(U_i-\widehat U_i)\sum_j
\frac{\partial\tau_i}{\partial\widehat U_j}\widehat U_j.
$$

For [Bohm momentum](bohm-boundary.md), replace $U_2-\widehat U_2$ by the actual
blended target jump
$p_s\,\operatorname{sign}(b\cdot n)c_s\widehat U_1+p_F U_2-\widehat U_2$.
The existing boundary assembly differentiates this target; the new helper
adds its product with $\partial\tau_2/\partial\widehat U_j$.

| Equation row | Added matrix | Added affine RHS |
|---|---|---|
| Element numerical flux | `Aul` | `S` |
| Interior trace continuity | `All` | `fh` |
| Bohm trace equation | `All` | `fh` |
| Element with prescribed trace | None after trace elimination | None |

The physics submodule computes tau values and derivatives from state and
geometry. The HDG helper only adds the supplied coefficient products to the
specified matrix and RHS entries; its caller chooses the ordinary or Bohm
jump and the destination equation. Interior faces call it for both element
balance and trace-flux continuity, which are different rows of the HDG system.

The advective term has degree zero in $\widehat U$. Each thermal term has
degree minus one, so its derivative contracted with $\widehat U$ equals the
negative thermal term, including above its temperature cap. Neutral diffusion
has its own density-dependent rates; no common homogeneity cancellation is
assumed. These identities apply away from the density floor.

ABS/MAX use active-branch derivatives without extra smoothing. At a spectral
MAX tie the velocity branch is selected. Temperature derivatives vanish below
or at their floor and above their cap; the uncapped derivative is used at the
cap. At the singular zero of the square-root spectral estimate, its derivative
is held at zero. This convention avoids an infinite tangent; it does not make
the coefficient differentiable there.
