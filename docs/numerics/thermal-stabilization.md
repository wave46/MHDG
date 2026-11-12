# Thermal stabilization

The NGammaTiTe family selects diagonal stabilization with `stab` in 2D and
TOR3D (other models use model-specific selectors):

| Value | Coefficients |
|---|---|
| 1 | Prescribed per-equation values from `tau(:)` |
| 2 | Automatic convection and diffusion coefficients |

Automatic stabilization combines the convective spectral radius,
perpendicular diffusion and parallel thermal conduction. For each energy equation,

$$
\tau_s=a+\frac{p}{h_e}D_{\perp,s}
+\frac{p}{h_e}|b\cdot n|\,\|b_{\mathrm{active}}\|
\frac{\kappa_{\parallel,s}(\widehat T_s,\widehat n_s)}{\widehat U_1}.
$$

Here $n$ is the outward unit normal. On 2D and TOR3D side faces, $p/h_e$
is the poloidal degree divided by the stored element size (the shortest vertex
edge for triangles). On TOR3D cross-section faces it is
$p_{\mathrm{tor}}/(R\Delta\varphi)$, where $\Delta\varphi$ is the toroidal
cell angle and $R$ is in solver coordinates. The magnetic vector
$b_{\mathrm{active}}$ contains the poloidal components in 2D and all three
components in TOR3D.
The wrapper passes $p/h_e$ as `length_factor`; it stays fixed during each
Newton assembly. Temperatures and transport coefficients
are in solver units. For the coupled four-equation plasma convection block,
the common scalar is its spectral radius:

$$
a=|b\cdot n|\max\left(\frac53|u_\parallel|,
                       |u_\parallel|+c_{\rm flux}\right),\qquad
u_\parallel=\frac{\widehat U_2}{\widehat U_1},\qquad
c_{\rm flux}=\sqrt{\frac53 M_{\rm ref}(\widehat T_i+\widehat T_e)}.
$$

The flux Jacobian has eigenvalues $(b\cdot n)$ times
$u_\parallel-c_{\rm flux}$, $u_\parallel$, $u_\parallel+c_{\rm flux}$
and $5u_\parallel/3$. Here $c_{\rm flux}$ is the adiabatic acoustic speed
of this convection block, with ratio of specific heats $\gamma=5/3$. The same scalar is added
to every diagonal row, including the neutral row as extra numerical damping.

Here $\kappa_{\parallel,s}$ is the parallel conductivity before local
free-streaming saturation, including the conductivity ceiling and selected
closed-surface correction. Classical stabilization uses the $T^{5/2}$ law;
physical conductivity uses `epn`. With `neoclassical_heat_conduction` enabled,
closed-core thermal stabilization uses the corrected, capped conductivity and
its density and temperature derivatives. Geometry stays fixed during Newton
assembly. See the [physical correction](../models/heat-flux-limiter.md#closed-field-lines-neoclassical-correction).

The free-streaming multiplier bounds the physical parallel conductive flux.
The conductivity penalty controls the element–trace temperature jump. Neutral
stabilization similarly uses the unlimited neutral diffusion coefficient.
The penalty can carry numerical heat flux larger than the physical limited
conduction; assess wall heat-flux accuracy using mesh resolution and the
size of the element–trace temperature jump.

The model-local `physics_stabilization` submodule owns coefficient evaluation;
HDG assembly owns jumps, matrix placement and boundary conditions. Values are
recomputed from the current trace state and geometry at each assembly. The
[free-streaming multiplier](../models/heat-flux-limiter.md#open-field-lines-free-streaming-limiter)
is held fixed in volume and face linearizations.

For the five-equation 2D plasma–neutral model, the tangent differentiates all
five diagonal coefficients with respect to $\widehat U$: the common advective
estimate, both thermal terms and the unlimited neutral diffusion coefficient.
Supplied plasma diffusion profiles remain fixed during this assembly.
These stabilization coefficients are independent of the element gradient $Q$.
TOR3D and other model families use frozen-coefficient tangents.

For an ordinary face, stabilization contributes
$\tau_i(\widehat U)(U_i-\widehat U_i)$. Its trace derivative is

$$
\frac{\partial}{\partial\widehat U_j}
  [\tau_i(\widehat U)(U_i-\widehat U_i)]
=-\tau_i\delta_{ij}
 +(U_i-\widehat U_i)\frac{\partial\tau_i}{\partial\widehat U_j}.
$$

The assembly adds $+\tau_i$ in the element-state column
and $-\tau_i$ in the trace-state column. The HDG coefficient helper adds the second
term above, with shape-function products and the face integration weight.
It also adds the matching absolute-Newton affine RHS:

$$
(U_i-\widehat U_i)\sum_j
\frac{\partial\tau_i}{\partial\widehat U_j}\widehat U_j.
$$

For [Bohm momentum](bohm-boundary.md), replace $U_2-\widehat U_2$ by the actual
blended target jump
$p_s\,\operatorname{sign}(b\cdot n)c_s\widehat U_1+p_F U_2-\widehat U_2$.
The boundary assembly differentiates this target; the coefficient helper
adds its product with $\partial\tau_2/\partial\widehat U_j$.

| Equation row | Added matrix | Added affine RHS |
|---|---|---|
| Element numerical flux | `Aul` | `S` |
| Interior trace continuity | `All` | `fh` |
| Bohm trace equation | `All` | `fh` |
| Element with prescribed trace | None after trace elimination | None |

The physics submodule computes tau values and derivatives from state and
geometry. The HDG helper adds the supplied coefficient products to the
specified matrix and RHS entries; its caller chooses the ordinary or Bohm
jump and the destination equation. Interior faces call it for both element
balance and trace-flux continuity, which are different rows of the HDG system.

The advective term has degree zero in $\widehat U$. Classical thermal terms
have degree minus one, so their derivatives contracted with $\widehat U$
equal the negative thermal terms, including above the conductivity ceiling.
The neoclassical conductivity contributes its additional density response
through the shared conductivity derivative. Neutral diffusion
is differentiated through its density-dependent rates. These identities apply
away from the density floor.

The characteristic-speed interpretation requires admissible temperatures.
For trial states the coefficient continues its square root with the absolute
value of the radicand. ABS/MAX use active-branch derivatives. At a spectral MAX
tie the velocity branch is selected. Temperature derivatives vanish below
or at their floor and above their cap; the uncapped derivative is used at the
cap. At the singular zero of the square-root spectral estimate, its derivative
is held at zero. The tangent remains finite at this singular state.
