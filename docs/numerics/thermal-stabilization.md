# Thermal stabilization

For the 2D temperature model, automatic diagonal stabilization (`stab=5`)
combines an advective estimate, perpendicular diffusion and parallel thermal
conduction. For each energy equation,

$$
\tau_s=a+\frac{p}{h_e}D_{\perp,s}
+\frac{p}{h_e}|b\cdot n|\,\|b_{1:2}\|
\frac{\kappa_{\parallel,s}(\min(T_s,T_{\max,s}))}{\widehat U_1}.
$$

Here $p$ is polynomial degree, $h_e$ is the element size used by the solver,
and $n$ is the outward unit normal. Temperatures and transport coefficients
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
