# Diffusion in HDG

For a scalar field $w$, the physical diffusive flux has the form
$\boldsymbol\Gamma_d=-D\mathsf P\nabla w$, with the isotropic or perpendicular
projection $\mathsf P$. Neutral diffusion uses $D_n$ and, when enabled, the
multiplier $\phi$. Turbulent diffusion $D_k$ contributes to plasma and k
transport; it does not replace $D_n$.

In the 2D NGammaTiTe family, $U$ is the element state, $\widehat U$ the
independent face trace, and $Q$ the element gradient variable. Constitutive
coefficients use

$$
D_n=D_n(U),\quad D_k=D_k(U)\quad\text{in the volume},
\qquad
D_n=D_n(\widehat U),\quad D_k=D_k(\widehat U)\quad\text{on faces}.
$$

The face convention includes interior, periodic, MPI partition and physical
boundary faces. Gradients remain element quantities. The neutral limiter
uses $(U,Q)$ in the volume and $(\widehat U,Q)$ on faces, including its
optional pressure-transport contribution.

$D_k$, $\phi$ and stabilization $\tau$ remain lagged in the nonlinear
linearization. Stabilization uses unlimited diffusion and acts on the jump
$U-\widehat U$. Neither convergence nor weak numerical-flux continuity
requires this jump to vanish. Recycling is described in the
[wall model](../models/wall-recycling.md).

Source: [diffusion coefficients](../../src/Models/NGammaTiTe/physics.f90),
[element assembly](../../src/HDG/hdg_ComputeJacobian.f90),
[boundary assembly](../../src/HDG/hdg_BC.f90).
