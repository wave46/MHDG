# Parallel heat-flux limiter

For ions and electrons, the conductive heat flux is

$$
q_{\parallel,s}=\mathrm{FL}_s q_{\mathrm{SH},s},\qquad
q_{\mathrm{SH},s}=-\kappa_{\parallel,s}\nabla_\parallel T_s,\qquad
\mathrm{FL}_s=\frac{1}{1+|q_{\mathrm{SH},s}|/(c_{\mathrm{fl},s}q_{\mathrm{cap},s})}.
$$

The free-streaming flux scales as $q_{\mathrm{fs},s}\propto nT_s^{3/2}$,
with the species mass factor included by the physics helpers.
`c_fli` and `c_fle` specify the fractions of free streaming used for the ion
and electron caps. `flux_limiter` enables both. This limits conductive heat
flux; the prescribed sheath energy flux is assembled separately.
`T_fluxlim_maxi` and `T_fluxlim_maxe` separately cap the temperature used
in the parallel conductivity.

The existing regularization uses $T_s\ge10^{-5}$ eV in the free-streaming
cap and $q_{\mathrm{cap},s}=\max(|q_{\mathrm{fs},s}|,10^{-10})$ in solver
flux units. Below the temperature floor, $q_{\mathrm{fs},s}$ still depends
linearly on density. Below the outer cap floor, $q_{\mathrm{cap},s}$ is
constant and its state derivative vanishes.

Write the Spitzer–Harm Taylor expression as $q^Q+q^U+q^0$, with
$q^Q=(\partial_Q q_{\mathrm{SH}})Q$, $q^U=(\partial_U q_{\mathrm{SH}})U$,
and $q^0=-(\partial_U q_{\mathrm{SH}})U_0$. All coefficients are evaluated
at $(U_0,Q_0)$. With the outer cap floor inactive, the limited expression is

$$
q_{\parallel,\mathrm{lin}}=\mathrm{FL}^2(q^Q+q^U+q^0)
+\frac{\mathrm{FL}^2 q_{\mathrm{SH}}r}{q_{\mathrm{fs}}}
(\partial_U q_{\mathrm{fs}})U,\qquad r=\frac{|q_{\mathrm{SH}}|}{c_{\mathrm{fl}}q_{\mathrm{fs}}}.
$$

For $h=\mathrm{FL} q_{\mathrm{SH}}$, the scalar derivatives are

$$
\frac{\partial h}{\partial q_{\mathrm{SH}}}=\mathrm{FL}^2,\qquad
\frac{\partial h}{\partial q_{\mathrm{fs}}}=
\begin{cases}
(\mathrm{FL} q_{\mathrm{SH}})(1-\mathrm{FL})/q_{\mathrm{fs}}, & |q_{\mathrm{fs}}|>10^{-10},\\
0, & |q_{\mathrm{fs}}|\le10^{-10}.
\end{cases}
$$

The affine linearization is $h_{\mathrm{lin}}=A_U U+A_Q Q+c_0$, where
$c_0=h_0-A_U U_0-A_Q Q_0$ and all coefficients are evaluated at the
linearization state $(U_0,Q_0)$. The assembled right-hand side stores
$-c_0$, the Jacobian contraction minus $h_0$. This recovers the actual
limited flux at that state, including when either floor is active.
When the outer cap floor is active, the helper adds
$\mathrm{FL}(\mathrm{FL}-1)q_{\mathrm{SH}}$ to this right-hand-side quantity; the
corresponding constant in the flux has the opposite sign.

In [physics.f90](../../src/Models/NGammaTiTe/physics.f90),
`compute_free_streaming_heat_flux_*` and their derivative helpers evaluate
the same floored coefficient. `compute_flux_limiter` supplies the value;
`compute_heat_flux_limiter_linearization` supplies that value, its cap
sensitivity and the affine correction. Volume, interior/exterior faces,
Bohm conditions and strong-Dirichlet contributions use these definitions.
State choices and deliberately lagged stabilization remain as specified
by their callers.
