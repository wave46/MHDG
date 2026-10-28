# Absolute and incremental Newton formulations

In `NUMER_LST`, `newton_formulation='absolute'` (default) or `'incremental'`
selects the linear unknown in the main `MHDG` driver. Incremental solves are
supported for 2D temperature builds. The manufactured-solution `Convergence`
driver requires `absolute`.

For the complete current state $X_k=(U_k,Q_k,\widehat U_k)$ and the assembled
tangent $J_k$, the equivalent systems are

$$
J_kX_{\mathrm{trial}}=J_kX_k-F(X_k),\qquad
J_k\Delta X=-F(X_k),\qquad X_{\mathrm{trial}}=X_k+\Delta X.
$$

The assembled tangent includes the selected coefficient derivatives and lagged
coefficients. Incremental form asks the linear solver to resolve the correction
directly. Residual reconstruction and local elimination involve cancellation;
roundoff and linear-solver tolerances can affect the numerical trajectory.

## Local elimination

Write $R_U,R_Q,R_{\widehat U}$ for the three current
[HDG residual blocks](hdg-residuals.md). The gradient residual already includes
the inverse gradient mass. The incremental local equations are

$$
\begin{aligned}
A_{uu}\Delta U+A_{uq}\Delta Q+A_{ul}\Delta\widehat U&=-R_U,\\
\Delta Q+M_Q^{-1}(A_{qu}\Delta U+A_{ql}\Delta\widehat U)&=-R_Q.
\end{aligned}
$$

`HDG_mapping` eliminates the element unknowns in favour of face unknowns.
Define $B=A_{uu}-A_{uq}M_Q^{-1}A_{qu}$ and
$C=A_{ul}-A_{uq}M_Q^{-1}A_{ql}$. Its incremental maps are

$$
\begin{aligned}
\Delta U&=-B^{-1}C\Delta\widehat U+U_0,
&U_0&=B^{-1}(-R_U+A_{uq}R_Q),\\
\Delta Q&=-M_Q^{-1}(A_{qu}\Delta U+A_{ql}\Delta\widehat U)-R_Q.
\end{aligned}
$$

The constant part of the second map is $Q_0=-M_Q^{-1}A_{qu}U_0-R_Q$, stored
as `L0`. A restart or mesh transfer can leave a nonzero gradient residual;
$R_Q$ carries this mismatch into the correction.

`HDG_assembly` inserts these maps into the face equations and sums neighbouring
contributions into the global trace system. With the map slopes stored as
`UU` and `LL`, its element contributions are

$$
K^e=A_{lu}\,\mathrm{UU}+A_{lq}\,\mathrm{LL}+A_{ll},\qquad
b^e=-R_{\widehat U}-A_{lu}U_0-A_{lq}Q_0.
$$

`element_conservation_residual`, `element_gradient_residual` and
`element_trace_residual` contain the current element
contributions $R_U^e,R_Q^e,R_{\widehat U}^e$, before neighbouring trace rows
are summed. Each array column holds one element's equation residual vector.
Diagnostic norms are calculated separately. Mapping takes
the first two to form its constant terms; assembly takes the third to form
the trace load. Calls without these optional arguments use the absolute
affine loads from the physical assembly arrays.

## Driver and boundary handling

Each [Newton iteration](../simulation/solver-flow.md) assembles at the current
state and evaluates any required residual. Incremental solves and optional
monitoring share this evaluation. After the selected stopping check, local
elimination and global assembly prepare the trace solve; local reconstruction
then recovers the element unknowns.

Unknown trace increments are solved globally. Prescribed traces have zero
increments; their boundary data are already eliminated into local loads.
An absolute solve reconstructs $X_{\mathrm{trial}}$; an incremental solve
reconstructs $\Delta X$. Damping gives, respectively,

$$
X_{k+1}=X_k+\lambda(X_{\mathrm{trial}}-X_k),\qquad
X_{k+1}=X_k+\lambda\Delta X.
$$

A failed solve or nonfinite trial restores the saved current state.

The periodic assembler retains separate unknowns on paired boundary faces.
Their current values can differ in an initial guess, after mesh transfer, or
during damped iteration. A copied row replaces one face column by its partner.
If its original condensed residual is $r_a$, replacing $\widehat U_a$ by
$\widehat U_b$ gives $r_b=r_a+K_{aa}^e(\widehat U_b-\widehat U_a)$.
Its incremental right-hand side is therefore
$-r_a-K_{aa}^e(\widehat U_b-\widehat U_a)$, using corresponding trace nodes.

The selected formulation is printed with numerical settings and saved as
`newton_formulation` in the HDF5 numerical-parameter group. It is independent of
`newton_convergence` and `residual_diagnostics_mode`.
