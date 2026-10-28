# Attributing a Newton response to equation terms

A temperature change at one point can depend on heating, transport and boundary
conditions throughout the coupled domain. Discrete-adjoint attribution measures
how each equation term contributes to that point's **linearized Newton change**.
The method uses the assembled residual and tangent, including their chosen
lagged coefficients. This guide describes the mathematics; it does not define
a production command or runtime option.

## The Newton response

Let $X=(U,Q,\widehat U)$ contain the element states, gradients and traces.
At the current state $X_k$, let $R=R(X_k)$ be the signed equation residual
and $J$ the assembled tangent. The undamped Newton correction satisfies

$$
J\Delta X=-R.
$$

This correction also exists when solving for the absolute trial state. With
the affine right-hand side $b=JX_k-R$, the two forms are

$$
JX_{\mathrm{trial}}=b,\qquad
\Delta X=X_{\mathrm{trial}}-X_k.
$$

After damping by $\lambda$, the updated state is $X_k+\lambda\Delta X$.
The attribution therefore applies to both
[Newton formulations](newton-formulation.md).

Choose one scalar observable $\mathcal O(X)$, such as electron temperature
at a wall point. Its first-order change is

$$
\delta\mathcal O_{\mathrm{lin}}=w^T\Delta X,
\qquad w=\lambda\nabla_X\mathcal O(X_k).
$$

The vector $w$ describes how each unknown affects the observable. Damping
enters once, through $w$. For a large correction, this first-order change can
differ from the finite difference $\mathcal O(X_k+\lambda\Delta X)-\mathcal O(X_k)$.

## Adjoint weights and contributions

Solve the transposed system

$$
J^Tz=w.
$$

The adjoint vector $z$ weights the equation residuals according to their effect
on the chosen observable. Although $w$ may be supported at one point, coupling
can give $z$ nonzero entries throughout the domain. Taking the scalar product
with the Newton equation gives

$$
w^T\Delta X=z^TJ\Delta X=-z^TR.
$$

Partition the complete residual into groups $R_c$, where $c$ labels a mechanism
such as conduction, exchange or sheath loss. Then

$$
R=\sum_cR_c,\qquad
a_c=-z^TR_c,\qquad
\boxed{\sum_ca_c=\delta\mathcal O_{\mathrm{lin}}.}
$$

One adjoint solve supplies weights for every group. Each $a_c$ has the
observable's units, for example eV for temperature. A negative heating
contribution means its residual drives a negative temperature correction
through this coupled tangent; the heating power itself can remain positive.
All groups use the same $J$. Changing a physical term and reconverging the
nonlinear model is a separate calculation.

These contractions use the signed equation vectors described in the
[HDG residual guide](hdg-residuals.md), before diagnostic norm reductions.

## Defining the observable

For trace electron temperature, let $N$ be the row of face interpolation
weights at the observation point. The interpolated trace density and
electron-energy variable are

$$
\widehat n=N\widehat U_1,\qquad
\widehat E_e=N\widehat U_4,\qquad
\widehat T_e=C_T\frac{\widehat E_e}{\widehat n},
$$

where $C_T$ converts the conserved-variable ratio to temperature in eV.
For $\mathcal O=\widehat T_e$, the nonzero entries of $w$ are

$$
\begin{aligned}
w_{\widehat U_1}&=-\lambda C_T
 \frac{\widehat E_e}{\widehat n^2}N^T,\\
w_{\widehat U_4}&=\lambda C_T\frac{1}{\widehat n}N^T.
\end{aligned}
$$

Both density and energy changes affect temperature. The same differentiation
defines observables inside an element, involving gradients, or averaged over
a region. Element and trace values are distinct; an observation on an element
interface must specify the element side or the trace. A prescribed trace has
no free correction when its prescribed data remain fixed.

## The condensed HDG adjoint

HDG eliminates element unknowns locally and assembles a global trace matrix
$K$. The adjoint can reuse that condensation. Write the residual as
$(R_U,R_Q,R_{\widehat U})$, using the inverse gradient mass $M_Q^{-1}$
(`iAqq`) in the gradient equation. The mixed tangent acts as

$$
J\Delta X=
\begin{pmatrix}
A_{uu}\Delta U+A_{uq}\Delta Q+A_{ul}\Delta\widehat U\\
\Delta Q+M_Q^{-1}(A_{qu}\Delta U+A_{ql}\Delta\widehat U)\\
A_{lu}\Delta U+A_{lq}\Delta Q+A_{ll}\Delta\widehat U
\end{pmatrix}.
$$

The subscripts indicate the equation row and unknown column; $l$ denotes
the trace. After eliminating $\Delta Q$, define the element blocks

$$
\begin{aligned}
B&=A_{uu}-A_{uq}M_Q^{-1}A_{qu},
& C&=A_{ul}-A_{uq}M_Q^{-1}A_{ql},\\
G&=A_{lu}-A_{lq}M_Q^{-1}A_{qu},
& D&=A_{ll}-A_{lq}M_Q^{-1}A_{ql}.
\end{aligned}
$$

| Block | Role after eliminating the gradient |
| --- | --- |
| $B$ | Element-state operator (`Mu` in mapping) |
| $C$ | Trace contribution to the element equations (`Ml` in mapping) |
| $G$ | Element-state contribution to the trace equations |
| $D$ | Trace contribution to the trace equations |

The remaining equations are

$$
\begin{aligned}
B\Delta U+C\Delta\widehat U&=-R_U+A_{uq}R_Q,\\
G\Delta U+D\Delta\widehat U&=-R_{\widehat U}+A_{lq}R_Q.
\end{aligned}
$$

Eliminating $\Delta U$ gives the element trace matrix

$$
K^e=D-GB^{-1}C.
$$

Here $B^{-1}$ denotes a local solve. This is the same matrix stored as `Kel`
by production assembly, equivalently
$A_{lu}\,\mathrm{UU}+A_{lq}\,\mathrm{LL}+A_{ll}$, where `UU` and `LL`
map trace corrections to element-state and gradient corrections. Assembly
combines the $K^e$ matrices into $K$ using the face orientations and global
trace mapping.

For an observable depending only on free traces, solve

$$
K^T\widehat z=w_{\widehat U}.
$$

Gather the global trace adjoint to each element as $\widehat z_e$. The local
weights that directly contract with its three residual blocks are

$$
\begin{aligned}
v_U^e&=B^{-T}G^T\widehat z_e,\\
v_Q^e&=A_{lq}^T\widehat z_e-A_{uq}^Tv_U^e,\\
v_{\widehat U}^e&=-\widehat z_e.
\end{aligned}
$$

These are the corresponding entries of $-z$; $B^{-T}$ is a local transposed
solve. Thus

$$
a_c=\sum_{e\;\mathrm{owned}}
\left[
(v_U^e)^TR_{U,c}^e+
(v_Q^e)^TR_{Q,c}^e+
(v_{\widehat U}^e)^TR_{\widehat U,c}^e
\right].
$$

Each physical element contributes once. Shared trace rows receive all their
element contributions, with the production orientation and ownership rules.
Volume or gradient observables supply additional local adjoint loads, which
also modify the condensed adjoint right-hand side; their global matrix is
still $K^T$.

## Physical grouping and verification

A transport group includes its volume weak term, element-face flux and
applicable trace or boundary term. Numerical stabilization contributes to both
element and trace equations. An exchange group includes its opposite ion and
electron energy entries. Group definitions follow the equations, so a single
wall-row contribution cannot stand in for an entire transport mechanism.

Check both residual completeness and response closure:

$$
\varepsilon_R=R-\sum_cR_c,\qquad
\varepsilon_{\mathcal O}=w^T\Delta X-\sum_ca_c.
$$

Report the projected remainder $-z^T\varepsilon_R$ alongside the contributions.
Independent solves $J\Delta X_c=-R_c$ for selected groups verify that
$w^T\Delta X_c=a_c$. Comparing the reconstructed full correction with the
actual Newton update checks the state ordering, damping and interpolation.
The result explains the selected linearized step, with its chosen grouping
and tangent, rather than a complete nonlinear trajectory.

For the general discrete-adjoint derivation, see
[dolfin-adjoint: differentiating functionals](https://www.dolfin-adjoint.org/en/latest/documentation/maths/3-gradients.html).
