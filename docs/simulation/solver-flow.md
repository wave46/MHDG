# Solver flow

This is the main `MHDG` driver flow. A Newton iteration changes the solution
at one fixed time endpoint; it does not advance physical time. The complete
iterate consists of element states $U$, gradients $Q$, and face traces
$\widehat U$.

```mermaid
flowchart TD
    A[Read inputs; load mesh and initial or restart state]
    B[Prepare fields, sources, parallel data and mesh matrices]
    C[Prepare next time-step trial and time-dependent inputs]
    D[Assemble element equations and boundary conditions at current iterate]
    E[Optional residual evaluation and reporting]
    F{Residual stopping selected and satisfied; admissible?}
    G[Eliminate element unknowns; assemble and solve global trace system]
    H[Reconstruct element states and gradients; damp complete update]
    I{Update stopping selected and satisfied; admissible?}
    J[Continue Newton or invoke configured recovery]
    K[Accept nonlinear solution]
    L[Record time-change measure; output; update solution history]
    M{Run complete?}
    N[Save final solution and release resources]
    A --> B --> C --> D --> E --> F
    F -- Yes --> K
    F -- No --> G --> H --> I
    I -- Yes --> K
    I -- No --> J
    J -- Continue or recover --> D
    J -- Retry physical step after rollback --> C
    K --> L --> M
    M -- No --> C
    M -- Yes --> N
```

## Assembly, solve and acceptance

Assembly evaluates coefficients and sources at the current iterate and forms
the element matrices and right-hand sides, including time-history terms when
active. Boundary assembly supplies the face equations. Element elimination
(static condensation) leaves a global system for unknown traces; prescribed
traces enter as boundary data. After its solution, local reconstruction recovers
$U$ and $Q$. The damping factor applies to all three fields together.

The default `newton_convergence='update'` checks the existing volume-update
measure after reconstruction and damping. Residual reporting is optional and
does not change that stopping rule. The experimental `residual` and `both`
choices can accept the assembled current state before another solve; `both`
also requires the preceding update to pass. The next assembly checks the state
produced by the previous update, at the same time endpoint. See
[residual evaluation and convergence](../numerics/hdg-residuals.md) for the
norms, controls and precise update measure.

## Recovery and outer continuation

Configured Newton recovery can restore a checkpoint or adapt the mesh, project
the solution/history and rebuild mesh-dependent solver data before continuing.
A failed physical step can instead restore the accepted state and retry a
shorter interval when enabled; exhausted recovery stops the run. Optional
prediction supplies an initial guess and falls back to the accepted state if
the predicted solve fails.

After nonlinear acceptance, the driver evaluates the separate time-change
measure. Transient runs advance history and time; legacy time relaxation can
stop at steady state. Pseudo-time or steady diffusion continuation can reduce
diffusion and solve again. Output and periodic mesh adaptation follow their
configured schedules. See [time stepping](time-stepping.md),
[configuration](configuration.md) and [restart/output](restart-and-output.md).

## Code entry points

| Responsibility | Location |
| --- | --- |
| Startup, outer loop, output and diffusion continuation | [`MHDG.f90`](../../src/MHDG.f90) |
| Time-step trials, Newton loop, damping and recovery | [`Main_utils.f90`](../../src/Utils/Main_utils.f90) |
| Element and boundary assembly | [`hdg_ComputeJacobian.f90`](../../src/HDG/hdg_ComputeJacobian.f90), [`hdg_BC.f90`](../../src/HDG/hdg_BC.f90) |
| Element elimination and global trace assembly | [`hdg_Mapping.f90`](../../src/HDG/hdg_Mapping.f90), [`hdg_Assembly.f90`](../../src/HDG/hdg_Assembly.f90) |

The separate manufactured-solution driver, `Convergence`, retains its own
update-based loop; it can report residuals without using them for acceptance.
