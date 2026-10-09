PROGRAM test_divergence_refinement
  USE Main_utils, ONLY: reset_newton_refinements, prepare_divergence_refinement, newton_refinement_required, &
       update_uiter_qiter_best, &
       physical_solution_is_admissible, update_best_newton_checkpoint, seed_divergence_checkpoint, &
       divergence_checkpoint_available, errNR_adapt, ir_adapt, &
       uiter_best, qiter_best, divergence_refinements, oscillation_refinements, prepare_oscillation_refinement, &
       project_u0_newmesh, Mesh_prec, ir
  USE globals
  USE MPI_OMP, ONLY: MPIvar
  USE reference_element, ONLY: create_reference_element, free_reference_element_pol
  USE adaptivity_indicator_module, ONLY: compute_error_oscillations
  USE mpi, ONLY: MPI_INIT, MPI_FINALIZE, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_COMM_WORLD
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  INTEGER :: ierr_test, j, h, e, node, eq, offset, expected_ranks, limit_case
  CHARACTER(LEN=32) :: rank_argument, mode_argument
  REAL*8 :: best_u(2), best_q(4), history(2,3), previous_coordinates(3,2), new_coordinates(4,2)
  REAL*8 :: expected(12,3), local_coordinate(2)
  LOGICAL :: recovered

  CALL MPI_INIT(ierr_test)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPIvar%glob_id, ierr_test)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPIvar%glob_size, ierr_test)
  CALL get_command_argument(1,rank_argument)
  READ(rank_argument,*) expected_ranks
  IF (MPIvar%glob_size /= expected_ranks) ERROR STOP 'test launcher does not match the linked MPI library'
  CALL get_command_argument(2,mode_argument)
  IF (mode_argument == 'input') THEN
    CALL read_input()
    WRITE(*,'(A,1X,I0)') 'DIV_REFINEMENT_LIMIT',adapt%max_divergence_refinements
    WRITE(*,'(A,1X,I0)') 'OSC_REFINEMENT_LIMIT',adapt%max_oscillation_refinements
    CALL MPI_FINALIZE(ierr_test)
    STOP
  ENDIF
  IF (adapt%max_divergence_refinements /= 4) ERROR STOP 'default divergence refinement limit'
  IF (adapt%max_oscillation_refinements /= 10) ERROR STOP 'default oscillation refinement limit'
  utils%printint = 0
  phys%neq = 2
  phys%idx_rhon_eq = -1
  Mesh%Nelems = 1
  Mesh%Nnodesperelem = 1
#ifdef PARALL
  ALLOCATE(Mesh%ghostElems(1))
  Mesh%ghostElems = 0
#endif
  ALLOCATE(sol%u(2),sol%q(4),sol%u_conv(2),sol%q_conv(4),sol%u0(2,3))
  sol%u = [1.d0,2.d0]
  sol%q = [3.d0,4.d0,5.d0,6.d0]
  history = RESHAPE([7.d0,8.d0,9.d0,10.d0,11.d0,12.d0],[2,3])
  sol%u0 = history
  time%it = 16
  time%ik = 16
  sol%Nt = 16
  time%t = 0.32d0
  time%dt = 0.02d0
  time%t_ME = 0.30d0
  time%dt_ME = 0.02d0
  switch%ME = .TRUE.
  phys%puff = 3.d21
  phys%feedback_integral_error = 123.d0
  phys%feedback_previous_error = 321.d0
  adapt%osc_check = -3.d0
  numer%tNR = 1.d-4
  numer%div = 1.d3
  numer%nrp = 12

  adapt%adaptivity = .FALSE.
  adapt%div_adapt = .TRUE.
  IF (newton_refinement_required(0.5d0*numer%tNR,.FALSE.,1)) ERROR STOP 'physical recovery enabled with adaptation off'
  IF (newton_refinement_required(0.1d0,.FALSE.,numer%nrp)) ERROR STOP 'NR-limit recovery enabled with adaptation off'
  IF (.NOT. newton_refinement_required(2.d0*numer%div,.TRUE.,1)) ERROR STOP 'inactive divergence stop changed'
  CALL reset_newton_refinements()
  IF (prepare_divergence_refinement()) ERROR STOP 'recovery enabled with adaptation off'
  adapt%adaptivity = .TRUE.
  adapt%div_adapt = .FALSE.
  IF (newton_refinement_required(0.5d0*numer%tNR,.FALSE.,1)) ERROR STOP 'physical recovery enabled with div_adapt off'
  IF (newton_refinement_required(0.1d0,.FALSE.,numer%nrp)) ERROR STOP 'NR-limit recovery enabled with div_adapt off'
  IF (prepare_divergence_refinement()) ERROR STOP 'recovery enabled with div_adapt off'

  adapt%div_adapt = .TRUE.
  IF (newton_refinement_required(0.5d0*numer%tNR,.TRUE.,1)) ERROR STOP 'admissible convergence requested refinement'
  IF (newton_refinement_required(0.1d0,.FALSE.,numer%nrp-1)) ERROR STOP 'transient negative iterate requested refinement'
  IF (newton_refinement_required(numer%tNR,.FALSE.,1)) ERROR STOP 'physical convergence boundary changed'
  IF (.NOT. newton_refinement_required(0.1d0,.FALSE.,numer%nrp)) ERROR STOP 'invalid NR-limit state did not request refinement'
  IF (newton_refinement_required(0.1d0,.TRUE.,numer%nrp)) ERROR STOP 'admissible NR-limit behavior changed'
  IF (.NOT. newton_refinement_required(2.d0*numer%div,.FALSE.,1)) ERROR STOP 'active divergence trigger changed'
  CALL reset_newton_refinements()
  best_u = [10.d0,20.d0]+MPIvar%glob_id
  best_q = [30.d0,40.d0,50.d0,60.d0]+MPIvar%glob_id
  CALL update_uiter_qiter_best(uiter_best,qiter_best,best_u,best_q)
  ! This selected best iterate need not satisfy osc_check. Each rank has
  ! distinct values, so recovery must use its own local checkpoint.
  DO limit_case = 1,3
    SELECT CASE(limit_case)
    CASE(1)
      adapt%max_divergence_refinements = 0
    CASE(2)
      adapt%max_divergence_refinements = 1
    CASE(3)
      adapt%max_divergence_refinements = 4
    END SELECT
    CALL reset_newton_refinements()
    CALL update_uiter_qiter_best(uiter_best,qiter_best,best_u,best_q)
    DO j = 1,adapt%max_divergence_refinements
      sol%u = 9.d30
      sol%q = 8.d30
      recovered = prepare_divergence_refinement()
      IF (.NOT. recovered) ERROR STOP 'valid divergence recovery rejected'
      IF (ANY(sol%u /= best_u) .OR. ANY(sol%q /= best_q)) ERROR STOP 'wrong recovery state'
      IF (ANY(sol%u_conv /= best_u) .OR. ANY(sol%q_conv /= best_q)) ERROR STOP 'wrong projection checkpoint'
      IF (divergence_refinements /= j) ERROR STOP 'retry accounting'
      CALL assert_physical_state_unchanged()
    ENDDO
    sol%u = 9.d30
    IF (prepare_divergence_refinement()) ERROR STOP 'unbounded divergence refinement'
    IF (ANY(sol%u /= 9.d30)) ERROR STOP 'exhausted recovery changed state'
  ENDDO
  adapt%max_divergence_refinements = 2

  ! A fresh timestep resets the budget AND the fallback solution.
  sol%u = [71.d0,72.d0]
  sol%q = [73.d0,74.d0,75.d0,76.d0]
  CALL reset_newton_refinements()
  sol%u = 9.d30
  IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'new timestep budget not reset'
  IF (ANY(sol%u /= [71.d0,72.d0])) ERROR STOP 'old timestep fallback retained'

  CALL check_refinement_budgets()
  CALL assert_physical_state_unchanged()
  DEALLOCATE(sol%u,sol%q,sol%u_conv,sol%q_conv,sol%u0,uiter_best,qiter_best)
  NULLIFY(uiter_best,qiter_best)
#ifdef PARALL
  DEALLOCATE(Mesh%ghostElems)
#endif
#if defined(NGAMMA) && defined(TEMPERATURE) && defined(NEUTRAL)
  CALL check_physical_checkpoint_selection()
  CALL check_oscillation_checkpoint_selection()
#endif

  ! Exercise production history projection on a split affine triangle. This
  ! performs interpolation and MPI gathering, without Gmsh or a plasma solve.
  phys%neq = 2
  phys%lscale = 0.001901d0
  time%tis = 3
  Mesh%ndim = 2
  Mesh%Nelems = 2
  Mesh%Nnodes = 4
  Mesh%Nnodesperelem = 3
  Mesh_prec%ndim = 2
  Mesh_prec%Nelems = 1
  Mesh_prec%Nnodes = 3
  Mesh_prec%Nnodesperelem = 3
  ALLOCATE(Mesh%X(4,2),Mesh%T(2,3),Mesh_prec%X(3,2),Mesh_prec%T(1,3),sol%u0(6,3))
  CALL create_reference_element(refElPol,2,1,verbose=0)
  Mesh_prec%X(:,1) = 2.d0+0.2d0*MPIvar%glob_id+0.05d0*(refElPol%coord2D(:,1)+1.d0)
  Mesh_prec%X(:,2) = 0.05d0*(refElPol%coord2D(:,2)+1.d0)
  Mesh_prec%T(1,:) = [1,2,3]
  Mesh%X(1:3,:) = Mesh_prec%X
  Mesh%X(4,:) = 0.5d0*(Mesh%X(1,:)+Mesh%X(2,:))
  Mesh%T(1,:) = [1,4,3]
  Mesh%T(2,:) = [4,2,3]
  Mesh_prec%X = Mesh_prec%X/phys%lscale
  Mesh%X = Mesh%X/phys%lscale
  previous_coordinates = Mesh_prec%X
  new_coordinates = Mesh%X
#ifdef PARALL
  Mesh_prec%Nel_glob = MPIvar%glob_size
  Mesh_prec%Nno_glob = 3*MPIvar%glob_size
  ALLOCATE(Mesh_prec%ghostElems(1),Mesh_prec%loc2glob_el(1),Mesh_prec%loc2glob_nodes(3))
  Mesh_prec%ghostElems = 0
  Mesh_prec%loc2glob_el = MPIvar%glob_id+1
  Mesh_prec%loc2glob_nodes = [1,2,3]+3*MPIvar%glob_id
#endif
  DO h=1,3
    DO node=1,3
      DO eq=1,2
        sol%u0((node-1)*2+eq,h) = history_field(Mesh_prec%X(node,:)*phys%lscale,eq,h)
      ENDDO
    ENDDO
    DO e=1,2
      DO node=1,3
        local_coordinate = Mesh%X(Mesh%T(e,node),:)*phys%lscale
        DO eq=1,2
          offset = ((e-1)*3+node-1)*2+eq
          expected(offset,h) = history_field(local_coordinate,eq,h)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  adapt%NR_adapt = .TRUE.
  adapt%freq_NR_adapt = 5
  ir = 5
  CALL project_u0_newmesh(preserve_history=.TRUE.)
  IF (.NOT. ALL(ieee_is_finite(sol%u0))) ERROR STOP 'nonfinite projected history'
  IF (MAXVAL(ABS(sol%u0-expected)) > 1.d-10) ERROR STOP 'time history not preserved on new mesh'
  IF (MAXVAL(ABS(Mesh%X-new_coordinates)) > 1.d-10 .OR. &
      MAXVAL(ABS(Mesh_prec%X-previous_coordinates)) > 1.d-10) ERROR STOP 'projection changed coordinate units'

  ! Existing callers omit preserve_history: keep their existing first-column
  ! projection behavior. The new all-column behavior is opt-in for recovery.
  DEALLOCATE(sol%u0)
  ALLOCATE(sol%u0(6,3))
  DO h=1,3
    DO node=1,3
      DO eq=1,2
        sol%u0((node-1)*2+eq,h) = history_field(Mesh_prec%X(node,:)*phys%lscale,eq,h)
      ENDDO
    ENDDO
  ENDDO
  adapt%NR_adapt = .FALSE.
  CALL project_u0_newmesh()
  IF (MAXVAL(ABS(sol%u0(:,1)-expected(:,1))) > 1.d-10 .OR. ANY(sol%u0(:,2:3) /= 0.d0)) &
    ERROR STOP 'default projection behavior changed'
  CALL free_reference_element_pol(refElPol)
  IF (MPIvar%glob_id == 0) WRITE(*,*) &
    'divergence refinement: checkpoint selection, retry bounds, MPI recovery and history projection PASS; ranks=', &
    MPIvar%glob_size
  CALL MPI_FINALIZE(ierr_test)

CONTAINS

  SUBROUTINE check_refinement_budgets()
    INTEGER :: attempt
    adapt%adaptivity = .TRUE.
    adapt%div_adapt = .TRUE.
    adapt%osc_adapt = .TRUE.
    adapt%max_divergence_refinements = 2
    adapt%max_oscillation_refinements = 2
    CALL reset_newton_refinements()
    ! Alternate both triggers, as can happen within one ME timestep. The
    ! post-projection checkpoint reseed and NR restart must retain both budgets.
    DO attempt = 1,2
       ir = 1
       IF (.NOT. prepare_oscillation_refinement()) ERROR STOP 'allowed oscillation refinement rejected'
       IF (oscillation_refinements /= attempt .OR. divergence_refinements /= attempt-1) &
            ERROR STOP 'oscillation refinement altered divergence budget'
       IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'alternating divergence refinement rejected'
       CALL seed_divergence_checkpoint()
       IF (oscillation_refinements /= attempt .OR. divergence_refinements /= attempt) &
            ERROR STOP 'remesh checkpoint reseed reset refinement budgets'
       CALL assert_physical_state_unchanged()
    ENDDO
    IF (prepare_oscillation_refinement()) ERROR STOP 'unbounded oscillation refinement'
    IF (prepare_divergence_refinement()) ERROR STOP 'alternating triggers bypassed divergence cap'
    IF (oscillation_refinements /= 2 .OR. divergence_refinements /= 2) &
         ERROR STOP 'exhausted requests changed counters'
    ! With divergence recovery off, oscillation recovery still has its own
    ! budget. Ten requests fit; the next is rejected without consuming a retry.
    adapt%div_adapt = .FALSE.
    adapt%max_oscillation_refinements = 10
    CALL reset_newton_refinements()
    IF (oscillation_refinements /= 0 .OR. divergence_refinements /= 0) &
         ERROR STOP 'next timestep did not reset both budgets'
    DO attempt = 1,10
       IF (.NOT. prepare_oscillation_refinement()) ERROR STOP 'new timestep oscillation budget not available'
    ENDDO
    IF (prepare_oscillation_refinement()) ERROR STOP 'default oscillation cap not enforced'
    IF (oscillation_refinements /= 10 .OR. divergence_refinements /= 0) &
         ERROR STOP 'oscillation-only budget accounting'
    adapt%max_oscillation_refinements = 0
    CALL reset_newton_refinements()
    IF (prepare_oscillation_refinement()) ERROR STOP 'zero oscillation cap allowed a retry'
    adapt%max_oscillation_refinements = 10
    adapt%osc_adapt = .FALSE.
    IF (prepare_oscillation_refinement()) ERROR STOP 'oscillation retry enabled with osc_adapt off'
    adapt%osc_adapt = .TRUE.
    adapt%adaptivity = .FALSE.
    IF (prepare_oscillation_refinement()) ERROR STOP 'oscillation retry enabled with adaptation off'
    IF (oscillation_refinements /= 0) ERROR STOP 'disabled oscillator consumed a retry'
    adapt%adaptivity = .TRUE.
    adapt%div_adapt = .TRUE.
    CALL assert_physical_state_unchanged()
  END SUBROUTINE check_refinement_budgets

  SUBROUTINE check_oscillation_checkpoint_selection()
    TYPE(Mesh_type) :: checkpoint_mesh
    REAL*8, ALLOCATABLE :: oscillations(:)
    REAL*8 :: minimum, maximum
    INTEGER :: n_osc, checkpoint_iteration, point
    LOGICAL :: allowed
    Mesh%Ndim = 2
    Mesh%Nelems = 1
    Mesh%Nnodes = 3
    Mesh%Nnodesperelem = 3
    phys%neq = 5
    phys%npv = 11
    phys%Mref = 12.d0
    phys%idx_rhon_eq = 5
    phys%idx_rhon_pv = 11
    utils%timing = .FALSE.
    adapt%quant_ind = 1
    adapt%thr_ind = 1.d-5
    adapt%osc_check = -3.d0
    ALLOCATE(adapt%n_quant_ind(1))
    adapt%n_quant_ind = [1]
#ifdef PARALL
    ALLOCATE(Mesh%ghostElems(1))
    Mesh%ghostElems = 0
#endif
    CALL create_reference_element(refElPol,2,1,verbose=0)
    ALLOCATE(sol%u(15),sol%q(30),sol%u_conv(15),sol%q_conv(30))
    DO point=1,3
       sol%u((point-1)*5+1:point*5) = [1.d0,0.2d0,2.d0,3.d0,0.01d0]
    ENDDO
    sol%q = 0.1d0
    sol%u_conv = 7.d0
    sol%q_conv = 8.d0
    checkpoint_iteration = 7
    checkpoint_mesh%Nelems = 77
    ! Smooth density satisfies osc_check, but negative energy on one rank
    ! must retain the previous checkpoint, its mesh and its iteration.
    IF (MPIvar%glob_id == MPIvar%glob_size-1) sol%u(4) = -1.d-8
    allowed = physical_solution_is_admissible(sol%u)
    CALL compute_error_oscillations(oscillations,minimum,maximum,n_osc,8,checkpoint_iteration, &
         checkpoint_mesh,checkpoint_allowed=allowed)
    IF (maximum > adapt%osc_check) ERROR STOP 'checkpoint gate fixture is not smooth'
    IF (checkpoint_iteration /= 7 .OR. checkpoint_mesh%Nelems /= 77) &
         ERROR STOP 'rejected checkpoint changed iteration or mesh'
    IF (ANY(sol%u_conv /= 7.d0) .OR. ANY(sol%q_conv /= 8.d0)) &
         ERROR STOP 'invalid smooth iterate overwrote oscillation checkpoint'
    ! Omitting the optional gate retains the existing checkpoint behavior.
    sol%u(4) = 3.d0
    CALL compute_error_oscillations(oscillations,minimum,maximum,n_osc,9,checkpoint_iteration,checkpoint_mesh)
    IF (checkpoint_iteration /= 9 .OR. checkpoint_mesh%Nelems /= Mesh%Nelems) &
         ERROR STOP 'default oscillation checkpoint behavior changed'
    IF (ANY(sol%u_conv /= sol%u) .OR. ANY(sol%q_conv /= sol%q)) &
         ERROR STOP 'default oscillation checkpoint copied wrong state'
    CALL free_mesh_loc(checkpoint_mesh)
    CALL free_reference_element_pol(refElPol)
    DEALLOCATE(oscillations,adapt%n_quant_ind,sol%u,sol%q,sol%u_conv,sol%q_conv)
#ifdef PARALL
    DEALLOCATE(Mesh%ghostElems)
#endif
  END SUBROUTINE check_oscillation_checkpoint_selection

  SUBROUTINE check_physical_checkpoint_selection()
    REAL*8 :: admissible_u(10), accepted(10), saved_error, rho, momentum, kinetic
    LOGICAL :: valid, updated
    INTEGER :: bad_rank, saved_iteration, retry_count
    phys%neq = 5
    phys%idx_rhon_eq = 5
    Mesh%Nelems = 2
    Mesh%Nnodesperelem = 1
#ifdef PARALL
    ALLOCATE(Mesh%ghostElems(2))
    Mesh%ghostElems = 0
#endif
    ALLOCATE(sol%u(10),sol%q(20),sol%u_conv(10),sol%q_conv(20))
    admissible_u = [1.d0,0.2d0,2.d0,3.d0,0.01d0, 1.d0,0.1d0,2.d0,3.d0,0.02d0]
    admissible_u = admissible_u*(1.d0+0.1d0*MPIvar%glob_id)
    sol%u = admissible_u
    sol%q = 0.1d0
    adapt%max_divergence_refinements = 3
    CALL reset_newton_refinements()
    IF (.NOT. divergence_checkpoint_available) ERROR STOP 'admissible initial fallback rejected'
    valid = physical_solution_is_admissible(sol%u)
    updated = update_best_newton_checkpoint(0.3d0,8,valid)
    IF (.NOT. updated) ERROR STOP 'admissible NR 8 not recorded'
    accepted = sol%u
    saved_error = errNR_adapt
    saved_iteration = ir_adapt

    ! Only one rank has negative density: every rank must retain NR 8.
    bad_rank = MPIvar%glob_size-1
    IF (MPIvar%glob_id == bad_rank) sol%u(1) = -1.d-4
    valid = physical_solution_is_admissible(sol%u)
    IF (valid) ERROR STOP 'rank-local negative density not globally rejected'
    updated = update_best_newton_checkpoint(0.095d0,10,valid)
    IF (updated .OR. ANY(uiter_best /= accepted)) ERROR STOP 'invalid lower-error iterate replaced recovery state'
    IF (errNR_adapt /= saved_error .OR. ir_adapt /= saved_iteration) ERROR STOP 'invalid iterate changed best score/iteration'
    IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'valid recovery checkpoint lost after bad iterate'
    IF (ANY(sol%u /= accepted)) ERROR STOP 'recovery did not restore admissible NR 8'

    ! A low residual must not complete a timestep with rank-local negative
    ! energy. It uses the same checkpoint and budget as the preceding divergence.
    sol%u = admissible_u
    IF (MPIvar%glob_id == bad_rank) sol%u(4) = -1.d-8
    valid = physical_solution_is_admissible(sol%u)
    IF (.NOT. newton_refinement_required(0.5d0*numer%tNR,valid,11)) &
         ERROR STOP 'inadmissible low-residual solution did not request refinement'
    IF (update_best_newton_checkpoint(0.5d0*numer%tNR,11,valid)) &
         ERROR STOP 'inadmissible low-residual solution replaced checkpoint'
    IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'physical convergence recovery rejected'
    IF (ANY(sol%u /= accepted) .OR. ANY(sol%u_conv /= accepted)) &
         ERROR STOP 'physical convergence recovery restored wrong checkpoint'
    IF (divergence_refinements /= 2) &
         ERROR STOP 'physical convergence recovery used a separate budget'

    ! The final NR iterate is still invalid with an unconverged residual.
    ! Recover before leaving the loop, using the remaining shared retry.
    sol%u = admissible_u
    IF (MPIvar%glob_id == bad_rank) sol%u(4) = -1.d-8
    valid = physical_solution_is_admissible(sol%u)
    IF (.NOT. newton_refinement_required(0.1d0,valid,numer%nrp)) &
         ERROR STOP 'inadmissible NR-limit solution did not request refinement'
    IF (update_best_newton_checkpoint(0.1d0,numer%nrp,valid)) &
         ERROR STOP 'inadmissible NR-limit solution replaced checkpoint'
    IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'physical NR-limit recovery rejected'
    IF (ANY(sol%u /= accepted) .OR. ANY(sol%u_conv /= accepted)) &
         ERROR STOP 'physical NR-limit recovery restored wrong checkpoint'
    IF (divergence_refinements /= adapt%max_divergence_refinements) &
         ERROR STOP 'physical NR-limit recovery used a separate budget'
    sol%u = admissible_u
    IF (MPIvar%glob_id == bad_rank) sol%u(4) = -1.d-8
    valid = physical_solution_is_admissible(sol%u)
    IF (.NOT. newton_refinement_required(0.5d0*numer%tNR,valid,1)) &
         ERROR STOP 'invalid convergence accepted after retry limit'
    IF (prepare_divergence_refinement()) ERROR STOP 'physical convergence bypassed retry limit'
    IF (.NOT. newton_refinement_required(0.1d0,valid,numer%nrp)) &
         ERROR STOP 'invalid NR-limit solution accepted after retry limit'
    IF (prepare_divergence_refinement()) ERROR STOP 'physical NR-limit recovery bypassed retry limit'

    sol%u = admissible_u
    IF (MPIvar%glob_id == bad_rank) sol%u(4) = -1.d-8
    IF (physical_solution_is_admissible(sol%u)) ERROR STOP 'negative electron energy accepted'
    sol%u = admissible_u
    IF (MPIvar%glob_id == bad_rank) sol%u(3) = 0.d0
    IF (physical_solution_is_admissible(sol%u)) ERROR STOP 'negative ion internal energy accepted'
    sol%u = admissible_u
    IF (MPIvar%glob_id == bad_rank) sol%u(5) = -1.d-8
    IF (physical_solution_is_admissible(sol%u)) ERROR STOP 'negative neutral density accepted'
    ! Zero electron/neutral energies and cancellation roundoff are allowed.
    sol%u = admissible_u
    rho = sol%u(1);momentum = sol%u(2)
    kinetic = 0.5d0*momentum*(momentum/rho)
    sol%u(3) = kinetic*(1.d0-4.d0*EPSILON(1.d0))
    sol%u(4:5) = 0.d0
    IF (.NOT. physical_solution_is_admissible(sol%u)) ERROR STOP 'roundoff thermal energy rejected'
#ifdef PARALL
    Mesh%ghostElems(2) = 1
    sol%u(6:10) = -1.d0
    IF (.NOT. physical_solution_is_admissible(sol%u)) ERROR STOP 'ghost values selected checkpoint validity'
    Mesh%ghostElems = 0
#endif

    ! Projection/reinitialization invalidates old-mesh buffers without spending
    ! another retry. A subsequent good iterate can seed this mesh's checkpoint.
    sol%u = admissible_u
    sol%u(1) = -1.d-4
    retry_count = divergence_refinements
    CALL seed_divergence_checkpoint()
    IF (divergence_checkpoint_available) ERROR STOP 'invalid projected fallback accepted'
    IF (prepare_divergence_refinement()) ERROR STOP 'old-mesh checkpoint used after bad projection'
    IF (divergence_refinements /= retry_count) ERROR STOP 'unavailable checkpoint spent a retry'
    sol%u = admissible_u
    valid = physical_solution_is_admissible(sol%u)
    IF (.NOT. update_best_newton_checkpoint(0.2d0,3,valid)) ERROR STOP 'good post-projection iterate not accepted'
    IF (.NOT. divergence_checkpoint_available) ERROR STOP 'checkpoint not restored after valid iterate'
    DEALLOCATE(sol%u,sol%q,sol%u_conv,sol%q_conv,uiter_best,qiter_best)
    NULLIFY(uiter_best,qiter_best)
    adapt%max_divergence_refinements = 2
#ifdef PARALL
    DEALLOCATE(Mesh%ghostElems)
#endif
  END SUBROUTINE check_physical_checkpoint_selection

  SUBROUTINE assert_physical_state_unchanged()
    IF (ANY(sol%u0 /= history)) ERROR STOP 'recovery changed time history'
    IF (time%it /= 16 .OR. time%ik /= 16 .OR. sol%Nt /= 16) ERROR STOP 'recovery advanced step'
    IF (time%t /= 0.32d0 .OR. time%dt /= 0.02d0 .OR. &
        time%t_ME /= 0.30d0 .OR. time%dt_ME /= 0.02d0) ERROR STOP 'recovery changed time or equilibrium time'
    IF (phys%puff /= 3.d21 .OR. phys%feedback_integral_error /= 123.d0 .OR. &
        phys%feedback_previous_error /= 321.d0) ERROR STOP 'recovery changed feedback'
  END SUBROUTINE assert_physical_state_unchanged

  REAL*8 FUNCTION history_field(x,eq,h) RESULT(value)
    REAL*8, INTENT(IN) :: x(2)
    INTEGER, INTENT(IN) :: eq,h
    value = 10.d0*h+eq+(0.5d0*h+eq)*x(1)-(h+2.d0*eq)*x(2)
  END FUNCTION history_field
END PROGRAM test_divergence_refinement
