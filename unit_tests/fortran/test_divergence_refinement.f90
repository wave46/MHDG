PROGRAM test_divergence_refinement
  USE Main_utils, ONLY: reset_divergence_refinement, prepare_divergence_refinement, &
       update_uiter_qiter_best, &
       uiter_best, qiter_best, divergence_refinements, &
       project_u0_newmesh, Mesh_prec, ir
  USE globals
  USE MPI_OMP, ONLY: MPIvar
  USE reference_element, ONLY: create_reference_element, free_reference_element_pol
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
    CALL MPI_FINALIZE(ierr_test)
    STOP
  ENDIF
  IF (adapt%max_divergence_refinements /= 2) ERROR STOP 'default divergence refinement limit'
  utils%printint = 0
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

  adapt%adaptivity = .FALSE.
  adapt%div_adapt = .TRUE.
  CALL reset_divergence_refinement()
  IF (prepare_divergence_refinement()) ERROR STOP 'recovery enabled with adaptation off'
  adapt%adaptivity = .TRUE.
  adapt%div_adapt = .FALSE.
  IF (prepare_divergence_refinement()) ERROR STOP 'recovery enabled with div_adapt off'

  adapt%div_adapt = .TRUE.
  CALL reset_divergence_refinement()
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
    CALL reset_divergence_refinement()
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
  CALL reset_divergence_refinement()
  sol%u = 9.d30
  IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'new timestep budget not reset'
  IF (ANY(sol%u /= [71.d0,72.d0])) ERROR STOP 'old timestep fallback retained'

  CALL assert_physical_state_unchanged()
  DEALLOCATE(sol%u,sol%q,sol%u_conv,sol%q_conv,sol%u0,uiter_best,qiter_best)
  NULLIFY(uiter_best,qiter_best)

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
