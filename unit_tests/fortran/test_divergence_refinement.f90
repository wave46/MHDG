PROGRAM test_divergence_refinement
  USE Main_utils, ONLY: reset_newton_refinements, prepare_divergence_refinement, newton_refinement_required, &
       update_uiter_qiter_best, &
       physical_solution_is_admissible, update_best_newton_checkpoint, seed_divergence_checkpoint, &
       divergence_checkpoint_available, errNR_adapt, ir_adapt, &
       uiter_best, qiter_best, divergence_refinements, oscillation_refinements, prepare_oscillation_refinement, &
       project_u0_newmesh, Mesh_prec, ir
  USE adaptivity_projection_module, ONLY: mesh_transfer_t
#ifdef PARALL
  USE Communications, ONLY: gather_mesh, gather_solution
#endif
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
  TYPE(mesh_transfer_t) :: transfer
#ifdef PARALL
  INTEGER, POINTER :: old_connectivity(:,:) => NULL()
  REAL*8, POINTER :: old_coordinates(:,:) => NULL()
#endif

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
    WRITE(*,'(A,1X,L1)') 'CHECK_NEUTRAL_POSITIVITY',adapt%check_neutral_positivity
    CALL MPI_FINALIZE(ierr_test)
    STOP
  ENDIF
  IF (adapt%max_divergence_refinements /= 4) ERROR STOP 'default divergence refinement limit'
  IF (adapt%max_oscillation_refinements /= 10) ERROR STOP 'default oscillation refinement limit'
  IF (.NOT. adapt%check_neutral_positivity) ERROR STOP 'default neutral positivity policy'
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
  CALL check_neutral_positivity_policy()
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
#ifdef PARALL
  CALL gather_mesh(Mesh_prec, old_connectivity, old_coordinates)
  CALL transfer%build(old_connectivity, old_coordinates, Mesh%T, Mesh%X)
  DEALLOCATE(old_connectivity, old_coordinates)
#else
  CALL transfer%build(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X)
#endif
  CALL project_u0_newmesh(preserve_history=.TRUE., transfer=transfer)
  CALL transfer%clear()
  IF (.NOT. ALL(ieee_is_finite(sol%u0))) ERROR STOP 'nonfinite projected history'
  IF (MAXVAL(ABS(sol%u0-expected)) > 1.d-10) ERROR STOP 'time history not preserved on new mesh'
  IF (MAXVAL(ABS(Mesh%X-new_coordinates)) > 1.d-10 .OR. &
      MAXVAL(ABS(Mesh_prec%X-previous_coordinates)) > 1.d-10) ERROR STOP 'projection changed coordinate units'

  ! Oscillation recovery preserves every column too, including when an NR
  ! frequency trigger coincides. Reuse the affine field on the original mesh.
  DEALLOCATE(sol%u0)
  ALLOCATE(sol%u0(6,3))
  DO h=1,3
    DO node=1,3
      DO eq=1,2
        sol%u0((node-1)*2+eq,h) = history_field(Mesh_prec%X(node,:)*phys%lscale,eq,h)
      ENDDO
    ENDDO
  ENDDO
  adapt%osc_adapt = .TRUE.
  adapt%osc_tol = -0.5d0
  CALL project_u0_newmesh(preserve_history=adapt%osc_adapt .AND. -0.4d0 .GT. adapt%osc_tol)
  IF (MAXVAL(ABS(sol%u0-expected)) > 1.d-10) ERROR STOP 'oscillation recovery lost time history'

  ! Callers omitting preserve_history keep their first-column projection behavior.
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
  CALL check_transfer_reuse()
  CALL check_mesh_copy_lifetime()
#ifdef PARALL
  CALL check_transfer_ownership()
#endif
  IF (MPIvar%glob_id == 0) WRITE(*,*) &
    'divergence refinement: checkpoint selection, retry bounds, MPI recovery and history projection PASS; ranks=', &
    MPIvar%glob_size
  CALL MPI_FINALIZE(ierr_test)

CONTAINS

  SUBROUTINE check_mesh_copy_lifetime()
    TYPE(Mesh_type) :: original, saved, restored
    REAL*8 :: coordinates(3,2)
    INTEGER :: trial

    coordinates(:,1) = [0.d0,1.d0,0.d0]+2*MPIvar%glob_id
    coordinates(:,2) = [0.d0,0.d0,1.d0]
    original%Ndim = 2
    original%Nelems = 1
    original%Nnodes = 3
    original%Nnodesperelem = 3
    original%Nnodesperface = 2
    original%Nfaces = 3
    original%Nextfaces = 3
    original%Nintfaces = 0
    original%elemType = 0
    original%ndir = 0
    original%ukf = 3
    original%Nnodes_toroidal = 0
    original%xmin = MINVAL(coordinates(:,1))
    original%xmax = MAXVAL(coordinates(:,1))
    original%ymin = 0.d0
    original%ymax = 1.d0
    original%puff_area = 0.d0
    original%pump_area = 0.d0
    original%core_area = 0.d0
    ALLOCATE(original%X(3,2),original%T(1,3),original%face_info(3,2))
    original%X = coordinates
    original%T(1,:) = [1,2,3]
    original%face_info(:,1) = 1
    original%face_info(:,2) = [1,2,3]
    ALLOCATE(original%flag_elems_rho(1),original%flag_elems_sc(1),original%minrho_elems(1), &
         original%sour_elems(1),original%diff_elems(1))
    original%flag_elems_rho = 1
    original%flag_elems_sc = 1
    original%minrho_elems = 2.d0
    original%sour_elems = 3.d0
    original%diff_elems = 4.d0
#ifdef PARALL
    original%nghostfaces = 0
    original%nghostelems = 0
    original%Nel_glob = MPIvar%glob_size
    original%Nno_glob = 3*MPIvar%glob_size
    original%Nfa_glob = 3*MPIvar%glob_size
    original%Ndir_glob = 0
    original%Ngho_glob = 0
    original%Nextfaces_glob = 3*MPIvar%glob_size
    original%Nintfaces_glob = 0
    original%Nextfaces_nogho = 3
    original%Nintfaces_nogho = 0
    ALLOCATE(original%loc2glob_nodes(3),original%fc2sd(1),original%pr2sd(1))
    original%loc2glob_nodes = [1,2,3]+3*MPIvar%glob_id
    original%fc2sd = 2
    original%pr2sd = MOD(MPIvar%glob_id+1,MPIvar%glob_size)+1
#endif
    CALL deep_copy_mesh_struct(original,saved)
    CALL free_mesh_loc(original)
    DO trial = 1, 2
       CALL deep_copy_mesh_struct(saved,restored)
       IF (ANY(restored%X /= coordinates) .OR. ANY(restored%T(1,:) /= [1,2,3])) &
            ERROR STOP 'mesh copy lost accepted coordinates/connectivity'
       IF (.NOT. ALLOCATED(restored%face_info)) ERROR STOP 'mesh copy lost boundary face information'
       IF (ANY(restored%face_info(:,2) /= [1,2,3])) ERROR STOP 'mesh copy changed boundary face information'
#ifdef PARALL
       IF (restored%Nextfaces_glob /= 3*MPIvar%glob_size .OR. restored%Nintfaces_glob /= 0 .OR. &
           restored%Nextfaces_nogho /= 3 .OR. restored%Nintfaces_nogho /= 0) ERROR STOP 'mesh copy lost MPI face counts'
       IF (ANY(restored%loc2glob_nodes /= [1,2,3]+3*MPIvar%glob_id) .OR. restored%fc2sd(1) /= 2 .OR. &
           restored%pr2sd(1) /= MOD(MPIvar%glob_id+1,MPIvar%glob_size)+1) ERROR STOP 'mesh copy lost MPI exchange maps'
       restored%loc2glob_nodes = -1
       restored%fc2sd = -1
#endif
       restored%X = -1.d0
       restored%T = -1
       restored%face_info = -1
       ALLOCATE(restored%Xg(1,2),restored%Xgf(1,2),restored%Xgb(1,2))
       CALL free_mesh_loc(restored)
       IF (ASSOCIATED(restored%X) .OR. ASSOCIATED(restored%T) .OR. ALLOCATED(restored%face_info) .OR. &
           ALLOCATED(restored%flag_elems_rho) .OR. ALLOCATED(restored%flag_elems_sc) .OR. &
           ALLOCATED(restored%minrho_elems) .OR. ALLOCATED(restored%sour_elems) .OR. &
           ALLOCATED(restored%diff_elems) .OR. ALLOCATED(restored%Xg) .OR. ALLOCATED(restored%Xgf) .OR. &
           ALLOCATED(restored%Xgb)) ERROR STOP 'mesh cleanup retained owned storage'
       CALL free_mesh_loc(restored)
    ENDDO
    CALL free_mesh_loc(saved)
    IF (MPIvar%glob_id == 0) WRITE(*,*) 'Mesh copy: independent restore, MPI face counts/maps and cleanup PASS'
  END SUBROUTINE check_mesh_copy_lifetime

  SUBROUTINE check_transfer_reuse()
    TYPE(mesh_transfer_t) :: mapping
    INTEGER :: degree, n, element, point, equation, component, old_node
    INTEGER, ALLOCATABLE :: old_t(:,:), new_t(:,:)
    REAL*8, ALLOCATABLE :: old_x(:,:), new_x(:,:), parent_points(:,:)
    REAL*8, POINTER :: old_u(:) => NULL(), old_q(:) => NULL(), new_u(:) => NULL(), new_q(:) => NULL()
    REAL*8 :: x, y, expected_value

    ! Curved quadratic geometry puts some targets outside the straight-sided
    ! candidate triangle. Distinct equation/gradient fields expose stride errors.
    DO degree = 1, 3
       CALL create_reference_element(refElPol,2,degree,verbose=0)
       n = refElPol%Nnodes2D
       ALLOCATE(old_t(1,n), new_t(2,n), old_x(n,2), new_x(2*n,2), parent_points(2*n,2))
       ALLOCATE(old_u(n*phys%neq),old_q(n*phys%neq*2))
       old_t(1,:) = [(point,point=1,n)]
       DO point = 1, n
          x = 0.5d0*(refElPol%coord2D(point,1)+1.d0)
          y = 0.5d0*(refElPol%coord2D(point,2)+1.d0)
          old_x(point,:) = curved_coordinates(x,y,degree)
          DO equation = 1, phys%neq
             old_u((point-1)*phys%neq+equation) = transfer_field(x,y,equation,degree)
             DO component = 1, 2
                old_q(((point-1)*phys%neq+equation-1)*2+component) = &
                     transfer_field(x,y,equation+component+3,degree)
             ENDDO
          ENDDO
          parent_points(point,:) = [0.5d0*x,y]
          parent_points(n+point,:) = [0.5d0+0.5d0*x-0.5d0*y,y]
       ENDDO
       DO element = 1, 2
          DO point = 1, n
             old_node = (element-1)*n+point
             new_t(element,point) = old_node
             new_x(old_node,:) = curved_coordinates(parent_points(old_node,1),parent_points(old_node,2),degree)
          ENDDO
       ENDDO
       CALL mapping%build(old_t, old_x, new_t, new_x)
       CALL mapping%apply_to_pointers(old_u, new_u, old_q, new_q)
       DO point = 1, 2*n
          x = parent_points(point,1)
          y = parent_points(point,2)
          DO equation = 1, phys%neq
             expected_value = transfer_field(x,y,equation,degree)
             IF (ABS(new_u((point-1)*phys%neq+equation)-expected_value) > 1.d-9) &
                  ERROR STOP 'cached solution polynomial projection'
             DO component = 1, 2
                expected_value = transfer_field(x,y,equation+component+3,degree)
                IF (ABS(new_q(((point-1)*phys%neq+equation-1)*2+component)-expected_value) > 1.d-9) &
                     ERROR STOP 'cached gradient polynomial projection'
             ENDDO
          ENDDO
       ENDDO
       ! Applying another field must use the same weights without retaining values.
       old_u = -2.d0*old_u+3.d0
       CALL mapping%apply(old_u, new_u)
       DO point = 1, 2*n
          DO equation = 1, phys%neq
             expected_value = -2.d0*transfer_field(parent_points(point,1),parent_points(point,2),equation,degree)+3.d0
             IF (ABS(new_u((point-1)*phys%neq+equation)-expected_value) > 1.d-9) &
                  ERROR STOP 'repeated transfer retained previous solution'
          ENDDO
       ENDDO
       ! A new mesh replaces the old map, including destination size and ordering.
       CALL mapping%build(old_t, old_x, old_t, old_x)
       CALL mapping%apply_to_pointers(old_u, new_u)
       IF (SIZE(new_u) /= SIZE(old_u)) ERROR STOP 'rebuilt transfer destination size'
       IF (MAXVAL(ABS(new_u-old_u)) > 1.d-9) ERROR STOP 'rebuilt transfer retained old geometry'
       CALL mapping%clear()
       DEALLOCATE(old_t,new_t,old_x,new_x,parent_points,old_u,old_q,new_u,new_q)
       CALL free_reference_element_pol(refElPol)
    ENDDO
    IF (MPIvar%glob_id == 0) WRITE(*,*) 'Mesh transfer: P1-P3, curved geometry, u/q, repeated apply and rebuild PASS'
  END SUBROUTINE check_transfer_reuse

  FUNCTION curved_coordinates(x,y,degree) RESULT(coordinates)
    REAL*8, INTENT(IN) :: x,y
    INTEGER, INTENT(IN) :: degree
    REAL*8 :: coordinates(2)

    coordinates = [x,y]
    IF (degree > 1) coordinates(2) = y-0.3d0*x*(1.d0-x-y)
  END FUNCTION curved_coordinates

  REAL*8 FUNCTION transfer_field(x,y,equation,degree) RESULT(value)
    REAL*8, INTENT(IN) :: x,y
    INTEGER, INTENT(IN) :: equation,degree

    value = equation+0.3d0*equation*x-0.7d0*y
    IF (degree >= 2) value = value+0.2d0*x*y-0.1d0*equation*x*x
    IF (degree == 3) value = value+0.1d0*equation*x*y*y
  END FUNCTION transfer_field

#ifdef PARALL
  SUBROUTINE check_transfer_ownership()
    TYPE(Mesh_type) :: source_mesh
    TYPE(mesh_transfer_t) :: mapping
    INTEGER, POINTER :: source_t(:,:) => NULL()
    REAL*8, POINTER :: source_x(:,:) => NULL(), source_u(:) => NULL(), source_q(:) => NULL()
    REAL*8, POINTER :: projected_u(:) => NULL(), projected_q(:) => NULL()
    REAL*8 :: local_u(12), local_q(24), destination_x(6,2), point_x(2), expected_value
    INTEGER :: destination_t(2,3), owner, point, equation, component, index

    CALL create_reference_element(refElPol,2,1,verbose=0)
    source_mesh%Ndim = 2
    source_mesh%Nelems = 2
    source_mesh%Nnodes = 6
    source_mesh%Nnodesperelem = 3
    source_mesh%Nel_glob = MPIvar%glob_size
    source_mesh%Nno_glob = 3*MPIvar%glob_size
    ALLOCATE(source_mesh%X(6,2),source_mesh%T(2,3),source_mesh%ghostElems(2))
    ALLOCATE(source_mesh%loc2glob_el(2),source_mesh%loc2glob_nodes(6))
    owner = MOD(MPIvar%glob_id+1,MPIvar%glob_size)
    source_mesh%ghostElems = [0,1]
    source_mesh%loc2glob_el = [MPIvar%glob_id+1,owner+1]
    source_mesh%loc2glob_nodes(:3) = [(point,point=1,3)]+3*MPIvar%glob_id
    source_mesh%loc2glob_nodes(4:) = [(point,point=1,3)]+3*owner
    source_mesh%T(1,:) = [1,2,3]
    source_mesh%T(2,:) = [4,5,6]
    destination_t = source_mesh%T
    ! Wrong ghost geometry and values must never enter the global source.
    source_mesh%X = -1.d6
    local_u = -1.d6
    local_q = -1.d6
    DO point = 1, 3
       point_x = 0.5d0*(refElPol%coord2D(point,:)+1.d0)
       source_mesh%X(point,:) = point_x+[2.d0*MPIvar%glob_id,0.d0]
       destination_x(point,:) = point_x+[2.d0*owner,0.d0]
       destination_x(point+3,:) = source_mesh%X(point,:)
       DO equation = 1, phys%neq
          index = (point-1)*phys%neq+equation
          local_u(index) = equation+SUM(source_mesh%X(point,:))
          local_q(2*index-1:2*index) = [2.d0,3.d0]*local_u(index)
       ENDDO
    ENDDO
    CALL gather_mesh(source_mesh,source_t,source_x)
    CALL mapping%build(source_t,source_x,destination_t,destination_x)
    CALL gather_solution(source_mesh,3,u_in=local_u,q_in=local_q,u_glob=source_u,q_glob=source_q)
    CALL mapping%apply_to_pointers(source_u,projected_u,source_q,projected_q)
    DO point = 1, 6
       DO equation = 1, phys%neq
          index = (point-1)*phys%neq+equation
          expected_value = equation+SUM(destination_x(point,:))
          IF (ABS(projected_u(index)-expected_value) > 1.d-9) ERROR STOP 'transfer source ownership'
          DO component = 1, 2
             IF (ABS(projected_q(2*index-2+component)-(component+1)*expected_value) > 1.d-9) THEN
                WRITE(*,*) 'Ghost gradient rank/point/equation/component/actual/expected: ', &
                     MPIvar%glob_id,point,equation,component,projected_q(2*index-2+component),(component+1)*expected_value
                ERROR STOP 'transfer ghost destination gradient'
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(source_t,source_x,source_u,source_q,projected_u,projected_q)
    DEALLOCATE(source_mesh%X,source_mesh%T,source_mesh%ghostElems,source_mesh%loc2glob_el,source_mesh%loc2glob_nodes)
    CALL free_reference_element_pol(refElPol)
    IF (MPIvar%glob_id == 0) WRITE(*,*) 'Mesh transfer: owned source and local/remote ghost destinations PASS'
  END SUBROUTINE check_transfer_ownership
#endif

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
    ! The relaxed policy must also allow a smooth iterate with negative raw
    ! neutrals to refresh the oscillation checkpoint, without changing its u.
    IF (MPIvar%glob_id == MPIvar%glob_size-1) sol%u(5) = -1.d-8
    adapt%check_neutral_positivity = .FALSE.
    allowed = physical_solution_is_admissible(sol%u)
    IF (.NOT. allowed) ERROR STOP 'relaxed neutral policy rejected smooth checkpoint'
    CALL compute_error_oscillations(oscillations,minimum,maximum,n_osc,10,checkpoint_iteration, &
         checkpoint_mesh,checkpoint_allowed=allowed)
    IF (checkpoint_iteration /= 10 .OR. ANY(sol%u_conv /= sol%u)) &
         ERROR STOP 'relaxed oscillation checkpoint lost raw neutral values'
    adapt%check_neutral_positivity = .TRUE.
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

  SUBROUTINE check_neutral_positivity_policy()
    REAL*8 :: accepted(5)
    LOGICAL :: valid
    INTEGER :: bad_rank
    phys%neq = 5
    phys%idx_rhon_eq = 5
    Mesh%Nelems = 1
    Mesh%Nnodesperelem = 1
#ifdef PARALL
    ALLOCATE(Mesh%ghostElems(1))
    Mesh%ghostElems = 0
#endif
    ALLOCATE(sol%u(5),sol%q(10),sol%u_conv(5),sol%q_conv(10))
    accepted = [1.d0,0.2d0,2.d0,3.d0,0.01d0]
    bad_rank = MPIvar%glob_size-1
    IF (MPIvar%glob_id == bad_rank) accepted(5) = -1.d-8
    sol%u = accepted
    sol%q = 0.1d0
    adapt%check_neutral_positivity = .TRUE.
    CALL reset_newton_refinements()
    IF (divergence_checkpoint_available) ERROR STOP 'strict policy seeded negative-neutral checkpoint'
    valid = physical_solution_is_admissible(sol%u)
    IF (.NOT. newton_refinement_required(0.5d0*numer%tNR,valid,1)) &
         ERROR STOP 'strict neutral convergence policy changed'

    adapt%check_neutral_positivity = .FALSE.
    ! This is also the post-projection reseeding path that formerly lost the
    ! fallback solely because neutrals were negative on one rank.
    CALL seed_divergence_checkpoint()
    IF (.NOT. divergence_checkpoint_available .OR. ANY(uiter_best /= accepted)) &
         ERROR STOP 'relaxed policy did not seed projected checkpoint'
    valid = physical_solution_is_admissible(sol%u)
    IF (.NOT. valid) ERROR STOP 'rank-local negative neutrals rejected under relaxed policy'
    IF (newton_refinement_required(0.5d0*numer%tNR,valid,2)) &
         ERROR STOP 'negative neutrals alone forced convergence refinement'
    IF (newton_refinement_required(0.1d0,valid,numer%nrp)) &
         ERROR STOP 'negative neutrals alone forced NR-limit refinement'
    IF (.NOT. update_best_newton_checkpoint(0.5d0*numer%tNR,2,valid)) &
         ERROR STOP 'relaxed policy rejected best Newton iterate'

    ! Relaxing neutrals must leave all plasma checks and genuine divergence
    ! recovery active. Reject worse physics even with a smaller residual.
    IF (MPIvar%glob_id == bad_rank) sol%u(1) = -1.d-4
    valid = physical_solution_is_admissible(sol%u)
    IF (valid) ERROR STOP 'relaxed neutral policy accepted negative plasma density'
    IF (update_best_newton_checkpoint(0.25d0*numer%tNR,3,valid)) &
         ERROR STOP 'relaxed neutral policy replaced checkpoint with invalid plasma'
    IF (ANY(uiter_best /= accepted)) ERROR STOP 'invalid plasma changed relaxed checkpoint'
    sol%u = accepted
    IF (MPIvar%glob_id == bad_rank) sol%u(4) = -1.d-8
    IF (physical_solution_is_admissible(sol%u)) ERROR STOP 'relaxed neutral policy accepted negative electron energy'
    sol%u = accepted
    IF (MPIvar%glob_id == bad_rank) sol%u(3) = 0.d0
    IF (physical_solution_is_admissible(sol%u)) ERROR STOP 'relaxed neutral policy accepted negative ion internal energy'
    IF (.NOT. newton_refinement_required(2.d0*numer%div,.TRUE.,3)) &
         ERROR STOP 'relaxed neutral policy disabled divergence trigger'
    IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'relaxed checkpoint unavailable for divergence'
    IF (ANY(sol%u /= accepted)) ERROR STOP 'recovery altered negative raw neutral values'

    ! A completed state with negative neutrals can seed the next timestep too.
    CALL reset_newton_refinements()
    IF (.NOT. divergence_checkpoint_available .OR. divergence_refinements /= 0) &
         ERROR STOP 'relaxed next-timestep checkpoint unavailable'
    IF (.NOT. prepare_divergence_refinement()) ERROR STOP 'relaxed next-timestep recovery unavailable'
    IF (ANY(sol%u /= accepted)) ERROR STOP 'next-timestep recovery changed raw neutrals'
    adapt%check_neutral_positivity = .TRUE.
    DEALLOCATE(sol%u,sol%q,sol%u_conv,sol%q_conv,uiter_best,qiter_best)
    NULLIFY(uiter_best,qiter_best)
#ifdef PARALL
    DEALLOCATE(Mesh%ghostElems)
#endif
    IF (MPIvar%glob_id == 0) WRITE(*,*) 'Neutral positivity policy: strict/relaxed convergence and checkpoints PASS'
  END SUBROUTINE check_neutral_positivity_policy

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
