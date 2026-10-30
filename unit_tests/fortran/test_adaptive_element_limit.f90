PROGRAM test_adaptive_element_limit
  USE globals
  USE MPI_OMP, ONLY: MPIvar
  USE mpi, ONLY: MPI_INIT, MPI_FINALIZE, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_COMM_WORLD, MPI_BCAST, MPI_INTEGER
  USE iso_c_binding, ONLY: c_size_t
  USE gmsh, ONLY: gmsh_t
  USE reference_element, ONLY: create_reference_element
  USE adaptivity_general_module, ONLY: adaptively_refine_mesh
  USE Main_utils, ONLY: adaptivity, mesh_name, order, count_adapt, restart_adapt, ir, Mesh_prec, &
       recover_newton_divergence, recover_newton_domain_failure, divergence_checkpoint_available, &
       divergence_refinements, uiter_best, qiter_best, update_newton_checkpoint_and_refine, &
       physical_iterate_valid, errNR, errNR_adapt, ir_check
  USE solver_status, ONLY: solver_status_t, solve_success, solve_retry
  IMPLICIT NONE
  TYPE(solver_status_t) :: status
  TYPE(gmsh_t) :: mesh_reader
  INTEGER(c_size_t), ALLOCATABLE :: tags(:), nodes(:)
  REAL*8 :: coordinates(4,2), old_time
  INTEGER :: ierr, e, node, candidate_count, old_it, old_ik, old_nt, candidate_unit
  LOGICAL :: recovered
  CHARACTER(LEN=20) :: mode

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIvar%glob_id,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,MPIvar%glob_size,ierr)
  CALL get_command_argument(1,mode)
  utils%timing = .FALSE.
  utils%printint = 0
  Mesh%ndim = 2
  Mesh%elemType = 0
  Mesh%Nelems = 2
  Mesh%Nnodes = 4
  Mesh%Nnodesperelem = 3
  Mesh%Nnodesperface = 2
  Mesh%Nfaces = 5
  Mesh%Nintfaces = 1
  Mesh%Nextfaces = 4
  Mesh%ndir = 0
  ALLOCATE(Mesh%X(4,2),Mesh%T(2,3))
  coordinates(:,1) = [1.d0,1.1d0,1.1d0,1.d0]
  coordinates(:,2) = [0.d0,0.d0,0.1d0,0.1d0]
  Mesh%X = coordinates
  Mesh%T(1,:) = [1,2,3]
  Mesh%T(2,:) = [1,3,4]
#ifdef PARALL
  Mesh%Nel_glob = 2
  Mesh%Nno_glob = 4
  ALLOCATE(Mesh%ghostElems(2),Mesh%loc2glob_el(2),Mesh%loc2glob_nodes(4))
  ! Each rank also holds the other triangle as a ghost.
  DO e=1,2
     Mesh%ghostElems(e) = MERGE(0,1,MOD(e-1,MPIvar%glob_size) == MPIvar%glob_id)
  ENDDO
  Mesh%loc2glob_el = [1,2]
  Mesh%loc2glob_nodes = [1,2,3,4]
#endif
  CALL create_reference_element(refElPol,2,1,verbose=0)
  phys%neq = 5
  phys%npv = 11
  phys%Mref = 1.d0
  phys%idx_rhon_eq = 5
  phys%idx_rhon_pv = 11
  ALLOCATE(sol%u(30),sol%q(60),sol%u_conv(30),sol%q_conv(60))
  DO node=1,6
     sol%u(5*node-4:5*node) = [1.d0,0.d0,3.d0,3.d0,0.1d0]
  ENDDO
  sol%q = 0.d0
  sol%u_conv = sol%u
  sol%q_conv = sol%q
  adapt%adaptivity = .TRUE.
  adapt%evaluator = 1
  adapt%shockcp_adapt = 1
  adapt%quant_ind = 1
  adapt%thr_ind = 1.d-8
  ALLOCATE(adapt%n_quant_ind(1))
  adapt%n_quant_ind = [1]
  adapt%geometry_path = 'square.geo'
  mesh_name = './res/cap_mesh'
  order = 1

  ! Exercise real generation with unlimited, exact and oversized candidate limits.
  IF (mode == 'generate') THEN
     IF (adapt%max_elements /= 0) ERROR STOP 'limit must be optional by default'
     CALL adaptively_refine_mesh(mesh_name,1,order,status)
     IF (status%outcome /= solve_success) ERROR STOP 'unlimited candidate rejected'
     IF (MPIvar%glob_id == 0) THEN
        CALL mesh_reader%initialize()
        CALL mesh_reader%option%setNumber('General.Verbosity',0.d0)
        CALL mesh_reader%open('./res/temp.msh')
        CALL mesh_reader%model%mesh%getElementsByType(2,tags,nodes)
        candidate_count = SIZE(tags)
        DEALLOCATE(tags,nodes)
        CALL mesh_reader%finalize()
     ENDIF
     CALL MPI_BCAST(candidate_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     IF (candidate_count <= 1) ERROR STOP 'fixture must exceed one triangle'
     adapt%max_elements = candidate_count
     CALL adaptively_refine_mesh(mesh_name,2,2,status)
     IF (status%outcome /= solve_success) ERROR STOP 'exact-cap high-order candidate rejected'
     IF (MPIvar%glob_id == 0) THEN
        CALL mesh_reader%initialize()
        CALL mesh_reader%option%setNumber('General.Verbosity',0.d0)
        CALL mesh_reader%open('./res/temp.msh')
        CALL mesh_reader%model%mesh%getElementsByType(9,tags,nodes)
        IF (SIZE(tags) /= candidate_count) ERROR STOP 'cap counted boundary elements or high-order nodes'
        DEALLOCATE(tags,nodes)
        CALL mesh_reader%finalize()
        OPEN(NEWUNIT=candidate_unit, FILE='./res/temp.msh', STATUS='replace', ACTION='write')
        WRITE(candidate_unit,'(A)') 'preserve previous candidate'
        CLOSE(candidate_unit)
     ENDIF
     adapt%max_elements = candidate_count-1
     CALL adaptively_refine_mesh(mesh_name,3,order,status)
     CALL assert_cap_rejection(status)
     ! Candidate generation must never install a mesh, including accepted candidates.
     IF (ANY(Mesh%X /= coordinates) .OR. Mesh%Nelems /= 2) ERROR STOP 'candidate changed active mesh'
  ELSE
     adapt%max_elements = 1
     phys%lscale = 0.01d0
     Mesh%X = coordinates/phys%lscale
     coordinates = Mesh%X
     ALLOCATE(phys%B(4,3))
     phys%B = 17.d0
     count_adapt = 0
     ir = 2
     restart_adapt = .TRUE.
     time%legacy_input = .TRUE.
     time%t = 0.3d0
     time%dt = 0.1d0
     time%it = 3
     time%ik = 3
     sol%Nt = 3
     old_time = time%t
     old_it = time%it
     old_ik = time%ik
     old_nt = sol%Nt
     IF (mode == 'periodic') THEN
        CALL adaptivity()
        ERROR STOP 'periodic cap failure was ignored'
     ENDIF
     CALL adaptivity(status=status)
     CALL assert_cap_rejection(status)
     IF (.NOT. restart_adapt) ERROR STOP 'rejected restart adaptation consumed'
     IF (ABS(time%t-old_time) > 1.d-14 .OR. time%it /= old_it .OR. &
         time%ik /= old_ik .OR. sol%Nt /= old_nt) ERROR STOP 'rejected adaptation advanced legacy clock'
     IF (MAXVAL(ABS(Mesh%X-coordinates)) > 1.d-12) ERROR STOP 'rejected mesh changed coordinate units'
     IF (MAXVAL(ABS(Mesh_prec%X-coordinates)) > 1.d-12) ERROR STOP 'previous mesh changed coordinate units'
     IF (Mesh%Nelems /= 2 .OR. count_adapt /= 0) ERROR STOP 'rejected mesh installed or counted'
     IF (.NOT. ASSOCIATED(phys%B)) ERROR STOP 'rejection freed current field'
     IF (ANY(phys%B /= 17.d0)) ERROR STOP 'rejection modified current field'

     ! Divergence and source-domain recovery must propagate the cap reason.
     restart_adapt = .FALSE.
     adapt%div_adapt = .TRUE.
     divergence_checkpoint_available = .TRUE.
     divergence_refinements = 0
     ALLOCATE(uiter_best, SOURCE=sol%u)
     ALLOCATE(qiter_best, SOURCE=sol%q)
     status = solver_status_t()
     recovered = recover_newton_divergence(status=status)
     IF (recovered) ERROR STOP 'oversized divergence recovery reported success'
     CALL assert_cap_rejection(status)
     status = solver_status_t(solve_retry,0,'Original source-domain failure')
     recovered = recover_newton_domain_failure(status,.FALSE.)
     IF (recovered) ERROR STOP 'oversized source-domain recovery reported success'
     CALL assert_cap_rejection(status)

     ! Ordinary in-Newton remeshing must also return before consuming the rejection.
     adapt%osc_adapt = .FALSE.
     adapt%NR_adapt = .TRUE.
     adapt%freq_NR_adapt = 1
     physical_iterate_valid = .FALSE.
     errNR = 1.d0
     errNR_adapt = 0.5d0
     ir_check = 7
     CALL update_newton_checkpoint_and_refine(status)
     CALL assert_cap_rejection(status)
     IF (ir /= 2 .OR. ir_check /= 7) ERROR STOP 'rejected Newton mesh reset iteration'
  ENDIF
  IF (MPIvar%glob_id == 0) WRITE(*,'(A)') 'ADAPTIVE_ELEMENT_LIMIT_PASS'
  CALL MPI_FINALIZE(ierr)
CONTAINS
  SUBROUTINE assert_cap_rejection(rejection)
    TYPE(solver_status_t), INTENT(IN) :: rejection
    IF (rejection%outcome /= solve_retry) ERROR STOP 'cap failure is not retryable'
    IF (INDEX(rejection%message,'Adaptive mesh exceeds max_elements: candidate=') /= 1 .OR. &
        INDEX(rejection%message,', limit=') == 0) ERROR STOP 'cap rejection reason lost or differs between ranks'
  END SUBROUTINE assert_cap_rejection
END PROGRAM test_adaptive_element_limit
