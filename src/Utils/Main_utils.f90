!*****************************************
! project: MHDG
! file: Main_utils.f90
! date: 25/11/2024
! Set of functions and subroutine used in MHDG.f90

MODULE Main_utils

  USE solver_status, ONLY: solver_status_t, solve_success, solve_retry, finite_solver_values, &
       agree_solver_status, stop_on_solver_failure
  USE adaptivity_indicator_module, ONLY: compute_error_oscillations
  USE input_files, ONLY: resolve_input_filename
  USE in_out
  USE GMSH_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, convert_gmsh_to_hdf5, gmsh_mesh2d_write, hdf5_save_mesh_struct
  USE reference_element
  USE preprocess
  USE MPI_OMP
  USE printutils
  USE debug
  USE initialization
  USE time_stepping, ONLY: next_physical_time, physical_time_is_complete, &
       schedule_half_step, remove_reached_target, ensure_time_history_capacity
  USE time_stepping, ONLY: update_timestep_growth, reset_timestep_growth
  USE newton_predictor_module, ONLY: newton_predictor_t
  USE adaptivity_common_module, ONLY: set_order_mesh
  USE adaptivity_projection_module, ONLY: mesh_transfer_t
  USE external_heating
#ifdef WITH_PETSC
   USE solve_petsc, only: matPETSC, InitPETSC, FinalizePETSC
#endif
  USE adaptivity_general_module
  USE Postprocess, only: computeL2ErrorAnalyticSol
#ifdef PARALL
  USE Communications
  USE domain_decomposition_module
#endif
  USE HDG_LimitingTechniques
  USE flux_surface_transport_data
  USE transport_models_1d
  USE residual_norm_module, ONLY: computeResidual
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

  IMPLICIT NONE

  TYPE :: newton_result_t
     TYPE(solver_status_t) :: status
     INTEGER :: iterations = 0
     REAL*8 :: first_error = 0.d0
     LOGICAL :: mesh_changed = .FALSE.
     LOGICAL :: iteration_limit_reached = .FALSE.
     LOGICAL :: predictor_used = .FALSE.
     LOGICAL :: predictor_fallback = .FALSE.
  END TYPE newton_result_t

  TYPE :: accepted_step_state_t
     PRIVATE
     TYPE(Time_type) :: clock
     INTEGER :: steps
     TYPE(Mesh_type), ALLOCATABLE :: original_mesh
     REAL*8, ALLOCATABLE :: u(:), q(:), trace(:), impurity_concentrations(:)
     REAL*8, ALLOCATABLE :: u0(:,:)
     REAL*8 :: puff, recycling, recycling_pump, zeff, line_density
     REAL*8 :: integral_error, previous_error, integral_error_xpr, previous_error_xpr
     TYPE(Feedback_derivative_history_type) :: derivative, derivative_xpr
   CONTAINS
     PROCEDURE :: capture => capture_accepted_step
     PROCEDURE :: restore => restore_accepted_step
     PROCEDURE :: restore_plasma => restore_accepted_plasma
     PROCEDURE :: capture_mesh => capture_accepted_mesh
     PROCEDURE :: restore_mesh => restore_accepted_mesh
     FINAL :: release_accepted_mesh
  END TYPE accepted_step_state_t

  INTEGER                      :: Np, Nel, Nfp, Nf, Ndim, Nthreads
  INTEGER                      :: it, ir, ir_check, it0, nts, nu, nut, nb_args, IERR, k, i,is, count_adapt = 0, ir_adapt = 0, checkpoint = 1, divergence_counter_adapt = 0, convergence_counter = 1, order
  INTEGER                      :: divergence_refinements = 0
  INTEGER                      :: oscillation_refinements = 0
  LOGICAL                      :: divergence_checkpoint_available = .FALSE.
  LOGICAL                      :: physical_iterate_valid = .TRUE.

  LOGICAL, ALLOCATABLE         :: mkelms(:)
  REAL*8                       :: dt, dt0, errNR, errNR_adapt, errlstime
  REAL*8, ALLOCATABLE          :: uiter(:), L2err(:)
  REAL*8, POINTER              :: uiter_best(:) => NULL(), qiter_best(:) => NULL()
  CHARACTER(LEN=1024)          :: mesh_name, mesh_input_name, save_name
  CHARACTER ( len = 255 )      :: gmsh_filename
  CHARACTER ( len = 255 )      :: h5_filename
  CHARACTER ( len = 50 )       :: count_adapt_char
  REAL*8                       :: cputtot, runttot
  INTEGER                      :: OMP_GET_MAX_THREADS
  INTEGER                      :: clock_rate,clock_start, clock_end
  REAL*8                       :: time_start, time_finish
  REAL*8, ALLOCATABLE          :: oscillations(:)
  REAL*8                       :: max_osc = -100., min_osc = -100.
  INTEGER,ALLOCATABLE          :: vector_nodes_unique(:,:)
  INTEGER                      :: N_n_vertex, n_osc
  INTEGER*8                    :: file_id
  TYPE(Mesh_type)              :: Mesh_prec
  LOGICAL                      :: restart_adapt
#ifdef PARALL
  TYPE(Mesh_type)              :: Mesh_glob
  REAL*8, POINTER              :: X_glob(:,:) => NULL(), B_glob(:,:) => NULL()
  REAL*8, POINTER              :: u_glob(:) => NULL(), q_glob(:) => NULL(), magnetic_flux_glob(:) => NULL(), u_tilde_glob(:) => NULL()
  INTEGER, POINTER             :: T_glob(:,:) => NULL(), Tb_glob(:,:) => NULL(), F_glob(:,:) => NULL(), N_glob(:,:) => NULL(), flipface_glob(:,:) => NULL(), Tlin_glob(:,:) => NULL(), &
                                  & boundaryFlag_glob(:) => NULL(), extFaces_Glob(:,:) => NULL(), intFaces_Glob(:,:) => NULL()
  CHARACTER(70)                :: npr, nid
#endif
#ifdef TOR3D
  INTEGER                      :: Neq, ntorloc,Nfl, Np1dPol, Np1dTor, Ng1dPol, Ng1dTor,
#endif
  TYPE(Reference_element_type) :: refElPol_prec
  TYPE(newton_predictor_t), PRIVATE :: predictor_history


CONTAINS

  SUBROUTINE capture_accepted_step(this)
    CLASS(accepted_step_state_t), INTENT(OUT) :: this

    ! Capture before selecting a trial endpoint or updating its controls.
    this%clock = time
    this%steps = sol%Nt
    this%u = sol%u
    this%q = sol%q
    this%trace = sol%u_tilde
    this%u0 = sol%u0
    IF (ALLOCATED(phys%impurity_concentrations)) this%impurity_concentrations = phys%impurity_concentrations
    this%puff = phys%puff
    this%recycling = phys%Re
    this%recycling_pump = phys%Re_pump
    this%zeff = phys%Zeff
    this%line_density = phys%n_li
    this%integral_error = phys%feedback_integral_error
    this%previous_error = phys%feedback_previous_error
    this%integral_error_xpr = phys%feedback_integral_error_xpr
    this%previous_error_xpr = phys%feedback_previous_error_xpr
    this%derivative = phys%feedback_derivative_history
    this%derivative_xpr = phys%feedback_derivative_history_xpr
  END SUBROUTINE capture_accepted_step

  SUBROUTINE restore_accepted_step(this)
    CLASS(accepted_step_state_t), INTENT(IN) :: this

    CALL this%restore_plasma()
    time = this%clock
    sol%Nt = this%steps
    IF (ALLOCATED(this%impurity_concentrations)) THEN
       phys%impurity_concentrations = this%impurity_concentrations
    ELSEIF (ALLOCATED(phys%impurity_concentrations)) THEN
       DEALLOCATE(phys%impurity_concentrations)
    ENDIF
    phys%puff = this%puff
    phys%Re = this%recycling
    phys%Re_pump = this%recycling_pump
    phys%Zeff = this%zeff
    phys%n_li = this%line_density
    phys%feedback_integral_error = this%integral_error
    phys%feedback_previous_error = this%previous_error
    phys%feedback_integral_error_xpr = this%integral_error_xpr
    phys%feedback_previous_error_xpr = this%previous_error_xpr
    phys%feedback_derivative_history = this%derivative
    phys%feedback_derivative_history_xpr = this%derivative_xpr
  END SUBROUTINE restore_accepted_step

  SUBROUTINE restore_accepted_plasma(this)
    CLASS(accepted_step_state_t), INTENT(IN) :: this

    ! Restore on the accepted mesh, without changing trial clocks or actuators.
    IF (ASSOCIATED(sol%u)) DEALLOCATE(sol%u)
    IF (ASSOCIATED(sol%q)) DEALLOCATE(sol%q)
    IF (ASSOCIATED(sol%u_tilde)) DEALLOCATE(sol%u_tilde)
    ALLOCATE(sol%u, SOURCE=this%u)
    ALLOCATE(sol%q, SOURCE=this%q)
    ALLOCATE(sol%u_tilde, SOURCE=this%trace)
    sol%u0 = this%u0
  END SUBROUTINE restore_accepted_plasma

  SUBROUTINE capture_accepted_mesh(this)
    CLASS(accepted_step_state_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%original_mesh)) RETURN
    ALLOCATE(this%original_mesh)
    CALL deep_copy_mesh_struct(Mesh,this%original_mesh)
  END SUBROUTINE capture_accepted_mesh

  SUBROUTINE restore_accepted_mesh(this, changed)
    CLASS(accepted_step_state_t), INTENT(IN) :: this
    LOGICAL, INTENT(OUT) :: changed

    changed = ALLOCATED(this%original_mesh)
    IF (.NOT. changed) RETURN
    CALL free_mesh_loc(Mesh)
    CALL deep_copy_mesh_struct(this%original_mesh,Mesh)
  END SUBROUTINE restore_accepted_mesh

  SUBROUTINE release_accepted_mesh(this)
    TYPE(accepted_step_state_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%original_mesh)) CALL free_mesh_loc(this%original_mesh)
  END SUBROUTINE release_accepted_mesh

  SUBROUTINE advance_time_step(result, finished)
    TYPE(newton_result_t), INTENT(OUT) :: result
    LOGICAL, INTENT(OUT) :: finished
    REAL*8 :: failed_endpoint, previous_proposal
    INTEGER :: rejections, status
    LOGICAL :: grew, clean_step
    CHARACTER(LEN=160) :: message

    rejections = 0
    DO
       CALL solve_time_step_trial(result, finished, failed_endpoint)
       IF (finished) RETURN
       IF (time%legacy_input) RETURN
       IF (result%status%outcome == solve_success) THEN
          CALL remove_reached_target()
          previous_proposal = time%dt_next
          clean_step = rejections == 0 .AND. .NOT. result%mesh_changed .AND. .NOT. result%predictor_fallback
          IF (adapt%adaptivity .AND. adapt%time_adapt) THEN
             IF (MOD(it,adapt%freq_t_adapt) == 0) clean_step = .FALSE.
          ENDIF
          CALL update_timestep_growth(result%iterations, result%first_error, &
               clean_step, grew)
          IF (grew .AND. MPIvar%glob_id == 0) WRITE(6,'(/,A,ES16.8,A,ES16.8,A,ES16.8,A,/)') &
               'TIMESTEP INCREASED at ', time%t*simpar%refval_time, ' s: ', &
               previous_proposal*simpar%refval_time, ' -> ', time%dt_next*simpar%refval_time, ' s'
          RETURN
       ENDIF
       IF (result%status%outcome /= solve_retry .OR. .NOT. time%dt_reduce_on_failure) RETURN

       rejections = rejections+1
       CALL schedule_half_step(failed_endpoint, status, message)
       IF (status /= 0) THEN
          IF (MPIvar%glob_id == 0) WRITE(6,*) 'Cannot retry timestep: ', TRIM(message), &
               '; last nonlinear failure: ', TRIM(result%status%message)
          result%status%message = message
          RETURN
       ENDIF
       CALL print_timestep_retry(failed_endpoint, rejections, result%status%message)
    ENDDO
  END SUBROUTINE advance_time_step

  SUBROUTINE print_timestep_retry(failed_endpoint, rejection, reason)
    REAL*8, INTENT(IN) :: failed_endpoint
    INTEGER, INTENT(IN) :: rejection
    CHARACTER(LEN=*), INTENT(IN) :: reason
    REAL*8 :: accepted_time, failed_time, retry_time

    IF (MPIvar%glob_id /= 0) RETURN
    accepted_time = time%t*simpar%refval_time
    failed_time = failed_endpoint*simpar%refval_time
    retry_time = accepted_time+(failed_time-accepted_time)/2.d0
    WRITE(6,'(/,A)') '***************************************************************'
    WRITE(6,'(A,I0,2A)') 'TIMESTEP REJECTED (', rejection, '): ', TRIM(reason)
    WRITE(6,'(A,ES16.8,A)') '  Restored accepted time: ', accepted_time, ' s'
    WRITE(6,'(A,ES16.8,A)') '  Failed endpoint:        ', failed_time, ' s'
    WRITE(6,'(A,ES16.8,A)') '  Retry endpoint:         ', retry_time, ' s'
    WRITE(6,'(A,ES16.8,A,ES16.8,A)') '  Timestep:               ', &
         failed_time-accepted_time, ' -> ', retry_time-accepted_time, ' s'
    WRITE(6,'(A,/)') '***************************************************************'
  END SUBROUTINE print_timestep_retry

  SUBROUTINE solve_time_step_trial(result, finished, attempted_endpoint)
    TYPE(newton_result_t), INTENT(OUT) :: result
    LOGICAL, INTENT(OUT) :: finished
    REAL*8, INTENT(OUT) :: attempted_endpoint
    TYPE(accepted_step_state_t) :: accepted
    LOGICAL :: accepted_restart_adapt

    IF (.NOT. time%legacy_input) THEN
       CALL accepted%capture()
       accepted_restart_adapt = restart_adapt
    ENDIF
    CALL prepare_time_step(finished)
    attempted_endpoint = time%t_target
    IF (finished) RETURN
    CALL update_uconv_qconv(sol%u,sol%q)
    CALL reset_newton_refinements()
    CALL update_time_dependent_inputs()
    IF (time%legacy_input) THEN
       time%t = time%t+time%dt
    ELSE
       time%t = time%t_target
    ENDIF
    time%it = time%it+1
    time%ik = time%ik+1
    CALL ensure_time_history_capacity(sol%Nt+1)
    sol%Nt = sol%Nt+1
    IF (MPIvar%glob_id == 0) THEN
       WRITE(6,'(" *",60("*"),"**")')
       WRITE(6,'(" *",20X,"Time iteration   = ",I5,16X," *")') time%it
       WRITE(6,'(" *",60("*"),"**")')
    ENDIF
    IF (time%legacy_input) THEN
       CALL solve_nonlinear_system(result)
    ELSE
       CALL solve_predicted_nonlinear_system(result,accepted)
       result%mesh_changed = ALLOCATED(accepted%original_mesh)
       IF (result%status%outcome == solve_retry) THEN
          CALL rollback_time_step(accepted)
          restart_adapt = accepted_restart_adapt
          CALL predictor_history%clear()
       ELSEIF (time%newton_predictor .AND. result%status%outcome == solve_success) THEN
          IF (result%mesh_changed .OR. result%predictor_fallback) THEN
             CALL predictor_history%clear()
          ELSE
             CALL predictor_history%remember(accepted%u,accepted%q,accepted%trace,accepted%clock%t)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE solve_time_step_trial

  SUBROUTINE solve_predicted_nonlinear_system(result, accepted)
    TYPE(newton_result_t), INTENT(OUT) :: result
    TYPE(accepted_step_state_t), INTENT(INOUT) :: accepted
    LOGICAL :: predicted
    INTEGER :: predicted_iterations, initial_checkpoint_iteration
    REAL*8 :: weight

    initial_checkpoint_iteration = ir_check
    CALL apply_newton_prediction(accepted,predictor_history,predicted,weight)
    CALL solve_nonlinear_system(result,accepted,predicted)
    IF (.NOT. predicted) RETURN
    result%predictor_used = .TRUE.
    IF (result%status%outcome /= solve_retry) RETURN

    IF (MPIvar%glob_id == 0) WRITE(6,*) 'NEWTON PREDICTOR FALLBACK at [s]: ', &
         time%t*simpar%refval_time, '; ', TRIM(result%status%message)
    predicted_iterations = result%iterations
    CALL reset_predicted_newton_solve(accepted,initial_checkpoint_iteration)
    CALL solve_nonlinear_system(result,accepted)
    result%iterations = result%iterations+predicted_iterations
    result%predictor_used = .TRUE.
    result%predictor_fallback = .TRUE.
  END SUBROUTINE solve_predicted_nonlinear_system

  SUBROUTINE apply_newton_prediction(accepted, history, predicted, weight)
    TYPE(accepted_step_state_t), INTENT(IN) :: accepted
    TYPE(newton_predictor_t), INTENT(IN) :: history
    LOGICAL, INTENT(OUT) :: predicted
    REAL*8, INTENT(OUT) :: weight
    TYPE(solver_status_t) :: status
    INTEGER :: reduction
    REAL*8 :: rate_change
    LOGICAL :: rates_available

    predicted = .FALSE.
    weight = 0.d0
    IF (.NOT. time%newton_predictor .OR. restart_adapt) RETURN
    CALL measure_predictor_rate_change(history,accepted%u,accepted%q,accepted%trace, &
         accepted%clock%t,rate_change,rates_available)
    IF (.NOT. rates_available) THEN
       weight = 0.d0
       RETURN
    ENDIF
    IF (rate_change > time%newton_predictor_rate_limit) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'NEWTON PREDICTOR SKIPPED: relative rate change ', &
            rate_change, ' exceeds limit ',time%newton_predictor_rate_limit
       weight = 0.d0
       RETURN
    ENDIF
    weight = history%weight(accepted%clock%t,time%dt)
    IF (weight <= 0.d0) RETURN
    DO reduction = 0, 3
       status = solver_status_t()
       IF (.NOT. history%apply(sol%u,sol%q,sol%u_tilde,weight)) &
            status = solver_status_t(solve_retry,0,'Predictor history is unavailable on this mesh')
       CALL agree_solver_status(status)
       IF (status%outcome /= solve_success) EXIT
       CALL check_newton_values(status)
       IF (status%outcome == solve_success) THEN
          predicted = physical_solution_is_admissible(sol%u)
          IF (predicted) EXIT
       ENDIF
       CALL accepted%restore_plasma()
       weight = weight/2.d0
    ENDDO
    IF (.NOT. predicted) THEN
       CALL accepted%restore_plasma()
       weight = 0.d0
    ELSE
       ! The first NR correction is measured from this guess, not from BDF u0.
       uiter = sol%u
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'NEWTON PREDICTOR weight: ', weight
    ENDIF
  END SUBROUTINE apply_newton_prediction

  SUBROUTINE measure_predictor_rate_change(history,u,q,trace,accepted_time,change,ready)
    TYPE(newton_predictor_t), INTENT(IN) :: history
    REAL*8, INTENT(IN) :: u(:),q(:),trace(:),accepted_time
    REAL*8, INTENT(OUT) :: change
    LOGICAL, INTENT(OUT) :: ready
    LOGICAL, ALLOCATABLE :: volume_owned(:),trace_owned(:)
    REAL*8 :: sums(3,phys%neq,3),denominator
    INTEGER :: equation,field
#ifdef PARALL
    INTEGER :: ierr_rates
#endif

    change = 0.d0
    CALL predictor_owned_points(SIZE(u)/phys%neq,SIZE(trace)/phys%neq,volume_owned,trace_owned)
    CALL history%rate_sums(u,q,trace,accepted_time,phys%neq,volume_owned,trace_owned,sums,ready)
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,ready,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr_rates)
#endif
    IF (.NOT. ready) RETURN
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,sums,SIZE(sums),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr_rates)
#endif
    ready = finite_solver_values(RESHAPE(sums,[SIZE(sums)]))
    IF (.NOT. ready) RETURN
    DO field = 1, 3
       DO equation = 1, phys%neq
          denominator = MAX(sums(1,equation,field),sums(2,equation,field))
          IF (denominator > 0.d0) change = MAX(change,SQRT(sums(3,equation,field)/denominator))
       ENDDO
    ENDDO
  END SUBROUTINE measure_predictor_rate_change

  SUBROUTINE predictor_owned_points(volume_points,trace_points,volume_owned,trace_owned)
    INTEGER, INTENT(IN) :: volume_points,trace_points
    LOGICAL, ALLOCATABLE, INTENT(OUT) :: volume_owned(:),trace_owned(:)
    INTEGER :: element,face,first,itor,tor_owned,nodes_per_element
#ifdef TOR3D
    INTEGER :: block_size,offset
#endif

    ALLOCATE(volume_owned(volume_points),trace_owned(trace_points))
    volume_owned = .FALSE.
    trace_owned = .FALSE.
    tor_owned = 1
    nodes_per_element = Mesh%Nnodesperelem
#ifdef TOR3D
    tor_owned = numer%ntor/MPIvar%ntor
    nodes_per_element = refElTor%Nnodes3D
#endif
    DO itor = 1, tor_owned
       DO element = 1, Mesh%Nelems
#ifdef PARALL
          IF (Mesh%ghostElems(element) /= 0) CYCLE
#endif
          first = ((itor-1)*Mesh%Nelems+element-1)*nodes_per_element+1
          volume_owned(first:first+nodes_per_element-1) = .TRUE.
       ENDDO
    ENDDO
#ifdef TOR3D
    ! Own each left poloidal face and lateral face once; omit the toroidal halo.
    block_size = Mesh%Nelems*refElPol%Nnodes2D+(Mesh%Nfaces-Mesh%Ndir)*refElTor%Nfl
    DO itor = 1, tor_owned
       offset = (itor-1)*block_size
       DO element = 1, Mesh%Nelems
#ifdef PARALL
          IF (Mesh%ghostElems(element) /= 0) CYCLE
#endif
          first = offset+(element-1)*refElPol%Nnodes2D+1
          trace_owned(first:first+refElPol%Nnodes2D-1) = .TRUE.
       ENDDO
       DO face = 1, Mesh%Nfaces-Mesh%Ndir
#ifdef PARALL
          IF (Mesh%ghostFaces(face) /= 0) CYCLE
#endif
          first = offset+Mesh%Nelems*refElPol%Nnodes2D+(face-1)*refElTor%Nfl+1
          trace_owned(first:first+refElTor%Nfl-1) = .TRUE.
       ENDDO
    ENDDO
#else
    DO face = 1, Mesh%Nfaces
#ifdef PARALL
       IF (Mesh%ghostFaces(face) /= 0) CYCLE
#endif
       first = (face-1)*Mesh%Nnodesperface+1
       trace_owned(first:first+Mesh%Nnodesperface-1) = .TRUE.
    ENDDO
#endif
  END SUBROUTINE predictor_owned_points

  SUBROUTINE reset_predicted_newton_solve(accepted,checkpoint_iteration)
    TYPE(accepted_step_state_t), INTENT(IN) :: accepted
    INTEGER, INTENT(IN) :: checkpoint_iteration

    ! Predicted solves stop before remeshing; keep target time, controls and budgets.
    CALL discard_trial_linear_system()
    CALL accepted%restore_plasma()
    CALL update_uconv_qconv(sol%u,sol%q)
    errNR_adapt = HUGE(1.d0)
    ir_adapt = 0
    ir_check = checkpoint_iteration
    divergence_counter_adapt = 0
    physical_iterate_valid = .TRUE.
    CALL seed_divergence_checkpoint()
    CALL free_mesh_loc(Mesh_prec)
    CALL deep_copy_mesh_struct(Mesh,Mesh_prec)
    IF (ASSOCIATED(sol%u_tilde0)) sol%u_tilde0 = 0.d0
    CALL reset_timestep_growth()
  END SUBROUTINE reset_predicted_newton_solve

  SUBROUTINE rollback_time_step(accepted)
    TYPE(accepted_step_state_t), INTENT(IN) :: accepted
    LOGICAL :: mesh_changed
    REAL*8 :: endpoint_time

    CALL discard_trial_linear_system()
    ! Free derived arrays before their mesh is replaced.
    IF (ALLOCATED(accepted%original_mesh)) CALL free_before_adaptivity()
    CALL accepted%restore_mesh(mesh_changed)
    CALL accepted%restore()
    endpoint_time = time%t_target
    time%t_target = time%t
    IF (mesh_changed) THEN
       CALL rebuild_mesh_solver_state(time%it)
    ELSE
       IF (magnetic_field_is_time_dependent()) CALL load_magnetic_field_Jtor(time%it)
       IF (switch%external_heating) CALL load_external_heating(physical_time=time%t*simpar%refval_time)
    ENDIF
    time%t_target = endpoint_time
    CALL update_uconv_qconv(sol%u,sol%q)
    CALL update_uiter()
    CALL reset_newton_refinements()
    ir = 1
    ir_check = 0
    CALL free_mesh_loc(Mesh_prec)
    CALL deep_copy_mesh_struct(Mesh,Mesh_prec)
    IF (ASSOCIATED(sol%u_tilde0)) DEALLOCATE(sol%u_tilde0)
    ALLOCATE(sol%u_tilde0(SIZE(sol%u_tilde)))
    sol%u_tilde0 = 0.d0
    IF (MPIvar%glob_id == 0) WRITE(6,*) 'Rejected timestep restored to accepted time [s]: ', &
         time%t*simpar%refval_time, ' mesh restored: ', mesh_changed
  END SUBROUTINE rollback_time_step

  SUBROUTINE discard_trial_linear_system()
#ifdef WITH_PASTIX
    USE solve_pastix, ONLY: terminate_mat_PASTIX
#endif
#ifdef WITH_PETSC
    USE solve_petsc, ONLY: terminate_PETSC
#endif
#ifdef WITH_PSBLAS
    USE solve_psblas, ONLY: terminate_PSBLAS
#endif
    IF (.NOT. matK%start) THEN
       SELECT CASE (lssolver%sollib)
       CASE (1)
#ifdef WITH_PASTIX
          CALL terminate_mat_PASTIX()
#endif
       CASE (2)
#ifdef WITH_PSBLAS
          CALL terminate_PSBLAS()
#endif
       CASE (3)
#ifdef WITH_PETSC
          CALL terminate_PETSC()
#endif
       END SELECT
    ENDIF
    matK%start = .TRUE.
    CALL free_mat()
  END SUBROUTINE discard_trial_linear_system

  SUBROUTINE handle_newton_result(result)
    TYPE(newton_result_t), INTENT(IN) :: result

    IF (result%status%outcome == solve_success) RETURN
    IF (result%iteration_limit_reached .AND. time%legacy_input .AND. &
        (switch%steady .OR. switch%psdtime .OR. (nb_args == 1 .AND. .NOT. switch%ME))) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) &
            'Warning: Newton iteration limit reached; retaining finite initialization/steady-state output. Residual: ', errNR
    ELSE
       CALL stop_on_solver_failure(result%status, 'Newton solve')
    ENDIF
  END SUBROUTINE handle_newton_result

  SUBROUTINE solve_nonlinear_system(result, accepted, predicted)
    TYPE(accepted_step_state_t), OPTIONAL, INTENT(INOUT) :: accepted
    TYPE(newton_result_t), INTENT(OUT) :: result
    LOGICAL, OPTIONAL, INTENT(IN) :: predicted
    LOGICAL :: use_prediction

    result = newton_result_t()
    use_prediction = .FALSE.
    IF (PRESENT(predicted)) use_prediction = predicted
    IF (.NOT. use_prediction) CALL update_uiter()
    errNR = HUGE(1.d0)
    ir = 1
    DO WHILE (ir <= numer%nrp)
       result%iterations = result%iterations+1
       CALL compute_newton_update(result%status,accepted)
       IF (result%status%outcome /= solve_success) RETURN
       CALL check_newton_convergence(result%status)
       IF (result%iterations == 1) result%first_error = errNR
       IF (result%status%outcome /= solve_success) RETURN

       IF (errNR < numer%tNR .AND. physical_iterate_valid) THEN
          CALL save_converged_newton_solution()
          RETURN
       ENDIF
       IF (newton_refinement_required(errNR,physical_iterate_valid,ir)) THEN
          IF (use_prediction) THEN
             result%status = solver_status_t(solve_retry,0,'Predicted Newton solve requires divergence recovery')
             RETURN
          ENDIF
          IF (recover_newton_divergence(accepted)) CYCLE
          result%status = solver_status_t(solve_retry, 0, 'Newton divergence or admissibility recovery exhausted')
          RETURN
       ENDIF
       IF (use_prediction .AND. newton_mesh_refinement_requested()) THEN
          result%status = solver_status_t(solve_retry,0,'Predicted Newton solve requires mesh refinement')
          RETURN
       ENDIF
       CALL update_newton_checkpoint_and_refine(result%status,accepted)
       IF (result%status%outcome /= solve_success) RETURN

       IF (MPIvar%glob_id == 0) THEN
          WRITE(6,*) '*********************************'
          WRITE(6,*) ' '
          WRITE(6,*) ' '
       ENDIF
       ir = ir+1
    ENDDO
    result%iteration_limit_reached = .TRUE.
    result%status = solver_status_t(solve_retry, 0, 'Newton iteration limit reached')
  END SUBROUTINE solve_nonlinear_system

  SUBROUTINE compute_newton_update(status, accepted)
    TYPE(accepted_step_state_t), OPTIONAL, INTENT(INOUT) :: accepted
    TYPE(solver_status_t), INTENT(OUT) :: status

    ! update nonconstant dumping factor
    CALL update_dumpnr()

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) "***** NR iteration: ", ir, "*****"
       WRITE (6, *) "NR dumping factor:  ",  numer%dumpnr
    ENDIF

    !update 1D transport stuff if used
    IF (switch%transport_1d) THEN
        CALL update_reduced_transport_profiles()
    ENDIF

    ! Compute Jacobian
    CALL HDG_computeJacobian()
    ! Set boundary conditions
    CALL hdg_BC()
    ! Compute elemental mapping
    CALL hdg_Mapping()
    ! Assembly the global matrix
    CALL hdg_Assembly()
    ! Solve linear system
#ifdef WITH_PETSC
    CALL solve_global_system(status, ir)
#else
    CALL solve_global_system(status)
#endif
    IF (status%outcome /= solve_success) RETURN
    ! Compute element-by-element solution
    CALL compute_element_solution()
    CALL check_newton_values(status)
    IF (status%outcome /= solve_success) RETURN
    IF (adapt%adaptivity .AND. restart_adapt) THEN
       CALL adaptivity(accepted=accepted)
       DEALLOCATE(uiter)
       ALLOCATE(uiter(SIZE(sol%u)))
       uiter = 0.
    ENDIF
  END SUBROUTINE compute_newton_update

  SUBROUTINE check_newton_values(status)
    TYPE(solver_status_t), INTENT(OUT) :: status

    status = solver_status_t()
    IF (.NOT. finite_solver_values(sol%u) .OR. .NOT. finite_solver_values(sol%q) .OR. &
        .NOT. finite_solver_values(sol%u_tilde)) &
         status = solver_status_t(solve_retry, 0, 'Nonfinite Newton solution or gradient')
    CALL agree_solver_status(status)
  END SUBROUTINE check_newton_values

  SUBROUTINE check_newton_convergence(status)
    TYPE(solver_status_t), INTENT(OUT) :: status

    status = solver_status_t()
    ! Compute error on oscillations, print max value of oscillation and save solution as check-point if oscillations are lower than threshold
    physical_iterate_valid = .TRUE.
    IF (adapt%adaptivity .AND. adapt%div_adapt) THEN
       physical_iterate_valid = physical_solution_is_admissible(sol%u)
       IF (.NOT. physical_iterate_valid .AND. MPIvar%glob_id .EQ. 0) &
            WRITE(*,*) 'NR iterate is not physically admissible; recovery checkpoint retained.'
    ENDIF
    IF (adapt%adaptivity) THEN
       CALL compute_error_oscillations(oscillations, min_osc, max_osc, n_osc, ir, ir_check, Mesh_prec, &
            checkpoint_allowed=physical_iterate_valid)
    ENDIF
    ! Save solution
    IF (switch%saveNR) THEN
       CALL setSolName(save_name, mesh_name, ir, .FALSE., .TRUE.)
       CALL HDF5_save_solution(save_name)
    END IF

    ! Check convergence of Newton-Raphson
    errNR = computeResidual(sol%u, uiter, 1.)
    errNR = errNR/numer%dumpnr

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (*, *)   "Error:                  ", errNR
#ifdef WITH_PETSC
       IF(lssolver%sollib .EQ. 3) THEN
          WRITE (*, *) "Relative Residue PETSc: ", matPETSC%residue
          WRITE (*,*)  "Number of Iterations:   ", matPETSC%its
          WRITE (*,*)  "Converged Reason:       ", matPETSC%convergedReason
       END IF
#endif
    ENDIF

    IF (.NOT. finite_solver_values([errNR])) &
         status = solver_status_t(solve_retry, 0, 'Nonfinite Newton residual')
    CALL agree_solver_status(status)
  END SUBROUTINE check_newton_convergence

  SUBROUTINE save_converged_newton_solution()
    ! The converged iterate is the completed solution checkpoint.
    WRITE(*,*) "Solution saved as checkpoint."
    CALL update_uconv_qconv(sol%u, sol%q)
    errNR_adapt = 1e10
    ir_adapt = 0
    ir_check = 1
    CALL free_mesh_loc(Mesh_prec)
    CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
  END SUBROUTINE save_converged_newton_solution

  LOGICAL FUNCTION recover_newton_divergence(accepted) RESULT(recovered)
    TYPE(accepted_step_state_t), OPTIONAL, INTENT(INOUT) :: accepted
    IF (.NOT. physical_iterate_valid .AND. MPIvar%glob_id .EQ. 0) THEN
       IF (errNR .LT. numer%tNR) THEN
          WRITE(*,*) 'NR residual converged, but solution is physically inadmissible; requesting refinement.'
       ELSEIF (ir .GE. numer%nrp .AND. errNR .LE. numer%div) THEN
          WRITE(*,*) 'NR iteration limit reached with physically inadmissible solution; requesting refinement.'
       ENDIF
    ENDIF
    recovered = prepare_divergence_refinement()
    IF (recovered) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (.NOT. physical_iterate_valid .AND. errNR .LE. numer%div) THEN
             WRITE(*,*) 'NR physical admissibility recovery: refining from best admissible iterate; retry ', &
                  divergence_refinements, ' of ', adapt%max_divergence_refinements
          ELSE
             WRITE(*,*) 'NR divergence: refining from best physically admissible iterate; retry ', &
                  divergence_refinements, ' of ', adapt%max_divergence_refinements
          ENDIF
       ENDIF
       ! Reload this attempt's equilibrium, without advancing the
       ! physical timestep or updating the puff/feedback controller.
       CALL adaptivity(preserve_history=.TRUE.,accepted=accepted)
       CALL update_uconv_qconv(sol%u, sol%q)
       IF (SIZE(uiter) .NE. SIZE(sol%u)) THEN
          DEALLOCATE(uiter)
          ALLOCATE(uiter(SIZE(sol%u)))
       ENDIF
       uiter = sol%u
       ir_check = 0
       CALL free_mesh_loc(Mesh_prec)
       CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
       ir = 1
       RETURN
    ENDIF
    IF (adapt%adaptivity .AND. adapt%div_adapt .AND. MPIvar%glob_id .EQ. 0) THEN
       IF (.NOT. divergence_checkpoint_available) THEN
          WRITE(*,*) 'NR refinement recovery stopped: no physically admissible checkpoint on this mesh.'
       ELSE
          WRITE(*,*) 'NR refinement recovery stopped: per-timestep refinement limit reached: ', &
               adapt%max_divergence_refinements
          WRITE(*,*) 'Timestep: ', time%it, ' NR iteration: ', ir, ' max oscillation: ', max_osc
          WRITE(*,*) 'Divergence refinements: ', divergence_refinements, ' oscillation refinements: ', &
               oscillation_refinements
       ENDIF
    ENDIF
  END FUNCTION recover_newton_divergence

  SUBROUTINE update_newton_checkpoint_and_refine(status, accepted)
    TYPE(accepted_step_state_t), OPTIONAL, INTENT(INOUT) :: accepted
    TYPE(solver_status_t), INTENT(OUT) :: status

    status = solver_status_t()
    uiter = sol%u
    !! ADAPTIVITY
    IF(update_best_newton_checkpoint(errNR,ir,physical_iterate_valid)) THEN
       ! if the NR is the lowest reached so far, then save it as best check-point
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(*,*) "Solution saved as last checkpoint."
       ENDIF

       divergence_counter_adapt = 0
    ELSEIF(errNR .GT. errNR_adapt) THEN
       divergence_counter_adapt = divergence_counter_adapt + 1
    ENDIF

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "ir_check: ", ir_check
    ENDIF

    ! Call adaptivity if one of the following conditions is respected
    IF (newton_mesh_refinement_requested()) THEN

       IF (adapt%osc_adapt .AND. max_osc .GT. adapt%osc_tol) THEN
          IF (.NOT. prepare_oscillation_refinement()) THEN
             IF (MPIvar%glob_id .EQ. 0) THEN
                WRITE(*,*) 'Oscillation refinement stopped: per-timestep refinement limit reached: ', &
                     adapt%max_oscillation_refinements
                WRITE(*,*) 'Timestep: ', time%it, ' NR iteration: ', ir, ' max oscillation: ', max_osc
                WRITE(*,*) 'Divergence refinements: ', divergence_refinements, ' oscillation refinements: ', &
                     oscillation_refinements
                FLUSH(6)
             ENDIF
             status = solver_status_t(solve_retry, 0, 'Oscillation refinement budget exhausted')
             RETURN
          ENDIF
          IF (MPIvar%glob_id .EQ. 0) WRITE(*,*) 'Oscillation refinement: retry ', &
               oscillation_refinements, ' of ', adapt%max_oscillation_refinements
       ENDIF

       CALL adaptivity(preserve_history=adapt%osc_adapt .AND. max_osc .GT. adapt%osc_tol,accepted=accepted)

       IF((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0)) THEN
          ! Update check-point solution to the one projected on the new mesh
          DEALLOCATE(sol%u_conv)
          ALLOCATE(sol%u_conv(SIZE(sol%u)))
          sol%u_conv = sol%u
          DEALLOCATE(sol%q_conv)
          ALLOCATE(sol%q_conv(SIZE(sol%q)))
          sol%q_conv = sol%q
       ENDIF

       ! update uiter to new mapped solution
       CALL update_uiter()

       IF(ir_check .NE. numer%nrp) THEN
          ir = ir_check
       ELSE
          ir = 0
       ENDIF

       CALL free_mesh_loc(Mesh_prec)
       CALL deep_copy_mesh_struct(Mesh, Mesh_prec)

    ENDIF
  END SUBROUTINE update_newton_checkpoint_and_refine

  LOGICAL FUNCTION newton_mesh_refinement_requested() RESULT(requested)
    requested = adapt%adaptivity .AND. ( &
         (adapt%osc_adapt .AND. max_osc > adapt%osc_tol) .OR. &
         (adapt%NR_adapt .AND. MOD(ir,adapt%freq_NR_adapt) == 0))
  END FUNCTION newton_mesh_refinement_requested

  SUBROUTINE exit_if_restart_complete()
    IF (time%legacy_input) THEN
       IF (.NOT. switch%ME .OR. time%it /= nts) RETURN
       IF (MPIvar%glob_id == 0) WRITE(6,*) &
            'ME final step already reached: ', time%it, ' (nts = ', nts, ')'
    ELSE
       IF (.NOT. physical_time_is_complete()) RETURN
       IF (MPIvar%glob_id == 0) WRITE(6,*) &
            'Physical final time already reached [s]: ', time%t*simpar%refval_time
    ENDIF

    ! Normal-run allocations and linear solver state are not initialized yet.
#ifdef WITH_PETSC
    IF (lssolver%sollib .EQ. 3) CALL FinalizePETSC()
#endif
    CALL MPI_finalize(IERR)
    STOP
  END SUBROUTINE exit_if_restart_complete

  SUBROUTINE prepare_time_step(finished)
    LOGICAL, INTENT(OUT) :: finished
    REAL*8, ALLOCATABLE :: frame_times(:)
    REAL*8 :: endpoint_time
    INTEGER :: status
    CHARACTER(LEN=256) :: message

    IF (time%legacy_input) THEN
       finished = it > nts
       RETURN
    ENDIF

    IF (time%mode == 'prescribed') THEN
       CALL magnetic_time_points(frame_times)
       CALL next_physical_time(endpoint_time, finished, status, message, frame_times)
    ELSE
       CALL next_physical_time(endpoint_time, finished, status, message)
    ENDIF
    IF (status /= 0) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'Cannot select physical timestep: ', TRIM(message)
       STOP 1
    ENDIF
    IF (finished) RETURN
    time%t_target = endpoint_time
    time%dt = endpoint_time-time%t
    time%t_ME = endpoint_time*simpar%refval_time
    time%dt_ME = time%dt*simpar%refval_time
  END SUBROUTINE prepare_time_step

  SUBROUTINE update_time_dependent_inputs()
    CALL update_magnetic_field_for_step()
    IF (switch%ME .OR. .NOT. time%legacy_input) CALL SetParticleSource()
    CALL update_heating_for_step()
  END SUBROUTINE update_time_dependent_inputs

  SUBROUTINE update_magnetic_field_for_step()
    IF (switch%ME .AND. time%legacy_input) THEN
       ! Initial/restart output represents the accepted frame. Sample the next
       ! endpoint and its controls only when an advance is actually attempted.
       CALL load_magnetic_field_Jtor(time%it + 1)
       IF (.NOT. ieee_is_finite(time%dt_ME) .OR. time%dt_ME <= 0.d0) THEN
          IF (MPIvar%glob_id == 0) WRITE(6,*) 'ME endpoint timestep must be positive: ', time%dt_ME
          STOP 1
       ENDIF
       IF (.NOT. ieee_is_finite(time%t_ME) .OR. &
           ABS(time%t_ME - (time%t*simpar%refval_time + time%dt_ME)) > &
           MAX(1.d-12,1.d-10*ABS(time%t_ME))) THEN
          IF (MPIvar%glob_id == 0) WRITE(6,*) &
               'ME frame time must match the timestep endpoint: ', time%t_ME, &
               time%t*simpar%refval_time + time%dt_ME
          STOP 1
       ENDIF
       time%dt = time%dt_ME/simpar%refval_time
    ELSEIF (.NOT. time%legacy_input) THEN
       IF (magnetic_field_is_time_dependent()) CALL load_magnetic_field_Jtor()
    ENDIF
  END SUBROUTINE update_magnetic_field_for_step

  SUBROUTINE update_heating_for_step()
    IF (.NOT. switch%external_heating) RETURN

    IF (.NOT. time%legacy_input) THEN
       CALL load_external_heating(physical_time=time%t_target*simpar%refval_time)
    ELSEIF (switch%ME) THEN
       CALL load_external_heating(physical_time=time%t_ME)
    ELSEIF (.NOT. switch%steady) THEN
       CALL load_external_heating(physical_time=(time%t+time%dt)*simpar%refval_time)
    ENDIF
  END SUBROUTINE update_heating_for_step

#ifdef PARALL
  SUBROUTINE domain_decomposition()
    IF(MPIvar%glob_size .GT. 1) THEN
       ! split the mesh
       
      CALL free_mesh_loc(Mesh_glob)
      CALL deep_copy_mesh_struct(Mesh, Mesh_glob)
      Mesh_glob%X = Mesh_glob%X*phys%lscale
      CALL split_mesh(MPIvar%glob_size, 3 , .FALSE.)
      CALL mesh_preprocess(ierr)

      WRITE (nid, *) MPIvar%glob_id + 1
      WRITE (npr, *) MPIvar%glob_size
      h5_filename = TRIM(ADJUSTL(mesh_name)) // '_' // TRIM(ADJUSTL(nid)) // '_' // TRIM(ADJUSTL(npr)) // '.h5'
      Mesh%X = Mesh%X*phys%lscale
      CALL HDF5_save_mesh_struct(Mesh, h5_filename)
      Mesh%X = Mesh%X/phys%lscale
      !CALL HDF5_save_mesh(h5_filename, Mesh%Ndim, mesh%Nelems, mesh%Nextfaces, mesh%Nnodes, mesh%Nnodesperelem, mesh%Nnodesperface, mesh%elemType, mesh%T, mesh%X, mesh%Tb, mesh%boundaryFlag)
       
    ELSE
       ALLOCATE(Mesh%ghostelems(Mesh%Nelems))
       ALLOCATE(Mesh%ghostfaces(Mesh%Nfaces))
       ALLOCATE(Mesh%ghostpro(Mesh%Nfaces))
       ALLOCATE(Mesh%ghostloc(Mesh%Nfaces))
       ALLOCATE(Mesh%loc2glob_el(Mesh%Nelems))
       ALLOCATE(Mesh%loc2glob_fa(Mesh%Nfaces))
       ALLOCATE(Mesh%loc2glob_nodes(Mesh%Nnodes))

       Mesh%ghostelems = 0
       Mesh%ghostfaces = 0
       Mesh%ghostpro = 0
       Mesh%ghostLoc = 0
       Mesh%loc2glob_el = [(i, i = 1, Mesh%Nelems)]
       Mesh%loc2glob_fa = [(i, i = 1, Mesh%Nfaces)]
       Mesh%loc2glob_nodes = [(i, i = 1, Mesh%Nnodes)]

       Mesh%Nel_glob = Mesh%Nelems
       Mesh%Nfa_glob = Mesh%Nfaces
       Mesh%Nno_glob = Mesh%Nnodes
       Mesh%Nintfaces_glob = Mesh%Nintfaces
       Mesh%Nextfaces_glob = Mesh%Nextfaces

       ! The initial global-mesh preprocess skips MPI boundary areas. Ghost
       ! ownership is now defined, so compute them for this one-rank partition.
       CALL computePuffArea()
       CALL computePumpArea()
       CALL computeCoreArea()
    ENDIF

  ENDSUBROUTINE domain_decomposition

  SUBROUTINE solution_decomposition()
    REAL*8, POINTER         :: u_glob(:), u_tilde_glob(:), q_glob(:)
    REAL*8, ALLOCATABLE     :: u_3d(:,:,:), u_tilde_3d(:,:,:), q_4D(:,:,:,:)
    REAL*8, ALLOCATABLE     :: u_3d_local(:,:,:), u_tilde_3d_local(:,:,:), q_4d_local(:,:,:,:)
    INTEGER                 :: i, index

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "*************************************************"
       WRITE(*,*) "           SOLUTION DECOMPOSITION                "
       WRITE(*,*) "*************************************************"
    ENDIF

    u_glob => sol%u
    u_tilde_glob => sol%u_tilde
    q_glob => sol%q

    ALLOCATE(u_3d(Mesh%Nel_glob, Mesh%Nnodesperelem, phys%neq))
    ALLOCATE(u_tilde_3d(Mesh%Nfa_glob, Mesh%Nnodesperface, phys%neq))
    ALLOCATE(q_4d(Mesh%Nel_glob, Mesh%Nnodesperelem, phys%neq, Mesh%Ndim))
    u_tilde_3d = 0.
    u_3d = 0.
    q_4D = 0.

    CALL reshape_transpose_permute(u_glob, u_3d, phys%neq, Mesh%Nel_glob, Mesh%Nnodesperelem)
    CALL reshape_transpose_permute(u_tilde_glob, u_tilde_3d, phys%neq, Mesh%Nfa_glob, Mesh%Nnodesperface)
    CALL reshape_transpose_permute_4D(q_glob, q_4D, Mesh%Ndim, phys%neq, Mesh%Nel_glob, Mesh%Nnodesperelem)

    ALLOCATE(u_3d_local(Mesh%Nelems, Mesh%Nnodesperelem, phys%neq))
    ALLOCATE(u_tilde_3d_local(Mesh%Nfaces, Mesh%Nnodesperface, phys%neq))
    ALLOCATE(q_4d_local(Mesh%Nelems, Mesh%Nnodesperelem, phys%neq, Mesh%Ndim))
    u_3d_local = 0.
    u_tilde_3d_local = 0.
    q_4d_local = 0.

    DO i = 1, Mesh%Nelems
      index = Mesh%loc2glob_el(i)
      u_3d_local(i,:,:) = u_3d(index,:,:)
      q_4d_local(i,:,:,:) = q_4D(index,:,:,:)
    ENDDO

    DO i = 1, Mesh%Nfaces
      u_tilde_3d_local(i,:,:) = u_tilde_3d(Mesh%loc2glob_fa(i),:,:)
    ENDDO

    DEALLOCATE(sol%u, sol%u_tilde, sol%q)

    ALLOCATE(sol%u(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq))
    ALLOCATE(sol%u_tilde(Mesh%Nfaces*Mesh%Nnodesperface*phys%neq))
    ALLOCATE(sol%q(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq*Mesh%Ndim))
    sol%u = 0
    sol%u_tilde = 0
    sol%q = 0

    CALL flatten_row_major(u_3d_local, sol%u, SIZE(u_3d_local,1), SIZE(u_3d_local,2), SIZE(u_3d_local,3))
    CALL flatten_row_major(u_tilde_3d_local, sol%u_tilde, SIZE(u_tilde_3d_local,1), SIZE(u_tilde_3d_local,2), SIZE(u_tilde_3d_local,3))
    CALL flatten_row_major_4D(q_4d_local, sol%q, SIZE(q_4d_local,1), SIZE(q_4d_local,2), SIZE(q_4d_local,3),SIZE(q_4d_local,4))


    DEALLOCATE(u_3d, u_tilde_3d, q_4d)
    DEALLOCATE(u_3d_local, u_tilde_3d_local, q_4d_local)
    NULLIFY(u_glob, u_tilde_glob, q_glob)

  ENDSUBROUTINE solution_decomposition
#endif

  SUBROUTINE update_reduced_transport_profiles()
    IF (.NOT. switch%transport_1d) RETURN

    CALL transport_model_1d%set_config(rho_edge=transport_model_input%rho_edge, rho_core=transport_model_input%rho_core, &
         rho_diffusion_model_max=transport_model_input%rho_diffusion_model_max, c_bohm_i=transport_model_input%c_bohm_i, &
         transport_region_policy=transport_model_input%transport_region_policy, &
         c_gyrobohm_i=transport_model_input%c_gyrobohm_i, c_bohm_e=transport_model_input%c_bohm_e, &
         c_gyrobohm_e=transport_model_input%c_gyrobohm_e, c_bohm_n=transport_model_input%c_bohm_n, &
         c_bohm_n_rho_slope=transport_model_input%c_bohm_n_rho_slope, &
         prandtl=transport_model_input%prandtl, pinch_model=transport_model_input%pinch_model, &
         pinch_equations=transport_model_input%pinch_equations, &
         c_pinch=transport_model_input%c_pinch, nu_th=transport_model_input%nu_th, &
         vpinch_const_phys=transport_model_input%vpinch_const_phys, rho_pinch_axis_width=transport_model_input%rho_pinch_axis_width, &
         rho_pinch_model_max=transport_model_input%rho_pinch_model_max, rho_pinch_edge_width=transport_model_input%rho_pinch_edge_width, &
         rho_blend_width=transport_model_input%rho_blend_width, diff_n_min_phys=transport_model_input%diff_n_min_phys, &
         diff_u_min_phys=transport_model_input%diff_u_min_phys, diff_e_min_phys=transport_model_input%diff_e_min_phys, &
         diff_ee_min_phys=transport_model_input%diff_ee_min_phys)
    CALL fs_transport%build_profiles(transport_model_1d%config%transport_region_policy)
    CALL transport_model_1d%update_from_flux_surfaces(fs_transport)
  END SUBROUTINE update_reduced_transport_profiles

  SUBROUTINE adaptivity(preserve_history, accepted)
#ifdef WITH_PASTIX
    USE solve_pastix, only: terminate_mat_PASTIX
#endif
#ifdef WITH_PETSC
    USE solve_petsc, only: terminate_PETSC, FinalizePETSC
#endif
#ifdef WITH_PSBLAS
    USE solve_psblas, only: terminate_PSBLAS
#endif
    TYPE(accepted_step_state_t), OPTIONAL, INTENT(INOUT) :: accepted
    LOGICAL, OPTIONAL, INTENT(IN) :: preserve_history
    TYPE(mesh_transfer_t) :: transfer
    INTEGER :: field_frame

    IF (PRESENT(accepted)) CALL accepted%capture_mesh()
    CALL reset_timestep_growth()
    CALL predictor_history%clear()
    field_frame = time%it

    ! Start timing
    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tps1)
       CALL system_CLOCK(timing%cks1, timing%clock_rate1)
    END IF

    IF(restart_adapt .AND. time%legacy_input) THEN
       time%t = time%t + time%dt
       time%it = time%it + 1
       time%ik = time%ik + 1
       sol%Nt = sol%Nt + 1
    ENDIF

    !! Start refining procedure
    count_adapt = count_adapt + 1

    ! Deep copy previous Mesh and reference element
    CALL free_mesh_loc(Mesh_prec)
    CALL deep_copy_mesh_struct(Mesh,Mesh_prec)
    CALL free_reference_element_pol(refElPol_prec)
    CALL deep_copy_refel_struct(refElPol,refElPol_prec)

    ! Rescale to dimensional values
    Mesh%X = Mesh%X*phys%lscale
    Mesh_prec%X = Mesh_prec%X*phys%lscale

    ! Free elmat, mat, magnetic field, Jtor, puff (to save some memory for the adaptivity)
    CALL free_before_adaptivity()
    CALL adaptively_refine_mesh(mesh_name,count_adapt,order)


#ifdef PARALL

    Mesh%X = Mesh%X/phys%lscale

    ! Save new mesh and solution
    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE (count_adapt_char, *) count_adapt
       CALL HDF5_save_mesh("./res/new_mesh_n" // TRIM(ADJUSTL(count_adapt_char)) // ".h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
       !CALL HDF5_save_solution("./res/projected_solution_n" // TRIM(ADJUSTL(count_adapt_char)))
    ENDIF

    ! Domain decomposition on the new mesh
    CALL split_mesh(MPIvar%glob_size, 3, .FALSE.)
    CALL mesh_preprocess(ierr)
    Mesh%X = Mesh%X*phys%lscale

    ! Restart communication
    CALL init_com()

    ! Ghost cells are removed from the mesh and the solution and then the result is gathered over the processes
    CALL gather_mesh(Mesh_prec,  T_glob = T_glob, X_glob = X_glob)
    CALL gather_solution(Mesh_in = Mesh_prec, Nnodesperelem = Mesh_prec%Nnodesperelem, u_in = sol%u_conv, q_in = sol%q_conv, u_glob = u_glob, q_glob = q_glob)
    ! Project the check-point solution to the new mesh
    CALL transfer%build(T_glob, X_glob, Mesh%T, Mesh%X)
    CALL transfer%apply_to_pointers(u_glob, sol%u, q_glob, sol%q)

    DEALLOCATE(X_glob,T_glob, u_glob, q_glob)
    NULLIFY(X_glob, T_glob, u_glob, q_glob)

#else
    Mesh%X = Mesh%X*phys%lscale
    Mesh%X = Mesh%X/phys%lscale

    ! Project the check-point solution to the new mesh
    CALL transfer%build(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X)
    CALL transfer%apply_to_pointers(sol%u_conv, sol%u, sol%q_conv, sol%q)

#endif
    IF(MPIvar%glob_id .EQ. 0) THEN
#ifndef PARALL
       WRITE(*,*) "Number of elements previous mesh: ", Mesh_prec%Nelems
       WRITE(*,*) "Number of elements current mesh:  ", Mesh%Nelems
#else
       WRITE(*,*) "Number of elements previous mesh: ", Mesh_prec%Nel_glob
       WRITE(*,*) "Number of elements current mesh:  ", Mesh%Nel_glob
#endif
    ENDIF

    ! Extract new trace solution
    CALL update_solution_arrays()

    ! Rescale back to adimensional values
    Mesh%X = Mesh%X/phys%lscale
    Mesh_prec%X = Mesh_prec%X/phys%lscale

    ! restart to first ever iteration
    matK%start = .TRUE.

    IF (lssolver%sollib .EQ. 1) THEN
#ifdef WITH_PASTIX
       CALL terminate_mat_PASTIX()
       ! MPI finalization
#endif
    ELSEIF (lssolver%sollib .EQ. 2) THEN
#ifdef WITH_PSBLAS
       CALL terminate_PSBLAS()
#endif
    ELSEIF (lssolver%sollib .EQ. 3) THEN
#ifdef WITH_PETSC
       CALL terminate_PETSC()
       CALL FinalizePETSC()
       CALL InitPETSC()
#endif
    ENDIF

    CALL rebuild_mesh_solver_state(field_frame)

    !! UPDATE VARIABLES
    errNR_adapt = 1e10
    CALL seed_divergence_checkpoint()
    IF(restart_adapt) THEN
       IF (time%legacy_input) THEN
          time%t = time%t - time%dt
          time%it = time%it - 1
          time%ik = time%ik - 1
          sol%Nt = sol%Nt - 1
       ENDIF
       restart_adapt = .FALSE.
    ENDIF

    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tpe1)
       CALL system_CLOCK(timing%cke1, timing%clock_rate1)
       timing%runtadapt = timing%runtadapt + (timing%cke1-timing%cks1)/REAL(timing%clock_rate1)
       timing%cputadapt = timing%cputadapt + timing%tpe1-timing%tps1
    END IF

    ! Newton remeshing transfers history; between-step callers instead reseed u0.
    IF (PRESENT(preserve_history)) CALL project_u0_newmesh(preserve_history, transfer)

  ENDSUBROUTINE adaptivity


  SUBROUTINE rebuild_mesh_solver_state(field_frame)
    INTEGER, INTENT(IN) :: field_frame

    ! reset parameters
    CALL set_parameters()

    ! Re-initialize magnetic field (the Mesh is needed)
    CALL initialize_magnetic_field()

    ! Re-load magnetic field and Jtor
    CALL load_magnetic_field_Jtor(field_frame)

    IF (switch%external_heating) THEN
       CALL initialize_external_heating()
       CALL load_external_heating(field_frame)
    ENDIF


    ! Re-Allocation and initialization of the elemental matrices
    CALL init_elmat()

    ! Re-Initialize shock capturing
    IF ((switch%shockcp .GT. 0) .OR. ((adapt%adaptivity) .AND. (adapt%shockcp_adapt .GT. 0)))  THEN
       CALL initializeShockCapturing()
    ENDIF

    ! Re-Compute first equation (definition of the gradient)
    CALL HDG_precalculatedfirstequation()

  END SUBROUTINE rebuild_mesh_solver_state

  SUBROUTINE free_before_adaptivity()
    IF(ASSOCIATED(phys%B)) DEALLOCATE(phys%B)
    IF(ASSOCIATED(phys%magnetic_flux)) DEALLOCATE(phys%magnetic_flux)
    IF(ASSOCIATED(phys%Bperturb)) DEALLOCATE(phys%Bperturb)
    IF(ASSOCIATED(phys%magnetic_psi)) DEALLOCATE(phys%magnetic_psi)
    IF(ASSOCIATED(phys%Jtor))     DEALLOCATE(phys%Jtor)
    IF(ASSOCIATED(phys%omega)) DEALLOCATE(phys%omega)
    IF(ASSOCIATED(phys%q_cyl)) DEALLOCATE(phys%q_cyl)
    IF(ASSOCIATED(phys%external_heating_ions)) DEALLOCATE(phys%external_heating_ions)
    IF(ASSOCIATED(phys%external_heating_electrons)) DEALLOCATE(phys%external_heating_electrons)

    NULLIFY(phys%B, phys%magnetic_flux, phys%Bperturb, phys%magnetic_psi, phys%Jtor, phys%omega, phys%q_cyl)
    NULLIFY(phys%external_heating_ions, phys%external_heating_electrons)
    IF (ALLOCATED(phys%diff_nn_Vol)) DEALLOCATE(phys%diff_nn_Vol)
    IF (ALLOCATED(phys%diff_nn_Fac)) DEALLOCATE(phys%diff_nn_Fac)
    IF (ALLOCATED(phys%diff_nn_Bou)) DEALLOCATE(phys%diff_nn_Bou)
    IF (ALLOCATED(phys%v_nn_Vol)) DEALLOCATE(phys%v_nn_Vol)
    IF (ALLOCATED(phys%v_nn_Fac)) DEALLOCATE(phys%v_nn_Fac)
    IF (ALLOCATED(phys%v_nn_Bou)) DEALLOCATE(phys%v_nn_Bou)
    IF (ALLOCATED(Mesh%Xg)) DEALLOCATE(Mesh%Xg)
    IF (ALLOCATED(Mesh%Xgf)) DEALLOCATE(Mesh%Xgf)
    IF (ALLOCATED(Mesh%Xgb)) DEALLOCATE(Mesh%Xgb)
    CALL free_el_mat()
    CALL free_mat()
  ENDSUBROUTINE free_before_adaptivity


  !************************************************
  ! Display results
  !************************************************
  SUBROUTINE displayResults()
    INTEGER              :: ieq
    REAL*8, ALLOCATABLE   :: uphy(:, :)
    REAL*8               :: Vmax(phys%npv), Vmin(phys%npv)
    nu = SIZE(sol%u)

    ALLOCATE (uphy(nu/phys%Neq, phys%npv))

    ! Compute physical variables
    CALL cons2phys(TRANSPOSE(RESHAPE(sol%u, (/phys%Neq, nu/phys%Neq/))), uphy)
    DO ieq = 1, phys%npv
       Vmax(ieq) = MAXVAL(uphy(:, ieq))
       Vmin(ieq) = MINVAL(uphy(:, ieq))
    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmax, phys%npv, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmin, phys%npv, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, '(" * ", 60("-"), "*")')
       WRITE (6, '(" * Time (adimensional) = ", E12.5, 27X, " *")') time%t
       WRITE (6, '(" * Dt (adimensional)        = ", E12.5, 27X, " *")') time%dt
       WRITE (6, '(" * ", 45("^"), 14X, " *")')
       WRITE (6, '(" * ", 10("_"), "      Minimum   ", 4X, "     Maximum   ", 14X, " *")')
       DO ieq = 1, phys%npv
          WRITE (6, '(" * ", A7, " -->", ES16.8, 3X, ES16.8, 13X, " *")') &
               & TRIM(phys%phyVarNam(ieq)), Vmin(ieq), Vmax(ieq)
       END DO
       WRITE (6, '(" * ", 60("-"), "*")')
       WRITE (6, '(" * Time residual  = ", 1X, 2(E16.8, 2X), 13X, " *")') sol%tres(sol%Nt)
       WRITE (6, *) '  '
       WRITE (6, *) '  '
       WRITE (6, *) '  '
    END IF
#ifdef PARALL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    DEALLOCATE (uphy)

  END SUBROUTINE displayResults

  !************************************************
  ! Set the name of the solution
  !************************************************
  SUBROUTINE setSolName(save_name, mesh_name, it, convNR, convT)
    CHARACTER(LEN=1024), INTENT(INOUT):: save_name
    CHARACTER(LEN=1024), INTENT(IN)   :: mesh_name
    INTEGER, INTENT(IN)                  :: it
    LOGICAL, INTENT(IN)                  :: convNR, convT
    CHARACTER(LEN=20)                 :: Num
    INTEGER                             :: l, i, output_step

    ! At the beginning, the save name is the mesh name..
    save_name = TRIM(ADJUSTL(mesh_name))

    ! Eliminate path info
    l = LEN(save_name)
    i = INDEX(save_name, '/', .TRUE.)
    save_name = save_name(i + 1:l)
    ! look for the P of the order of the polynomial (P4 for example)
    i = INDEX(save_name, 'P', .TRUE.)
    WRITE (Num, "(i10)") refElPol%nDeg
    save_name(i+1:i+1) = TRIM(ADJUSTL(Num))


#ifdef TOR3D
    ! Add the number of toroidal elements
    WRITE (Num, "(i10)") numer%ntor
    Num = TRIM(ADJUSTL(Num))
    save_name = TRIM(ADJUSTL(save_name))//"_Ntor"//Num

    ! Add the poloidal interpolation in the toroidal direction
    WRITE (Num, "(i10)") numer%ptor
    Num = TRIM(ADJUSTL(Num))
    save_name = TRIM(ADJUSTL(save_name))//"Ptor"//Num

#endif
#ifndef KEQUATION
    ! Diffusion
    WRITE (Num, "(E10.3)") phys%diff_n*simpar%refval_diffusion
#else
    WRITE (Num, "(E10.3)") (phys%diff_n+phys%diff_k_min)*simpar%refval_diffusion
#endif
    save_name = TRIM(ADJUSTL(save_name))//"_DPe"//TRIM(ADJUSTL(Num))
#ifdef TEMPERATURE
    WRITE (Num, "(E10.3)") phys%diff_pari
    save_name = TRIM(ADJUSTL(save_name))//"_DPai"//TRIM(ADJUSTL(Num))

    WRITE (Num, "(E10.3)") phys%diff_pare
    save_name = TRIM(ADJUSTL(save_name))//"_DPae"//TRIM(ADJUSTL(Num))
#endif
    ! Complete the save name: if not converged the NR, I put NR + the iteration number
    IF (.NOT. convNR) THEN
       WRITE (Num, "(I0.4)") it
       save_name = TRIM(save_name)//"_NR"//TRIM(Num)
    END IF

    ! Complete the save name: if not converged the time scheme, I put the iteration number
    IF (.NOT. convT .OR. (.NOT. time%legacy_input .AND. convNR)) THEN
       output_step = it
       IF (time%legacy_input) output_step = it/utils%freqsave
       WRITE (Num, "(I0.4)") output_step
       save_name = TRIM(save_name)//"_"//TRIM(Num)
    END IF

    IF (switch%decoup) THEN
       save_name = TRIM(ADJUSTL(save_name))//'_UNCP'
    ENDIF

    ! Add "Sol_"
#ifdef TOR3D
    save_name = 'Sol3D_'//TRIM(ADJUSTL(save_name))
#else
    save_name = 'Sol2D_'//TRIM(ADJUSTL(save_name))
#endif
    ! Add save_folder
    save_name = TRIM(ADJUSTL(input%save_folder))//save_name
  END SUBROUTINE setSolName

  SUBROUTINE compute_dt(errlstime)
    REAL*8, INTENT(in) :: errlstime

    IF ((errlstime*time%dt) < 1e-3) THEN
       WRITE (6, *) "******** Changing time step ***********"
       time%dt = time%dt*2.
    ENDIF
  END SUBROUTINE compute_dt


#ifdef TOR3D
  !**********************************************
  ! Definition of the toroidal discretization
  !**********************************************
  SUBROUTINE define_toroidal_discretization
    INTEGER :: i, ntorloc, itor, itorg, nnodes_toroidal
    INTEGER :: ind(numer%ptor + 1)
    REAL*8  :: tdiv(numer%ntor + 1), tel(numer%ptor + 1), htor

    ! Toroidal discretization
    htor = numer%tmax/numer%ntor
    tdiv = 0.
    DO i = 1, numer%ntor
       tdiv(i + 1) = i*htor
    END DO
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif

    nnodes_toroidal = ntorloc + 1 + (numer%ptor - 1)*ntorloc
    Mesh%Nnodes_toroidal = nnodes_toroidal
    ALLOCATE (Mesh%toroidal(nnodes_toroidal))

    DO itor = 1, ntorloc
#ifdef PARALL
       itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
       IF (itorg == numer%ntor + 1) itorg = 1
#else
       itorg = itor
#endif
       tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
       ind = (itor - 1)*numer%ptor + (/(i, i=1, numer%ptor + 1)/)
       Mesh%toroidal(ind) = tel
    END DO
  END SUBROUTINE define_toroidal_discretization
#endif

  SUBROUTINE update_solution_arrays()
    DEALLOCATE(sol%u_tilde)
    DEALLOCATE(sol%u_tilde0)
    ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
    ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
    CALL extractFaceSolution()
    sol%u_tilde0 = sol%u_tilde

    IF(restart_adapt) THEN
       DEALLOCATE(sol%u0)
       ALLOCATE (sol%u0(SIZE(sol%u), time%tis))
       sol%u0 = 0.
       sol%u0(:, 1) = sol%u
    ENDIF

    IF(.NOT. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))) THEN
       ! Update check-point solution to the one projected on the new mesh
       DEALLOCATE(sol%u_conv)
       ALLOCATE(sol%u_conv(SIZE(sol%u)))
       sol%u_conv = sol%u
       DEALLOCATE(sol%q_conv)
       ALLOCATE(sol%q_conv(SIZE(sol%q)))
       sol%q_conv = sol%q
    ENDIF

  ENDSUBROUTINE update_solution_arrays

  SUBROUTINE update_solution()
    IF (time%tis .GT. 1 .AND. it .LT. time%tis) THEN
       DO is = it, 1, -1
          sol%u0(:, is + 1) = sol%u0(:, is)
       END DO
    ELSEIF (time%tis .GT. 1) THEN
       DO is = time%tis, 2, -1
          sol%u0(:, is) = sol%u0(:, is - 1)
       END DO
    END IF
    sol%u0(:, 1) = sol%u
  ENDSUBROUTINE update_solution

  SUBROUTINE update_uiter()
    IF(SIZE(uiter) .NE. SIZE(sol%u0(:,1))) THEN
       DEALLOCATE(uiter)
       ALLOCATE(uiter(SIZE(sol%u0(:,1))))
    ENDIF
    uiter = sol%u0(:, 1)
  ENDSUBROUTINE update_uiter

  SUBROUTINE update_uiter_qiter_best(uiter_best, qiter_best, u, q)
    REAL*8, POINTER, INTENT(INOUT) :: uiter_best(:), qiter_best(:)
    REAL*8, INTENT(IN)           :: u(:), q(:)

    IF(ASSOCIATED(uiter_best)) THEN
       IF(SIZE(uiter_best) .NE. SIZE(u)) THEN
          DEALLOCATE(uiter_best)
          DEALLOCATE(qiter_best)
          ALLOCATE(uiter_best(SIZE(u)))
          ALLOCATE(qiter_best(SIZE(q)))
       ENDIF
    ELSE
       ALLOCATE(uiter_best(SIZE(u)))
       ALLOCATE(qiter_best(SIZE(q)))
    ENDIF
    ! A new best iterate must refresh the checkpoint even on the same mesh.
    uiter_best = u
    qiter_best = q
  ENDSUBROUTINE update_uiter_qiter_best

  LOGICAL FUNCTION physical_solution_is_admissible(u) RESULT(valid)
    REAL*8, INTENT(IN) :: u(:)
    INTEGER :: point, element, offset, nodes_per_element
#ifdef PARALL
    INTEGER :: ierr_valid
#endif
#ifdef TEMPERATURE
    REAL*8 :: kinetic, thermal, roundoff
#endif

    valid = .TRUE.
    nodes_per_element = Mesh%Nnodesperelem
#ifdef TOR3D
    nodes_per_element = refElTor%Nnodes3D
#endif
    DO point = 1, SIZE(u)/phys%neq
#ifdef PARALL
       ! Evaluate owned values once; stale ghost copies do not select a checkpoint.
       element = MOD((point-1)/nodes_per_element,Mesh%Nelems)+1
       IF (Mesh%ghostElems(element) .NE. 0) CYCLE
#endif
       offset = (point-1)*phys%neq
       IF (.NOT. ALL(ieee_is_finite(u(offset+1:offset+phys%neq)))) THEN
          valid = .FALSE.
          EXIT
       ENDIF
#ifdef NGAMMA
       IF (.NOT. (u(offset+1) .GT. 0.d0)) THEN
          valid = .FALSE.
          EXIT
       ENDIF
#endif
#ifdef TEMPERATURE
       IF (phys%neq .GE. 4) THEN
          ! Raw thermal energies, before cons2phys clips density/temperature.
          kinetic = 0.5d0*u(offset+2)*(u(offset+2)/u(offset+1))
          thermal = u(offset+3)-kinetic
          roundoff = 64.d0*EPSILON(1.d0)*MAX(ABS(u(offset+3)),ABS(kinetic))
          IF (.NOT. ieee_is_finite(kinetic) .OR. .NOT. (thermal .GE. -roundoff) .OR. &
              .NOT. (u(offset+4) .GE. 0.d0)) THEN
             valid = .FALSE.
             EXIT
          ENDIF
       ENDIF
#endif
#ifdef NEUTRAL
       IF (adapt%check_neutral_positivity .AND. phys%idx_rhon_eq .GT. 0 .AND. phys%idx_rhon_eq .LE. phys%neq) THEN
          IF (.NOT. (u(offset+phys%idx_rhon_eq) .GE. 0.d0)) THEN
             valid = .FALSE.
             EXIT
          ENDIF
       ENDIF
#endif
    ENDDO
#ifdef PARALL
    ! Physical admissibility is local, unlike the already global NR error.
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,valid,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr_valid)
#endif
  END FUNCTION physical_solution_is_admissible

  SUBROUTINE seed_divergence_checkpoint()
    IF (.NOT. (adapt%adaptivity .AND. adapt%div_adapt)) RETURN
    divergence_checkpoint_available = physical_solution_is_admissible(sol%u)
    physical_iterate_valid = divergence_checkpoint_available
    errNR_adapt = HUGE(1.d0)
    ir_adapt = 0
    IF (divergence_checkpoint_available) THEN
       CALL update_uiter_qiter_best(uiter_best,qiter_best,sol%u,sol%q)
    ELSEIF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) 'Recovery checkpoint unavailable: initial/projected state is not physically admissible.'
    ENDIF
  END SUBROUTINE seed_divergence_checkpoint

  LOGICAL FUNCTION update_best_newton_checkpoint(error,iteration,admissible) RESULT(updated)
    REAL*8, INTENT(IN) :: error
    INTEGER, INTENT(IN) :: iteration
    LOGICAL, INTENT(IN) :: admissible
    updated = .FALSE.
    IF (.NOT. (admissible .AND. error .LT. errNR_adapt)) RETURN
    errNR_adapt = error
    ir_adapt = iteration
    CALL update_uiter_qiter_best(uiter_best,qiter_best,sol%u,sol%q)
    IF (adapt%adaptivity .AND. adapt%div_adapt) divergence_checkpoint_available = .TRUE.
    updated = .TRUE.
  END FUNCTION update_best_newton_checkpoint

  SUBROUTINE reset_newton_refinements()
    divergence_refinements = 0
    oscillation_refinements = 0
    divergence_checkpoint_available = .FALSE.
    ! Seed this timestep's fallback; reseeding after remeshing retains the budget.
    CALL seed_divergence_checkpoint()
  ENDSUBROUTINE reset_newton_refinements

  LOGICAL FUNCTION prepare_oscillation_refinement() RESULT(prepared)
    prepared = .FALSE.
    IF (.NOT. (adapt%adaptivity .AND. adapt%osc_adapt)) RETURN
    IF (oscillation_refinements .GE. adapt%max_oscillation_refinements) RETURN
    oscillation_refinements = oscillation_refinements + 1
    prepared = .TRUE.
  END FUNCTION prepare_oscillation_refinement

  LOGICAL FUNCTION newton_refinement_required(error,admissible,iteration) RESULT(required)
    REAL*8, INTENT(IN) :: error
    LOGICAL, INTENT(IN) :: admissible
    INTEGER, INTENT(IN) :: iteration
    ! Preserve the existing divergence stop when recovery is disabled. Under
    ! active recovery, an invalid state retries at convergence or the NR limit.
    required = error .GT. numer%div
    IF (adapt%adaptivity .AND. adapt%div_adapt) THEN
       required = required .OR. (.NOT. admissible .AND. &
            (error .LT. numer%tNR .OR. iteration .GE. numer%nrp))
    ENDIF
  END FUNCTION newton_refinement_required

  LOGICAL FUNCTION prepare_divergence_refinement() RESULT(prepared)
    prepared = .FALSE.
    IF (.NOT. (adapt%adaptivity .AND. adapt%div_adapt)) RETURN
    IF (divergence_refinements .GE. adapt%max_divergence_refinements) RETURN
    IF (.NOT. divergence_checkpoint_available) RETURN

    ! The global NR error and shared retry count select the same recovery path
    ! on every rank. Initialization and adaptation maintain the local buffers.
    divergence_refinements = divergence_refinements + 1
    ! Estimate/refine from the best physically admissible iterate.
    ! Time history and feedback state stay intact for the same timestep.
    sol%u = uiter_best
    sol%q = qiter_best
    CALL update_uconv_qconv(sol%u, sol%q)
    prepared = .TRUE.
  END FUNCTION prepare_divergence_refinement


  SUBROUTINE set_parameters()
    order = refElPol%nDeg
#ifdef TOR3D
    Ndim = 3                                               ! Number of dimensions
#ifdef PARALL
    ntorloc = numer%ntor/MPIvar%ntor
#else
    ntorloc = numer%ntor
#endif
    Nel = Mesh%Nelems*ntorloc                             ! Number of 3D elements
    Np = refElPol%Nnodes2D*refElTor%Nnodes1D             ! Number of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*refElTor%Nnodes1D             ! Number of nodes in the lateral faces
    Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                                     ! Number of faces in the 2D mesh
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) + phys%Neq*refElPol%Nnodes2D*Mesh%Nelems ! Size of utilde
    ELSE
       nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
    ENDIF
#else
    nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
#endif
    Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
    nu = phys%Neq*Nel*Np
#else
    Ndim = Mesh%ndim
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfp = refElPol%Nfacenodes*Nf
    nut = Mesh%Nfaces*Mesh%Nnodesperface*phys%Neq
    nu = Mesh%Nelems*Mesh%Nnodesperelem*phys%Neq
#endif
  ENDSUBROUTINE set_parameters

  SUBROUTINE initialize_solu0_uiter_uconv()
    CALL predictor_history%clear()
    ALLOCATE(uiter(SIZE(sol%u)))
    ALLOCATE(uiter_best(SIZE(sol%u)))
    ALLOCATE(qiter_best(SIZE(sol%q)))
    ALLOCATE(sol%u0(SIZE(sol%u), time%tis))
    ALLOCATE (sol%u_conv(SIZE(sol%u)))
    ALLOCATE (sol%q_conv(SIZE(sol%q)))
    sol%u0 = 0.
    sol%u0(:, 1) = sol%u
    uiter = sol%u
    uiter_best = sol%u
    qiter_best = sol%q
    sol%u_conv = sol%u
    sol%q_conv = sol%q
  ENDSUBROUTINE initialize_solu0_uiter_uconv

  SUBROUTINE update_u0(u)
    REAL*8, INTENT(IN)      :: u(:)

    IF(SIZE(sol%u0,1) .NE. SIZE(u)) THEN
       DEALLOCATE(sol%u0)
       ALLOCATE (sol%u0(SIZE(u), time%tis))
    ENDIF
    sol%u0 = 0.
    sol%u0(:, 1) = u
  ENDSUBROUTINE update_u0

  SUBROUTINE update_uconv_qconv(u_conv, q_conv)
    REAL*8, INTENT(IN)        :: u_conv(:), q_conv(:)
    IF(SIZE(sol%u_conv) .NE. SIZE(u_conv)) THEN
       DEALLOCATE(sol%u_conv)
       ALLOCATE(sol%u_conv(SIZE(u_conv)))
    ENDIF
    IF(SIZE(sol%q_conv) .NE. SIZE(q_conv)) THEN
       DEALLOCATE(sol%q_conv)
       ALLOCATE(sol%q_conv(SIZE(q_conv)))
    ENDIF
    sol%u_conv = u_conv
    sol%q_conv = q_conv
  ENDSUBROUTINE update_uconv_qconv

  SUBROUTINE update_dumpnr()
    numer%dumpnr = numer%dumpnr_min+(numer%dumpnr_max-numer%dumpnr_min)/2.*(1+TANH((ir-numer%dumpnr_n0)/numer%dumpnr_width))
  ENDSUBROUTINE update_dumpnr

   SUBROUTINE update_delta_te(te,rho_poloidal)
      REAL*8,INTENT(IN) :: te(:), rho_poloidal(:)
      INTEGER :: i, count_0_8, count_1
      REAL*8 :: threshold, Te_0_8, Te_1

      threshold = 1.0E-3 
      count_0_8 = 0
      count_1 = 0
      Te_0_8 = 0.0
      Te_1 = 0.0

      DO i = 1, SIZE(rho_poloidal)
         IF (ABS(rho_poloidal(i) - 1.0) < threshold) THEN
            count_1 = count_1 + 1
            Te_1 = Te_1 + te(i)
         ELSE IF (ABS(rho_poloidal(i) - 0.8) < threshold) THEN
            count_0_8 = count_0_8 + 1
            Te_0_8 = Te_0_8 + te(i)
         END IF
      END DO

#ifdef PARALL
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, count_1, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, count_0_8, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, Te_1, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, Te_0_8, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      IF (count_1 > 0) THEN
         Te_1 = Te_1 / REAL(count_1)
      ELSE
         Te_1 = 0.
      END IF

      IF (count_0_8 > 0) THEN
         Te_0_8 = Te_0_8 / REAL(count_0_8)
      ELSE
         Te_0_8 = 0.
      END IF

      phys%delta_te = (Te_0_8-Te_1)/Te_1
   
ENDSUBROUTINE update_delta_te   


  SUBROUTINE initialize_solution()
    IF (nb_args .EQ. 2) THEN
       ! restart simulation: load solution from file (the name is given in argument)
       CALL HDF5_load_solution(save_name)
       ALLOCATE(sol%u_tilde0(SIZE(sol%u_tilde)))
    ELSE
       CALL init_sol()
    END IF
  ENDSUBROUTINE initialize_solution

  SUBROUTINE load_mesh()
    CHARACTER(LEN=1024) :: ordered_mesh_name, mesh_filename
    INTEGER             :: mesh_flip

    IF(switch%readMeshFromSol) THEN
       CALL HDF5_load_mesh_from_solution(save_name)
    ELSE

       ! Load the mesh file from gmsh input or h5
       
      mesh_flip = 1
      IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) mesh_flip = 0
      CALL resolve_input_filename(mesh_input_name, '.msh', mesh_filename, ierr)
      IF (ierr /= 0) THEN
         IF (MPIvar%glob_id == 0) WRITE(*,*) 'Mesh input file not found: ', TRIM(mesh_filename)
         STOP 1
      ENDIF
      CALL load_gmsh_mesh(mesh_filename, mesh_flip)
      CALL create_reference_element(refElPol, 2, verbose = 1)

      IF((switch%set_2d_order) .AND. (refElPol%nDeg .NE. switch%order_2d) ) THEN
         CALL set_order_mesh(mesh_filename, adapt%geometry_path, switch%order_2d, &
              ordered_mesh_name)
         CALL free_mesh()
         CALL load_gmsh_mesh(ordered_mesh_name, mesh_flip)
         CALL free_reference_element_pol(refElPol)
      ENDIF

      IF(switch%gmsh2h5) THEN
         IF(MPIvar%glob_id .EQ. 0) THEN
           h5_filename = TRIM(ADJUSTL(mesh_name)) // '.h5'
           CALL convert_gmsh_to_hdf5(h5_filename, SIZE(Mesh%X,2), SIZE(Mesh%T,1), SIZE(Mesh%Tb,1), SIZE(Mesh%X,1), SIZE(Mesh%T,2), refElPol%nDeg+1, 0, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
         ENDIF
      ENDIF


      IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
         Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
      END IF
      IF (adapt%shockcp_adapt .GT. 0) THEN
         gmsh_filename = mesh_filename
         CALL copy_file(gmsh_filename,"./res/temp.msh")
      ENDIF
       
    ENDIF
  ENDSUBROUTINE load_mesh

  SUBROUTINE check_input_arguments()
    INTEGER :: n, ext_length, i
    CHARACTER(LEN=4), PARAMETER :: mesh_extensions(2) = [CHARACTER(LEN=4) :: '.msh', '.h5']

    ! Check the number of input arguments
    nb_args = iargc()

    IF (nb_args .LT. 1) THEN
       PRINT *, " Error: the mesh name is needed"
       STOP
    END IF

    IF (nb_args .GT. 2) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) 'Error: mesh-to-mesh restart projection is unsupported.'
          WRITE (6, *) 'Usage: MHDG mesh [restart_solution]'
       ENDIF
       CALL MPI_finalize(IERR)
       STOP 1
    ENDIF

    CALL getarg(1, mesh_name)
    mesh_name = ADJUSTL(TRIM(mesh_name))
    mesh_input_name = mesh_name
    n = LEN_TRIM(mesh_name)
    DO i = 1, SIZE(mesh_extensions)
       ext_length = LEN_TRIM(mesh_extensions(i))
       IF (n > ext_length) THEN
          IF (mesh_name(n-ext_length+1:n) == TRIM(mesh_extensions(i))) THEN
             mesh_name = mesh_name(:n-ext_length)
             EXIT
          ENDIF
       ENDIF
    ENDDO

    IF (nb_args .GT. 1) THEN
       CALL getarg(2, save_name)
       save_name = ADJUSTL(TRIM(save_name))
       IF(MPIvar%glob_id .EQ. 0) THEN
         PRINT *, " Restart simulation with solution: ", save_name
       ENDIF
    END IF

  ENDSUBROUTINE check_input_arguments

  SUBROUTINE check_for_NaNs()
    DO i = 1, SIZE(sol%u)
       IF(ISNAN(sol%u(i))) THEN
          WRITE (6, *) "NaN detected. STOPPING."
          STOP
       END IF
    END DO
  ENDSUBROUTINE check_for_NaNs

  SUBROUTINE project_u0_newmesh(preserve_history, transfer)
    LOGICAL, OPTIONAL, INTENT(IN) :: preserve_history
    TYPE(mesh_transfer_t), OPTIONAL, INTENT(IN) :: transfer
    TYPE(mesh_transfer_t) :: local_transfer
    REAL*8, ALLOCATABLE :: previous_history(:,:)
    LOGICAL :: retain_history

    retain_history = .NOT. time%legacy_input .OR. switch%ME
    IF (PRESENT(preserve_history)) retain_history = retain_history .OR. preserve_history

    IF(retain_history .OR. .NOT. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))) THEN
       previous_history = sol%u0
    ELSE
       ALLOCATE(previous_history(SIZE(sol%u_conv),SIZE(sol%u0,2)))
       previous_history = 0.d0
       previous_history(:,1) = sol%u_conv
    ENDIF

    DEALLOCATE(sol%u0)
    ALLOCATE(sol%u0(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq,time%tis))
    sol%u0 = 0.d0

    IF (PRESENT(transfer)) THEN
       CALL project_saved_history(transfer, previous_history, retain_history)
    ELSE
       ! Standalone callers build once and reuse across existing history columns.
       Mesh%X = Mesh%X*phys%lscale
       Mesh_prec%X = Mesh_prec%X*phys%lscale
#ifdef PARALL
       CALL gather_mesh(Mesh_in=Mesh_prec, T_glob=T_glob, X_glob=X_glob)
       CALL local_transfer%build(T_glob, X_glob, Mesh%T, Mesh%X)
       DEALLOCATE(T_glob, X_glob)
       NULLIFY(T_glob, X_glob)
#else
       CALL local_transfer%build(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X)
#endif
       Mesh%X = Mesh%X/phys%lscale
       Mesh_prec%X = Mesh_prec%X/phys%lscale
       CALL project_saved_history(local_transfer, previous_history, retain_history)
    ENDIF
  END SUBROUTINE project_u0_newmesh

  SUBROUTINE project_saved_history(transfer, previous_history, retain_history)
    TYPE(mesh_transfer_t), INTENT(IN) :: transfer
    REAL*8, INTENT(IN) :: previous_history(:,:)
    LOGICAL, INTENT(IN) :: retain_history
    INTEGER :: history_id, n_history

    n_history = 1
    IF (retain_history) n_history = SIZE(previous_history,2)
    DO history_id = 1, n_history
#ifdef PARALL
       CALL gather_solution(Mesh_in=Mesh_prec, Nnodesperelem=Mesh_prec%Nnodesperelem, &
            u_in=previous_history(:,history_id), u_glob=u_glob)
       CALL transfer%apply(u_glob, sol%u0(:,history_id))
       DEALLOCATE(u_glob)
       NULLIFY(u_glob)
#else
       CALL transfer%apply(previous_history(:,history_id), sol%u0(:,history_id))
#endif
    ENDDO
  END SUBROUTINE project_saved_history

  SUBROUTINE reduce_diffusion()
    ! Update the diffusion, the elemental matrices and the solution
    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) "************************************************"
       WRITE (6, *) "Reducing diffusion: ", phys%diff_n*switch%diffred*simpar%refval_diffusion

#ifdef KEQUATION
       WRITE (6, *) "K diffusion min: ", phys%diff_k_min*switch%diffred*simpar%refval_diffusion
#endif
       WRITE (6, *) "************************************************"
    END IF
    phys%diff_n = phys%diff_n*switch%diffred
    phys%diff_u = phys%diff_u*switch%diffred
#ifdef TEMPERATURE
    phys%diff_e = phys%diff_e*switch%diffred
    phys%diff_ee = phys%diff_ee*switch%diffred
#endif
#ifdef VORTICITY
    phys%diff_vort = phys%diff_vort*switch%diffred
    phys%diff_pot = phys%diff_pot*switch%diffred
#endif
#ifdef KEQUATION
    phys%diff_k_min = phys%diff_k_min*switch%diffred
    phys%diff_k_max = phys%diff_k_max!*switch%diffred
#endif
  ENDSUBROUTINE reduce_diffusion

  SUBROUTINE print_timing_infos()
    ! Code timing
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%timing) THEN
          cputtot = 1e-8
          runttot = 1e-8
          cputtot = timing%cputpre + timing%cputjac + timing%cputbcd + timing%cputmap + timing%cputass + timing%cputglb + timing%cputsol + timing%cputadapt
          runttot = timing%runtpre + timing%runtjac + timing%runtbcd + timing%runtmap + timing%runtass + timing%runtglb + timing%runtsol + timing%runtadapt
#ifdef PARALL
          cputtot = cputtot + timing%cputcom
          runttot = runttot + timing%runtcom
#endif
          WRITE(6, *) " "
          WRITE(6, *) " "
          WRITE(6, *) " "
          WRITE(6, '(" *", 90("*"), "**")')
          WRITE(6, '(" *", 36X, "CODE TIMING ( Nthreads = ",i2,")", 26X, " *")') Nthreads
          WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
          WRITE(6, '(" *", 2X,  "Precal. matr     : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
               &timing%cputpre,timing%cputpre/cputtot*100,timing%runtpre,timing%runtpre/runttot*100,timing%cputpre/timing%runtpre/Nthreads
          WRITE(6, '(" *", 2X,  "Jacobian         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputjac,timing%cputjac/cputtot*100,timing%runtjac,timing%runtjac/runttot*100,timing%cputjac/timing%runtjac/Nthreads
          WRITE(6, '(" *", 2X,  "Mapping          : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputmap,timing%cputmap/cputtot*100,timing%runtmap,timing%runtmap/runttot*100,timing%cputmap/timing%runtmap/Nthreads
          WRITE(6, '(" *", 2X,  "Boundary cond.   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputbcd,timing%cputbcd/cputtot*100,timing%runtbcd,timing%runtbcd/runttot*100,timing%cputbcd/timing%runtbcd/Nthreads
          WRITE(6, '(" *", 2X,  "Assembly         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputass,timing%cputass/cputtot*100,timing%runtass,timing%runtass/runttot*100,timing%cputass/timing%runtass/Nthreads
          WRITE(6, '(" *", 2X,  "Solve glob. syst.: ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputglb,timing%cputglb/cputtot*100,timing%runtglb,timing%runtglb/runttot*100,timing%cputglb/timing%runtglb/Nthreads
          WRITE(6, '(" *", 2X,  "Element solution : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputsol,timing%cputsol/cputtot*100,timing%runtsol,timing%runtsol/runttot*100,timing%cputsol/timing%runtsol/Nthreads
          IF(adapt%adaptivity .EQV. .TRUE.) THEN
             WRITE(6, '(" *", 2X,  "Adaptivity       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
                  &timing%cputadapt,timing%cputadapt/cputtot*100,timing%runtadapt,timing%runtadapt/runttot*100,timing%cputadapt/timing%runtadapt/Nthreads
          ENDIF
#ifdef PARALL
          WRITE(6, '(" *", 2X,  "Communications   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputcom,timing%cputcom/cputtot*100,timing%runtcom,timing%runtcom/runttot*100,timing%cputcom/timing%runtcom/Nthreads
#endif
          WRITE(6, '(" *", 2X, "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")')   &
               cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          WRITE(6, '(" *", 90("*"), "**")')
          WRITE(6, *) " "
          WRITE(6, *) " "
          WRITE(6, *) " "
       END IF

       IF (lssolver%timing) THEN
          cputtot = timing%clstime1 + timing%clstime2 + timing%clstime3 + timing%clstime4 + timing%clstime5 + timing%clstime6
          runttot = timing%rlstime1 + timing%rlstime2 + timing%rlstime3 + timing%rlstime4 + timing%rlstime5 + timing%rlstime6
          WRITE(6, '(" *", 90("*"), "**")')
          IF (lssolver%sollib .EQ. 1) THEN
             WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PASTIX ( Nthreads = ",i2,")", 14X, " *")') Nthreads
          ELSE IF (lssolver%sollib .EQ. 2) THEN
             WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PSBLAS ( Nthreads = ",i2,")", 14X, " *")') Nthreads
          ELSE IF (lssolver%sollib .EQ. 3) THEN
             WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PETSc ( Nthreads = ",i2,")", 14X, " *")') Nthreads
          ENDIF
          WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
          IF (lssolver%sollib .EQ. 1) THEN
             WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
                  &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
             WRITE(6, '(" *", 2X,  "Check mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
             WRITE(6, '(" *", 2X,  "Anal. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
             WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
             WRITE(6, '(" *", 2X,  "LU decomp.       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
             WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime6,timing%clstime6/cputtot*100,timing%rlstime6,timing%rlstime6/runttot*100,timing%clstime6/timing%rlstime6/Nthreads
             WRITE(6, '(" *", 2X,  "Total time      : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",7X,F5.1 , 10X, " *")') &
                  &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          ELSEIF (lssolver%sollib .EQ. 2) THEN
             WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
                  &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
             WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
             WRITE(6, '(" *", 2X,  "Build prec       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
             WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
             WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
             WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') &
                  &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          ELSEIF (lssolver%sollib .EQ. 3) THEN
             WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
                  &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
             WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
             WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
             WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
                  &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
             WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') &
                  &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
          ENDIF
          WRITE (6, '(" *", 90("*"), "**")')
          WRITE (6, *) " "
          WRITE (6, *) " "
          WRITE (6, *) " "
       END IF
    END IF
  ENDSUBROUTINE print_timing_infos

  SUBROUTINE free_main()
#ifdef WITH_PASTIX
    USE solve_pastix, only: terminate_mat_PASTIX
#endif
#ifdef WITH_PETSC
    USE solve_petsc, only: terminate_PETSC, FinalizePETSC
#endif
#ifdef WITH_PSBLAS
    USE solve_psblas, only: terminate_PSBLAS
#endif

    IF(ALLOCATED(uiter)) DEALLOCATE (uiter)
    IF(ASSOCIATED(uiter_best)) DEALLOCATE(uiter_best)
    IF(ASSOCIATED(qiter_best)) DEALLOCATE(qiter_best)
    NULLIFY(uiter_best)
    NULLIFY(qiter_best)
    IF(ALLOCATED(mkelms)) DEALLOCATE (mkelms)
    IF(ALLOCATED(sol%u0)) DEALLOCATE (sol%u0)
    IF(ALLOCATED(oscillations)) DEALLOCATE(oscillations)

    CALL free_splines(splines)

    IF (lssolver%sollib .EQ. 1) THEN
#ifdef WITH_PASTIX
       CALL terminate_mat_PASTIX()
       ! MPI finalization
       CALL MPI_finalize(IERR)
#endif
    ELSEIF (lssolver%sollib .EQ. 2) THEN
#ifdef WITH_PSBLAS
       CALL terminate_PSBLAS()
#endif
    ELSEIF (lssolver%sollib .EQ. 3) THEN
#ifdef WITH_PETSC
       CALL terminate_PETSC()
       CALL FinalizePETSC()
       CALL MPI_finalize(IERR)
#endif
    ENDIF

    CALL free_all
    CALL free_mesh
#ifdef PARALL
    CALL free_mesh_loc(Mesh_glob)
#endif
    CALL free_mesh_loc(Mesh_prec)
    CALL free_reference_element
    CALL free_reference_element_pol(refElPol_prec)

  ENDSUBROUTINE free_main


ENDMODULE Main_utils
