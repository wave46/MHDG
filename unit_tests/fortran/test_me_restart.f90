PROGRAM test_me_restart
  USE initialization, ONLY: init_sim, next_physical_time, ensure_time_history_capacity
  USE in_out, ONLY: HDF5_load_solution, HDF5_load_mesh_from_solution, save_feedback_state, load_feedback_state
  USE GMSH_io_module, ONLY: load_gmsh_mesh
  USE Main_utils, ONLY: check_input_arguments, mesh_name, mesh_input_name, save_name, setSolName
  USE Main_utils, ONLY: prepare_physical_time_step, initialize_solu0_uiter_uconv
  USE globals
  USE MPI_OMP, ONLY: MPIvar, init_MPI_OMP
#ifdef PARALL
  USE Communications, ONLY: gather_nodal_values
#endif
  USE mpi, ONLY: MPI_FINALIZE
  USE interpolation, ONLY: binarysearch, interpolate
  USE external_heating, ONLY: initialize_external_heating, load_external_heating, close_heating_input
  USE magnetic_field, ONLY: equilibrium_input_filename, equilibrium_grid_dimensions, load_stored_toroidal_current
  USE magnetic_field, ONLY: open_control_series, load_puff_from_file, load_zeff, load_impurity_concentration
  USE magnetic_field, ONLY: initialize_magnetic_field, load_magnetic_field_nodes, sample_magnetic_input, close_magnetic_input
  USE magnetic_field, ONLY: magnetic_field_is_time_dependent, SetParticleSource, initialize_puff, clear_control_curves
  USE magnetic_field, ONLY: adjust_puff_feedback, adjust_recycling_feedback, adjust_impurity_concentration_feedback
  USE equilibrium_series, ONLY: equilibrium_sample_t
  USE HDF5_io_module, ONLY: HDF5_open, HDF5_close, HDF5_dataset_shape, HDF5_create, &
       HDF5_group_create, HDF5_group_open, HDF5_group_close, HDF5_real_saving, &
       HDF5_integer_saving, HDF5_array1D_saving
  USE HDF5, ONLY: HID_T, HSIZE_T, H5Lexists_f
  IMPLICIT NONE
  INTEGER :: nts
  REAL :: dt, restored_time
  CHARACTER(LEN=1000) :: checkpoint
  REAL :: knots(3), values(2,2), expected_heating
  INTEGER(HID_T) :: field_file
  INTEGER(HSIZE_T), ALLOCATABLE :: dimensions(:)
  INTEGER :: ierr, ip, jp, configured_nts
  LOGICAL :: trajectory
  TYPE(equilibrium_sample_t) :: field_sample
  REAL*8 :: current(1)
  CHARACTER(LEN=1000) :: field_filename
  CHARACTER(LEN=1024) :: output_name, expected_output, reference_mesh_name
  INTEGER :: sample, control_length
  REAL*8 :: initial_puff, control_times(3)
  CHARACTER(LEN=1000) :: control_filename
  CHARACTER(LEN=32) :: mode
  CHARACTER(LEN=1000) :: saved_field_path, saved_jtor_path
  LOGICAL :: saved_me, has_reference_current

  CALL get_command_argument(1, mode)
  IF (mode == 'physical-restart') THEN
     CALL check_physical_restart()
     STOP
  ENDIF
  IF (mode == 'time-controls') THEN
     CALL check_physical_controls()
     WRITE(*,'(A)') 'PHYSICAL_CONTROLS_PASS'
     CALL MPI_FINALIZE(ierr)
     STOP
  ENDIF
  IF (mode == 'time-sampling') THEN
     CALL check_physical_sampling()
     WRITE(*,'(A)') 'PHYSICAL_SAMPLING_PASS'
     STOP
  ENDIF
  IF (mode == 'time-clock') THEN
     CALL check_physical_clock()
     WRITE(*,'(A)') 'PHYSICAL_CLOCK_PASS'
     STOP
  ENDIF
  IF (mode == 'time-input' .OR. mode == 'time-input-init') THEN
     MPIvar%glob_id = 0
     MPIvar%glob_size = 1
     CALL read_input()
     IF (mode == 'time-input-init') THEN
        CALL init_sim(nts, dt)
        WRITE(*,'(A,2(1X,ES24.16),2(1X,I0))') 'TIME_INIT_STATE', &
             time%t*simpar%refval_time, time%dt*simpar%refval_time, sol%Nt, SIZE(sol%time)
     ENDIF
     WRITE(*,'(A,2(1X,A),4(1X,ES24.16),3(1X,I0))') 'TIME_INPUT_STATE', &
          TRIM(time%mode), TRIM(time%units), time%dt0, time%t_start, time%t_final, &
          simpar%refval_time, MERGE(1,0,time%dt_reduce_on_failure), MERGE(1,0,time%dt_grow_on_easy), time%nts
     STOP
  ENDIF
  IF (mode == 'nodal') THEN
     CALL init_MPI_OMP()
     CALL check_nodal_input()
     CALL MPI_FINALIZE(ierr)
     STOP
  ENDIF

  ! Exact control endpoints and grid corners must leave index+1 in bounds.
  knots = [0.d0, 0.5d0, 0.56d0]
  IF (binarysearch(3,knots,knots(1),1.d-12) /= 1) ERROR STOP 'first interpolation endpoint'
  IF (binarysearch(3,knots,knots(3),1.d-12) /= 2) ERROR STOP 'last interpolation endpoint'
  values = RESHAPE([1.d0,2.d0,3.d0,4.d0],[2,2])
  IF (interpolate(2,knots(1:2),2,knots(1:2),values,0.5d0,0.5d0,1.d-12) /= 4.d0) &
       ERROR STOP 'last grid corner'

  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  CALL check_input_arguments()
  CALL read_input()
  CALL init_sim(nts, dt)

  ! Actual control readers must interpolate identically with inferred lengths.
  initial_puff = phys%puff
  control_filename = 'controls'
  IF (INDEX(mesh_input_name, '.msh') > 0) control_filename = 'controls.h5'
  CALL open_control_series(control_filename, 'target_density', field_file, control_length)
  CALL HDF5_close(field_file)
  control_times = [0.d0, 0.25d0, 0.56d0]
  IF (ALLOCATED(phys%impurity_concentrations)) DEALLOCATE(phys%impurity_concentrations)
  ALLOCATE(phys%impurity_concentrations(1))
  DO sample = 1, SIZE(control_times)
     time%t_ME = control_times(sample)
     CALL load_puff_from_file(control_filename, input%puff_dimension)
     CALL load_zeff(control_filename, input%zeff_dimension)
     CALL load_impurity_concentration(control_filename, input%impurity_concentration_dimension)
     IF (ABS(phys%puff-(10.d0+20.d0*time%t_ME)) > 1.d-12 .OR. &
         ABS(phys%Zeff-(2.d0+3.d0*time%t_ME)) > 1.d-12 .OR. &
         ABS(phys%impurity_concentrations(1)-(0.1d0+0.01d0*time%t_ME)) > 1.d-12) &
          ERROR STOP 'control interpolation with inferred lengths'
  ENDDO
  phys%puff = initial_puff
  time%t_ME = 0.d0

  ! Full filenames and legacy stems must load identical physical mesh data.
  CALL load_gmsh_mesh(mesh_input_name, 0)
  CALL check_mesh()
  DEALLOCATE(Mesh%X, Mesh%T, Mesh%Tb, Mesh%boundaryFlag, Mesh%face_info)
  CALL HDF5_load_mesh_from_solution(save_name)
  CALL check_mesh()
  DEALLOCATE(Mesh%X, Mesh%T, Mesh%Tb, Mesh%boundaryFlag)
  refElPol%nDeg = 1
  CALL setSolName(output_name, mesh_name, 0, .TRUE., .FALSE.)
  reference_mesh_name = 'inputs/unit_mesh_P1'
  CALL setSolName(expected_output, reference_mesh_name, 0, .TRUE., .FALSE.)
  IF (output_name /= expected_output) ERROR STOP 'output name changed with full input filename'
  CALL check_output_names()

  ! Legacy numbered heating retains the magnetic frame schedule without timestamps.
  Mesh%Nnodes = 1
  ALLOCATE(Mesh%X(1,2))
  Mesh%X = 0.5d0/phys%lscale
  input%external_heating_from_grid = .TRUE.
  input%external_heating_path = 'heat_0000.h5'
  IF (switch%ME) input%external_heating_path = 'heat'
  configured_nts = time%nts
  time%nts = 1
  CALL initialize_external_heating()
  CALL load_external_heating()
  IF (ANY(phys%external_heating_ions /= 0.d0) .OR. ANY(phys%external_heating_electrons /= 0.d0)) &
       ERROR STOP 'initial heating frame'
  IF (switch%ME) THEN
     CALL load_external_heating(1)
     expected_heating = 7.d0*simpar%refval_time/simpar%refval_specenergydens/simpar%refval_mass
     IF (ABS(phys%external_heating_ions(1)-expected_heating) > 1.d-12*ABS(expected_heating) .OR. &
         ABS(phys%external_heating_electrons(1)-2.d0*expected_heating) > 1.d-12*ABS(expected_heating)) &
          ERROR STOP 'endpoint heating frame'
  ENDIF
  CALL close_heating_input()
  ! A new trajectory has its own grid and nonuniform timestamps, even with fixed B.
  input%external_heating_path = 'heat_trajectory.h5'
  DO sample = 1, 3
     control_times = [0.d0, 0.08d0, 0.02d0]
     CALL load_external_heating(physical_time=control_times(sample))
     expected_heating = (1.d0+87.5d0*control_times(sample))* &
          simpar%refval_time/simpar%refval_specenergydens/simpar%refval_mass
     IF (ABS(phys%external_heating_ions(1)-expected_heating) > 1.d-12*ABS(expected_heating) .OR. &
         ABS(phys%external_heating_electrons(1)-2*expected_heating) > 1.d-12*ABS(expected_heating)) &
          ERROR STOP 'independent heating interpolation and normalization'
  ENDDO
  ! Rebuild on a changed mesh, reusing dimensional cached frames without rescaling them.
  DEALLOCATE(Mesh%X,phys%external_heating_ions,phys%external_heating_electrons)
  Mesh%Nnodes = 2
  ALLOCATE(Mesh%X(2,2))
  Mesh%X(1,:) = [0.1d0,0.2d0]/phys%lscale
  Mesh%X(2,:) = [0.3d0,0.4d0]/phys%lscale
  CALL initialize_external_heating()
  CALL load_external_heating(physical_time=0.08d0)
  DO sample = 1, 2
     expected_heating = (SUM(Mesh%X(sample,:))*phys%lscale+7.d0)* &
          simpar%refval_time/simpar%refval_specenergydens/simpar%refval_mass
     IF (ABS(phys%external_heating_ions(sample)-expected_heating) > 1.d-12*ABS(expected_heating)) &
          ERROR STOP 'heating cached-frame projection on changed mesh'
  ENDDO
  CALL close_heating_input()
  DEALLOCATE(Mesh%X,phys%external_heating_ions,phys%external_heating_electrons)
  Mesh%Nnodes = 1

  ! Read the same non-square grid without relying on configured dimensions.
  field_filename = equilibrium_input_filename(input%field_path, 0)
  CALL HDF5_open(field_filename, field_file, ierr)
  IF (ierr /= 0) ERROR STOP 'open equilibrium fixture'
  CALL H5Lexists_f(field_file, 'equilibrium', trajectory, ierr)
  IF (.NOT. trajectory) THEN
    CALL equilibrium_grid_dimensions(field_file, &
         [CHARACTER(LEN=8) :: 'r2D', 'z2D', 'flux2D', 'Br2D', 'Bz2D', 'Bphi2D'], ip, jp)
    IF (ip /= 3 .OR. jp /= 4) ERROR STOP 'equilibrium grid orientation'
    CALL HDF5_dataset_shape(field_file, 'time', dimensions, ierr)
    IF (ierr /= 0) ERROR STOP 'control-vector shape query'
    IF (SIZE(dimensions) /= 1) ERROR STOP 'control-vector rank'
    IF (dimensions(1) /= 1) ERROR STOP 'control-vector length'
  ENDIF
  CALL HDF5_close(field_file)
  ALLOCATE(Mesh%X(1,2))
  Mesh%X = 0.5d0/phys%lscale
  time%dt_ME = 0.d0
  time%t_ME = 0.d0
  ! This reader probe has two frames; restore the restart's full nts afterwards.
  CALL sample_magnetic_input(0, field_sample)
  IF (ANY(SHAPE(field_sample%psi) /= [3,4]) .OR. field_sample%psi(3,4) /= 0.5d0) &
       ERROR STOP 'provider magnetic initial frame/orientation'
  CALL load_stored_toroidal_current(current, 0)
  IF (ABS(current(1)-3.5d0) > 1.d-12) ERROR STOP 'initial current interpolation'
  IF (switch%ME) THEN
     field_filename = equilibrium_input_filename(input%field_path, 1)
     time%dt_ME = 0.08d0
     time%t_ME = 0.08d0
     CALL sample_magnetic_input(1, field_sample)
     IF (ABS(field_sample%time-0.08d0) > 1.d-12 .OR. field_sample%frame_left /= 1) &
          ERROR STOP 'provider magnetic endpoint'
     CALL load_stored_toroidal_current(current, 1)
     IF (ABS(current(1)-13.5d0) > 1.d-12) ERROR STOP 'endpoint current interpolation'
  ENDIF
  CALL close_magnetic_input()
  ! Missing embedded reference current is optional only when requested by the caller.
  saved_field_path = input%field_path
  saved_jtor_path = input%jtor_path
  saved_me = switch%ME
  input%field_path = 'field_0000.h5'
  input%jtor_path = ''
  switch%ME = .FALSE.
  CALL load_stored_toroidal_current(current, available=has_reference_current)
  IF (has_reference_current .OR. ANY(current /= 0.d0)) ERROR STOP 'missing reference current availability'
  CALL close_magnetic_input()
  input%field_path = saved_field_path
  input%jtor_path = saved_jtor_path
  switch%ME = saved_me
  time%nts = configured_nts
  DEALLOCATE(Mesh%X)

  ! Dimensions of the synthetic triangle checkpoint, after the one-node probes.
  Mesh%Nnodes = 3
  Mesh%Nelems = 1
  Mesh%Nfaces = 3
  Mesh%Nnodesperface = 2
  refElPol%Nnodes2D = 3
  refElPol%Nfaces = 3
  refElPol%Nfacenodes = 2
  phys%feedback_integral_error = -11.d0
  phys%feedback_previous_error = -12.d0
  checkpoint = save_name
  CALL HDF5_load_solution(checkpoint)
  restored_time = 0.d0
  IF (sol%Nt > 0) restored_time = sol%time(sol%Nt)
  IF (time%it /= time%ik .OR. time%it /= sol%Nt) STOP 2
  WRITE (*,'(A,2(1X,I0),5(1X,ES24.16))') 'ME_RESTART_STATE', &
       time%it, sol%Nt, time%t, restored_time, phys%puff, &
       phys%feedback_integral_error, phys%feedback_previous_error

CONTAINS
  SUBROUTINE check_output_names()
    CHARACTER(LEN=1024) :: base, name, previous
    INTEGER :: saved_frequency
    LOGICAL :: saved_legacy, saved_decoupled

    saved_frequency = utils%freqsave
    saved_legacy = time%legacy_input
    saved_decoupled = switch%decoup
    switch%decoup = .FALSE.
    utils%freqsave = 5
    time%legacy_input = .TRUE.
    CALL setSolName(base, mesh_name, 0, .TRUE., .TRUE.)
    CALL setSolName(name, mesh_name, 6, .TRUE., .FALSE.)
    IF (name /= TRIM(base)//'_0001') ERROR STOP 'legacy output sequence changed'
    CALL setSolName(name, mesh_name, 50000, .TRUE., .FALSE.)
    IF (name /= TRIM(base)//'_10000') ERROR STOP 'legacy output index above 9999'

    time%legacy_input = .FALSE.
    CALL setSolName(name, mesh_name, 0, .TRUE., .FALSE.)
    IF (name /= TRIM(base)//'_0000') ERROR STOP 'physical initial output'
    CALL setSolName(previous, mesh_name, 6, .TRUE., .FALSE.)
    IF (previous /= TRIM(base)//'_0006') ERROR STOP 'physical accepted-step output index'
    CALL setSolName(name, mesh_name, 7, .TRUE., .FALSE.)
    IF (name /= TRIM(base)//'_0007' .OR. name == previous) ERROR STOP 'physical substep collision'
    CALL setSolName(previous, mesh_name, 7, .TRUE., .TRUE.)
    IF (previous /= name) ERROR STOP 'physical final output lost accepted index'
    CALL setSolName(name, mesh_name, 10000, .TRUE., .TRUE.)
    IF (name /= TRIM(base)//'_10000') ERROR STOP 'physical output index above 9999'

    CALL setSolName(name, mesh_name, 2, .FALSE., .TRUE.)
    IF (name /= TRIM(base)//'_NR0002') ERROR STOP 'Newton diagnostics gained a physical step index'

    utils%freqsave = saved_frequency
    time%legacy_input = saved_legacy
    switch%decoup = saved_decoupled
  END SUBROUTINE check_output_names

  SUBROUTINE check_physical_restart()
    TYPE(equilibrium_sample_t) :: sample
    LOGICAL :: finished, moving
    REAL*8 :: accepted_time, accepted_target

    MPIvar%glob_id = 0
    MPIvar%glob_size = 1
    CALL read_input()
    CALL init_sim(nts, dt)
    Mesh%Nnodes = 3
    Mesh%Nelems = 1
    Mesh%Nfaces = 3
    Mesh%Nnodesperface = 2
    refElPol%Nnodes2D = 3
    refElPol%Nfaces = 3
    refElPol%Nfacenodes = 2
    checkpoint = 'checkpoint.h5'
    CALL HDF5_load_solution(checkpoint)
    CALL initialize_puff()
    CALL initialize_solu0_uiter_uconv()
    IF (ANY(sol%u0(:,1) /= sol%u)) ERROR STOP 'restart BDF1 accepted state'
    IF (time%it /= time%ik .OR. sol%Nt /= time%it) ERROR STOP 'physical checkpoint counters'
    accepted_time = time%t*simpar%refval_time
    accepted_target = time%t_target*simpar%refval_time
    moving = magnetic_field_is_time_dependent()
    CALL sample_magnetic_input(sample=sample)
    CALL prepare_physical_time_step(finished)
    IF (time%t*simpar%refval_time /= accepted_time) ERROR STOP 'restart preparation advanced accepted clock'
    WRITE(*,'(A,2(1X,I0),4(1X,ES24.16),1X,L1,4(1X,ES24.16))') 'PHYSICAL_RESTART_STATE', &
         time%it, sol%Nt, accepted_time, accepted_target, time%t_target*simpar%refval_time, &
         sample%fraction, finished, phys%puff, phys%feedback_integral_error, &
         phys%feedback_previous_error, REAL(SIZE(sol%time),8)
    CALL close_magnetic_input()
  END SUBROUTINE check_physical_restart

  SUBROUTINE check_physical_controls()
    INTEGER :: selector, trial, repetition
    REAL*8 :: interval, expected_puff, expected_recycling, expected_impurity, first_result(7), result(7)

    CALL init_MPI_OMP()
    CALL read_input()
    CALL init_sim(nts, dt)
    switch%ME = .FALSE.
    switch%impurity_radiation = .TRUE.
    ALLOCATE(phys%impurity_concentrations(2))
    phys%impurity_concentrations = [0.02d0, 0.03d0]
    phys%puff = 8.d0
    phys%Zeff = 1.7d0
    switch%target_variable = 0
    input%puff_path = ''
    input%zeff_path = ''
    input%impurity_concentration_path = ''
    time%t_target = 0.25d0/simpar%refval_time
    time%t_ME = 99.d0
    CALL initialize_puff()
    CALL SetParticleSource()
    IF (phys%puff /= 8.d0 .OR. phys%Zeff /= 1.7d0 .OR. &
        ANY(phys%impurity_concentrations /= [0.02d0,0.03d0])) ERROR STOP 'empty control paths changed constants'

    input%puff_path = 'controls.h5'
    input%impurity_concentration_path = 'controls.h5'
    CALL initialize_puff()
    IF (ABS(phys%puff-15.d0) > 1.d-12 .OR. phys%Zeff /= 1.7d0 .OR. &
        ABS(phys%impurity_concentrations(1)-0.1025d0) > 1.d-12 .OR. &
        phys%impurity_concentrations(2) /= 0.03d0) ERROR STOP 'fixed-field sources or fixed Zeff'
    input%zeff_path = 'controls.h5'
    time%t_target = 0.56d0/simpar%refval_time
    CALL SetParticleSource()
    IF (ABS(phys%puff-21.2d0) > 1.d-12 .OR. ABS(phys%Zeff-3.68d0) > 1.d-12 .OR. &
        ABS(phys%impurity_concentrations(1)-0.1056d0) > 1.d-12) ERROR STOP 'scheduled control endpoint'

    CALL load_gmsh_mesh('inputs/unit_mesh_P1.msh', 0)
    refElPol%Nnodes2D = 3
    ALLOCATE(sol%u(3*phys%Neq))
    sol%u = 0.d0
#ifdef PARALL
    ALLOCATE(Mesh%ghostElems(1))
    Mesh%ghostElems = 0
#endif
    time%it = 5
    time%dt_ME = 123.d0
    input%target_density_path = 'controls.h5'
    input%target_density_xpr_path = 'controls.h5'
    input%puff_path = ''
    input%zeff_path = ''
    input%impurity_concentration_path = ''
    phys%feedback_propotional_gain = 0.1d0
    phys%feedback_integral_gain = 0.2d0
    phys%feedback_derivative_gain = 0.3d0
    phys%feedback_propotional_gain_xpr = 0.001d0
    phys%feedback_integral_gain_xpr = 0.002d0
    phys%feedback_derivative_gain_xpr = 0.003d0
    DO selector = 1, 3
       switch%target_variable = selector
       DO trial = 1, 2
          interval = 0.02d0*trial
          time%dt = interval/simpar%refval_time
          time%t = time%t_target-time%dt
          DO repetition = 1, 2
             phys%puff = 10.d0
             phys%Re = 0.5d0
             phys%impurity_concentrations(1) = 0.02d0
             phys%feedback_integral_error = 3.d0
             phys%feedback_previous_error = 0.2d0
             phys%feedback_integral_error_xpr = 4.d0
             phys%feedback_previous_error_xpr = 0.3d0
             phys%feedback_derivative_history = Feedback_derivative_history_type( &
                  0.5d0, 0.3d0, time%t*simpar%refval_time-interval, .TRUE.)
             phys%feedback_derivative_history_xpr = Feedback_derivative_history_type( &
                  0.5d0, 0.2d0, time%t*simpar%refval_time-interval, .TRUE.)
             CALL initialize_puff()
             IF (phys%puff /= 10.d0 .OR. phys%Re /= 0.5d0 .OR. &
                 phys%impurity_concentrations(1) /= 0.02d0 .OR. &
                 phys%feedback_integral_error /= 3.d0 .OR. phys%feedback_integral_error_xpr /= 4.d0 .OR. &
                 phys%feedback_previous_error /= 0.2d0 .OR. phys%feedback_previous_error_xpr /= 0.3d0 .OR. &
                 phys%feedback_derivative_history%target /= 0.5d0 .OR. &
                 phys%feedback_derivative_history_xpr%target /= 0.5d0) &
                  ERROR STOP 'startup advanced restored feedback state'
             CALL SetParticleSource()
             expected_puff = 10.d0
             expected_recycling = 0.5d0
             expected_impurity = 0.02d0
             IF (selector == 2) THEN
                expected_recycling = 0.5d0 + 0.65d0 + 0.09d0/interval
             ELSE
                expected_puff = 10.d0 + 0.65d0 + 0.09d0/interval
             ENDIF
             IF (selector == 3) expected_impurity = 0.02d0 + 0.0085d0 + 0.0006d0/interval
             IF (ABS(phys%puff-expected_puff) > 1.d-12 .OR. &
                 ABS(phys%Re-expected_recycling) > 1.d-12 .OR. &
                 ABS(phys%impurity_concentrations(1)-expected_impurity) > 1.d-12 .OR. &
                 ABS(phys%feedback_integral_error-(3.d0+0.5d0*interval)) > 1.d-12) &
                  ERROR STOP 'controller did not use the physical numerical interval'
             IF (selector == 3) THEN
                IF (ABS(phys%feedback_integral_error_xpr-(4.d0+0.5d0*interval)) > 1.d-12) &
                     ERROR STOP 'XPR physical numerical interval'
             ENDIF
             result = [phys%puff, phys%Re, phys%impurity_concentrations(1), &
                  phys%feedback_integral_error, phys%feedback_previous_error, &
                  phys%feedback_integral_error_xpr, phys%feedback_previous_error_xpr]
             IF (repetition == 1) first_result = result
             IF (ANY(result /= first_result)) ERROR STOP 'controller retry from the same accepted state'
          ENDDO
       ENDDO
    ENDDO
    CALL check_feedback_timing()
    CALL clear_control_curves()
  END SUBROUTINE check_physical_controls

  REAL*8 FUNCTION apply_feedback_test(selector, target, measured) RESULT(command)
    INTEGER, INTENT(IN) :: selector
    REAL*8, INTENT(IN) :: target, measured

    SELECT CASE (selector)
    CASE (1)
       CALL adjust_puff_feedback(target, measured)
       command = phys%puff
    CASE (2)
       CALL adjust_recycling_feedback(target, measured)
       command = phys%Re
    CASE (3)
       CALL adjust_impurity_concentration_feedback(target, measured)
       command = phys%impurity_concentrations(1)
    END SELECT
  END FUNCTION apply_feedback_test

  SUBROUTINE check_feedback_timing()
    REAL*8, PARAMETER :: intervals(9) = [0.005d0,0.004d0,0.003d0,0.002d0,0.001d0, &
         0.0005d0,0.00025d0,0.000325d0,0.0004225d0]
    REAL*8, PARAMETER :: retries(4) = [0.001d0,0.0005d0,0.00025d0,0.000325d0]
    TYPE(Feedback_derivative_history_type) :: saved_history(2)
    REAL*8 :: command, previous_command, target, measured, saved_values(6), expected(6), actual(6)
    INTEGER :: selector, step, attempt
    INTEGER(HID_T) :: file_id, simulation_group, physics_group

    simpar%refval_time = 1.d0
    time%legacy_input = .FALSE.
    phys%feedback_propotional_gain = 0.d0
    phys%feedback_integral_gain = 0.d0
    phys%feedback_derivative_gain = 0.0002d0
    phys%feedback_propotional_gain_xpr = 0.d0
    phys%feedback_integral_gain_xpr = 0.d0
    phys%feedback_derivative_gain_xpr = 0.0002d0
    DO selector = 1, 3
       time%t = 0.4d0
       time%it = 0
       phys%puff = 10.d0
       phys%Re = 0.8d0
       phys%impurity_concentrations(1) = 0.02d0
       phys%feedback_integral_error = 0.d0
       phys%feedback_integral_error_xpr = 0.d0
       IF (selector == 3) THEN
          phys%feedback_derivative_history_xpr = Feedback_derivative_history_type()
       ELSE
          phys%feedback_derivative_history = Feedback_derivative_history_type()
       ENDIF
       SELECT CASE (selector)
       CASE (1)
          command = phys%puff
       CASE (2)
          command = phys%Re
       CASE (3)
          command = phys%impurity_concentrations(1)
       END SELECT
       DO step = 1, SIZE(intervals)
          time%dt = intervals(step)
          time%t_target = time%t+time%dt
          target = 100.d0+2.d0*time%t_target
          measured = 50.d0+3.d0*time%t
          previous_command = command
          command = apply_feedback_test(selector, target, measured)
          IF (step == 1) THEN
             IF (command /= previous_command) ERROR STOP 'first feedback derivative must be zero'
          ELSE
             IF (ABS(command-previous_command+0.0002d0) > 1.d-12) &
                  ERROR STOP 'timestep change created a derivative jump'
          ENDIF
          time%t = time%t_target
          time%it = time%it+1
       ENDDO
    ENDDO

    ! Tentative evaluations are replayed from the same accepted controller state.
    saved_values = [phys%puff,phys%impurity_concentrations(1),phys%feedback_integral_error, &
         phys%feedback_integral_error_xpr,phys%feedback_previous_error,phys%feedback_previous_error_xpr]
    saved_history = [phys%feedback_derivative_history,phys%feedback_derivative_history_xpr]
    DO attempt = 1, SIZE(retries)
       phys%puff = saved_values(1)
       phys%impurity_concentrations(1) = saved_values(2)
       phys%feedback_integral_error = saved_values(3)
       phys%feedback_integral_error_xpr = saved_values(4)
       phys%feedback_previous_error = saved_values(5)
       phys%feedback_previous_error_xpr = saved_values(6)
       phys%feedback_derivative_history = saved_history(1)
       phys%feedback_derivative_history_xpr = saved_history(2)
       time%dt = retries(attempt)
       time%t_target = time%t+time%dt
       target = 100.d0+2.d0*time%t_target
       measured = 50.d0+3.d0*time%t
       command = apply_feedback_test(1, target, measured)
       command = apply_feedback_test(3, target, measured)
       IF (ABS(phys%puff-saved_values(1)+0.0002d0) > 1.d-12 .OR. &
           ABS(command-saved_values(2)+0.0002d0) > 1.d-12 .OR. &
           ABS(phys%feedback_integral_error-saved_values(3)-(target-measured)*time%dt) > 1.d-12 .OR. &
           ABS(phys%feedback_integral_error_xpr-saved_values(4)-(target-measured)*time%dt) > 1.d-12) &
            ERROR STOP 'feedback retry accumulated rejected state'
    ENDDO
    time%t = time%t_target
    time%it = time%it+1

    ! Round-trip actual feedback persistence, then compare the next update.
    switch%target_variable = 3
    CALL HDF5_create('feedback_state.h5', file_id, ierr)
    CALL HDF5_group_create('simulation_parameters', file_id, simulation_group, ierr)
    CALL HDF5_group_create('physics', simulation_group, physics_group, ierr)
    CALL HDF5_real_saving(physics_group, phys%puff, 'puff')
    CALL HDF5_real_saving(physics_group, phys%Re, 'recycling')
    CALL HDF5_integer_saving(physics_group, SIZE(phys%impurity_concentrations), 'n_impurities')
    CALL HDF5_array1D_saving(physics_group, phys%impurity_concentrations, SIZE(phys%impurity_concentrations), &
         'impurity_concentrations')
    CALL save_feedback_state(physics_group)
    CALL HDF5_group_close(physics_group, ierr)
    CALL HDF5_group_close(simulation_group, ierr)
    CALL HDF5_close(file_id)
    time%dt = 1.3d0*time%dt
    time%t_target = time%t+time%dt
    target = 100.d0+2.d0*time%t_target
    measured = 50.d0+3.d0*time%t
    command = apply_feedback_test(1, target, measured)
    command = apply_feedback_test(3, target, measured)
    expected = [phys%puff,phys%impurity_concentrations(1),phys%feedback_integral_error, &
         phys%feedback_integral_error_xpr,phys%feedback_previous_error,phys%feedback_previous_error_xpr]
    phys%puff = 0.d0
    phys%impurity_concentrations(1) = 0.d0
    phys%feedback_integral_error = 0.d0
    phys%feedback_integral_error_xpr = 0.d0
    phys%feedback_previous_error = 0.d0
    phys%feedback_previous_error_xpr = 0.d0
    phys%feedback_derivative_history = Feedback_derivative_history_type()
    phys%feedback_derivative_history_xpr = Feedback_derivative_history_type()
    CALL HDF5_open('feedback_state.h5', file_id, ierr)
    CALL HDF5_group_open(file_id, 'simulation_parameters', simulation_group, ierr)
    CALL load_feedback_state(simulation_group)
    CALL HDF5_group_close(simulation_group, ierr)
    CALL HDF5_close(file_id)
    command = apply_feedback_test(1, target, measured)
    command = apply_feedback_test(3, target, measured)
    actual = [phys%puff,phys%impurity_concentrations(1),phys%feedback_integral_error, &
         phys%feedback_integral_error_xpr,phys%feedback_previous_error,phys%feedback_previous_error_xpr]
    IF (ANY(ABS(actual-expected) > 1.d-12)) ERROR STOP 'feedback restart changed the next update'

    WRITE(*,'(A)') 'FEEDBACK_TIMING_PASS'
  END SUBROUTINE check_feedback_timing

  SUBROUTINE check_physical_sampling()
    TYPE(equilibrium_sample_t) :: sample
    CHARACTER(LEN=16) :: expected_kind
    LOGICAL :: moving
    REAL*8 :: reference_current(2), fraction, scale, expected, r, z
    INTEGER :: node

    MPIvar%glob_id = 0
    MPIvar%glob_size = 1
    CALL read_input()
    CALL get_command_argument(2, expected_kind)
    moving = expected_kind == 'moving'
    time%t = 0.02d0/simpar%refval_time
    time%t_target = 0.04d0/simpar%refval_time
    time%dt = time%t_target-time%t
    time%dt_ME = 0.02d0
    time%t_ME = 0.04d0
    time%it = 99
    time%nts = 9999
    IF (magnetic_field_is_time_dependent() .NEQV. moving) ERROR STOP 'input-frame magnetic evolution inference'
    CALL sample_magnetic_input(9999, sample)
    fraction = MERGE(0.5d0, 0.d0, moving)
    IF (ABS(sample%time-0.04d0) > 1.d-14 .OR. ABS(sample%fraction-fraction) > 1.d-14) &
         ERROR STOP 'physical sampling used numerical step/frame index'
    IF (ANY(ABS(sample%bphi-(1.d0+fraction)) > 1.d-14)) ERROR STOP 'physical magnetic interpolation'

    Mesh%Nnodes = 1
    ALLOCATE(Mesh%X(1,2))
    Mesh%X = 0.5d0/phys%lscale
    CALL load_stored_toroidal_current(reference_current(:1), 9999)
    IF (ABS(reference_current(1)-(3.5d0+10*fraction)) > 1.d-12) ERROR STOP 'physical current interpolation'

    time%t_target = 0.06d0/simpar%refval_time
    fraction = MERGE(0.75d0, 0.d0, moving)
    CALL sample_magnetic_input(9999, sample)
    IF (ANY(ABS(sample%bphi-(1.d0+fraction)) > 1.d-14)) ERROR STOP 'next physical magnetic sample'
    input%external_heating_from_grid = .TRUE.
    scale = simpar%refval_time/simpar%refval_specenergydens/simpar%refval_mass
    IF (moving) THEN
       input%external_heating_path = 'heat'
       CALL initialize_external_heating()
       CALL load_external_heating(9999)
       expected = 7.d0*0.06d0/0.08d0*scale
       IF (ABS(phys%external_heating_ions(1)-expected) > 1.d-12*ABS(expected)) &
            ERROR STOP 'numbered heating used numerical count as frame count'
       CALL close_heating_input()
       DEALLOCATE(phys%external_heating_ions, phys%external_heating_electrons)
    ENDIF
    input%external_heating_path = 'heat_trajectory.h5'
    CALL initialize_external_heating()
    CALL load_external_heating(9999)
    expected = (1.d0+87.5d0*0.06d0)*scale
    IF (ABS(phys%external_heating_ions(1)-expected) > 1.d-12*ABS(expected)) &
         ERROR STOP 'heating did not use its own clock at the physical target'

    DEALLOCATE(Mesh%X, phys%external_heating_ions, phys%external_heating_electrons)
    Mesh%Nnodes = 2
    ALLOCATE(Mesh%X(2,2))
    Mesh%X(1,:) = [0.1d0,0.2d0]/phys%lscale
    Mesh%X(2,:) = [0.3d0,0.4d0]/phys%lscale
    CALL load_stored_toroidal_current(reference_current, 9999)
    CALL initialize_external_heating()
    CALL load_external_heating(9999)
    DO node = 1, 2
       r = Mesh%X(node,1)*phys%lscale
       z = Mesh%X(node,2)*phys%lscale
       IF (ABS(reference_current(node)-(1+2*r+3*z+10*fraction)) > 1.d-12) &
            ERROR STOP 'current representation on changed nodes'
       expected = (r+z+87.5d0*0.06d0)*scale
       IF (ABS(phys%external_heating_ions(node)-expected) > 1.d-12*ABS(expected) .OR. &
           ABS(phys%external_heating_electrons(node)-2*expected) > 1.d-12*ABS(expected)) &
            ERROR STOP 'heating representation on changed nodes'
    ENDDO
    CALL close_heating_input()
    CALL close_magnetic_input()
  END SUBROUTINE check_physical_sampling

  SUBROUTINE check_physical_clock()
    REAL*8 :: target, accepted, frames(4), expected(3)
    INTEGER :: step, status
    LOGICAL :: finished
    CHARACTER(LEN=256) :: message

    simpar%refval_time = 0.01d0
    time%legacy_input = .FALSE.
    time%mode = 'fixed'
    time%t = 0.012d0/simpar%refval_time
    time%dt0 = 0.04d0/simpar%refval_time
    time%t_final = 0.105d0/simpar%refval_time
    time%it = 47
    expected = [0.052d0, 0.092d0, 0.105d0]
    DO step = 1, SIZE(expected)
       accepted = time%t
       CALL prepare_physical_time_step(finished)
       target = time%t_target
       IF (finished) ERROR STOP 'fixed physical endpoint selection'
       IF (ABS(target*simpar%refval_time-expected(step)) > 1.d-15) ERROR STOP 'fixed final interval'
       IF (ABS(time%dt-(target-accepted)) > 1.d-15 .OR. time%dt0 /= 4.d0) ERROR STOP 'fixed interval preparation'
       IF (time%t /= accepted .OR. time%it /= 47) ERROR STOP 'endpoint selection mutated accepted state'
       time%t = target
    ENDDO
    CALL next_physical_time(target, finished, status, message)
    IF (status /= 0 .OR. .NOT. finished) ERROR STOP 'physical final time stopping'

    time%mode = 'prescribed'
    frames = [0.d0, 0.02d0, 0.07d0, 0.11d0]
    time%t = 0.012d0/simpar%refval_time
    time%t_final = 0.10d0/simpar%refval_time
    expected = [0.02d0, 0.07d0, 0.10d0]
    DO step = 1, SIZE(expected)
       CALL next_physical_time(target, finished, status, message, frames)
       IF (status /= 0 .OR. finished) ERROR STOP 'prescribed physical endpoint selection'
       IF (ABS(target*simpar%refval_time-expected(step)) > 1.d-15) ERROR STOP 'nonuniform/final prescribed interval'
       time%t = target
    ENDDO
    time%t = 0.035d0/simpar%refval_time
    CALL next_physical_time(target, finished, status, message, frames)
    IF (status /= 0 .OR. ABS(target*simpar%refval_time-0.07d0) > 1.d-15) &
         ERROR STOP 'intermediate time did not target next supplied frame'
    time%t_final = 0.12d0/simpar%refval_time
    CALL next_physical_time(target, finished, status, message, frames)
    IF (status == 0) ERROR STOP 'prescribed final time exceeds input coverage'
    time%t = frames(4)/simpar%refval_time
    time%t_final = (frames(4)+EPSILON(frames(4)))/simpar%refval_time
    CALL next_physical_time(target, finished, status, message, frames)
    IF (status /= 0) ERROR STOP 'roundoff at last supplied timestamp'
    IF (.NOT. finished .AND. target /= time%t_final) ERROR STOP 'last supplied timestamp final remainder'

    CALL ensure_time_history_capacity(3)
    sol%Nt = 3
    sol%time(:3) = [1.d0,2.d0,3.d0]
    sol%tres(:3) = [4.d0,5.d0,6.d0]
    CALL ensure_time_history_capacity(33)
    CALL ensure_time_history_capacity(34)
    IF (SIZE(sol%time) < 34 .OR. SIZE(sol%tres) /= SIZE(sol%time)) ERROR STOP 'history capacity'
    IF (ANY(sol%time(:3) /= [1.d0,2.d0,3.d0]) .OR. ANY(sol%tres(:3) /= [4.d0,5.d0,6.d0])) &
         ERROR STOP 'history growth lost accepted entries'
    IF (ANY(sol%time(4:) /= 0.d0) .OR. ANY(sol%tres(4:) /= 0.d0) .OR. sol%Nt /= 3) &
         ERROR STOP 'history growth modified accepted count or unfilled entries'
  END SUBROUTINE check_physical_clock

  SUBROUTINE check_nodal_input()
    INTEGER :: i, frame, last_frame
    INTEGER, ALLOCATABLE :: node_ids(:)
    REAL*8 :: node, expected_flux, expected_field(3)

    CALL get_command_argument(2, input%field_path)
    CALL get_command_argument(3, mode)
    switch%ME = mode == 'moving'
    switch%ohmicsrc = .FALSE.
    switch%RMP = .FALSE.
    switch%Ripple = .FALSE.
    simpar%refval_charge = 1.d0
    simpar%refval_mass = 1.d0
    simpar%refval_time = 1.d0
    Mesh%Nnodes = 6
#ifdef PARALL
    Mesh%Nno_glob = 6
    IF (MPIvar%glob_size > 1) Mesh%Nnodes = 3
#endif
    ALLOCATE(Mesh%X(Mesh%Nnodes,2), node_ids(Mesh%Nnodes))
    node_ids = [(i, i=1,Mesh%Nnodes)]
#ifdef PARALL
    IF (MPIvar%glob_size > 1) THEN
       ! Deliberately reordered, overlapping partitions, including shared nodes.
       node_ids = [MPIvar%glob_id+3, MPIvar%glob_id+1, MPIvar%glob_id+2]
    ENDIF
    ALLOCATE(Mesh%loc2glob_nodes(Mesh%Nnodes))
    Mesh%loc2glob_nodes = node_ids
#endif
    Mesh%X(:,1) = 2.d0 + 0.1d0*node_ids
    Mesh%X(:,2) = 0.2d0*node_ids
    CALL initialize_magnetic_field()
    last_frame = 0
    IF (switch%ME) last_frame = 1
    DO frame = 0, last_frame
       CALL load_magnetic_field_nodes(frame)
       DO i = 1, Mesh%Nnodes
          node = REAL(node_ids(i), KIND=8)
          expected_field = [node+0.1d0, 2.d0*node+0.2d0, 3.d0*node+0.3d0] + 10.d0*frame
          expected_flux = 10.d0-node+100.d0*frame
          IF (MAXVAL(ABS(phys%B(i,:)-expected_field)) > 1.d-12 .OR. &
              ABS(phys%magnetic_flux(i)-expected_flux) > 1.d-12 .OR. &
              ABS(phys%magnetic_psi(i)-(6.d0-node)/16.d0) > 1.d-12 .OR. &
              ABS(phys%omega(i)-SQRT(SUM(expected_field**2))) > 1.d-12) &
               ERROR STOP 'nodal values differ from global mesh ordering'
       ENDDO
       IF (ABS(phys%Flux2Dmin-(4.d0+100.d0*frame)) > 1.d-12 .OR. &
           ABS(phys%Flux2Dmax-(9.d0+100.d0*frame)) > 1.d-12 .OR. &
           ABS(phys%r_axis-2.6d0) > 1.d-12 .OR. ABS(phys%z_axis-1.2d0) > 1.d-12) &
            ERROR STOP 'global nodal flux range or magnetic axis'
       IF (switch%ME) THEN
          IF (ABS(time%t_ME-0.08d0*frame) > 1.d-12 .OR. &
              ABS(time%dt_ME-0.08d0*frame) > 1.d-12) ERROR STOP 'nodal endpoint time'
       ENDIF
    ENDDO
#ifdef PARALL
    CALL check_nodal_gather()
#endif
    IF (MPIvar%glob_id == 0) WRITE(*,'(A,1X,I0)') 'GLOBAL_NODAL_INPUT_PASS', MPIvar%glob_size
  END SUBROUTINE check_nodal_input

#ifdef PARALL
  SUBROUTINE check_nodal_gather()
    TYPE(Mesh_type) :: partition
    REAL*8, ALLOCATABLE :: local_values(:)
    REAL*8, POINTER :: gathered(:)
    REAL*8 :: expected(6)
    INTEGER :: i

    expected = [(REAL(i-3, KIND=8), i=1,6)]
    partition%Nno_glob = 6
    ALLOCATE(partition%T(2,3), partition%ghostElems(2))
    IF (MPIvar%glob_size == 1) THEN
       ALLOCATE(partition%loc2glob_nodes(6))
       partition%loc2glob_nodes = [1,2,3,4,5,6]
       partition%T(1,:) = [1,2,3]
       partition%T(2,:) = [4,5,6]
       partition%ghostElems = 0
    ELSE
       ! Four overlapping partitions cover all six nodes; node 4 is ghost-only locally.
       ALLOCATE(partition%loc2glob_nodes(4))
       partition%loc2glob_nodes(:3) = MPIvar%glob_id + [1,2,3]
       partition%loc2glob_nodes(4) = 6
       IF (MPIvar%glob_id >= 2) partition%loc2glob_nodes(4) = 1
       partition%T(1,:) = [1,2,3]
       partition%T(2,:) = [2,3,4]
       partition%ghostElems = [0,1]
    ENDIF
    local_values = expected(partition%loc2glob_nodes)
    IF (MPIvar%glob_size > 1) local_values(4) = 999.d0

    CALL gather_nodal_values(partition, local_values, gathered, allgather=.TRUE.)
    IF (ANY(gathered /= expected)) ERROR STOP 'shared nodal samples were summed or ghosts counted'
    DEALLOCATE(gathered)
    CALL gather_nodal_values(partition, local_values, gathered, allgather=.FALSE.)
    IF (MPIvar%glob_id == 0) THEN
       IF (ANY(gathered /= expected)) ERROR STOP 'root-only nodal gathering differs'
    ENDIF
    DEALLOCATE(gathered, local_values, partition%T, partition%ghostElems, partition%loc2glob_nodes)
  END SUBROUTINE check_nodal_gather
#endif

  SUBROUTINE check_mesh()
    IF (Mesh%Nnodes /= 3 .OR. Mesh%Nelems /= 1) ERROR STOP 'mesh dimensions'
    IF (ABS(MINVAL(Mesh%X(:,1))*phys%lscale-2.d0) > 1.d-12 .OR. &
        ABS(MAXVAL(Mesh%X(:,1))*phys%lscale-3.d0) > 1.d-12 .OR. &
        ABS(MAXVAL(Mesh%X(:,2))*phys%lscale-1.d0) > 1.d-12) ERROR STOP 'mesh coordinates'
    IF (MINVAL(Mesh%T) /= 1 .OR. MAXVAL(Mesh%T) /= 3) ERROR STOP 'mesh connectivity'
  END SUBROUTINE check_mesh
END PROGRAM test_me_restart
