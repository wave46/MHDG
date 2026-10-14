PROGRAM test_me_restart
  USE initialization, ONLY: init_sim
  USE in_out, ONLY: HDF5_load_solution, HDF5_load_mesh_from_solution
  USE GMSH_io_module, ONLY: load_gmsh_mesh
  USE Main_utils, ONLY: check_input_arguments, mesh_name, mesh_input_name, save_name, setSolName
  USE globals
  USE MPI_OMP, ONLY: MPIvar, init_MPI_OMP
  USE mpi, ONLY: MPI_FINALIZE
  USE interpolation, ONLY: binarysearch, interpolate
  USE external_heating, ONLY: initialize_external_heating, load_external_heating
  USE magnetic_field, ONLY: equilibrium_input_filename, equilibrium_grid_dimensions, load_stored_toroidal_current
  USE magnetic_field, ONLY: open_control_series, load_puff_from_file, load_zeff, load_impurity_concentration
  USE magnetic_field, ONLY: initialize_magnetic_field, load_magnetic_field_nodes, sample_magnetic_input, close_magnetic_input
  USE equilibrium_series, ONLY: equilibrium_sample_t
  USE HDF5_io_module, ONLY: HDF5_open, HDF5_close, HDF5_dataset_shape
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

  CALL get_command_argument(1, mode)
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

  ! The legacy heating reader follows the same initial/endpoint frame selection.
  Mesh%Nnodes = 1
  ALLOCATE(Mesh%X(1,2))
  Mesh%X = 0.5d0/phys%lscale
  input%external_heating_from_grid = .TRUE.
  input%external_heating_path = 'heat_0000.h5'
  IF (switch%ME) input%external_heating_path = 'heat'
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
  DEALLOCATE(Mesh%X,phys%external_heating_ions,phys%external_heating_electrons)

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
  configured_nts = time%nts
  time%nts = 1
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
  SUBROUTINE check_nodal_input()
    INTEGER :: i, frame, last_frame
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
    Mesh%Nno_glob = 6
    Mesh%Nnodes = 6
    IF (MPIvar%glob_size > 1) Mesh%Nnodes = 3
    ALLOCATE(Mesh%X(Mesh%Nnodes,2), Mesh%loc2glob_nodes(Mesh%Nnodes))
    IF (MPIvar%glob_size == 1) THEN
       Mesh%loc2glob_nodes = [1,2,3,4,5,6]
    ELSE
       ! Deliberately reordered, overlapping partitions, including shared nodes.
       Mesh%loc2glob_nodes = [MPIvar%glob_id+3, MPIvar%glob_id+1, MPIvar%glob_id+2]
    ENDIF
    Mesh%X(:,1) = 2.d0 + 0.1d0*Mesh%loc2glob_nodes
    Mesh%X(:,2) = 0.2d0*Mesh%loc2glob_nodes
    CALL initialize_magnetic_field()
    last_frame = 0
    IF (switch%ME) last_frame = 1
    DO frame = 0, last_frame
       CALL load_magnetic_field_nodes(frame)
       DO i = 1, Mesh%Nnodes
          node = REAL(Mesh%loc2glob_nodes(i), KIND=8)
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
    IF (MPIvar%glob_id == 0) WRITE(*,'(A,1X,I0)') 'GLOBAL_NODAL_INPUT_PASS', MPIvar%glob_size
  END SUBROUTINE check_nodal_input

  SUBROUTINE check_mesh()
    IF (Mesh%Nnodes /= 3 .OR. Mesh%Nelems /= 1) ERROR STOP 'mesh dimensions'
    IF (ABS(MINVAL(Mesh%X(:,1))*phys%lscale-2.d0) > 1.d-12 .OR. &
        ABS(MAXVAL(Mesh%X(:,1))*phys%lscale-3.d0) > 1.d-12 .OR. &
        ABS(MAXVAL(Mesh%X(:,2))*phys%lscale-1.d0) > 1.d-12) ERROR STOP 'mesh coordinates'
    IF (MINVAL(Mesh%T) /= 1 .OR. MAXVAL(Mesh%T) /= 3) ERROR STOP 'mesh connectivity'
  END SUBROUTINE check_mesh
END PROGRAM test_me_restart
