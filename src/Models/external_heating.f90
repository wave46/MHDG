MODULE external_heating
  USE prec_const, ONLY: float
  USE globals
  USE MPI_OMP, ONLY: MPIvar
  USE interpolation, ONLY: interpolate
  USE input_files, ONLY: resolve_input_filename, numbered_stem, numbered_filename
  USE heating_series, ONLY: heating_series_t, heating_sample_t
  USE magnetic_field, ONLY: magnetic_time_points
  IMPLICIT NONE
  TYPE(heating_series_t), SAVE, PRIVATE :: heating_input
  LOGICAL, SAVE, PRIVATE :: input_open = .FALSE., time_dependent = .FALSE., mesh_ready = .FALSE.
  REAL(float), SAVE, PRIVATE :: applied_time = 0
CONTAINS

  SUBROUTINE initialize_external_heating
    INTEGER :: nnodes
#ifdef TOR3D
    nnodes = Mesh%Nnodes*Mesh%Nnodes_toroidal
#else
    nnodes = Mesh%Nnodes
#endif
    ALLOCATE(phys%external_heating_ions(nnodes), phys%external_heating_electrons(nnodes))
    phys%external_heating_ions = 0
    phys%external_heating_electrons = 0
    mesh_ready = .FALSE.
  END SUBROUTINE initialize_external_heating

  SUBROUTINE load_external_heating(frame, physical_time)
    INTEGER, OPTIONAL, INTENT(IN) :: frame
    REAL(float), OPTIONAL, INTENT(IN) :: physical_time
    IF (input%external_heating_from_grid) THEN
       CALL load_external_heating_from_grid(frame, physical_time)
    ELSE
       CALL load_external_heating_nodes()
    ENDIF
  END SUBROUTINE load_external_heating

  SUBROUTINE load_external_heating_from_grid(frame, physical_time)
    INTEGER, OPTIONAL, INTENT(IN) :: frame
    REAL(float), OPTIONAL, INTENT(IN) :: physical_time
    TYPE(heating_sample_t) :: power
    REAL(float), ALLOCATABLE :: r(:), z(:), times(:)
    REAL(float) :: target_time, x, y
    INTEGER :: ierr, i, ind
#ifdef TOR3D
    INTEGER :: j
#endif
    CHARACTER(LEN=1000) :: message

    CALL open_heating_input()
    target_time = time%t*simpar%refval_time
    IF (PRESENT(physical_time)) THEN
       target_time = physical_time
    ELSEIF (PRESENT(frame) .AND. switch%ME) THEN
       CALL magnetic_time_points(times)
       IF (frame < 0 .OR. frame >= SIZE(times)) &
            CALL check_heating_input(-1, 'Requested heating endpoint is outside the magnetic schedule')
       target_time = times(frame+1)
    ENDIF
    IF (mesh_ready) THEN
       IF (.NOT. time_dependent .OR. target_time == applied_time) RETURN
    ENDIF
    CALL heating_input%sample_at(target_time, power, ierr, message)
    CALL check_heating_input(ierr, message)
    CALL heating_input%grid(r, z)
    r = r/phys%lscale
    z = z/phys%lscale
    power%ions = power%ions*simpar%refval_time/simpar%refval_specenergydens/simpar%refval_mass
    power%electrons = power%electrons*simpar%refval_time/simpar%refval_specenergydens/simpar%refval_mass
    DO i = 1, Mesh%Nnodes
       x = Mesh%X(i,1)
       y = Mesh%X(i,2)
       ind = i
#ifdef TOR3D
       DO j = 1, Mesh%Nnodes_toroidal
          ind = (j-1)*Mesh%Nnodes+i
#endif
          phys%external_heating_ions(ind) = interpolate(SIZE(z), z, SIZE(r), r, power%ions, y, x, 1.d-12)
          phys%external_heating_electrons(ind) = interpolate(SIZE(z), z, SIZE(r), r, power%electrons, y, x, 1.d-12)
#ifdef TOR3D
       ENDDO
#endif
    ENDDO
    applied_time = target_time
    mesh_ready = .TRUE.
  END SUBROUTINE load_external_heating_from_grid

  SUBROUTINE open_heating_input()
    INTEGER :: ierr, n
    REAL(float), ALLOCATABLE :: times(:)
    CHARACTER(LEN=1000) :: filename, stem, message
    LOGICAL :: numbered, full_filename

    IF (input_open) RETURN
    n = LEN_TRIM(input%external_heating_path)
    full_filename = .FALSE.
    IF (n >= 3) full_filename = input%external_heating_path(n-2:n) == '.h5'
    CALL resolve_input_filename(input%external_heating_path, '.h5', filename, ierr)
    IF (ierr /= 0 .AND. switch%ME .AND. .NOT. full_filename) THEN
       CALL numbered_stem(input%external_heating_path, stem, ierr, message)
       CALL check_heating_input(ierr, message)
       filename = numbered_filename(stem, 0)
    ENDIF
    CALL heating_input%open(filename, ierr, message)
    CALL check_heating_input(ierr, message)
    numbered = .FALSE.
    n = LEN_TRIM(filename)
    IF (switch%ME .AND. n > 8) numbered = filename(n-7:n) == '_0000.h5'
    IF (numbered .AND. .NOT. heating_input%is_trajectory()) THEN
       CALL magnetic_time_points(times)
       CALL heating_input%open_numbered(input%external_heating_path, times(:time%nts+1), ierr, message)
       CALL check_heating_input(ierr, message)
    ENDIF
    CALL heating_input%time_points(times)
    time_dependent = SIZE(times) > 1
    input_open = .TRUE.
    IF (MPIvar%glob_id == 0) WRITE(6,*) 'Heating provider opened: ', TRIM(input%external_heating_path)
  END SUBROUTINE open_heating_input

  SUBROUTINE check_heating_input(ierr, message)
    INTEGER, INTENT(IN) :: ierr
    CHARACTER(LEN=*), INTENT(IN) :: message
    IF (ierr == 0) RETURN
    IF (MPIvar%glob_id == 0) WRITE(6,*) 'Heating input error: ', TRIM(message)
    STOP 1
  END SUBROUTINE check_heating_input

  SUBROUTINE close_heating_input()
    CALL heating_input%clear()
    input_open = .FALSE.
    time_dependent = .FALSE.
    mesh_ready = .FALSE.
    applied_time = 0
  END SUBROUTINE close_heating_input

  SUBROUTINE load_external_heating_nodes
    WRITE(*,*) 'Nodal external heating is not implemented'
    STOP 1
  END SUBROUTINE load_external_heating_nodes
END MODULE external_heating
