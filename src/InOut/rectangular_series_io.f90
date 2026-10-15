MODULE rectangular_series_io
  USE prec_const, ONLY: float
  USE HDF5
  USE HDF5_io_module, ONLY: HDF5_open, HDF5_close, HDF5_dataset_shape, &
       HDF5_array2D_reading, HDF5_frame_reading, HDF5_real_vector_reading
  USE input_files, ONLY: resolve_input_filename
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_rectangular_grid, read_matrix
  PUBLIC :: locate_series_time

  TYPE, PUBLIC :: rectangular_series_t
     PRIVATE
     INTEGER(HID_T) :: file_id = -1
     LOGICAL :: trajectory = .FALSE.
     CHARACTER(LEN=32) :: label = 'Series'
     REAL(float), ALLOCATABLE :: r(:), z(:), times(:)
     INTEGER :: cache_frame(2) = 0
   CONTAINS
     PROCEDURE :: open => rectangular_open
     PROCEDURE :: inspect_component => rectangular_inspect_component
     PROCEDURE :: read_component => rectangular_read_component
     PROCEDURE :: bracket => rectangular_bracket
     PROCEDURE :: reserve_frame => rectangular_reserve_frame
     PROCEDURE :: mark_frame => rectangular_mark_frame
     PROCEDURE :: grid => rectangular_grid
     PROCEDURE :: same_grid => rectangular_same_grid
     PROCEDURE :: time_points => rectangular_time_points
     PROCEDURE :: set_times => rectangular_set_times
     PROCEDURE :: frame_count => rectangular_frame_count
     PROCEDURE :: is_trajectory => rectangular_is_trajectory
     PROCEDURE :: file => rectangular_file
     PROCEDURE :: clear => rectangular_clear
     FINAL :: rectangular_finalize
  END TYPE rectangular_series_t
CONTAINS

  SUBROUTINE rectangular_open(this, path, group, ierr, message, label)
    CLASS(rectangular_series_t), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: path, group
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label
    CHARACTER(LEN=1000) :: filename

    CALL this%clear()
    IF (PRESENT(label)) this%label = label
    message = TRIM(this%label)//' input file not found: '//TRIM(path)
    CALL resolve_input_filename(path, '.h5', filename, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot open '//TRIM(this%label)//' input: '//TRIM(filename)
    CALL HDF5_open(filename, this%file_id, ierr)
    IF (ierr /= 0) THEN
       this%file_id = -1
       RETURN
    ENDIF
    message = 'Cannot identify '//TRIM(this%label)//' input layout'
    CALL H5Lexists_f(this%file_id, group, this%trajectory, ierr)
    IF (ierr == 0) CALL read_series_coordinates(this, ierr, message)
    IF (ierr /= 0) CALL this%clear()
  END SUBROUTINE rectangular_open

  SUBROUTINE read_series_coordinates(this, ierr, message)
    CLASS(rectangular_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    IF (this%trajectory) THEN
       message = 'Cannot read '//TRIM(this%label)//' time vector'
       CALL HDF5_real_vector_reading(this%file_id, 'time', this%times, ierr)
       IF (ierr /= 0) RETURN
       message = 'Cannot read '//TRIM(this%label)//' R grid'
       CALL HDF5_real_vector_reading(this%file_id, 'grid/r', this%r, ierr)
       IF (ierr /= 0) RETURN
       message = 'Cannot read '//TRIM(this%label)//' Z grid'
       CALL HDF5_real_vector_reading(this%file_id, 'grid/z', this%z, ierr)
       IF (ierr /= 0) RETURN
    ELSE
       CALL read_rectangular_grid(this%file_id, this%r, this%z, ierr, message)
       IF (ierr /= 0) RETURN
       this%times = [0._float]
    ENDIF
    ierr = -1
    message = TRIM(this%label)//' grid axes must be finite, increasing and contain at least two points'
    IF (.NOT. valid_axis(this%r) .OR. .NOT. valid_axis(this%z)) RETURN
    message = TRIM(this%label)//' times must be finite and strictly increasing'
    IF (.NOT. valid_times(this%times)) RETURN
    ierr = 0
    message = ''
  END SUBROUTINE read_series_coordinates

  SUBROUTINE rectangular_inspect_component(this, name, exists, ierr, message)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(OUT) :: exists
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER(HSIZE_T), ALLOCATABLE :: dimensions(:)
    INTEGER(HSIZE_T) :: expected(3)
    INTEGER :: rank

    message = 'Cannot read '//TRIM(this%label)//' dataset: '//TRIM(name)
    CALL H5Lexists_f(this%file_id, name, exists, ierr)
    IF (ierr /= 0) RETURN
    IF (.NOT. exists) THEN
       message = ''
       RETURN
    ENDIF
    CALL HDF5_dataset_shape(this%file_id, name, dimensions, ierr)
    IF (ierr /= 0) RETURN
    IF (this%trajectory) THEN
       expected = INT([SIZE(this%r), SIZE(this%z), SIZE(this%times)], HSIZE_T)
       rank = 3
    ELSE
       expected(:2) = INT([SIZE(this%z), SIZE(this%r)], HSIZE_T)
       rank = 2
    ENDIF
    ierr = -1
    message = TRIM(this%label)//' dataset shape differs from grid/time: '//TRIM(name)
    IF (SIZE(dimensions) /= rank) RETURN
    IF (ANY(dimensions /= expected(:rank))) RETURN
    ierr = 0
    message = ''
  END SUBROUTINE rectangular_inspect_component

  SUBROUTINE rectangular_read_component(this, name, frame, values, ierr)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: frame
    REAL(float), ALLOCATABLE, INTENT(OUT) :: values(:,:)
    INTEGER, INTENT(OUT) :: ierr

    CALL read_grid_component(this%file_id, this%trajectory, name, frame, SIZE(this%r), SIZE(this%z), values, ierr)
  END SUBROUTINE rectangular_read_component

  SUBROUTINE rectangular_bracket(this, t, left, right, actual_time, weight, ierr, message)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    REAL(float), INTENT(IN) :: t
    INTEGER, INTENT(OUT) :: left, right, ierr
    REAL(float), INTENT(OUT) :: actual_time, weight
    CHARACTER(LEN=*), INTENT(OUT) :: message

    ierr = -1
    message = TRIM(this%label)//' series is not open'
    IF (this%file_id < 0) RETURN
    CALL locate_series_time(this%times, t, left, right, actual_time, weight, ierr, message)
  END SUBROUTINE rectangular_bracket

  SUBROUTINE rectangular_reserve_frame(this, frame, keep_frame, slot, cached)
    CLASS(rectangular_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: frame, keep_frame
    INTEGER, INTENT(OUT) :: slot
    LOGICAL, INTENT(OUT) :: cached

    cached = .TRUE.
    DO slot = 1, 2
       IF (this%cache_frame(slot) == frame) RETURN
    ENDDO
    cached = .FALSE.
    slot = 1
    IF (this%cache_frame(slot) == keep_frame) slot = 2
    this%cache_frame(slot) = 0
  END SUBROUTINE rectangular_reserve_frame

  SUBROUTINE rectangular_mark_frame(this, slot, frame)
    CLASS(rectangular_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: slot, frame
    this%cache_frame(slot) = frame
  END SUBROUTINE rectangular_mark_frame

  SUBROUTINE rectangular_grid(this, r, z)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: r(:), z(:)
    IF (ALLOCATED(this%r)) r = this%r
    IF (ALLOCATED(this%z)) z = this%z
  END SUBROUTINE rectangular_grid

  LOGICAL FUNCTION rectangular_same_grid(this, other)
    CLASS(rectangular_series_t), INTENT(IN) :: this, other
    rectangular_same_grid = .FALSE.
    IF (.NOT. ALLOCATED(this%r) .OR. .NOT. ALLOCATED(other%r)) RETURN
    IF (SIZE(this%r) /= SIZE(other%r) .OR. SIZE(this%z) /= SIZE(other%z)) RETURN
    rectangular_same_grid = ALL(this%r == other%r) .AND. ALL(this%z == other%z)
  END FUNCTION rectangular_same_grid

  SUBROUTINE rectangular_time_points(this, times)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: times(:)
    IF (ALLOCATED(this%times)) times = this%times
  END SUBROUTINE rectangular_time_points

  SUBROUTINE rectangular_set_times(this, times, ierr, message)
    CLASS(rectangular_series_t), INTENT(INOUT) :: this
    REAL(float), INTENT(IN) :: times(:)
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    ierr = -1
    message = TRIM(this%label)//' times must be finite and strictly increasing'
    IF (.NOT. valid_times(times)) RETURN
    this%times = times
    this%cache_frame = 0
    ierr = 0
    message = ''
  END SUBROUTINE rectangular_set_times

  INTEGER FUNCTION rectangular_frame_count(this)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    rectangular_frame_count = 0
    IF (ALLOCATED(this%times)) rectangular_frame_count = SIZE(this%times)
  END FUNCTION rectangular_frame_count

  LOGICAL FUNCTION rectangular_is_trajectory(this)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    rectangular_is_trajectory = this%trajectory
  END FUNCTION rectangular_is_trajectory

  FUNCTION rectangular_file(this) RESULT(file_id)
    CLASS(rectangular_series_t), INTENT(IN) :: this
    INTEGER(HID_T) :: file_id
    ! Borrowed handle for format-specific metadata; this object owns closing it.
    file_id = this%file_id
  END FUNCTION rectangular_file

  SUBROUTINE rectangular_clear(this)
    CLASS(rectangular_series_t), INTENT(INOUT) :: this
    IF (this%file_id >= 0) CALL HDF5_close(this%file_id)
    this%file_id = -1
    this%trajectory = .FALSE.
    this%label = 'Series'
    this%cache_frame = 0
    IF (ALLOCATED(this%r)) DEALLOCATE(this%r)
    IF (ALLOCATED(this%z)) DEALLOCATE(this%z)
    IF (ALLOCATED(this%times)) DEALLOCATE(this%times)
  END SUBROUTINE rectangular_clear

  SUBROUTINE rectangular_finalize(this)
    TYPE(rectangular_series_t), INTENT(INOUT) :: this
    CALL this%clear()
  END SUBROUTINE rectangular_finalize

  SUBROUTINE read_rectangular_grid(file_id, r, z, ierr, message)
    INTEGER(HID_T), INTENT(IN) :: file_id
    REAL(float), ALLOCATABLE, INTENT(OUT) :: r(:), z(:)
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    REAL(float), ALLOCATABLE :: r2d(:,:), z2d(:,:)

    message = 'Cannot read rectangular R grid'
    CALL read_matrix(file_id, 'r2D', r2d, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot read rectangular Z grid'
    CALL read_matrix(file_id, 'z2D', z2d, ierr)
    IF (ierr /= 0) RETURN
    ierr = -1
    message = 'R,Z grids have different shapes'
    IF (ANY(SHAPE(r2d) /= SHAPE(z2d))) RETURN
    r = r2d(1,:)
    z = z2d(:,1)
    message = 'Input requires a rectangular R,Z grid with finite increasing axes'
    IF (.NOT. valid_axis(r) .OR. .NOT. valid_axis(z)) RETURN
    IF (ANY(r2d /= SPREAD(r, 1, SIZE(z))) .OR. ANY(z2d /= SPREAD(z, 2, SIZE(r)))) RETURN
    ierr = 0
    message = ''
  END SUBROUTINE read_rectangular_grid

  LOGICAL FUNCTION valid_axis(axis)
    REAL(float), INTENT(IN) :: axis(:)
    valid_axis = .FALSE.
    IF (SIZE(axis) < 2) RETURN
    IF (.NOT. ALL(ieee_is_finite(axis))) RETURN
    valid_axis = ALL(axis(2:) > axis(:SIZE(axis)-1))
  END FUNCTION valid_axis

  SUBROUTINE read_matrix(file_id, dataset_name, values, ierr)
    INTEGER(HID_T), INTENT(IN) :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: dataset_name
    REAL(float), ALLOCATABLE, INTENT(OUT) :: values(:,:)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER(HSIZE_T), ALLOCATABLE :: dimensions(:)

    CALL HDF5_dataset_shape(file_id, dataset_name, dimensions, ierr)
    IF (ierr /= 0) RETURN
    ierr = -1
    IF (SIZE(dimensions) /= 2) RETURN
    IF (ANY(dimensions < 1) .OR. ANY(dimensions > INT(HUGE(1), HSIZE_T))) RETURN
    ALLOCATE(values(INT(dimensions(1)), INT(dimensions(2))))
    CALL HDF5_array2D_reading(file_id, values, dataset_name, ierr)
  END SUBROUTINE read_matrix

  LOGICAL FUNCTION valid_times(times)
    REAL(float), INTENT(IN) :: times(:)
    valid_times = .FALSE.
    IF (SIZE(times) < 1) RETURN
    IF (.NOT. ALL(ieee_is_finite(times))) RETURN
    valid_times = ALL(times(2:) > times(:SIZE(times)-1))
  END FUNCTION valid_times

  SUBROUTINE read_grid_component(file_id, trajectory, name, frame, nr, nz, values, ierr)
    INTEGER(HID_T), INTENT(IN) :: file_id
    LOGICAL, INTENT(IN) :: trajectory
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: frame, nr, nz
    REAL(float), ALLOCATABLE, INTENT(OUT) :: values(:,:)
    INTEGER, INTENT(OUT) :: ierr
    REAL(float), ALLOCATABLE :: plane(:,:)

    IF (trajectory) THEN
       ALLOCATE(plane(nr,nz))
       CALL HDF5_frame_reading(file_id, name, frame, plane, ierr)
       IF (ierr /= 0) RETURN
       ! Python (time,z,r) reverses to Fortran (r,z,time); solver planes are (z,r).
       values = TRANSPOSE(plane)
    ELSE
       CALL read_matrix(file_id, name, values, ierr)
       IF (ierr /= 0) RETURN
    ENDIF
    IF (.NOT. ALL(ieee_is_finite(values))) ierr = -1
  END SUBROUTINE read_grid_component

  SUBROUTINE locate_series_time(times, t, left, right, actual_time, weight, ierr, message)
    REAL(float), INTENT(IN) :: times(:), t
    INTEGER, INTENT(OUT) :: left, right, ierr
    REAL(float), INTENT(OUT) :: actual_time, weight
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: middle, n
    REAL(float) :: tolerance

    ierr = -1
    message = 'Requested series time must be finite'
    IF (.NOT. ieee_is_finite(t)) RETURN
    n = SIZE(times)
    left = 1
    right = 1
    actual_time = t
    weight = 0
    IF (n > 1) THEN
       tolerance = 64*EPSILON(t)*MAX(1._float, ABS(t), ABS(times(1)), ABS(times(n)))
       message = 'Requested time is outside trajectory coverage'
       IF (t < times(1)-tolerance .OR. t > times(n)+tolerance) RETURN
       actual_time = MAX(times(1), MIN(t, times(n)))
       right = n
       DO WHILE (right-left > 1)
          middle = (left+right)/2
          IF (times(middle) <= actual_time) THEN
             left = middle
          ELSE
             right = middle
          ENDIF
       ENDDO
       IF (ABS(actual_time-times(left)) <= tolerance) THEN
          right = left
       ELSEIF (ABS(actual_time-times(right)) <= tolerance) THEN
          left = right
       ELSE
          weight = (actual_time-times(left))/(times(right)-times(left))
       ENDIF
    ENDIF
    ierr = 0
    message = ''
  END SUBROUTINE locate_series_time

END MODULE rectangular_series_io
