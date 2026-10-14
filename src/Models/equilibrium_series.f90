MODULE equilibrium_series
  USE prec_const, ONLY: float
  USE HDF5
  USE HDF5_io_module, ONLY: HDF5_open, HDF5_close, HDF5_dataset_shape, &
       HDF5_real_vector_reading, HDF5_real_reading, HDF5_array2D_reading, HDF5_frame_reading
  USE input_files, ONLY: resolve_input_filename
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: equilibrium_sample_t
     REAL(float), ALLOCATABLE :: psi(:,:), bphi(:,:), br(:,:), bz(:,:), jtor(:,:)
     REAL(float) :: psi_sep = 0, time = 0, fraction = 0
     INTEGER :: frame_left = 0, frame_right = 0 ! Zero-based input positions.
  END TYPE equilibrium_sample_t

  TYPE, PUBLIC :: equilibrium_series_t
     PRIVATE
     INTEGER(HID_T) :: file_id = -1
     LOGICAL :: trajectory = .FALSE., has_poloidal = .FALSE., has_jtor = .FALSE.
     REAL(float), ALLOCATABLE :: r(:), z(:), times(:), psi_sep(:)
     CHARACTER(LEN=32) :: datasets(5)
     TYPE(equilibrium_sample_t) :: cache(2)
     INTEGER :: cache_frame(2) = 0
   CONTAINS
     PROCEDURE :: open => series_open
     PROCEDURE :: sample_at => series_sample_at
     PROCEDURE :: grid => series_grid
     PROCEDURE :: time_points => series_time_points
     PROCEDURE :: clear => series_clear
     FINAL :: series_finalize
  END TYPE equilibrium_series_t

CONTAINS

  SUBROUTINE series_open(this, path, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: path
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    CHARACTER(LEN=1000) :: filename

    CALL this%clear()
    message = 'Equilibrium input file not found: '//TRIM(path)
    CALL resolve_input_filename(path, '.h5', filename, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot open equilibrium input: '//TRIM(filename)
    CALL HDF5_open(filename, this%file_id, ierr)
    IF (ierr /= 0) THEN
       this%file_id = -1
       RETURN
    ENDIF

    CALL read_equilibrium_metadata(this, ierr, message)
    IF (ierr == 0) CALL validate_coordinates(this, ierr, message)
    IF (ierr == 0) CALL validate_field_datasets(this, ierr, message)
    IF (ierr /= 0) CALL this%clear()
  END SUBROUTINE series_open

  SUBROUTINE read_equilibrium_metadata(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    message = 'Cannot identify equilibrium input layout'
    CALL H5Lexists_f(this%file_id, 'equilibrium', this%trajectory, ierr)
    IF (ierr /= 0) RETURN
    IF (this%trajectory) THEN
       CALL read_trajectory_metadata(this, ierr, message)
    ELSE
       CALL read_legacy_metadata(this, ierr, message)
    ENDIF
  END SUBROUTINE read_equilibrium_metadata

  SUBROUTINE read_trajectory_metadata(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    this%datasets = [CHARACTER(LEN=32) :: 'equilibrium/psi', 'equilibrium/Bphi', &
         'equilibrium/Br', 'equilibrium/Bz', 'equilibrium/Jtor']
    message = 'Cannot read trajectory time vector'
    CALL HDF5_real_vector_reading(this%file_id, 'time', this%times, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot read trajectory R grid'
    CALL HDF5_real_vector_reading(this%file_id, 'grid/r', this%r, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot read trajectory Z grid'
    CALL HDF5_real_vector_reading(this%file_id, 'grid/z', this%z, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot read trajectory psiSep vector'
    CALL HDF5_real_vector_reading(this%file_id, 'equilibrium/psiSep', this%psi_sep, ierr)
  END SUBROUTINE read_trajectory_metadata

  SUBROUTINE read_legacy_metadata(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    REAL(float), ALLOCATABLE :: r2d(:,:), z2d(:,:)

    this%datasets = [CHARACTER(LEN=32) :: 'flux2D', 'Bphi2D', 'Br2D', 'Bz2D', 'Jtor']
    message = 'Cannot read stationary R grid'
    CALL read_matrix(this%file_id, 'r2D', r2d, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot read stationary Z grid'
    CALL read_matrix(this%file_id, 'z2D', z2d, ierr)
    IF (ierr /= 0) RETURN
    ierr = -1
    message = 'Stationary R,Z grids have different shapes'
    IF (ANY(SHAPE(r2d) /= SHAPE(z2d))) RETURN
    this%r = r2d(1,:)
    this%z = z2d(:,1)
    message = 'Stationary equilibrium requires a shared rectilinear R,Z grid'
    IF (ANY(r2d /= SPREAD(this%r, 1, SIZE(this%z))) .OR. &
        ANY(z2d /= SPREAD(this%z, 2, SIZE(this%r)))) RETURN
    this%times = [0._float]
    message = 'Cannot read stationary psiSep'
    ALLOCATE(this%psi_sep(1))
    CALL HDF5_real_reading(this%file_id, this%psi_sep(1), 'psiSep', ierr)
  END SUBROUTINE read_legacy_metadata

  SUBROUTINE validate_coordinates(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    ierr = -1
    message = 'Equilibrium grid axes must be finite, increasing and contain at least two points'
    IF (.NOT. valid_axis(this%r) .OR. .NOT. valid_axis(this%z)) RETURN
    message = 'Equilibrium times must be finite and strictly increasing'
    IF (.NOT. ALL(ieee_is_finite(this%times))) RETURN
    IF (ANY(this%times(2:) <= this%times(:SIZE(this%times)-1))) RETURN
    message = 'psiSep must be finite and have one value per equilibrium'
    IF (SIZE(this%psi_sep) /= SIZE(this%times)) RETURN
    IF (.NOT. ALL(ieee_is_finite(this%psi_sep))) RETURN
    ierr = 0
    message = ''
  END SUBROUTINE validate_coordinates

  SUBROUTINE validate_field_datasets(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER(HSIZE_T), ALLOCATABLE :: dimensions(:)
    INTEGER(HSIZE_T) :: expected(3)
    INTEGER :: i, rank
    LOGICAL :: exists(5)

    IF (this%trajectory) THEN
       expected = INT([SIZE(this%r), SIZE(this%z), SIZE(this%times)], HSIZE_T)
       rank = 3
    ELSE
       expected(:2) = INT([SIZE(this%z), SIZE(this%r)], HSIZE_T)
       rank = 2
    ENDIF
    DO i = 1, SIZE(this%datasets)
       message = 'Cannot read equilibrium dataset: '//TRIM(this%datasets(i))
       CALL H5Lexists_f(this%file_id, TRIM(this%datasets(i)), exists(i), ierr)
       IF (ierr /= 0) RETURN
       IF (.NOT. exists(i)) THEN
          IF (i > 2) CYCLE
          ierr = -1
          RETURN
       ENDIF
       CALL HDF5_dataset_shape(this%file_id, TRIM(this%datasets(i)), dimensions, ierr)
       IF (ierr /= 0) RETURN
       ierr = -1
       message = 'Equilibrium dataset shape differs from grid/time: '//TRIM(this%datasets(i))
       IF (SIZE(dimensions) /= rank) RETURN
       IF (ANY(dimensions /= expected(:rank))) RETURN
    ENDDO
    message = 'Br and Bz must either both be present or both be absent'
    ierr = -1
    IF (exists(3) .NEQV. exists(4)) RETURN
    this%has_poloidal = exists(3)
    this%has_jtor = exists(5)
    ierr = 0
    message = ''
  END SUBROUTINE validate_field_datasets

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

  SUBROUTINE read_component(this, name, frame, values, ierr)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: frame
    REAL(float), ALLOCATABLE, INTENT(OUT) :: values(:,:)
    INTEGER, INTENT(OUT) :: ierr
    REAL(float), ALLOCATABLE :: plane(:,:)

    IF (this%trajectory) THEN
       ALLOCATE(plane(SIZE(this%r), SIZE(this%z)))
       CALL HDF5_frame_reading(this%file_id, name, frame, plane, ierr)
       IF (ierr /= 0) RETURN
       ! Python trajectory (time,z,r) reverses to Fortran (r,z,time).
       values = TRANSPOSE(plane)
    ELSE
       CALL read_matrix(this%file_id, name, values, ierr)
       IF (ierr /= 0) RETURN
    ENDIF
    IF (.NOT. ALL(ieee_is_finite(values))) ierr = -1
  END SUBROUTINE read_component

  SUBROUTINE cached_frame(this, frame, keep_frame, slot, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: frame, keep_frame
    INTEGER, INTENT(OUT) :: slot, ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    ierr = 0
    message = ''
    DO slot = 1, 2
       IF (this%cache_frame(slot) == frame) RETURN
    ENDDO
    slot = 1
    IF (this%cache_frame(slot) == keep_frame) slot = 2
    this%cache_frame(slot) = 0
    CALL read_equilibrium_frame(this, frame, this%cache(slot), ierr, message)
    IF (ierr == 0) this%cache_frame(slot) = frame
  END SUBROUTINE cached_frame

  SUBROUTINE read_equilibrium_frame(this, frame, sample, ierr, message)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: frame
    TYPE(equilibrium_sample_t), INTENT(OUT) :: sample
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: i

    ierr = 0
    message = ''
    DO i = 1, SIZE(this%datasets)
       IF (i == 3 .OR. i == 4) THEN
          IF (.NOT. this%has_poloidal) CYCLE
       ELSEIF (i == 5) THEN
          IF (.NOT. this%has_jtor) CYCLE
       ENDIF
       SELECT CASE (i)
       CASE (1)
          CALL read_component(this, TRIM(this%datasets(i)), frame, sample%psi, ierr)
       CASE (2)
          CALL read_component(this, TRIM(this%datasets(i)), frame, sample%bphi, ierr)
       CASE (3)
          CALL read_component(this, TRIM(this%datasets(i)), frame, sample%br, ierr)
       CASE (4)
          CALL read_component(this, TRIM(this%datasets(i)), frame, sample%bz, ierr)
       CASE (5)
          CALL read_component(this, TRIM(this%datasets(i)), frame, sample%jtor, ierr)
       END SELECT
       IF (ierr /= 0) THEN
          message = 'Cannot read finite equilibrium frame: '//TRIM(this%datasets(i))
          RETURN
       ENDIF
    ENDDO
    sample%psi_sep = this%psi_sep(frame)
  END SUBROUTINE read_equilibrium_frame

  SUBROUTINE series_sample_at(this, t, sample, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    REAL(float), INTENT(IN) :: t
    TYPE(equilibrium_sample_t), INTENT(OUT) :: sample
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: left, right, a, b
    REAL(float) :: actual_time, weight

    CALL locate_time_bracket(this, t, left, right, actual_time, weight, ierr, message)
    IF (ierr /= 0) RETURN
    CALL cached_frame(this, left, right, a, ierr, message)
    IF (ierr /= 0) RETURN
    CALL cached_frame(this, right, left, b, ierr, message)
    IF (ierr /= 0) RETURN
    IF (left == right) THEN
       sample = this%cache(a)
    ELSE
       CALL interpolate_frames(this%cache(a), this%cache(b), weight, sample)
    ENDIF
    sample%time = actual_time
    sample%frame_left = left-1
    sample%frame_right = right-1
    sample%fraction = weight
  END SUBROUTINE series_sample_at

  SUBROUTINE locate_time_bracket(this, t, left, right, actual_time, weight, ierr, message)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    REAL(float), INTENT(IN) :: t
    INTEGER, INTENT(OUT) :: left, right, ierr
    REAL(float), INTENT(OUT) :: actual_time, weight
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: middle, n
    REAL(float) :: tolerance

    ierr = -1
    message = 'Equilibrium series is not open'
    IF (this%file_id < 0) RETURN
    message = 'Requested equilibrium time must be finite'
    IF (.NOT. ieee_is_finite(t)) RETURN
    n = SIZE(this%times)
    left = 1
    right = 1
    actual_time = t
    weight = 0
    IF (this%trajectory .AND. n > 1) THEN
       tolerance = 64*EPSILON(t)*MAX(1._float, ABS(t), ABS(this%times(1)), ABS(this%times(n)))
       message = 'Requested equilibrium time is outside trajectory coverage'
       IF (t < this%times(1)-tolerance .OR. t > this%times(n)+tolerance) RETURN
       actual_time = MAX(this%times(1), MIN(t, this%times(n)))
       right = n
       DO WHILE (right-left > 1)
          middle = (left+right)/2
          IF (this%times(middle) <= actual_time) THEN
             left = middle
          ELSE
             right = middle
          ENDIF
       ENDDO
       IF (ABS(actual_time-this%times(left)) <= tolerance) THEN
          right = left
       ELSEIF (ABS(actual_time-this%times(right)) <= tolerance) THEN
          left = right
       ELSE
          weight = (actual_time-this%times(left))/(this%times(right)-this%times(left))
       ENDIF
    ENDIF
    ierr = 0
    message = ''
  END SUBROUTINE locate_time_bracket

  SUBROUTINE interpolate_frames(left, right, weight, sample)
    TYPE(equilibrium_sample_t), INTENT(IN) :: left, right
    REAL(float), INTENT(IN) :: weight
    TYPE(equilibrium_sample_t), INTENT(OUT) :: sample

    sample = left
    sample%psi = (1-weight)*left%psi + weight*right%psi
    sample%bphi = (1-weight)*left%bphi + weight*right%bphi
    sample%psi_sep = (1-weight)*left%psi_sep + weight*right%psi_sep
    IF (ALLOCATED(left%br)) THEN
       sample%br = (1-weight)*left%br + weight*right%br
       sample%bz = (1-weight)*left%bz + weight*right%bz
    ENDIF
    IF (ALLOCATED(left%jtor)) sample%jtor = (1-weight)*left%jtor + weight*right%jtor
  END SUBROUTINE interpolate_frames

  SUBROUTINE series_grid(this, r, z)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: r(:), z(:)
    IF (ALLOCATED(this%r)) r = this%r
    IF (ALLOCATED(this%z)) z = this%z
  END SUBROUTINE series_grid

  SUBROUTINE series_time_points(this, times)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: times(:)
    IF (ALLOCATED(this%times)) times = this%times
  END SUBROUTINE series_time_points

  SUBROUTINE series_clear(this)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    IF (this%file_id >= 0) CALL HDF5_close(this%file_id)
    this%file_id = -1
    this%trajectory = .FALSE.
    this%has_poloidal = .FALSE.
    this%has_jtor = .FALSE.
    this%cache_frame = 0
    this%cache = equilibrium_sample_t()
    IF (ALLOCATED(this%r)) DEALLOCATE(this%r)
    IF (ALLOCATED(this%z)) DEALLOCATE(this%z)
    IF (ALLOCATED(this%times)) DEALLOCATE(this%times)
    IF (ALLOCATED(this%psi_sep)) DEALLOCATE(this%psi_sep)
  END SUBROUTINE series_clear

  SUBROUTINE series_finalize(this)
    TYPE(equilibrium_series_t), INTENT(INOUT) :: this
    CALL this%clear()
  END SUBROUTINE series_finalize

END MODULE equilibrium_series
