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
     LOGICAL :: trajectory = .FALSE., numbered = .FALSE.
     LOGICAL :: has_poloidal = .FALSE., has_jtor = .FALSE.
     CHARACTER(LEN=1000) :: field_path = '', jtor_path = '', jtor_stem = ''
     REAL(float), ALLOCATABLE :: r(:), z(:), times(:), psi_sep(:)
     REAL(float), ALLOCATABLE :: jtor_r(:), jtor_z(:)
     CHARACTER(LEN=32) :: datasets(5)
     TYPE(equilibrium_sample_t) :: cache(2)
     INTEGER :: cache_frame(2) = 0
   CONTAINS
     PROCEDURE :: open => series_open
     PROCEDURE :: open_numbered => series_open_numbered
     PROCEDURE :: sample_at => series_sample_at
     PROCEDURE :: grid => series_grid
     PROCEDURE :: current_grid => series_current_grid
     PROCEDURE :: time_points => series_time_points
     PROCEDURE :: clear => series_clear
     FINAL :: series_finalize
  END TYPE equilibrium_series_t

CONTAINS

  SUBROUTINE series_open(this, path, ierr, message, jtor_path)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: path
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: jtor_path
    CHARACTER(LEN=1000) :: filename

    CALL this%clear()
    IF (PRESENT(jtor_path)) this%jtor_path = TRIM(ADJUSTL(jtor_path))
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
    IF (ierr == 0 .AND. LEN_TRIM(this%jtor_path) > 0) CALL read_current_metadata(this, ierr, message)
    IF (ierr /= 0) CALL this%clear()
  END SUBROUTINE series_open

  SUBROUTINE series_open_numbered(this, path, last_frame, ierr, message, jtor_path)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: last_frame
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: jtor_path
    CHARACTER(LEN=1000) :: field_stem, current_stem, filename, current_filename

    CALL this%clear()
    ierr = -1
    message = 'Numbered equilibrium requires frames 0000 through last_frame (0..9999)'
    IF (last_frame < 0 .OR. last_frame > 9999) RETURN
    CALL numbered_stem(path, field_stem, ierr, message)
    IF (ierr /= 0) RETURN
    current_stem = ''
    current_filename = ''
    IF (PRESENT(jtor_path)) THEN
       IF (LEN_TRIM(jtor_path) > 0) THEN
          CALL numbered_stem(jtor_path, current_stem, ierr, message)
          IF (ierr /= 0) RETURN
          current_filename = numbered_filename(current_stem, 0)
       ENDIF
    ENDIF
    filename = numbered_filename(field_stem, 0)
    CALL this%open(filename, ierr, message, current_filename)
    IF (ierr == 0) THEN
       ierr = -1
       message = 'Numbered equilibrium frames must use the legacy rectangular-grid layout'
       IF (.NOT. this%trajectory) ierr = 0
    ENDIF
    IF (ierr == 0) THEN
       this%numbered = .TRUE.
       this%field_path = field_stem
       this%jtor_stem = current_stem
       CALL read_numbered_timeline(this, last_frame, ierr, message)
    ENDIF
    IF (ierr == 0) CALL validate_coordinates(this, ierr, message)
    IF (ierr /= 0) CALL this%clear()
  END SUBROUTINE series_open_numbered

  SUBROUTINE numbered_stem(path, stem, ierr, message)
    CHARACTER(LEN=*), INTENT(IN) :: path
    CHARACTER(LEN=*), INTENT(OUT) :: stem, message
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: n

    stem = TRIM(ADJUSTL(path))
    n = LEN_TRIM(stem)
    ierr = -1
    message = 'Numbered input requires a stem or an initial _0000.h5 filename: '//TRIM(path)
    IF (n == 0) RETURN
    IF (n >= 3) THEN
       IF (stem(n-2:n) == '.h5') THEN
          IF (n <= 8) RETURN
          IF (stem(n-7:n) /= '_0000.h5') RETURN
          stem = stem(:n-8)
       ENDIF
    ENDIF
    ierr = 0
    message = ''
  END SUBROUTINE numbered_stem

  FUNCTION numbered_filename(stem, frame) RESULT(filename)
    CHARACTER(LEN=*), INTENT(IN) :: stem
    INTEGER, INTENT(IN) :: frame
    CHARACTER(LEN=1000) :: filename
    CHARACTER(LEN=4) :: index
    WRITE(index, '(I4.4)') frame
    filename = TRIM(stem)//'_'//index//'.h5'
  END FUNCTION numbered_filename

  SUBROUTINE read_numbered_timeline(this, last_frame, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: last_frame
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: frame
    REAL(float) :: current_time
    CHARACTER(LEN=1000) :: filename

    DEALLOCATE(this%times, this%psi_sep)
    ALLOCATE(this%times(last_frame+1), this%psi_sep(last_frame+1))
    DO frame = 0, last_frame
       filename = numbered_filename(this%field_path, frame)
       message = 'Cannot read numbered equilibrium metadata: '//TRIM(filename)
       CALL read_numbered_metadata(filename, this%times(frame+1), ierr, this%psi_sep(frame+1))
       IF (ierr /= 0) RETURN
       IF (LEN_TRIM(this%jtor_stem) == 0) CYCLE
       filename = numbered_filename(this%jtor_stem, frame)
       message = 'Cannot read numbered Jtor timestamp: '//TRIM(filename)
       CALL read_numbered_metadata(filename, current_time, ierr)
       IF (ierr /= 0) RETURN
       ierr = -1
       message = 'Paired equilibrium/Jtor timestamps differ: '//TRIM(filename)
       IF (.NOT. ieee_is_finite(current_time)) RETURN
       IF (ABS(current_time-this%times(frame+1)) > &
            64*EPSILON(current_time)*MAX(1._float, ABS(current_time), ABS(this%times(frame+1)))) RETURN
    ENDDO
    ierr = 0
    message = ''
  END SUBROUTINE read_numbered_timeline

  SUBROUTINE read_numbered_metadata(path, time, ierr, psi_sep)
    CHARACTER(LEN=*), INTENT(IN) :: path
    REAL(float), INTENT(OUT) :: time
    INTEGER, INTENT(OUT) :: ierr
    REAL(float), OPTIONAL, INTENT(OUT) :: psi_sep
    INTEGER(HID_T) :: file_id
    CHARACTER(LEN=1000) :: filename

    CALL resolve_input_filename(path, '.h5', filename, ierr)
    IF (ierr /= 0) RETURN
    CALL HDF5_open(filename, file_id, ierr)
    IF (ierr /= 0) RETURN
    CALL HDF5_real_reading(file_id, time, 'time', ierr)
    IF (ierr == 0 .AND. PRESENT(psi_sep)) CALL HDF5_real_reading(file_id, psi_sep, 'psiSep', ierr)
    CALL HDF5_close(file_id)
  END SUBROUTINE read_numbered_metadata

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

    this%datasets = [CHARACTER(LEN=32) :: 'flux2D', 'Bphi2D', 'Br2D', 'Bz2D', 'Jtor']
    CALL read_rectangular_grid(this%file_id, this%r, this%z, ierr, message)
    IF (ierr /= 0) RETURN
    this%times = [0._float]
    message = 'Cannot read stationary psiSep'
    ALLOCATE(this%psi_sep(1))
    CALL HDF5_real_reading(this%file_id, this%psi_sep(1), 'psiSep', ierr)
  END SUBROUTINE read_legacy_metadata

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
    message = 'Equilibrium requires a rectangular R,Z grid with finite increasing axes'
    IF (.NOT. valid_axis(r) .OR. .NOT. valid_axis(z)) RETURN
    IF (ANY(r2d /= SPREAD(r, 1, SIZE(z))) .OR. ANY(z2d /= SPREAD(z, 2, SIZE(r)))) RETURN
    ierr = 0
    message = ''
  END SUBROUTINE read_rectangular_grid

  SUBROUTINE read_current_metadata(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER(HID_T) :: file_id
    INTEGER(HSIZE_T), ALLOCATABLE :: dimensions(:)
    CHARACTER(LEN=1000) :: filename

    ierr = -1
    message = 'Trajectory Jtor must be stored in equilibrium/Jtor, not a separate file'
    IF (this%trajectory) RETURN
    CALL resolve_input_filename(this%jtor_path, '.h5', filename, ierr)
    message = 'Cannot open separate Jtor input: '//TRIM(this%jtor_path)
    IF (ierr /= 0) RETURN
    CALL HDF5_open(filename, file_id, ierr)
    IF (ierr /= 0) RETURN
    CALL read_rectangular_grid(file_id, this%jtor_r, this%jtor_z, ierr, message)
    IF (ierr == 0) THEN
       message = 'Cannot read separate Jtor dataset'
       CALL HDF5_dataset_shape(file_id, 'Jtor', dimensions, ierr)
    ENDIF
    CALL HDF5_close(file_id)
    IF (ierr /= 0) RETURN
    ierr = -1
    message = 'Jtor dataset shape differs from its R,Z grid'
    IF (SIZE(dimensions) /= 2) RETURN
    IF (ANY(dimensions /= INT([SIZE(this%jtor_z), SIZE(this%jtor_r)], HSIZE_T))) RETURN
    this%jtor_path = filename
    this%has_jtor = .TRUE.
    ierr = 0
    message = ''
  END SUBROUTINE read_current_metadata

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

    exists = .FALSE.
    IF (this%trajectory) THEN
       expected = INT([SIZE(this%r), SIZE(this%z), SIZE(this%times)], HSIZE_T)
       rank = 3
    ELSE
       expected(:2) = INT([SIZE(this%z), SIZE(this%r)], HSIZE_T)
       rank = 2
    ENDIF
    DO i = 1, SIZE(this%datasets)
       IF (i == 5 .AND. LEN_TRIM(this%jtor_path) > 0) CYCLE
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
    CALL load_equilibrium_frame(this, frame, this%cache(slot), ierr, message)
    IF (ierr == 0) this%cache_frame(slot) = frame
  END SUBROUTINE cached_frame

  SUBROUTINE load_equilibrium_frame(this, frame, sample, ierr, message)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: frame
    TYPE(equilibrium_sample_t), INTENT(OUT) :: sample
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    IF (this%numbered) THEN
       CALL load_numbered_frame(this, frame, sample, ierr, message)
    ELSE
       CALL read_equilibrium_frame(this, frame, sample, ierr, message)
       IF (ierr == 0 .AND. LEN_TRIM(this%jtor_path) > 0) &
            CALL read_current_frame(this%jtor_path, sample%jtor, ierr, message)
    ENDIF
  END SUBROUTINE load_equilibrium_frame

  SUBROUTINE load_numbered_frame(this, frame, sample, ierr, message)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: frame
    TYPE(equilibrium_sample_t), INTENT(OUT) :: sample
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    TYPE(equilibrium_series_t) :: input
    CHARACTER(LEN=1000) :: filename, current_filename, detail

    filename = numbered_filename(this%field_path, frame-1)
    current_filename = ''
    IF (LEN_TRIM(this%jtor_stem) > 0) current_filename = numbered_filename(this%jtor_stem, frame-1)
    CALL input%open(filename, ierr, message, current_filename)
    IF (ierr == 0) CALL validate_numbered_frame(this, input, ierr, message)
    IF (ierr == 0) THEN
       CALL read_equilibrium_frame(input, 1, sample, ierr, message)
       IF (ierr == 0 .AND. LEN_TRIM(current_filename) > 0) &
            CALL read_current_frame(input%jtor_path, sample%jtor, ierr, message)
    ENDIF
    CALL input%clear()
    IF (ierr /= 0) THEN
       detail = message
       message = TRIM(filename)//': '//TRIM(detail)
    ENDIF
  END SUBROUTINE load_numbered_frame

  SUBROUTINE validate_numbered_frame(this, input, ierr, message)
    CLASS(equilibrium_series_t), INTENT(IN) :: this, input
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    ierr = -1
    message = 'Numbered frames must retain their layout and optional fields'
    IF (input%trajectory) RETURN
    IF (this%has_poloidal .NEQV. input%has_poloidal) RETURN
    IF (this%has_jtor .NEQV. input%has_jtor) RETURN
    message = 'Numbered magnetic R,Z grid differs from frame 0000'
    IF (.NOT. same_axis(this%r, input%r) .OR. .NOT. same_axis(this%z, input%z)) RETURN
    IF (LEN_TRIM(this%jtor_stem) > 0) THEN
       message = 'Numbered Jtor R,Z grid differs from frame 0000'
       IF (.NOT. same_axis(this%jtor_r, input%jtor_r) .OR. &
            .NOT. same_axis(this%jtor_z, input%jtor_z)) RETURN
    ENDIF
    ierr = 0
    message = ''
  END SUBROUTINE validate_numbered_frame

  LOGICAL FUNCTION same_axis(left, right)
    REAL(float), INTENT(IN) :: left(:), right(:)
    same_axis = .FALSE.
    IF (SIZE(left) /= SIZE(right)) RETURN
    same_axis = ALL(left == right)
  END FUNCTION same_axis

  SUBROUTINE read_current_frame(path, values, ierr, message)
    CHARACTER(LEN=*), INTENT(IN) :: path
    REAL(float), ALLOCATABLE, INTENT(OUT) :: values(:,:)
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER(HID_T) :: file_id

    message = 'Cannot read finite separate Jtor frame: '//TRIM(path)
    CALL HDF5_open(path, file_id, ierr)
    IF (ierr /= 0) RETURN
    CALL read_matrix(file_id, 'Jtor', values, ierr)
    CALL HDF5_close(file_id)
    IF (ierr /= 0) RETURN
    IF (.NOT. ALL(ieee_is_finite(values))) ierr = -1
    IF (ierr == 0) message = ''
  END SUBROUTINE read_current_frame

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
          IF (.NOT. this%has_jtor .OR. LEN_TRIM(this%jtor_path) > 0) CYCLE
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
    IF (n > 1) THEN
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

  SUBROUTINE series_current_grid(this, r, z)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: r(:), z(:)
    IF (.NOT. this%has_jtor) RETURN
    IF (ALLOCATED(this%jtor_r)) THEN
       r = this%jtor_r
       z = this%jtor_z
    ELSE
       r = this%r
       z = this%z
    ENDIF
  END SUBROUTINE series_current_grid

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
    this%numbered = .FALSE.
    this%field_path = ''
    this%jtor_path = ''
    this%jtor_stem = ''
    this%has_poloidal = .FALSE.
    this%has_jtor = .FALSE.
    this%cache_frame = 0
    this%cache = equilibrium_sample_t()
    IF (ALLOCATED(this%r)) DEALLOCATE(this%r)
    IF (ALLOCATED(this%z)) DEALLOCATE(this%z)
    IF (ALLOCATED(this%times)) DEALLOCATE(this%times)
    IF (ALLOCATED(this%psi_sep)) DEALLOCATE(this%psi_sep)
    IF (ALLOCATED(this%jtor_r)) DEALLOCATE(this%jtor_r)
    IF (ALLOCATED(this%jtor_z)) DEALLOCATE(this%jtor_z)
  END SUBROUTINE series_clear

  SUBROUTINE series_finalize(this)
    TYPE(equilibrium_series_t), INTENT(INOUT) :: this
    CALL this%clear()
  END SUBROUTINE series_finalize

END MODULE equilibrium_series
