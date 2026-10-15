MODULE equilibrium_series
  USE prec_const, ONLY: float
  USE HDF5
  USE HDF5_io_module, ONLY: HDF5_open, HDF5_close, HDF5_dataset_shape, &
       HDF5_real_vector_reading, HDF5_real_reading
  USE input_files, ONLY: resolve_input_filename, numbered_stem, numbered_filename
  USE rectangular_series_io, ONLY: rectangular_series_t, read_rectangular_grid, read_matrix
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: equilibrium_sample_t
     REAL(float), ALLOCATABLE :: psi(:,:), bphi(:,:), br(:,:), bz(:,:), jtor(:,:)
     REAL(float), ALLOCATABLE :: poloidal_scale
     REAL(float) :: psi_sep = 0, time = 0, fraction = 0
     INTEGER :: frame_left = 0, frame_right = 0 ! Zero-based input positions.
  END TYPE equilibrium_sample_t

  TYPE, PUBLIC :: equilibrium_series_t
     PRIVATE
     TYPE(rectangular_series_t) :: series
     LOGICAL :: numbered = .FALSE.
     LOGICAL :: has_poloidal = .FALSE., has_jtor = .FALSE.
     CHARACTER(LEN=1000) :: field_path = '', jtor_path = '', jtor_stem = ''
     REAL(float), ALLOCATABLE :: psi_sep(:)
     REAL(float), ALLOCATABLE :: poloidal_scale
     REAL(float), ALLOCATABLE :: jtor_r(:), jtor_z(:)
     CHARACTER(LEN=32) :: datasets(5)
     TYPE(equilibrium_sample_t) :: cache(2)
   CONTAINS
     PROCEDURE :: open => series_open
     PROCEDURE :: open_numbered => series_open_numbered
     PROCEDURE :: sample_at => series_sample_at
     PROCEDURE :: grid => series_grid
     PROCEDURE :: current_grid => series_current_grid
     PROCEDURE :: is_trajectory => series_is_trajectory
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

    CALL this%clear()
    IF (PRESENT(jtor_path)) this%jtor_path = TRIM(ADJUSTL(jtor_path))
    CALL this%series%open(path, 'equilibrium', ierr, message, 'Equilibrium')
    IF (ierr /= 0) RETURN
    CALL read_equilibrium_metadata(this, ierr, message)
    IF (ierr == 0) CALL read_poloidal_scale(this, ierr, message)
    IF (ierr == 0) CALL validate_psi_sep(this, ierr, message)
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
       IF (.NOT. this%series%is_trajectory()) ierr = 0
    ENDIF
    IF (ierr == 0) THEN
       this%numbered = .TRUE.
       this%field_path = field_stem
       this%jtor_stem = current_stem
       CALL read_numbered_timeline(this, last_frame, ierr, message)
    ENDIF
    IF (ierr == 0) CALL validate_psi_sep(this, ierr, message)
    IF (ierr /= 0) CALL this%clear()
  END SUBROUTINE series_open_numbered



  SUBROUTINE read_numbered_timeline(this, last_frame, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: last_frame
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: frame
    REAL(float) :: current_time
    REAL(float), ALLOCATABLE :: times(:)
    CHARACTER(LEN=1000) :: filename

    DEALLOCATE(this%psi_sep)
    ALLOCATE(times(last_frame+1), this%psi_sep(last_frame+1))
    DO frame = 0, last_frame
       filename = numbered_filename(this%field_path, frame)
       message = 'Cannot read numbered equilibrium metadata: '//TRIM(filename)
       CALL read_numbered_metadata(filename, times(frame+1), ierr, this%psi_sep(frame+1))
       IF (ierr /= 0) RETURN
       IF (LEN_TRIM(this%jtor_stem) == 0) CYCLE
       filename = numbered_filename(this%jtor_stem, frame)
       message = 'Cannot read numbered Jtor timestamp: '//TRIM(filename)
       CALL read_numbered_metadata(filename, current_time, ierr)
       IF (ierr /= 0) RETURN
       ierr = -1
       message = 'Paired equilibrium/Jtor timestamps differ: '//TRIM(filename)
       IF (.NOT. ieee_is_finite(current_time)) RETURN
       IF (ABS(current_time-times(frame+1)) > &
            64*EPSILON(current_time)*MAX(1._float, ABS(current_time), ABS(times(frame+1)))) RETURN
    ENDDO
    CALL this%series%set_times(times, ierr, message)
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

    IF (this%series%is_trajectory()) THEN
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
    message = 'Cannot read trajectory psiSep vector'
    CALL HDF5_real_vector_reading(this%series%file(), 'equilibrium/psiSep', this%psi_sep, ierr)
  END SUBROUTINE read_trajectory_metadata

  SUBROUTINE read_legacy_metadata(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    this%datasets = [CHARACTER(LEN=32) :: 'flux2D', 'Bphi2D', 'Br2D', 'Bz2D', 'Jtor']
    message = 'Cannot read stationary psiSep'
    ALLOCATE(this%psi_sep(1))
    CALL HDF5_real_reading(this%series%file(), this%psi_sep(1), 'psiSep', ierr)
  END SUBROUTINE read_legacy_metadata

  SUBROUTINE read_poloidal_scale(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    CHARACTER(LEN=32) :: name
    LOGICAL :: exists

    name = 'poloidal_scale'
    IF (this%series%is_trajectory()) name = 'equilibrium/poloidal_scale'
    message = 'Cannot inspect poloidal_scale'
    CALL H5Lexists_f(this%series%file(), TRIM(name), exists, ierr)
    IF (ierr /= 0) RETURN
    IF (.NOT. exists) THEN
       message = ''
       RETURN
    ENDIF
    message = 'poloidal_scale must be a finite nonzero scalar'
    ALLOCATE(this%poloidal_scale)
    CALL HDF5_real_reading(this%series%file(), this%poloidal_scale, TRIM(name), ierr)
    IF (ierr /= 0) RETURN
    ierr = -1
    IF (.NOT. ieee_is_finite(this%poloidal_scale)) RETURN
    IF (this%poloidal_scale == 0._float) RETURN
    ierr = 0
    message = ''
  END SUBROUTINE read_poloidal_scale

  SUBROUTINE read_current_metadata(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER(HID_T) :: file_id
    INTEGER(HSIZE_T), ALLOCATABLE :: dimensions(:)
    CHARACTER(LEN=1000) :: filename

    ierr = -1
    message = 'Trajectory Jtor must be stored in equilibrium/Jtor, not a separate file'
    IF (this%series%is_trajectory()) RETURN
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

  SUBROUTINE validate_psi_sep(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    ierr = -1
    message = 'psiSep must be finite and have one value per equilibrium'
    IF (SIZE(this%psi_sep) /= this%series%frame_count()) RETURN
    IF (.NOT. ALL(ieee_is_finite(this%psi_sep))) RETURN
    ierr = 0
    message = ''
  END SUBROUTINE validate_psi_sep

  SUBROUTINE validate_field_datasets(this, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: i
    LOGICAL :: exists(5)

    exists = .FALSE.
    DO i = 1, SIZE(this%datasets)
       IF (i == 5 .AND. LEN_TRIM(this%jtor_path) > 0) CYCLE
       CALL this%series%inspect_component(TRIM(this%datasets(i)), exists(i), ierr, message)
       IF (ierr /= 0) RETURN
       IF (.NOT. exists(i) .AND. i <= 2) THEN
          ierr = -1
          message = 'Cannot read equilibrium dataset: '//TRIM(this%datasets(i))
          RETURN
       ENDIF
    ENDDO
    message = 'Br and Bz must either both be present or both be absent'
    ierr = -1
    IF (exists(3) .NEQV. exists(4)) RETURN
    this%has_poloidal = exists(3)
    this%has_jtor = exists(5)
    ierr = 0
    message = ''
  END SUBROUTINE validate_field_datasets

  SUBROUTINE cached_frame(this, frame, keep_frame, slot, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: frame, keep_frame
    INTEGER, INTENT(OUT) :: slot, ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    LOGICAL :: cached

    ierr = 0
    message = ''
    CALL this%series%reserve_frame(frame, keep_frame, slot, cached)
    IF (cached) RETURN
    CALL load_equilibrium_frame(this, frame, this%cache(slot), ierr, message)
    IF (ierr == 0) CALL this%series%mark_frame(slot, frame)
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
    IF (input%series%is_trajectory()) RETURN
    IF (this%has_poloidal .NEQV. input%has_poloidal) RETURN
    IF (this%has_jtor .NEQV. input%has_jtor) RETURN
    message = 'Numbered poloidal_scale differs from frame 0000'
    IF (ALLOCATED(this%poloidal_scale) .NEQV. ALLOCATED(input%poloidal_scale)) RETURN
    IF (ALLOCATED(this%poloidal_scale)) THEN
       IF (this%poloidal_scale /= input%poloidal_scale) RETURN
    ENDIF
    message = 'Numbered magnetic R,Z grid differs from frame 0000'
    IF (.NOT. this%series%same_grid(input%series)) RETURN
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
          CALL this%series%read_component(TRIM(this%datasets(i)), frame, sample%psi, ierr)
       CASE (2)
          CALL this%series%read_component(TRIM(this%datasets(i)), frame, sample%bphi, ierr)
       CASE (3)
          CALL this%series%read_component(TRIM(this%datasets(i)), frame, sample%br, ierr)
       CASE (4)
          CALL this%series%read_component(TRIM(this%datasets(i)), frame, sample%bz, ierr)
       CASE (5)
          CALL this%series%read_component(TRIM(this%datasets(i)), frame, sample%jtor, ierr)
       END SELECT
       IF (ierr /= 0) THEN
          message = 'Cannot read finite equilibrium frame: '//TRIM(this%datasets(i))
          RETURN
       ENDIF
    ENDDO
    sample%psi_sep = this%psi_sep(frame)
    IF (ALLOCATED(this%poloidal_scale)) sample%poloidal_scale = this%poloidal_scale
  END SUBROUTINE read_equilibrium_frame

  SUBROUTINE series_sample_at(this, t, sample, ierr, message)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    REAL(float), INTENT(IN) :: t
    TYPE(equilibrium_sample_t), INTENT(OUT) :: sample
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: left, right, a, b
    REAL(float) :: actual_time, weight

    CALL this%series%bracket(t, left, right, actual_time, weight, ierr, message)
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
    CALL this%series%grid(r, z)
  END SUBROUTINE series_grid

  SUBROUTINE series_current_grid(this, r, z)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: r(:), z(:)
    IF (.NOT. this%has_jtor) RETURN
    IF (ALLOCATED(this%jtor_r)) THEN
       r = this%jtor_r
       z = this%jtor_z
    ELSE
       CALL this%series%grid(r, z)
    ENDIF
  END SUBROUTINE series_current_grid

  LOGICAL FUNCTION series_is_trajectory(this)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    series_is_trajectory = this%series%is_trajectory()
  END FUNCTION series_is_trajectory

  SUBROUTINE series_time_points(this, times)
    CLASS(equilibrium_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: times(:)
    CALL this%series%time_points(times)
  END SUBROUTINE series_time_points

  SUBROUTINE series_clear(this)
    CLASS(equilibrium_series_t), INTENT(INOUT) :: this
    CALL this%series%clear()
    this%numbered = .FALSE.
    this%field_path = ''
    this%jtor_path = ''
    this%jtor_stem = ''
    this%has_poloidal = .FALSE.
    this%has_jtor = .FALSE.
    this%cache = equilibrium_sample_t()
    IF (ALLOCATED(this%psi_sep)) DEALLOCATE(this%psi_sep)
    IF (ALLOCATED(this%poloidal_scale)) DEALLOCATE(this%poloidal_scale)
    IF (ALLOCATED(this%jtor_r)) DEALLOCATE(this%jtor_r)
    IF (ALLOCATED(this%jtor_z)) DEALLOCATE(this%jtor_z)
  END SUBROUTINE series_clear

  SUBROUTINE series_finalize(this)
    TYPE(equilibrium_series_t), INTENT(INOUT) :: this
    CALL this%clear()
  END SUBROUTINE series_finalize

END MODULE equilibrium_series
