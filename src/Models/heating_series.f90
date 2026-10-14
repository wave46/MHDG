MODULE heating_series
  USE prec_const, ONLY: float
  USE rectangular_series_io, ONLY: rectangular_series_t
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC :: heating_sample_t
     REAL(float), ALLOCATABLE :: ions(:,:), electrons(:,:)
     REAL(float) :: time = 0, fraction = 0
     INTEGER :: frame_left = 0, frame_right = 0
  END TYPE heating_sample_t

  TYPE, PUBLIC :: heating_series_t
     PRIVATE
     TYPE(rectangular_series_t) :: series
     CHARACTER(LEN=32) :: datasets(2)
     TYPE(heating_sample_t) :: cache(2)
   CONTAINS
     PROCEDURE :: open => heating_open
     PROCEDURE :: sample_at => heating_sample_at
     PROCEDURE :: grid => heating_grid
     PROCEDURE :: time_points => heating_time_points
     PROCEDURE :: clear => heating_clear
     FINAL :: heating_finalize
  END TYPE heating_series_t

CONTAINS

  SUBROUTINE heating_open(this, path, ierr, message)
    CLASS(heating_series_t), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: path
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    CALL this%clear()
    CALL this%series%open(path, 'heating', ierr, message, 'Heating')
    IF (ierr /= 0) RETURN
    IF (this%series%is_trajectory()) THEN
       this%datasets = [CHARACTER(LEN=32) :: 'heating/ions', 'heating/electrons']
    ELSE
       this%datasets = [CHARACTER(LEN=32) :: 'heating2D_ions', 'heating2D_electrons']
    ENDIF
    CALL validate_heating_datasets(this, ierr, message)
    IF (ierr /= 0) CALL this%clear()
  END SUBROUTINE heating_open

  SUBROUTINE validate_heating_datasets(this, ierr, message)
    CLASS(heating_series_t), INTENT(IN) :: this
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: i
    LOGICAL :: exists

    DO i = 1, SIZE(this%datasets)
       CALL this%series%inspect_component(TRIM(this%datasets(i)), exists, ierr, message)
       IF (ierr /= 0) RETURN
       IF (.NOT. exists) THEN
          ierr = -1
          message = 'Cannot read heating dataset: '//TRIM(this%datasets(i))
          RETURN
       ENDIF
    ENDDO
  END SUBROUTINE validate_heating_datasets

  SUBROUTINE heating_sample_at(this, t, sample, ierr, message)
    CLASS(heating_series_t), INTENT(INOUT) :: this
    REAL(float), INTENT(IN) :: t
    TYPE(heating_sample_t), INTENT(OUT) :: sample
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    INTEGER :: left, right, a, b
    REAL(float) :: actual_time, weight

    CALL this%series%bracket(t, left, right, actual_time, weight, ierr, message)
    IF (ierr /= 0) RETURN
    CALL cached_heating_frame(this, left, right, a, ierr, message)
    IF (ierr /= 0) RETURN
    CALL cached_heating_frame(this, right, left, b, ierr, message)
    IF (ierr /= 0) RETURN
    sample = this%cache(a)
    IF (left /= right) THEN
       sample%ions = (1-weight)*this%cache(a)%ions + weight*this%cache(b)%ions
       sample%electrons = (1-weight)*this%cache(a)%electrons + weight*this%cache(b)%electrons
    ENDIF
    sample%time = actual_time
    sample%frame_left = left-1
    sample%frame_right = right-1
    sample%fraction = weight
  END SUBROUTINE heating_sample_at

  SUBROUTINE cached_heating_frame(this, frame, keep_frame, slot, ierr, message)
    CLASS(heating_series_t), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: frame, keep_frame
    INTEGER, INTENT(OUT) :: slot, ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    LOGICAL :: cached

    ierr = 0
    message = ''
    CALL this%series%reserve_frame(frame, keep_frame, slot, cached)
    IF (cached) RETURN
    CALL read_heating_frame(this, frame, this%cache(slot), ierr, message)
    IF (ierr == 0) CALL this%series%mark_frame(slot, frame)
  END SUBROUTINE cached_heating_frame

  SUBROUTINE read_heating_frame(this, frame, sample, ierr, message)
    CLASS(heating_series_t), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: frame
    TYPE(heating_sample_t), INTENT(OUT) :: sample
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message

    message = 'Cannot read finite heating frame: '//TRIM(this%datasets(1))
    CALL this%series%read_component(this%datasets(1), frame, sample%ions, ierr)
    IF (ierr /= 0) RETURN
    message = 'Cannot read finite heating frame: '//TRIM(this%datasets(2))
    CALL this%series%read_component(this%datasets(2), frame, sample%electrons, ierr)
    IF (ierr == 0) message = ''
  END SUBROUTINE read_heating_frame

  SUBROUTINE heating_grid(this, r, z)
    CLASS(heating_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: r(:), z(:)
    CALL this%series%grid(r, z)
  END SUBROUTINE heating_grid

  SUBROUTINE heating_time_points(this, times)
    CLASS(heating_series_t), INTENT(IN) :: this
    REAL(float), ALLOCATABLE, INTENT(OUT) :: times(:)
    CALL this%series%time_points(times)
  END SUBROUTINE heating_time_points

  SUBROUTINE heating_clear(this)
    CLASS(heating_series_t), INTENT(INOUT) :: this
    CALL this%series%clear()
    this%cache = heating_sample_t()
  END SUBROUTINE heating_clear

  SUBROUTINE heating_finalize(this)
    TYPE(heating_series_t), INTENT(INOUT) :: this
    CALL this%clear()
  END SUBROUTINE heating_finalize

END MODULE heating_series
