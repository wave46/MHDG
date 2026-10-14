MODULE input_files
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: resolve_input_filename, numbered_stem, numbered_filename
CONTAINS
  SUBROUTINE resolve_input_filename(path, extension, filename, ierr)
    CHARACTER(LEN=*), INTENT(IN) :: path, extension
    CHARACTER(LEN=*), INTENT(OUT) :: filename
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: n, ext_length
    LOGICAL :: exists

    ierr = -1
    filename = TRIM(ADJUSTL(path))
    n = LEN_TRIM(filename)
    IF (n == 0) RETURN
    INQUIRE(FILE=TRIM(filename), EXIST=exists)
    IF (exists) THEN
       ierr = 0
       RETURN
    ENDIF
    ext_length = LEN_TRIM(extension)
    IF (ext_length == 0) RETURN
    IF (n >= ext_length) THEN
       IF (filename(n-ext_length+1:n) == TRIM(extension)) RETURN
    ENDIF
    filename = TRIM(ADJUSTL(path))//TRIM(extension)
    INQUIRE(FILE=TRIM(filename), EXIST=exists)
    IF (exists) ierr = 0
  END SUBROUTINE resolve_input_filename
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
END MODULE input_files
