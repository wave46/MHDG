MODULE input_files
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: resolve_input_filename
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
END MODULE input_files
