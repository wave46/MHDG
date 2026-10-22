PROGRAM test_neutral_wall_input
  USE input_reader, ONLY: read_input
  USE globals, ONLY: phys
  USE MPI_OMP, ONLY: MPIvar
  IMPLICIT NONE
  CHARACTER(LEN=16) :: argument

  CALL get_command_argument(1, argument)
  IF (argument == '--build-kind') THEN
#ifdef KEQUATION
     WRITE (*,'(A)') 'k'
#else
     WRITE (*,'(A)') 'non-k'
#endif
     STOP
  ENDIF

  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  CALL read_input()
  WRITE (*,'(A,2(1X,ES24.16))') 'NEUTRAL_WALL_ALBEDOS', &
       phys%Re_n, phys%Re_n_pump
END PROGRAM test_neutral_wall_input
