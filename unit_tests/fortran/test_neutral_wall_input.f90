PROGRAM test_neutral_wall_input
  USE input_reader, ONLY: read_input
  USE globals, ONLY: phys, numer
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
  WRITE (*,'(A,1X,L1)') 'RECYCLE_PLASMA_DIFFUSION', phys%recycle_plasma_diffusion
  WRITE (*,'(A,1X,I0)') 'STABILIZATION_TYPE', numer%stab
  WRITE (*,'(A,5(1X,ES24.16))') 'STABILIZATION_TAU', numer%tau(1:5)
END PROGRAM test_neutral_wall_input
