PROGRAM test_neutral_wall_input
  USE input_reader, ONLY: read_input
  USE hdg_residual_module, ONLY: hdg_residual
  USE globals, ONLY: phys, numer, utils
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
  WRITE (*,'(A,1X,L1)') 'BOHM_MOMENTUM_ALWAYS_INCLUDE_DIFFUSION', phys%bohm_momentum_always_include_diffusion
  WRITE (*,'(A,1X,A)') 'NEWTON_FORMULATION',TRIM(numer%newton_formulation)
  WRITE (*,'(A,1X,A)') 'NEWTON_CONVERGENCE',TRIM(numer%newton_convergence)
  WRITE (*,'(A,1X,A)') 'RESIDUAL_REPORT_MODE',TRIM(hdg_residual%get_mode())
END PROGRAM test_neutral_wall_input
