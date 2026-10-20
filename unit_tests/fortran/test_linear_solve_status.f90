PROGRAM test_linear_solve_status
  USE solver_status
  USE globals
  USE Main_utils, ONLY: check_newton_values
  USE MPI_OMP, ONLY: init_MPI_OMP, MPIvar
  USE mpi, ONLY: MPI_FINALIZE
#ifdef WITH_PASTIX
  USE solve_pastix, ONLY: matPASTIX, init_mat_PASTIX, check_mat_PASTIX, anal_mat_PASTIX, &
       LU_mat_pastix, solve_mat_PASTIX, terminate_mat_PASTIX
#endif
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan, ieee_positive_inf
  IMPLICIT NONE
  TYPE(solver_status_t) :: status
  INTEGER :: ierr, code, pass, reason

  CALL init_MPI_OMP()
  status = solver_status_t()
  IF (MPIvar%glob_id == MPIvar%glob_size-1) status = solver_status_t(solve_retry, -3, 'Iteration limit')
  CALL agree_solver_status(status)
  IF (status%outcome /= solve_retry .OR. status%backend_code /= -3) ERROR STOP 'rank-local failure lost'
  IF (TRIM(status%message) /= 'Iteration limit') ERROR STOP 'failure details lost'
  status = solver_status_t(solve_retry, -9, 'Nonfinite result')
  IF (MPIvar%glob_id == MPIvar%glob_size-1) status = solver_status_t(solve_fatal, 7, 'Bad parameter')
  CALL agree_solver_status(status)
  IF (status%outcome /= solve_fatal .OR. status%backend_code /= 7) ERROR STOP 'fatal failure priority'

  DO reason = -9, -2
     status = petsc_convergence_status(reason)
     IF (reason == -7 .OR. reason == -8) THEN
        IF (status%outcome /= solve_fatal) ERROR STOP 'PETSc incompatible method marked retryable'
     ELSE
        IF (status%outcome /= solve_retry) ERROR STOP 'PETSc numerical failure not retryable'
     ENDIF
  ENDDO
  status = petsc_convergence_status(-11)
  IF (status%outcome /= solve_fatal) ERROR STOP 'unclassified PC failure marked retryable'
  status = petsc_convergence_status(0)
  IF (status%outcome == solve_success) ERROR STOP 'unfinished PETSc solve accepted'
  status = petsc_convergence_status(2)
  IF (status%outcome /= solve_success) ERROR STOP 'converged PETSc solve rejected'

  ALLOCATE(sol%u(2),sol%q(4),sol%u_tilde(2))
  sol%u = 1.d0
  sol%q = 2.d0
  sol%u_tilde = [11.d0,12.d0]
  CALL check_newton_values(status)
  IF (status%outcome /= solve_success) ERROR STOP 'finite Newton iterate rejected'
  IF (MPIvar%glob_id == MPIvar%glob_size-1) sol%q(1) = ieee_value(0.d0,ieee_positive_inf)
  CALL check_newton_values(status)
  IF (status%outcome /= solve_retry) ERROR STOP 'nonfinite gradient not rejected collectively'
  sol%q = 2.d0
  CALL check_newton_values(status)
  IF (status%outcome /= solve_success) ERROR STOP 'Newton status not reset after failure'
  DEALLOCATE(sol%u,sol%q,sol%u_tilde)

#ifdef WITH_PASTIX
  utils%timing = .FALSE.
  lssolver%timing = .FALSE.
  lssolver%sollib = 1
  switch%ckeramp = .FALSE.
  ALLOCATE(sol%u_tilde(2))
  sol%u_tilde = [11.d0,12.d0]
  DO pass = 1, 2
     CALL assemble_diagonal_system()
     CALL init_mat_PASTIX(matPASTIX)
     matK%start = .FALSE.
     CALL check_mat_PASTIX(matPASTIX)
     CALL anal_mat_PASTIX(matPASTIX,code)
     IF (code /= 0) ERROR STOP 'PaStiX analysis failed'
     CALL LU_mat_pastix(matPASTIX,code)
     IF (code /= 0) ERROR STOP 'PaStiX factorization failed'
     CALL solve_mat_PASTIX(matPASTIX,code)
     IF (code /= 0 .OR. MAXVAL(ABS(matPASTIX%rhs-[1.d0,2.d0])) > 1.d-12) &
          ERROR STOP 'PaStiX result changed'

     ! One bad rank must return on every rank, retain the trace, and clear the
     ! existing backend and its borrowed matrix storage before reinitialization.
     IF (MPIvar%glob_id == MPIvar%glob_size-1) THEN
        IF (pass == 1) matK%vals(1) = ieee_value(0.d0,ieee_quiet_nan)
        IF (pass == 2) rhs%vals(1) = ieee_value(0.d0,ieee_positive_inf)
     ENDIF
#ifdef WITH_PETSC
     CALL solve_global_system(status,1)
#else
     CALL solve_global_system(status)
#endif
     IF (status%outcome /= solve_retry) ERROR STOP 'nonfinite linear input accepted'
     IF (ANY(sol%u_tilde /= [11.d0,12.d0])) ERROR STOP 'failed solve changed trace'
     IF (.NOT. matK%start .OR. ASSOCIATED(matK%vals) .OR. ASSOCIATED(rhs%vals)) &
          ERROR STOP 'failed solve retained matrix storage'
     IF (ASSOCIATED(matPASTIX%spm) .OR. ASSOCIATED(matPASTIX%pastix_data)) ERROR STOP 'failed backend retained'
  ENDDO
  ! Cleanup is safe after the failed solve already released the backend.
  CALL terminate_mat_PASTIX()
  DEALLOCATE(sol%u_tilde)
#endif
  IF (MPIvar%glob_id == 0) WRITE(*,*) 'Linear solve status: collective failures, PETSc reasons, PaStiX reinitialization PASS'
  CALL MPI_FINALIZE(ierr)

CONTAINS

  SUBROUTINE assemble_diagonal_system()
    matK%start = .TRUE.
    matK%n = 2
    matK%nnz = 2
    ALLOCATE(matK%rowptr(3),matK%cols(2),matK%vals(2),matK%loc2glob(2),rhs%vals(2))
    matK%rowptr = [1,2,3]
    matK%cols = [1,2]+2*MPIvar%glob_id
    matK%loc2glob = matK%cols
    matK%vals = 2.d0
    rhs%vals = [2.d0,4.d0]
  END SUBROUTINE assemble_diagonal_system
END PROGRAM test_linear_solve_status
