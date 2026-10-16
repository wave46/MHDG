MODULE solver_status
  USE MPI_OMP, ONLY: MPIvar
  USE mpi, ONLY: MPI_ALLREDUCE, MPI_BCAST, MPI_ABORT, MPI_COMM_WORLD, MPI_INTEGER, MPI_CHARACTER, MPI_MAX, MPI_MIN
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: solver_status_t, solve_success, solve_retry, solve_fatal
  PUBLIC :: agree_solver_status, stop_on_solver_failure, finite_solver_values, petsc_convergence_status

  INTEGER, PARAMETER :: solve_success = 0, solve_retry = 1, solve_fatal = 2

  TYPE :: solver_status_t
     INTEGER :: outcome = solve_success
     INTEGER :: backend_code = 0
     CHARACTER(LEN=160) :: message = ''
  END TYPE solver_status_t

CONTAINS

  LOGICAL FUNCTION finite_solver_values(values)
    REAL*8, INTENT(IN) :: values(:)

    finite_solver_values = ALL(ieee_is_finite(values))
  END FUNCTION finite_solver_values

  FUNCTION petsc_convergence_status(reason) RESULT(status)
    INTEGER, INTENT(IN) :: reason
    TYPE(solver_status_t) :: status

    IF (reason > 0) RETURN
    status%backend_code = reason
    status%message = 'PETSc KSP did not converge'
    ! Stable KSPConvergedReason values: NULL, ITS, DTOL, BREAKDOWN, BREAKDOWN_BICG, NANORINF.
    SELECT CASE (reason)
    CASE (-2,-3,-4,-5,-6,-9)
       status%outcome = solve_retry
    CASE DEFAULT
       ! Includes incompatible methods/preconditioners and unclassified PC failures.
       status%outcome = solve_fatal
    END SELECT
  END FUNCTION petsc_convergence_status

  SUBROUTINE agree_solver_status(status)
    TYPE(solver_status_t), INTENT(INOUT) :: status
    INTEGER :: outcome, owner, local_owner, ierr

    CALL MPI_ALLREDUCE(status%outcome, outcome, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    IF (outcome == solve_success) RETURN
    local_owner = HUGE(1)
    IF (status%outcome == outcome) local_owner = MPIvar%glob_id
    CALL MPI_ALLREDUCE(local_owner, owner, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    status%outcome = outcome
    CALL MPI_BCAST(status%backend_code, 1, MPI_INTEGER, owner, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(status%message, LEN(status%message), MPI_CHARACTER, owner, MPI_COMM_WORLD, ierr)
  END SUBROUTINE agree_solver_status

  SUBROUTINE stop_on_solver_failure(status)
    TYPE(solver_status_t), INTENT(IN) :: status
    INTEGER :: ierr

    IF (status%outcome == solve_success) RETURN
    IF (MPIvar%glob_id == 0) THEN
       WRITE(*,*) 'Linear solve failed: ', TRIM(status%message)
       WRITE(*,*) 'Failure outcome, backend code: ', status%outcome, status%backend_code
       FLUSH(6)
    ENDIF
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    ERROR STOP 'Linear solve failed'
  END SUBROUTINE stop_on_solver_failure

END MODULE solver_status
