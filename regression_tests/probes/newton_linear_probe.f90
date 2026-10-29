! Test-only observation of the assembled mixed equations before relaxation.
MODULE newton_linear_probe
  USE hdg_residual_module, ONLY: hdg_residual_t, hdg_residual_norms_t, &
       hdg_residual_monitor_t, evaluate_hdg_residual
  USE MPI_OMP, ONLY: MPIvar
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: capture_newton_linear_residual, report_newton_linear_residual
  TYPE(hdg_residual_monitor_t), SAVE :: monitor
  TYPE(hdg_residual_norms_t), SAVE :: initial
CONTAINS
  SUBROUTINE capture_newton_linear_residual()
    TYPE(hdg_residual_t) :: residual
    CALL monitor%begin_solve()
    CALL evaluate_hdg_residual(residual)
    CALL monitor%measure(residual,initial)
  END SUBROUTINE

  SUBROUTINE report_newton_linear_residual()
    TYPE(hdg_residual_t) :: residual
    TYPE(hdg_residual_norms_t) :: trial
    REAL*8 :: relative
    ! No reassembly: this measures the original linear equations at the trial.
    CALL evaluate_hdg_residual(residual)
    CALL monitor%measure(residual,trial)
    relative=trial%total_rms/MAX(initial%total_rms,TINY(1.d0))
    IF(MPIvar%glob_id==0)WRITE(*,'(A,ES24.16,A,ES24.16,A,ES24.16)') &
         'Newton linear check: initial=',initial%total_rms, &
         '; trial=',trial%total_rms,'; relative=',relative
  END SUBROUTINE
END MODULE
