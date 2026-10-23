PROGRAM test_newton_damping
  USE globals
  USE Main_utils, ONLY: finish_newton_update
  USE element_solution_module, ONLY: compute_element_solution
  USE global_system_solver, ONLY: solve_global_system
  USE solver_status, ONLY: solver_status_t, solve_success, solve_retry
  USE MPI_OMP, ONLY: init_MPI_OMP, MPIvar
  USE mpi, ONLY: MPI_FINALIZE
#ifdef WITH_PASTIX
  USE solve_pastix, ONLY: terminate_mat_PASTIX
#endif
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan
  IMPLICIT NONE
  REAL*8 :: old_u(2), old_q(4), old_trace(6), trial_u(2), trial_q(4), trial_trace(6)
  REAL*8 :: factors(3), alpha, double_damped(2)
  TYPE(solver_status_t) :: status
  INTEGER :: i,j,pass,ierr

  CALL init_MPI_OMP()
  utils%timing=.FALSE.
  phys%neq=2
  Mesh%Ndim=2
  Mesh%Nelems=1
  Mesh%Nnodesperelem=1
  Mesh%Nnodesperface=1
  refElPol%Nfacenodes=1
  refElPol%Nfaces=3
  ALLOCATE(Mesh%F(1,3))
  Mesh%F(1,:)=[1,2,3]
  ALLOCATE(sol%u(2),sol%q(4),sol%u_tilde(6))
  ALLOCATE(elMat%UU(2,6,1),elMat%U0(2,1),elMat%LL(4,6,1),elMat%L0(4,1))
  DO j=1,6
     DO i=1,2
        elMat%UU(i,j,1)=0.07d0*(i+2*j)
     ENDDO
     DO i=1,4
        elMat%LL(i,j,1)=0.03d0*(2*i-j)
     ENDDO
  ENDDO
  ! Deliberately nonzero local residual: the old state does not lie on the
  ! newly assembled affine relation. This distinguishes full-state damping
  ! from trace-only damping and from damping local variables twice.
  elMat%U0(:,1)=[-2.d0,3.d0]
  elMat%L0(:,1)=[1.d0,-4.d0,2.d0,-1.d0]
  old_u=[1.d0,2.d0]+MPIvar%glob_id
  old_q=[-1.d0,3.d0,2.d0,-4.d0]+MPIvar%glob_id
  old_trace=[(0.2d0*i,i=1,6)]+MPIvar%glob_id
  trial_trace=old_trace+[2.d0,-1.d0,3.d0,1.d0,0.d0,0.d0]
  ! Last face represents prescribed values, already imposed in the snapshot.
  sol%u_tilde=trial_trace
  CALL compute_element_solution()
  trial_u=sol%u
  trial_q=sol%q
  factors=[0.d0,0.25d0,1.d0]
  DO pass=1,SIZE(factors)
     alpha=factors(pass)
     sol%u=trial_u
     sol%q=trial_q
     sol%u_tilde=trial_trace
     status=solver_status_t()
     CALL finish_newton_update(status,old_u,old_q,old_trace,alpha)
     IF (status%outcome/=solve_success) ERROR STOP 'finite Newton trial rejected'
     CALL close(sol%u,(1.d0-alpha)*old_u+alpha*trial_u)
     CALL close(sol%q,(1.d0-alpha)*old_q+alpha*trial_q)
     CALL close(sol%u_tilde,(1.d0-alpha)*old_trace+alpha*trial_trace)
     IF (ANY(sol%u_tilde(5:6)/=old_trace(5:6))) ERROR STOP 'prescribed trace relaxed'
     IF (alpha==1.d0) THEN
        IF (ANY(sol%u/=trial_u).OR.ANY(sol%q/=trial_q).OR.ANY(sol%u_tilde/=trial_trace)) &
             ERROR STOP 'unit damping changed the full trial'
     ENDIF
  ENDDO
  alpha=0.25d0
  double_damped=(1.d0-alpha)*old_u+alpha*(MATMUL(elMat%UU(:,:,1), &
       (1.d0-alpha)*old_trace+alpha*trial_trace)+elMat%U0(:,1))
  IF (MAXVAL(ABS(double_damped-((1.d0-alpha)*old_u+alpha*trial_u)))<0.1d0) &
       ERROR STOP 'fixture fails to expose double damping'

  sol%u=trial_u
  sol%q=trial_q
  sol%u_tilde=trial_trace
  IF (MPIvar%glob_id==MPIvar%glob_size-1) sol%q(1)=ieee_value(0.d0,ieee_quiet_nan)
  status=solver_status_t()
  CALL finish_newton_update(status,old_u,old_q,old_trace,alpha)
  IF (status%outcome/=solve_retry) ERROR STOP 'rank-local nonfinite trial accepted'
  CALL unchanged()
#ifdef WITH_PASTIX
  CALL check_trace_solve()
#endif
  sol%u=trial_u
  sol%q=trial_q
  sol%u_tilde=trial_trace
  status=solver_status_t(solve_retry,17,'mock linear failure')
  CALL finish_newton_update(status,old_u,old_q,old_trace,alpha)
  IF (status%outcome/=solve_retry) ERROR STOP 'rollback silently accepted a failed linear solve'
  IF (status%backend_code/=17) ERROR STOP 'failed linear status lost'
  CALL unchanged()
  IF (MPIvar%glob_id==0) WRITE(*,*) 'Full-state damping: condensed trial, alpha 0/0.25/1, prescribed trace, collective rollback PASS'
  CALL MPI_FINALIZE(ierr)
CONTAINS

#ifdef WITH_PASTIX
  SUBROUTINE check_trace_solve()
    INTEGER :: pass,j
    Mesh%Nfaces=3
    Mesh%ukf=2
    Mesh%Ndir=1
#ifdef PARALL
    Mesh%nghostfaces=0
    ALLOCATE(Mesh%ghostfaces(3),Mesh%fc2sd(0),Mesh%fc2rv(0),Mesh%pr2sd(0),Mesh%pr2rv(0))
    Mesh%ghostfaces=0
#endif
    lssolver%sollib=1
    lssolver%timing=.FALSE.
    switch%ckeramp=.FALSE.
    DO pass=1,2
       alpha=MERGE(0.25d0,1.d0,pass==1)
       numer%dumpnr=alpha
       matK%start=.TRUE.
       matK%n=4
       matK%nnz=4
       ALLOCATE(matK%rowptr(5),matK%cols(4),matK%vals(4),matK%loc2glob(4),rhs%vals(4))
       matK%rowptr=[1,2,3,4,5]
       matK%cols=[(j,j=1,4)]+4*MPIvar%glob_id
       matK%loc2glob=matK%cols
       matK%vals=2.d0
       rhs%vals=2.d0*trial_trace(1:4)
       sol%u=old_u
       sol%q=old_q
       sol%u_tilde=old_trace
#ifdef WITH_PETSC
       CALL solve_global_system(status,1)
#else
       CALL solve_global_system(status)
#endif
       IF (status%outcome/=solve_success) ERROR STOP 'full trace solve failed'
       ! The linear trial is independent of the Newton damping factor.
       CALL close(sol%u_tilde,trial_trace)
       CALL compute_element_solution()
       CALL finish_newton_update(status,old_u,old_q,old_trace,alpha)
       CALL close(sol%u,(1.d0-alpha)*old_u+alpha*trial_u)
       CALL close(sol%q,(1.d0-alpha)*old_q+alpha*trial_q)
       CALL close(sol%u_tilde,(1.d0-alpha)*old_trace+alpha*trial_trace)
       IF (status%outcome/=solve_success) ERROR STOP 'trace solver damping check failed'
       IF (ANY(sol%u_tilde(5:6)/=old_trace(5:6))) ERROR STOP 'linear solve overwrote prescribed trace'
       CALL terminate_mat_PASTIX()
    ENDDO
  END SUBROUTINE check_trace_solve
#endif
  SUBROUTINE close(actual,expected)
    REAL*8,INTENT(IN)::actual(:),expected(:)
    IF (MAXVAL(ABS(actual-expected))>2.d-13) ERROR STOP 'inconsistent state damping'
  END SUBROUTINE close
  SUBROUTINE unchanged()
    IF (ANY(sol%u/=old_u).OR.ANY(sol%q/=old_q).OR.ANY(sol%u_tilde/=old_trace)) &
         ERROR STOP 'failed trial changed saved state'
  END SUBROUTINE unchanged
END PROGRAM test_newton_damping
