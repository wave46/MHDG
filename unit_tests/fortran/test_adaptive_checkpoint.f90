PROGRAM test_adaptive_checkpoint
  ! Run from lib/ with make check-adaptive-checkpoint; no plasma solve or MPI
  ! initialization is needed to exercise these production checkpoint routines.
  USE Main_utils, ONLY: update_uiter_qiter_best, update_uconv_qconv
  USE globals, ONLY: sol
  IMPLICIT NONE
  REAL*8, POINTER :: best_u(:) => NULL(), best_q(:) => NULL()
  REAL*8 :: u(2), q(4), resized_u(3), resized_q(6)

  u = [1.d0,2.d0]
  q = [3.d0,4.d0,5.d0,6.d0]
  CALL update_uiter_qiter_best(best_u,best_q,u,q)
  IF (ANY(best_u /= u) .OR. ANY(best_q /= q)) ERROR STOP 'initial checkpoint'

  ! A better Newton iterate on the SAME mesh must replace the saved values.
  u = [7.d0,8.d0]
  q = [9.d0,10.d0,11.d0,12.d0]
  CALL update_uiter_qiter_best(best_u,best_q,u,q)
  IF (ANY(best_u /= u) .OR. ANY(best_q /= q)) ERROR STOP 'stale same-size checkpoint'
  u = -1.d0
  q = -2.d0
  IF (ANY(best_u /= [7.d0,8.d0]) .OR. ANY(best_q /= [9.d0,10.d0,11.d0,12.d0])) &
    ERROR STOP 'checkpoint must own a copy'

  ! Mesh adaptation changes both buffer sizes. Subsequent same-size updates
  ! must still refresh; the old implementation froze after this allocation.
  resized_u = [13.d0,14.d0,15.d0]
  resized_q = [16.d0,17.d0,18.d0,19.d0,20.d0,21.d0]
  CALL update_uiter_qiter_best(best_u,best_q,resized_u,resized_q)
  IF (SIZE(best_u) /= 3 .OR. SIZE(best_q) /= 6) ERROR STOP 'checkpoint resize'
  resized_u = resized_u+10.d0
  resized_q = resized_q+20.d0
  CALL update_uiter_qiter_best(best_u,best_q,resized_u,resized_q)
  IF (ANY(best_u /= resized_u) .OR. ANY(best_q /= resized_q)) ERROR STOP 'stale remeshed checkpoint'

  ! Convergence saves the actual completed solution, which may differ from
  ! the previously recorded best iterate; neither copy aliases the other.
  ALLOCATE(sol%u(2),sol%q(4),sol%u_conv(3),sol%q_conv(6))
  sol%u = u
  sol%q = q
  CALL update_uconv_qconv(best_u,best_q)
  CALL update_uconv_qconv(sol%u,sol%q)
  IF (SIZE(sol%u_conv) /= 2 .OR. SIZE(sol%q_conv) /= 4) ERROR STOP 'projection checkpoint resize'
  IF (ANY(sol%u_conv /= u) .OR. ANY(sol%q_conv /= q)) ERROR STOP 'completed-state checkpoint'
  IF (ANY(best_u /= resized_u) .OR. ANY(best_q /= resized_q)) ERROR STOP 'best-iterate checkpoint changed'
  sol%u = -3.d0
  sol%q = -4.d0
  IF (ANY(sol%u_conv /= u) .OR. ANY(sol%q_conv /= q)) ERROR STOP 'completed checkpoint aliases solution'
  DEALLOCATE(best_u,best_q,sol%u,sol%q,sol%u_conv,sol%q_conv)
  WRITE(*,*) 'adaptive checkpoint: repeated updates, remeshing, independent copies and completed-state refresh PASS'
END PROGRAM test_adaptive_checkpoint
