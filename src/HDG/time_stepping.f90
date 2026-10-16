MODULE time_stepping
  USE globals, ONLY: time, simpar, sol
  USE rectangular_series_io, ONLY: locate_series_time
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: next_physical_time, physical_time_is_complete, schedule_half_step
  PUBLIC :: remove_reached_target, ensure_time_history_capacity

CONTAINS

  SUBROUTINE next_physical_time(target, finished, ierr, message, frame_times)
    REAL*8, INTENT(OUT) :: target
    LOGICAL, INTENT(OUT) :: finished
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    REAL*8, OPTIONAL, INTENT(IN) :: frame_times(:) ! Physical seconds, as stored by the provider.
    REAL*8 :: tolerance

    target = time%t
    finished = .FALSE.
    ierr = -1
    message = 'Accepted time exceeds t_final'
    tolerance = 64*EPSILON(time%t)*MAX(1.d0, ABS(time%t), ABS(time%t_final))
    IF (time%t > time%t_final+tolerance) RETURN
    IF (physical_time_is_complete()) THEN
       finished = .TRUE.
       ierr = 0
       message = ''
       RETURN
    ENDIF

    SELECT CASE (time%mode)
    CASE ('fixed')
       target = MIN(time%t_final, time%t+time%dt0)
    CASE ('adaptive')
       target = MIN(time%t_final, time%t+time%dt_next)
    CASE ('prescribed')
       CALL next_prescribed_time(target, ierr, message, frame_times)
       IF (ierr /= 0) RETURN
    CASE DEFAULT
       message = 'Unknown physical time mode'
       RETURN
    END SELECT
    IF (ALLOCATED(time%pending_times)) THEN
       IF (SIZE(time%pending_times) > 0) target = MIN(time%t_final, time%pending_times(SIZE(time%pending_times)))
    ENDIF
    ierr = -1
    message = 'Timestep cannot advance physical time'
    IF (target <= time%t) RETURN
    message = 'Selected timestep is smaller than dt_min'
    IF (time%dt_reduce_on_failure .AND. target /= time%t_final) THEN
       IF (target-time%t < time%dt_min*(1.d0-64*EPSILON(target))) RETURN
    ENDIF
    ierr = 0
    message = ''
  END SUBROUTINE next_physical_time

  SUBROUTINE schedule_half_step(failed_endpoint, ierr, message)
    REAL*8, INTENT(IN) :: failed_endpoint
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    REAL*8 :: reduced_dt, midpoint
    INTEGER :: pending

    ! Called after rollback: time%t and the pending list belong to the accepted state.
    ierr = -1
    reduced_dt = (failed_endpoint-time%t)/2.d0
    message = 'Halving would reduce the timestep below dt_min'
    IF (reduced_dt < time%dt_min*(1.d0-64*EPSILON(reduced_dt))) RETURN
    midpoint = time%t+reduced_dt
    message = 'Cannot bisect timestep at floating-point precision'
    IF (midpoint <= time%t .OR. midpoint >= failed_endpoint) RETURN

    IF (time%mode == 'prescribed') THEN
       pending = 0
       IF (ALLOCATED(time%pending_times)) pending = SIZE(time%pending_times)
       IF (pending == 0) THEN
          time%pending_times = [failed_endpoint, midpoint]
       ELSE
          ! A failed pending target remains on the stack while its midpoint is tried.
          time%pending_times = [time%pending_times, midpoint]
       ENDIF
    ELSE
       time%dt_next = reduced_dt
    ENDIF
    ierr = 0
    message = ''
  END SUBROUTINE schedule_half_step

  SUBROUTINE remove_reached_target()
    INTEGER :: pending

    IF (.NOT. ALLOCATED(time%pending_times)) RETURN
    pending = SIZE(time%pending_times)
    IF (pending <= 1 .OR. physical_time_is_complete()) THEN
       DEALLOCATE(time%pending_times)
    ELSE
       time%pending_times = time%pending_times(:pending-1)
    ENDIF
  END SUBROUTINE remove_reached_target

  LOGICAL FUNCTION physical_time_is_complete()
    REAL*8 :: tolerance

    tolerance = 64*EPSILON(time%t)*MAX(1.d0, ABS(time%t), ABS(time%t_final))
    physical_time_is_complete = time%t >= time%t_final-tolerance
  END FUNCTION physical_time_is_complete

  SUBROUTINE next_prescribed_time(target, ierr, message, frame_times)
    REAL*8, INTENT(OUT) :: target
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(LEN=*), INTENT(OUT) :: message
    REAL*8, OPTIONAL, INTENT(IN) :: frame_times(:)
    REAL*8 :: actual_time, weight
    INTEGER :: left, right

    target = time%t
    ierr = -1
    message = 'Prescribed mode requires at least two equilibrium timestamps'
    IF (.NOT. PRESENT(frame_times)) RETURN
    IF (SIZE(frame_times) < 2) RETURN
    CALL locate_series_time(frame_times, time%t_final*simpar%refval_time, &
         left, right, actual_time, weight, ierr, message)
    IF (ierr /= 0) RETURN
    CALL locate_series_time(frame_times, time%t*simpar%refval_time, &
         left, right, actual_time, weight, ierr, message)
    IF (ierr /= 0) RETURN
    target = time%t_final
    IF (left < SIZE(frame_times)) target = MIN(target, frame_times(left+1)/simpar%refval_time)
  END SUBROUTINE next_prescribed_time

  SUBROUTINE ensure_time_history_capacity(required)
    INTEGER, INTENT(IN) :: required
    INTEGER :: capacity, previous
    REAL*8, ALLOCATABLE :: times(:), residuals(:)

    previous = 0
    IF (ALLOCATED(sol%time)) previous = SIZE(sol%time)
    IF (required <= previous) RETURN
    capacity = MAX(required, 16, 2*previous)
    ALLOCATE(times(capacity), residuals(capacity))
    times = 0.d0
    residuals = 0.d0
    IF (previous > 0) THEN
       times(:previous) = sol%time
       residuals(:previous) = sol%tres
    ENDIF
    CALL MOVE_ALLOC(times, sol%time)
    CALL MOVE_ALLOC(residuals, sol%tres)
  END SUBROUTINE ensure_time_history_capacity

END MODULE time_stepping
