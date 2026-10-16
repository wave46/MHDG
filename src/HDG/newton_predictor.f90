MODULE newton_predictor_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: newton_predictor_t

  TYPE :: newton_predictor_t
     PRIVATE
     REAL*8 :: previous_time = 0.d0
     REAL*8, ALLOCATABLE :: u(:), q(:), trace(:)
     REAL*8, ALLOCATABLE :: u_rate(:), q_rate(:), trace_rate(:)
   CONTAINS
     PROCEDURE :: remember => remember_predictor_state
     PROCEDURE :: clear => clear_predictor_state
     PROCEDURE :: weight => predictor_weight
     PROCEDURE :: apply => extrapolate_newton_guess
     PROCEDURE :: rate_sums => predictor_rate_sums
  END TYPE newton_predictor_t

CONTAINS

  SUBROUTINE remember_predictor_state(this, u, q, trace, accepted_time)
    CLASS(newton_predictor_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: u(:), q(:), trace(:), accepted_time
    REAL*8 :: elapsed

    IF (ALLOCATED(this%u)) THEN
       elapsed = accepted_time-this%previous_time
       IF (matching_history(this,u,q,trace) .AND. elapsed > 0.d0) THEN
          this%u_rate = (u-this%u)/elapsed
          this%q_rate = (q-this%q)/elapsed
          this%trace_rate = (trace-this%trace)/elapsed
       ELSE
          CALL this%clear()
       ENDIF
    ENDIF

    this%u = u
    this%q = q
    this%trace = trace
    this%previous_time = accepted_time
  END SUBROUTINE remember_predictor_state

  SUBROUTINE clear_predictor_state(this)
    CLASS(newton_predictor_t), INTENT(OUT) :: this
    ! INTENT(OUT) releases allocatable history and restores default components.
  END SUBROUTINE clear_predictor_state

  REAL*8 FUNCTION predictor_weight(this, accepted_time, next_dt) RESULT(weight)
    CLASS(newton_predictor_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: accepted_time, next_dt
    REAL*8 :: previous_dt

    weight = 0.d0
    IF (.NOT. ALLOCATED(this%u)) RETURN
    previous_dt = accepted_time-this%previous_time
    IF (previous_dt <= 0.d0 .OR. next_dt <= 0.d0) RETURN
    weight = MIN(1.d0, next_dt/previous_dt)
  END FUNCTION predictor_weight

  LOGICAL FUNCTION extrapolate_newton_guess(this, u, q, trace, weight) RESULT(applied)
    CLASS(newton_predictor_t), INTENT(IN) :: this
    REAL*8, INTENT(INOUT) :: u(:), q(:), trace(:)
    REAL*8, INTENT(IN) :: weight

    applied = .FALSE.
    IF (.NOT. matching_history(this,u,q,trace) .OR. weight <= 0.d0) RETURN
    u = u+weight*(u-this%u)
    q = q+weight*(q-this%q)
    trace = trace+weight*(trace-this%trace)
    applied = .TRUE.
  END FUNCTION extrapolate_newton_guess

  LOGICAL FUNCTION matching_history(this,u,q,trace) RESULT(matches)
    CLASS(newton_predictor_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: u(:),q(:),trace(:)

    matches = .FALSE.
    IF (.NOT. ALLOCATED(this%u)) RETURN
    matches = SIZE(u) == SIZE(this%u) .AND. SIZE(q) == SIZE(this%q) .AND. SIZE(trace) == SIZE(this%trace)
  END FUNCTION matching_history

  SUBROUTINE predictor_rate_sums(this,u,q,trace,accepted_time,neq,volume_owned,trace_owned,sums,ready)
    CLASS(newton_predictor_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: u(:),q(:),trace(:),accepted_time
    INTEGER, INTENT(IN) :: neq
    LOGICAL, INTENT(IN) :: volume_owned(:),trace_owned(:)
    REAL*8, INTENT(OUT) :: sums(3,neq,3)
    LOGICAL, INTENT(OUT) :: ready
    REAL*8 :: elapsed

    sums = 0.d0
    ready = .FALSE.
    IF (.NOT. ALLOCATED(this%u_rate)) RETURN
    IF (.NOT. matching_history(this,u,q,trace)) RETURN
    elapsed = accepted_time-this%previous_time
    IF (elapsed <= 0.d0) RETURN
    CALL sum_field_rates(u,this%u,this%u_rate,elapsed,neq,1,volume_owned,sums(:,:,1))
    CALL sum_field_rates(q,this%q,this%q_rate,elapsed,neq,SIZE(q)/SIZE(u),volume_owned,sums(:,:,2))
    CALL sum_field_rates(trace,this%trace,this%trace_rate,elapsed,neq,1,trace_owned,sums(:,:,3))
    ready = .TRUE.
  END SUBROUTINE predictor_rate_sums

  SUBROUTINE sum_field_rates(current,previous,old_rate,elapsed,neq,components,owned,sums)
    REAL*8, INTENT(IN) :: current(:),previous(:),old_rate(:),elapsed
    INTEGER, INTENT(IN) :: neq,components
    LOGICAL, INTENT(IN) :: owned(:)
    REAL*8, INTENT(OUT) :: sums(3,neq)
    INTEGER :: point,equation,component,index
    REAL*8 :: rate

    sums = 0.d0
    DO point = 1, SIZE(owned)
       IF (.NOT. owned(point)) CYCLE
       DO equation = 1, neq
          DO component = 1, components
             index = ((point-1)*neq+equation-1)*components+component
             rate = (current(index)-previous(index))/elapsed
             sums(1,equation) = sums(1,equation)+old_rate(index)**2
             sums(2,equation) = sums(2,equation)+rate**2
             sums(3,equation) = sums(3,equation)+(rate-old_rate(index))**2
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE sum_field_rates

END MODULE newton_predictor_module
