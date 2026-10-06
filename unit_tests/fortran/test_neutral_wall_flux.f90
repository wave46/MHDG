PROGRAM test_neutral_wall_flux
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys
  USE neutral_flux_limiter, ONLY: neutral_tn_source_ti, neutral_tn_source_fixed
  USE physics, ONLY: neutral_rt, neutral_thermal_speed, compute_neutral_thermal_speed, &
       compute_dneutral_thermal_speed_dU, compute_neutral_wall_flux, &
       compute_dneutral_wall_flux_dU, compute_neutral_flux_cap_speed, &
       initialize_neutral_flux_limiter_config
  IMPLICIT NONE

  phys%Mref = 1.25d0
  phys%idx_rhon_eq = 5
  neutral_rt%transport_ti_floor = 2.d-8

  CALL test_normalization()
  CALL test_state_derivatives()
  CALL test_limiter_independence()
  CALL test_inactive_and_zero_density()
  WRITE (*,'(A)') 'neutral thermal speed and wall flux checks: PASS'

CONTAINS

  SUBROUTINE test_normalization()
    REAL*8 :: U(6), cn, Fw, pi

    CALL assert_close(neutral_thermal_speed(2.25d0,4.d0),3.d0, &
         'scalar thermal speed retains Mref',1.d-14)
    CALL assert_close(neutral_thermal_speed(2.25d0,0.d0),0.d0, &
         'scalar primitive does not floor temperature',1.d-14)
    CALL make_state(0.6d0,0.4d0,U)
    CALL compute_neutral_thermal_speed(U,cn)
    CALL assert_close(cn,SQRT(phys%Mref*0.6d0), &
         'Ti wall thermal speed',1.d-14)
    pi = ACOS(-1.d0)
    CALL compute_neutral_wall_flux(U,0.d0,Fw)
    CALL assert_close(Fw,U(5)*cn/SQRT(2.d0*pi), &
         'one-sided Maxwellian flux at zero albedo',1.d-14)
    CALL compute_neutral_wall_flux(U,0.99d0,Fw)
    CALL assert_close(Fw,(1.d0-0.99d0)*U(5)*cn/SQRT(2.d0*pi), &
         'net atomic absorption at 0.99 albedo',1.d-14)
    CALL make_state(-neutral_rt%transport_ti_floor,0.d0,U)
    CALL compute_neutral_thermal_speed(U,cn)
    CALL assert_close(cn,SQRT(phys%Mref*neutral_rt%transport_ti_floor), &
         'wall speed reuses the transport temperature floor',1.d-14)
  END SUBROUTINE test_normalization

  SUBROUTINE test_state_derivatives()
    REAL*8 :: temperatures(6), albedos(3), scales(3), densities(2)
    REAL*8 :: U(6), V(6), jac(6), dcn(6)
    REAL*8 :: Fw, scaled_flux, cn, scaled_cn, velocity
    INTEGER :: it, ir, is, id, neq

    temperatures = (/0.6d0,0.5d0,1.d0,1.03d0,2.d0,-1.d0/)
    temperatures(2:) = temperatures(2:)*neutral_rt%transport_ti_floor
    albedos = (/0.d0,0.99d0,1.d0/)
    scales = (/1.d-3,1.d0,1.d3/)
    ! The small-density case would activate cons2phys's density floor.
    densities = (/1.d0,1.d-23/)
    DO neq = 5,6
      DO it = 1,SIZE(temperatures)
        velocity = 0.d0
        IF (it == 1) velocity = 0.4d0
        CALL make_state(temperatures(it),velocity,U)
        DO id = 1,SIZE(densities)
          V = U*densities(id)
          CALL compute_neutral_thermal_speed(V(:neq),cn)
          CALL compute_dneutral_thermal_speed_dU(V(:neq),dcn(:neq))
          CALL assert_close(DOT_PRODUCT(dcn(:neq),V(:neq)),0.d0, &
               'thermal-speed Euler identity',5.d-13,cn)
          DO ir = 1,SIZE(albedos)
            CALL compute_neutral_wall_flux(V(:neq),albedos(ir),Fw)
            CALL compute_dneutral_wall_flux_dU(V(:neq),albedos(ir),jac(:neq))
            CALL assert_close(DOT_PRODUCT(jac(:neq),V(:neq)),Fw, &
                 'wall-flux Euler identity',5.d-13)
            IF (Fw < 0.d0) ERROR STOP 'absorption must be outward-positive'
            IF (jac(4) /= 0.d0) ERROR STOP 'wall flux depends on electron energy'
            IF (neq == 6) THEN
              IF (jac(6) /= 0.d0) ERROR STOP 'wall flux depends on neutral momentum'
            ENDIF
            CALL check_finite_differences(V(:neq),albedos(ir),Fw,jac(:neq),cn,dcn(:neq))
            DO is = 1,SIZE(scales)
              CALL compute_neutral_wall_flux(scales(is)*V(:neq), &
                   albedos(ir),scaled_flux)
              CALL compute_neutral_thermal_speed(scales(is)*V(:neq),scaled_cn)
              CALL assert_close(scaled_flux,scales(is)*Fw, &
                   'wall-flux common scaling',5.d-13)
              CALL assert_close(scaled_cn,cn,'thermal-speed common scaling',5.d-13)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE test_state_derivatives

  SUBROUTINE check_finite_differences(U,albedo,Fw,jac,cn,dcn)
    REAL*8, INTENT(IN) :: U(:), albedo, Fw, jac(:), cn, dcn(:)
    REAL*8 :: V(SIZE(U)), steps(3), h, state_scale
    REAL*8 :: fp, fm, cp, cm, derivative_scale
    INTEGER :: j, istep

    steps = (/1.d-5,3.d-6,1.d-6/)
    DO j = 1,SIZE(U)
      state_scale = ABS(U(j))
      IF (state_scale == 0.d0) state_scale = U(1)*cn
      DO istep = 1,SIZE(steps)
        h = steps(istep)*state_scale
        V = U
        V(j) = U(j)+h
        CALL compute_neutral_wall_flux(V,albedo,fp)
        CALL compute_neutral_thermal_speed(V,cp)
        V(j) = U(j)-h
        CALL compute_neutral_wall_flux(V,albedo,fm)
        CALL compute_neutral_thermal_speed(V,cm)
        derivative_scale = MAX(ABS(jac(j)),ABS(Fw)/state_scale)
        CALL assert_close((fp-fm)/(2.d0*h),jac(j), &
             'wall-flux finite-difference derivative',2.d-7,derivative_scale)
        derivative_scale = MAX(ABS(dcn(j)),cn/state_scale)
        CALL assert_close((cp-cm)/(2.d0*h),dcn(j), &
             'thermal-speed finite-difference derivative',2.d-7,derivative_scale)
      ENDDO
    ENDDO
  END SUBROUTINE check_finite_differences

  SUBROUTINE test_limiter_independence()
    REAL*8 :: U(6), V(6), baseline, changed, jac(6), changed_jac(6), cn

    CALL make_state(0.6d0,0.4d0,U)
    phys%neutral_flux_limiter_tn_source_id = neutral_tn_source_ti
    phys%neutral_flux_limiter_tn = 0.d0
    phys%neutral_flux_limiter_eps = 0.d0
    phys%neutral_flux_limiter_fs_fraction = 1.d0
    phys%neutral_flux_limiter_fs_flux_min = 0.d0
    phys%diff_nn_min = 1.d0
    CALL initialize_neutral_flux_limiter_config()
    CALL compute_neutral_flux_cap_speed(U,cn)
    CALL assert_close(cn,SQRT(phys%Mref*0.6d0), &
         'Ti limiter cap speed retains Mref',1.d-14)
    CALL compute_neutral_wall_flux(U,0.99d0,baseline)
    CALL compute_dneutral_wall_flux_dU(U,0.99d0,jac)
    phys%neutral_flux_limiter_tn_source_id = neutral_tn_source_fixed
    phys%neutral_flux_limiter_tn = 9.d0
    CALL initialize_neutral_flux_limiter_config()
    CALL compute_neutral_flux_cap_speed(U,cn)
    CALL assert_close(cn,SQRT(phys%Mref*9.d0), &
         'fixed-Tn limiter cap speed retains Mref',1.d-14)
    ! The fixed-temperature cap must not evaluate a conservative Ti ratio.
    V = 0.d0
    CALL compute_neutral_flux_cap_speed(V,cn)
    CALL assert_close(cn,SQRT(phys%Mref*9.d0), &
         'fixed-Tn cap ignores the plasma state',1.d-14)
    CALL compute_neutral_wall_flux(U,0.99d0,changed)
    CALL compute_dneutral_wall_flux_dU(U,0.99d0,changed_jac)
    IF (changed /= baseline .OR. ANY(changed_jac /= jac)) &
         ERROR STOP 'limiter temperature selection changes wall absorption'
  END SUBROUTINE test_limiter_independence

  SUBROUTINE test_inactive_and_zero_density()
    REAL*8 :: U(6), Fw, jac(6), cn

    ! Inactive wall physics must skip even the singular temperature evaluation.
    U = 0.d0
    CALL compute_neutral_wall_flux(U,1.d0,Fw)
    CALL compute_dneutral_wall_flux_dU(U,1.d0,jac)
    IF (Fw /= 0.d0 .OR. ANY(jac /= 0.d0)) ERROR STOP 'inactive wall flux is nonzero'
    CALL make_state(0.6d0,0.4d0,U)
    U(5) = 0.d0
    CALL compute_neutral_thermal_speed(U,cn)
    CALL compute_neutral_wall_flux(U,0.d0,Fw)
    CALL compute_dneutral_wall_flux_dU(U,0.d0,jac)
    IF (Fw /= 0.d0) ERROR STOP 'zero neutral density has nonzero absorption'
    CALL assert_close(jac(5),cn/SQRT(2.d0*ACOS(-1.d0)), &
         'zero-density flux still has density derivative',1.d-14)
  END SUBROUTINE test_inactive_and_zero_density

  SUBROUTINE make_state(ti,velocity,U)
    REAL*8, INTENT(IN) :: ti, velocity
    REAL*8, INTENT(OUT) :: U(6)

    U = (/2.d0,0.d0,0.d0,0.8d0,0.03d0,-0.004d0/)
    U(2) = U(1)*velocity
    U(3) = U(1)*(1.5d0*phys%Mref*ti+0.5d0*velocity**2)
  END SUBROUTINE make_state

  SUBROUTINE assert_close(actual,expected,label,rtol,scale)
    REAL*8, INTENT(IN) :: actual, expected, rtol
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL*8, INTENT(IN), OPTIONAL :: scale
    REAL*8 :: normalization

    normalization = MAX(ABS(actual),ABS(expected),TINY(1.d0))
    IF (PRESENT(scale)) normalization = MAX(normalization,ABS(scale))
    IF (.NOT. ieee_is_finite(actual) .OR. &
         ABS(actual-expected) > rtol*normalization) THEN
      WRITE (*,'(A,2(1X,ES24.16))') label,actual,expected
      ERROR STOP 'neutral wall kernel assertion failed'
    ENDIF
  END SUBROUTINE assert_close

END PROGRAM test_neutral_wall_flux
