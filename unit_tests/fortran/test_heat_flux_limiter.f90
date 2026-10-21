PROGRAM test_heat_flux_limiter
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys, simpar, switch
  USE physics, ONLY: compute_free_streaming_heat_flux_ions, &
       compute_free_streaming_heat_flux_electrons, &
       compute_dfree_streaming_heat_flux_ions_dU, &
       compute_dfree_streaming_heat_flux_electrons_dU, &
       compute_spitzer_harm_flux_ions, compute_spitzer_harm_flux_electrons, &
       computeVi, computeVe, compute_dV_dUi, compute_dV_dUe, &
       computeAlphai, computeAlphae, compute_dAlpha_dUi, compute_dAlpha_dUe, &
       compute_flux_limiter, compute_heat_flux_limiter_linearization, heat_flux_cap_floor
  IMPLICIT NONE
  REAL*8 :: U(5), Q(2,5), b(3), V(5), dV(5,5), da(5), dqfs(5), ju(5), jq(2,5)
  REAL*8 :: up(5), qp(2,5), temperatures(4), densities(4)
  REAL*8 :: alpha, coefficient, qfs, x, phi, gain, correction, affine, realized
  REAL*8 :: cap, qscale, scale, step, plus, minus, numerical, analytic, error, worst
  INTEGER :: species, regime, direction, energy, j, checks, reconstructions
  REAL*8, PARAMETER :: c_fl = 0.2d0, perturbation = 2.d-7, tolerance = 3.d-6

  simpar%refval_temperature = 1.d0
  phys%Mref = 1.3d0
  phys%epn = 2.5d0
  phys%diff_pari = 0.4d0
  phys%diff_pare = 0.7d0
  phys%T_fluxlim_maxi = 100.d0
  phys%T_fluxlim_maxe = 100.d0
  switch%testcase = 50
  b = [0.3d0,-0.4d0,SQRT(0.75d0)]
  ! Ordinary, temperature floor only, both floors, outer cap floor only.
  temperatures = [2.d0,0.5d-5,0.5d-5,2.d0]
  densities = [1.d0,1.d0,1.d-5,1.d-14]
  checks = 0
  reconstructions = 0
  worst = 0.d0
  DO species = 1,2
    energy = species+2
    DO regime = 1,4
      U(1) = densities(regime)
      U(2) = 0.2d0*U(1)
      U(3) = U(1)*(1.5d0*phys%Mref*temperatures(regime)+0.02d0)
      U(4) = U(1)*1.5d0*phys%Mref*temperatures(regime)
      U(5) = 0.1d0*U(1)
      CALL coefficients(U,species,V,dV,alpha,da,coefficient,qfs,dqfs)
      cap = MAX(ABS(qfs),heat_flux_cap_floor)
      qscale = c_fl*cap/(coefficient*alpha*V(energy)*b(1))
      DO direction = -1,1
        Q = 0.d0
        Q(1,energy) = DBLE(direction)*0.7d0*qscale
        x = coefficient*alpha*DOT_PRODUCT(MATMUL(Q,V),b(1:2))
        CALL compute_heat_flux_limiter_linearization(qfs,x,c_fl,DOT_PRODUCT(dqfs,U),phi,gain,correction)
        ju = phi**2*coefficient*(DOT_PRODUCT(MATMUL(Q,V),b(1:2))*da + &
             alpha*MATMUL(TRANSPOSE(MATMUL(Q,dV)),b(1:2))) + gain*dqfs
        DO j = 1,5
          jq(:,j) = phi**2*coefficient*alpha*V(j)*b(1:2)
        ENDDO
        affine = phi**2*coefficient*alpha*DOT_PRODUCT( &
             MATMUL(TRANSPOSE(MATMUL(Q,dV)),b(1:2)),U)+correction
        realized = DOT_PRODUCT(ju,U)+SUM(jq*Q)-affine
        IF (.NOT. ieee_is_finite(realized) .OR. &
             ABS(realized-phi*x) > 2.d-12*c_fl*cap) ERROR STOP 'limited heat-flux reconstruction'
        reconstructions = reconstructions+1
        DO j = 1,4
          scale = U(1)
          IF (j == 2) scale = 0.2d0*U(1)
          IF (j == 3 .OR. j == 4) scale = U(1)*1.5d0*phys%Mref*temperatures(regime)
          step = perturbation*scale
          up = U
          up(j) = U(j)+step
          CALL value(up,Q,species,plus)
          up(j) = U(j)-step
          CALL value(up,Q,species,minus)
          numerical = (plus-minus)/(2.d0*perturbation)
          analytic = ju(j)*scale
          CALL compare(numerical,analytic,c_fl*cap)
        ENDDO
        qp = Q
        qp(1,energy) = Q(1,energy)+perturbation*qscale
        CALL value(U,qp,species,plus)
        qp(1,energy) = Q(1,energy)-perturbation*qscale
        CALL value(U,qp,species,minus)
        CALL compare((plus-minus)/(2.d0*perturbation),jq(1,energy)*qscale,c_fl*cap)
      ENDDO
    ENDDO
  ENDDO
  WRITE(*,*) 'Heat limiter: PASS; derivative checks, reconstructions, worst error:',checks,reconstructions,worst

CONTAINS

  SUBROUTINE coefficients(state,which,vec,dvec,a,derivative_a,c,fs,derivative_fs)
    REAL*8, INTENT(IN) :: state(5)
    INTEGER, INTENT(IN) :: which
    REAL*8, INTENT(OUT) :: vec(5),dvec(5,5),a,derivative_a(5),c,fs,derivative_fs(5)
    IF (which == 1) THEN
      CALL computeVi(state,vec)
      CALL compute_dV_dUi(state,dvec)
      a = computeAlphai(state)
      CALL compute_dAlpha_dUi(state,derivative_a)
      CALL compute_free_streaming_heat_flux_ions(state,fs)
      CALL compute_dfree_streaming_heat_flux_ions_dU(state,derivative_fs)
      c = phys%diff_pari
    ELSE
      CALL computeVe(state,vec)
      CALL compute_dV_dUe(state,dvec)
      a = computeAlphae(state)
      CALL compute_dAlpha_dUe(state,derivative_a)
      CALL compute_free_streaming_heat_flux_electrons(state,fs)
      CALL compute_dfree_streaming_heat_flux_electrons_dU(state,derivative_fs)
      c = phys%diff_pare
    ENDIF
    c = c*(2.d0/(3.d0*phys%Mref))**(1.d0+phys%epn)
  END SUBROUTINE coefficients

  SUBROUTINE value(state,gradient,which,result)
    REAL*8, INTENT(IN) :: state(5),gradient(2,5)
    INTEGER, INTENT(IN) :: which
    REAL*8, INTENT(OUT) :: result
    REAL*8 :: fs,sh,limiter
    IF (which == 1) THEN
      CALL compute_spitzer_harm_flux_ions(state,gradient,b,sh)
      CALL compute_free_streaming_heat_flux_ions(state,fs)
    ELSE
      CALL compute_spitzer_harm_flux_electrons(state,gradient,b,sh)
      CALL compute_free_streaming_heat_flux_electrons(state,fs)
    ENDIF
    CALL compute_flux_limiter(fs,sh,c_fl,limiter)
    result = -limiter*sh ! Assembly uses the opposite of the outward conductive flux.
  END SUBROUTINE value

  SUBROUTINE compare(numeric,exact,flux_scale)
    REAL*8, INTENT(IN) :: numeric,exact,flux_scale
    error = ABS(numeric-exact)/MAX(flux_scale,ABS(numeric),ABS(exact))
    IF (.NOT. ieee_is_finite(error) .OR. error > tolerance) THEN
      WRITE(*,*) 'species, regime, direction, component:',species,regime,direction,j
      WRITE(*,*) 'finite difference, analytic, scaled error:',numeric,exact,error
      ERROR STOP 'limited heat-flux derivative'
    ENDIF
    worst = MAX(worst,error)
    checks = checks+1
  END SUBROUTINE compare

END PROGRAM test_heat_flux_limiter
