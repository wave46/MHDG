PROGRAM test_ohmic_source
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys, simpar
  USE physics, ONLY: compute_Sohmic, compute_dSohmic_dU, initialize_ohmic_runtime_constants
  IMPLICIT NONE
  REAL*8 :: U(5), derivative(5), perturbed(5), density, temperature, energy_scale
  REAL*8 :: source, expected, plus, minus, step, numeric, exact, cap
  REAL*8 :: references(2), masses(2), densities(4), temperatures(8), invalid_densities(3), invalid_energies(3), worst
  INTEGER :: normalization, i, j, k, checks, invalid_n, invalid_e

  references = [1.d0,1000.d0]
  masses = [1.3d0,7.5d0]
  densities = [1.d-24,1.d-8,0.1502404d0,2.d0]
  invalid_densities = [0.d0,-1.d-24,-0.15d0]
  invalid_energies = [-1.d0,0.d0,1.d0]
  temperatures = [-0.224887d0,0.d0,0.1d0,0.999d0,1.d0,1.001d0,10.d0,1000.d0]
  phys%Zeff = 1.7d0
  phys%Pohmic = 0.23d0
  checks = 0
  worst = 0.d0
  DO normalization = 1,2
    simpar%refval_temperature = references(normalization)
    phys%Mref = masses(normalization)
    CALL initialize_ohmic_runtime_constants()
    energy_scale = 1.5d0*phys%Mref/simpar%refval_temperature
    cap = phys%Zeff*phys%Pohmic*simpar%refval_temperature**1.5d0
    DO i = 1,SIZE(densities)
      density = densities(i)
      DO j = 1,SIZE(temperatures)
        temperature = temperatures(j)
        U = 0.d0
        U(1) = density
        U(4) = density*energy_scale*temperature
        CALL compute_Sohmic(U,source)
        CALL compute_dSohmic_dU(U,derivative)
        expected = cap/MAX(temperature,1.d0)**1.5d0
        CALL compare(source,expected,expected,2.d-13)
        IF (.NOT. ALL(ieee_is_finite(derivative))) ERROR STOP 'nonfinite Ohmic Jacobian'
        IF (temperature < 1.d0 .AND. ANY(derivative /= 0.d0)) ERROR STOP 'frozen resistivity derivative'
        ! Centered differences away from the kink, including density dependence.
        IF (temperature == 1.d0) CYCLE
        DO k = 1,5
          step = 1.d-6*density
          IF (k == 4) step = step*energy_scale*MAX(ABS(temperature),1.d0)
          perturbed = U
          perturbed(k) = U(k)+step
          CALL compute_Sohmic(perturbed,plus)
          perturbed(k) = U(k)-step
          CALL compute_Sohmic(perturbed,minus)
          numeric = (plus-minus)/2.d0
          exact = derivative(k)*step
          CALL compare(numeric,exact,source*1.d-6,3.d-7)
        ENDDO
      ENDDO
      ! The source is continuous at the floor; use the frozen-side Jacobian there.
      U(4) = density*energy_scale
      CALL compute_Sohmic(U,source)
      CALL compare(source,cap,cap,2.d-13)
      U(4) = density*energy_scale*(1.d0-1.d-8)
      CALL compute_Sohmic(U,minus)
      U(4) = density*energy_scale*(1.d0+1.d-8)
      CALL compute_Sohmic(U,plus)
      CALL compare(minus,cap,cap,2.d-13)
      CALL compare(plus,cap,cap,2.d-8)
    ENDDO
    ! Reproduce the diagnosed negative-energy quadrature state directly.
    U = 0.d0
    U(1) = 0.1502404d0
    U(4) = -0.0126624d0
    CALL compute_Sohmic(U,source)
    CALL compute_dSohmic_dU(U,derivative)
    CALL compare(source,cap,cap,2.d-13)
    IF (ANY(derivative /= 0.d0)) ERROR STOP 'negative-energy Ohmic Jacobian'
    ! Nonpositive density has zero source and Jacobian for either energy sign.
    DO invalid_n = 1,SIZE(invalid_densities)
      DO invalid_e = 1,SIZE(invalid_energies)
        U = 0.d0
        U(1) = invalid_densities(invalid_n)
        U(4) = invalid_energies(invalid_e)*energy_scale
        CALL compute_Sohmic(U,source)
        CALL compute_dSohmic_dU(U,derivative)
        CALL compare(source,0.d0,cap,0.d0)
        IF (ANY(derivative /= 0.d0)) ERROR STOP 'nonpositive-density Ohmic Jacobian'
      ENDDO
    ENDDO
    ! At fixed positive energy, approach zero density from the valid side.
    U = 0.d0
    U(4) = energy_scale
    DO i = 1,2
      U(1) = densities(i)
      CALL compute_Sohmic(U,source)
      expected = cap*U(1)**1.5d0
      CALL compare(source,expected,expected,2.d-13)
    ENDDO
    ! Cache only normalization factors: Zeff and Pohmic may change at runtime.
    U(1) = 1.d0
    U(4) = -energy_scale
    phys%Zeff = 2.d0*phys%Zeff
    phys%Pohmic = 3.d0*phys%Pohmic
    CALL compute_Sohmic(U,source)
    CALL compare(source,6.d0*cap,6.d0*cap,2.d-13)
    U(4) = 10.d0*energy_scale
    CALL compute_Sohmic(U,source)
    CALL compare(source,6.d0*cap/10.d0**1.5d0,cap,2.d-13)
    phys%Zeff = phys%Zeff/2.d0
    phys%Pohmic = phys%Pohmic/3.d0
  ENDDO
  WRITE(*,*) 'Ohmic source: PASS; 1 eV floor, nonpositive density, normalization, continuity, Jacobian checks:',checks

CONTAINS

  SUBROUTINE compare(actual,target,scale,tolerance)
    REAL*8, INTENT(IN) :: actual,target,scale,tolerance
    REAL*8 :: error
    error = ABS(actual-target)/MAX(scale,ABS(actual),ABS(target),TINY(1.d0))
    IF (.NOT. ieee_is_finite(actual) .OR. error > tolerance) THEN
      WRITE(*,*) 'normalization, density, temperature, component:',normalization,i,j,k
      WRITE(*,*) 'actual, target, scaled error:',actual,target,error
      ERROR STOP 'Ohmic source check'
    ENDIF
    worst = MAX(worst,error)
    checks = checks+1
  END SUBROUTINE compare
END PROGRAM test_ohmic_source
