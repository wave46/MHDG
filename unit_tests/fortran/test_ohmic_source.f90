PROGRAM test_ohmic_source
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys, simpar
  USE physics, ONLY: compute_Sohmic, compute_dSohmic_dU, initialize_ohmic_runtime_constants
  IMPLICIT NONE
  REAL*8 :: U(5), derivative(5), perturbed(5), density, temperature, energy_scale
  REAL*8 :: source, expected, plus, minus, step, numeric, exact, reference_source
  REAL*8 :: references(2), masses(2), densities(4), temperatures(14), invalid_densities(3), invalid_energies(3)
  REAL*8 :: boundary, endpoint_derivative(5), left_derivative(5), right_derivative(5)
  INTEGER :: normalization, i, j, k, checks, invalid_n, invalid_e, endpoint

  references = [1.d0,1000.d0]
  masses = [1.3d0,7.5d0]
  densities = [1.d-24,1.d-8,0.1502404d0,2.d0]
  invalid_densities = [0.d0,-1.d-24,-0.15d0]
  invalid_energies = [-1.d0,0.d0,1.d0]
  temperatures = [-0.224887d0,0.d0,0.05d0,0.099d0,0.1d0,0.1001d0, &
       0.2d0,0.55d0,0.9d0,0.999d0,1.d0,1.001d0,10.d0,1000.d0]
  phys%Zeff = 1.7d0
  phys%Pohmic = 0.23d0
  checks = 0
  DO normalization = 1,2
    simpar%refval_temperature = references(normalization)
    phys%Mref = masses(normalization)
    CALL initialize_ohmic_runtime_constants()
    energy_scale = 1.5d0*phys%Mref/simpar%refval_temperature
    reference_source = phys%Zeff*phys%Pohmic*simpar%refval_temperature**1.5d0
    DO i = 1,SIZE(densities)
      density = densities(i)
      DO j = 1,SIZE(temperatures)
        temperature = temperatures(j)
        U = 0.d0
        U(1) = density
        U(4) = density*energy_scale*temperature
        CALL compute_Sohmic(U,source)
        CALL compute_dSohmic_dU(U,derivative)
        expected = expected_source(temperature)
        CALL compare(source,expected,reference_source,2.d-13)
        IF (.NOT. ALL(ieee_is_finite(derivative))) ERROR STOP 'nonfinite Ohmic Jacobian'
        IF (temperature < 0.1d0 .AND. ANY(derivative /= 0.d0)) ERROR STOP 'inactive Ohmic Jacobian'
        ! Temperature-based gating keeps the degree-zero source affine assembly valid.
        CALL compare(DOT_PRODUCT(derivative,U),0.d0,reference_source,2.d-13)
        IF (temperature == 0.1d0 .OR. temperature == 1.d0) CYCLE
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
          CALL compare(numeric,exact,reference_source*1.d-6,3.d-7)
        ENDDO
      ENDDO
      ! Both source and Jacobian approach the same endpoint values from each side.
      DO endpoint = 1,2
        boundary = 0.1d0
        IF (endpoint == 2) boundary = 1.d0
        U(4) = density*energy_scale*boundary
        CALL compute_dSohmic_dU(U,endpoint_derivative)
        U(4) = density*energy_scale*(boundary-1.d-8)
        CALL compute_Sohmic(U,minus)
        CALL compute_dSohmic_dU(U,left_derivative)
        U(4) = density*energy_scale*(boundary+1.d-8)
        CALL compute_Sohmic(U,plus)
        CALL compute_dSohmic_dU(U,right_derivative)
        expected = expected_source(boundary)
        CALL compare(minus,expected,reference_source,2.d-8)
        CALL compare(plus,expected,reference_source,2.d-8)
        DO k = 1,5
          step = density
          IF (k == 4) step = step*energy_scale
          CALL compare(left_derivative(k)*step,endpoint_derivative(k)*step,reference_source,1.d-5)
          CALL compare(right_derivative(k)*step,endpoint_derivative(k)*step,reference_source,1.d-5)
        ENDDO
      ENDDO
    ENDDO
    ! Reproduce the diagnosed negative-energy quadrature state directly.
    U = 0.d0
    U(1) = 0.1502404d0
    U(4) = -0.0126624d0
    CALL compute_Sohmic(U,source)
    CALL compute_dSohmic_dU(U,derivative)
    CALL compare(source,0.d0,reference_source,0.d0)
    IF (ANY(derivative /= 0.d0)) ERROR STOP 'negative-energy Ohmic Jacobian'
    DO invalid_n = 1,SIZE(invalid_densities)
      DO invalid_e = 1,SIZE(invalid_energies)
        U = 0.d0
        U(1) = invalid_densities(invalid_n)
        U(4) = invalid_energies(invalid_e)*energy_scale
        CALL compute_Sohmic(U,source)
        CALL compute_dSohmic_dU(U,derivative)
        CALL compare(source,0.d0,reference_source,0.d0)
        IF (ANY(derivative /= 0.d0)) ERROR STOP 'nonpositive-density Ohmic Jacobian'
      ENDDO
    ENDDO
    ! At fixed positive energy, approach zero density from the valid side.
    U = 0.d0
    U(4) = energy_scale
    DO i = 1,2
      U(1) = densities(i)
      CALL compute_Sohmic(U,source)
      expected = reference_source*U(1)**1.5d0
      CALL compare(source,expected,expected,2.d-13)
    ENDDO
    ! Cache normalization only: Zeff/Pohmic changes must affect transition and hot branches.
    U(1) = 1.d0
    phys%Zeff = 2.d0*phys%Zeff
    phys%Pohmic = 3.d0*phys%Pohmic
    DO j = 1,2
      temperature = 0.55d0
      IF (j == 2) temperature = 10.d0
      U(4) = temperature*energy_scale
      CALL compute_Sohmic(U,source)
      CALL compare(source,6.d0*expected_source(temperature),6.d0*reference_source,2.d-13)
    ENDDO
    phys%Zeff = phys%Zeff/2.d0
    phys%Pohmic = phys%Pohmic/3.d0
  ENDDO
  WRITE(*,*) 'Ohmic source: PASS; 0.1-1 eV gate, normalization, C1 endpoints, Jacobian and affine checks:',checks

CONTAINS

  REAL*8 FUNCTION expected_source(te) RESULT(value)
    REAL*8, INTENT(IN) :: te
    REAL*8 :: x, weight
    value = 0.d0
    IF (te <= 0.1d0) RETURN
    weight = 1.d0
    IF (te < 1.d0) THEN
      x = (te-0.1d0)/0.9d0
      weight = 3*x**2-2*x**3
    ENDIF
    value = reference_source*weight/te**1.5d0
  END FUNCTION expected_source

  SUBROUTINE compare(actual,target,scale,tolerance)
    REAL*8, INTENT(IN) :: actual,target,scale,tolerance
    REAL*8 :: error
    error = ABS(actual-target)/MAX(scale,ABS(actual),ABS(target),TINY(1.d0))
    IF (.NOT. ieee_is_finite(actual) .OR. error > tolerance) THEN
      WRITE(*,*) 'normalization, density, temperature, component:',normalization,i,j,k
      WRITE(*,*) 'actual, target, scaled error:',actual,target,error
      ERROR STOP 'Ohmic source check'
    ENDIF
    checks = checks+1
  END SUBROUTINE compare
END PROGRAM test_ohmic_source
