PROGRAM test_smooth_transition
  USE prec_const, ONLY: float
  USE smooth_transition, ONLY: smoothstep01, smoothstep01_derivative
  IMPLICIT NONE
  REAL(float) :: s(7), expected(7), expected_derivative(7), interior(5), numeric
  REAL(float), PARAMETER :: h = 1.e-6_float
  INTEGER :: i

  s = [-2.0_float,0.0_float,0.25_float,0.5_float,0.75_float,1.0_float,2.0_float]
  expected = [0.0_float,0.0_float,0.15625_float,0.5_float,0.84375_float,1.0_float,1.0_float]
  expected_derivative = [0.0_float,0.0_float,1.125_float,1.5_float,1.125_float,0.0_float,0.0_float]
  IF (ANY(ABS(smoothstep01(s)-expected) > EPSILON(h))) ERROR STOP 'smoothstep values/clipping'
  IF (ANY(ABS(smoothstep01_derivative(s)-expected_derivative) > EPSILON(h))) &
       ERROR STOP 'smoothstep derivatives/endpoints'
  interior = [0.1_float,0.25_float,0.5_float,0.75_float,0.9_float]
  DO i = 1,SIZE(interior)
    numeric = (smoothstep01(interior(i)+h)-smoothstep01(interior(i)-h))/(2*h)
    IF (ABS(numeric-smoothstep01_derivative(interior(i))) > 1.e-9_float) &
         ERROR STOP 'smoothstep finite-difference derivative'
  ENDDO
  WRITE(*,*) 'Smooth transition: PASS; clipped values, endpoint slopes, elemental arrays and finite differences'
END PROGRAM test_smooth_transition
