MODULE smooth_transition
  USE prec_const, ONLY: float
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: smoothstep01, smoothstep01_derivative

CONTAINS

  PURE ELEMENTAL FUNCTION smoothstep01(s) RESULT(weight)
    REAL(float), INTENT(IN) :: s
    REAL(float) :: weight, sc

    sc = MAX(MIN(s, 1.0_float), 0.0_float)
    weight = sc*sc*(3.0_float - 2.0_float*sc)
  END FUNCTION smoothstep01

  PURE ELEMENTAL FUNCTION smoothstep01_derivative(s) RESULT(derivative)
    REAL(float), INTENT(IN) :: s
    REAL(float) :: derivative, sc

    sc = MAX(MIN(s, 1.0_float), 0.0_float)
    derivative = 6.0_float*sc*(1.0_float - sc)
  END FUNCTION smoothstep01_derivative

END MODULE smooth_transition
