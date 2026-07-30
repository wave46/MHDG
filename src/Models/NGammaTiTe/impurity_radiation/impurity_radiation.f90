MODULE impurity_radiation_model
  USE globals, ONLY: phys, simpar
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_impurity_radiation_model
  PUBLIC :: adimensionalize_impurity_radiation_model
  PUBLIC :: compute_impurity_cooling

  INTEGER, PARAMETER :: cooling_coeff_count = 17
  REAL*8 :: log_energy_rate_scale = 0.d0

CONTAINS

  SUBROUTINE init_impurity_radiation_model()
    INTEGER :: i

    IF (phys%n_impurities < 1) THEN
       WRITE(6,*) 'Impurity radiation requires at least one configured entry'
       STOP
    END IF

    IF (ALLOCATED(phys%alpha_cooling_factor_impurities)) DEALLOCATE(phys%alpha_cooling_factor_impurities)
    ALLOCATE(phys%alpha_cooling_factor_impurities(cooling_coeff_count, phys%n_impurities))

    DO i = 1, phys%n_impurities
       CALL load_impurity_cooling_coefficients(phys%impurity_names(i), &
            &phys%alpha_cooling_factor_impurities(:,i))
    END DO
  END SUBROUTINE init_impurity_radiation_model

  SUBROUTINE load_impurity_cooling_coefficients(impurity_name, coefficients)
    CHARACTER(LEN=*), INTENT(IN) :: impurity_name
    REAL*8, INTENT(OUT) :: coefficients(:)

    SELECT CASE (TRIM(ADJUSTL(impurity_name)))
    CASE ('N')
       coefficients = (/-2.49348163e+01,  9.52628451e+00, -4.39511346e+00,  2.00446916e+00,&
       -2.27166819e+00,  9.95587194e-01,  1.11209167e+00, -1.13089140e+00,&
        2.22902131e-01,  1.37491696e-01, -9.74320916e-02,  2.88337222e-02,&
       -5.03419501e-03,  5.54252119e-04, -3.80153983e-05,  1.49061146e-06,&
       -2.56122449e-08/)
    CASE ('W')
       coefficients = (/-3.29798537e+01,  3.38045169e+01, -3.81202398e+01,  2.47450333e+01,&
       -9.04921504e+00,  1.91417658e+00, -2.32405139e-01,  1.46022765e-02,&
       -2.46177645e-04, -1.89159754e-05,  7.73790432e-07,  0.d0,&
        0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    CASE DEFAULT
       WRITE(6,*) 'Cooling factor not defined for impurity ', TRIM(ADJUSTL(impurity_name))
       STOP
    END SELECT
  END SUBROUTINE load_impurity_cooling_coefficients

  SUBROUTINE adimensionalize_impurity_radiation_model()
    REAL*8 :: log_temperature_shift
    INTEGER :: i

    log_temperature_shift = LOG(simpar%refval_temperature)
    log_energy_rate_scale = LOG(simpar%refval_density*simpar%refval_time* &
         &simpar%refval_charge/simpar%refval_mass*simpar%refval_time**2/simpar%refval_length**2)

    DO i = 1, phys%n_impurities
       CALL shift_log_polynomial(phys%alpha_cooling_factor_impurities(:,i), &
            &log_temperature_shift, log_energy_rate_scale)
    END DO
  END SUBROUTINE adimensionalize_impurity_radiation_model

  SUBROUTINE compute_impurity_cooling(te, cooling, dcooling_dte)
    REAL*8, INTENT(IN) :: te
    REAL*8, INTENT(OUT) :: cooling
    REAL*8, INTENT(OUT), OPTIONAL :: dcooling_dte
    REAL*8 :: entry_cooling, entry_derivative
    INTEGER :: i

    cooling = 0.d0
    IF (PRESENT(dcooling_dte)) dcooling_dte = 0.d0

    DO i = 1, phys%n_impurities
       CALL evaluate_cooling_polynomial(te, phys%alpha_cooling_factor_impurities(:,i), &
            &entry_cooling, entry_derivative)
       cooling = cooling + phys%impurity_concentrations(i)*entry_cooling
       IF (PRESENT(dcooling_dte)) THEN
          dcooling_dte = dcooling_dte + phys%impurity_concentrations(i)*entry_derivative
       END IF
    END DO
  END SUBROUTINE compute_impurity_cooling

  SUBROUTINE evaluate_cooling_polynomial(te, coefficients, cooling, dcooling_dte)
    REAL*8, INTENT(IN) :: te, coefficients(:)
    REAL*8, INTENT(OUT) :: cooling, dcooling_dte
    REAL*8 :: te_min, te_max, te_fit, log_cooling, log_slope

    te_min = 1.d-1/simpar%refval_temperature
    te_max = 3.d3/simpar%refval_temperature
    te_fit = MIN(MAX(te, te_min), te_max)

    CALL evaluate_log_polynomial(te_fit, coefficients, log_cooling, log_slope)
    log_cooling = log_cooling + log_slope*(LOG(te) - LOG(te_fit))

    IF (log_cooling > 6.d0*LOG(10.d0) + log_energy_rate_scale) THEN
       WRITE(6,*) 'Impurity cooling evaluation exceeded its valid range at Te = ', te
       STOP
    END IF

    IF (log_cooling < -100.d0 + log_energy_rate_scale) THEN
       cooling = 0.d0
       dcooling_dte = 0.d0
    ELSE
       cooling = EXP(log_cooling)/1.d6
       dcooling_dte = cooling*log_slope/te
    END IF
  END SUBROUTINE evaluate_cooling_polynomial

  SUBROUTINE evaluate_log_polynomial(x, coefficients, value, derivative)
    REAL*8, INTENT(IN) :: x, coefficients(:)
    REAL*8, INTENT(OUT) :: value, derivative
    REAL*8 :: log_x
    INTEGER :: i

    log_x = LOG(x)
    value = 0.d0
    derivative = 0.d0
    DO i = 1, SIZE(coefficients)
       value = value + coefficients(i)*log_x**(i-1)
    END DO
    DO i = 2, SIZE(coefficients)
       derivative = derivative + (i-1)*coefficients(i)*log_x**(i-2)
    END DO
  END SUBROUTINE evaluate_log_polynomial

  SUBROUTINE shift_log_polynomial(coefficients, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: coefficients(:)
    REAL*8, INTENT(IN) :: shift_x, log_scale
    REAL*8 :: shifted(SIZE(coefficients)), shift_factor
    INTEGER :: i, j, shift_power

    shifted = 0.d0
    DO i = 1, SIZE(coefficients)
       DO j = 1, i
          shift_power = i-j
          shift_factor = 1.d0
          IF (shift_power > 0) shift_factor = shift_x**shift_power
          shifted(j) = shifted(j) + coefficients(i)*binomial_coefficient(i-1,j-1)*shift_factor
       END DO
    END DO
    shifted(1) = shifted(1) + log_scale
    coefficients = shifted
  END SUBROUTINE shift_log_polynomial

  REAL*8 FUNCTION binomial_coefficient(n, k)
    INTEGER, INTENT(IN) :: n, k
    INTEGER :: i, kk

    IF (k < 0 .OR. k > n) THEN
       binomial_coefficient = 0.d0
       RETURN
    END IF
    IF (k == 0 .OR. k == n) THEN
       binomial_coefficient = 1.d0
       RETURN
    END IF

    kk = MIN(k, n-k)
    binomial_coefficient = 1.d0
    DO i = 1, kk
       binomial_coefficient = binomial_coefficient*DBLE(n-kk+i)/DBLE(i)
    END DO
  END FUNCTION binomial_coefficient

END MODULE impurity_radiation_model
