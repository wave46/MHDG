MODULE impurity_radiation_model
  USE globals, ONLY: phys, simpar
  USE globals, ONLY: input
  USE MPI_OMP, ONLY: MPIvar
  USE shift_logpoly, ONLY: shift_logpoly_1d_17
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_impurity_radiation_input
  PUBLIC :: init_impurity_radiation_model
  PUBLIC :: adimensionalize_impurity_radiation_model
  PUBLIC :: compute_impurity_weighted_cooling
  PUBLIC :: compute_impurity_weighted_dcooling_dU

  INTEGER, PARAMETER :: cooling_coeff_count = 17
  INTEGER, PARAMETER :: max_impurity_entries = 64

CONTAINS

  SUBROUTINE read_impurity_radiation_input()
    INTEGER :: uimpurity, ios, i
    INTEGER :: n_impurities
    CHARACTER(LEN=20) :: impurity_names(max_impurity_entries)
    REAL*8 :: impurity_concentrations(max_impurity_entries)
    NAMELIST /IMPURITY_RADIATION_LST/ n_impurities, impurity_names, impurity_concentrations

    IF (LEN_TRIM(input%impurity_model_path) == 0) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'impurity_model_path must be set when impurity_radiation is enabled'
       STOP
    END IF

    n_impurities = 0
    impurity_names = ''
    impurity_concentrations = 0.d0

    uimpurity = 102
    OPEN(uimpurity, file=TRIM(ADJUSTL(input%impurity_model_path)), status='old', iostat=ios)
    IF (ios /= 0) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'Could not open impurity radiation settings file: ', TRIM(ADJUSTL(input%impurity_model_path))
       STOP
    END IF

    READ(uimpurity, nml=IMPURITY_RADIATION_LST, iostat=ios)
    CLOSE(uimpurity)
    IF (ios /= 0) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'Could not read IMPURITY_RADIATION_LST from file: ', TRIM(ADJUSTL(input%impurity_model_path))
       STOP
    END IF

    IF (n_impurities < 1) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'n_impurities must be at least 1 for impurity radiation'
       STOP
    END IF
    IF (n_impurities > max_impurity_entries) THEN
       IF (MPIvar%glob_id == 0) WRITE(6,*) 'n_impurities exceeds maximum supported entries: ', n_impurities, max_impurity_entries
       STOP
    END IF

    DO i = 1, n_impurities
       impurity_names(i) = TRIM(ADJUSTL(impurity_names(i)))
       SELECT CASE (TRIM(impurity_names(i)))
       CASE ('N', 'W')
       CASE DEFAULT
          IF (MPIvar%glob_id == 0) WRITE(6,*) 'Unsupported impurity radiation species: ', TRIM(impurity_names(i))
          STOP
       END SELECT

       IF (impurity_concentrations(i) < 0.d0) THEN
          IF (MPIvar%glob_id == 0) WRITE(6,*) 'Impurity concentration must be non-negative for species ', TRIM(impurity_names(i)), ': ', impurity_concentrations(i)
          STOP
       END IF
    END DO

    IF (ALLOCATED(phys%impurity_names)) DEALLOCATE(phys%impurity_names)
    IF (ALLOCATED(phys%impurity_concentrations)) DEALLOCATE(phys%impurity_concentrations)

    phys%n_impurities = n_impurities
    ALLOCATE(phys%impurity_names(n_impurities))
    ALLOCATE(phys%impurity_concentrations(n_impurities))
    phys%impurity_names = impurity_names(1:n_impurities)
    phys%impurity_concentrations = impurity_concentrations(1:n_impurities)
  END SUBROUTINE read_impurity_radiation_input

  SUBROUTINE init_impurity_radiation_model()
    INTEGER :: i

    IF (phys%n_impurities < 1) THEN
       WRITE(6,*) 'Impurity radiation requires n_impurities >= 1'
       STOP
    END IF
    IF (.NOT. ALLOCATED(phys%impurity_names) .OR. &
       &.NOT. ALLOCATED(phys%impurity_concentrations)) THEN
       WRITE(6,*) 'Impurity radiation mixture was not initialized from input'
       STOP
    END IF

    IF (ALLOCATED(phys%alpha_cooling_factor_impurities)) DEALLOCATE(phys%alpha_cooling_factor_impurities)
    ALLOCATE(phys%alpha_cooling_factor_impurities(cooling_coeff_count, phys%n_impurities))

    DO i = 1, phys%n_impurities
       CALL load_impurity_cooling_coefficients(phys%impurity_names(i), phys%alpha_cooling_factor_impurities(:, i))
    END DO

  END SUBROUTINE init_impurity_radiation_model

  SUBROUTINE load_impurity_cooling_coefficients(impurity_name, coeffs)
    CHARACTER(LEN=*), INTENT(IN) :: impurity_name
    REAL*8, INTENT(OUT) :: coeffs(:)

    SELECT CASE (TRIM(ADJUSTL(impurity_name)))
    CASE ('N')
       coeffs = (/-2.49348163e+01,  9.52628451e+00, -4.39511346e+00,  2.00446916e+00,&
       -2.27166819e+00,  9.95587194e-01,  1.11209167e+00, -1.13089140e+00,&
        2.22902131e-01,  1.37491696e-01, -9.74320916e-02,  2.88337222e-02,&
       -5.03419501e-03,  5.54252119e-04, -3.80153983e-05,  1.49061146e-06,&
       -2.56122449e-08/)
    CASE ('W')
       coeffs = (/-3.29798537e+01,  3.38045169e+01, -3.81202398e+01,  2.47450333e+01,&
       -9.04921504e+00,  1.91417658e+00, -2.32405139e-01,  1.46022765e-02,&
       -2.46177645e-04, -1.89159754e-05,  7.73790432e-07,  0.d0,&
       0.d0,  0.d0, 0.d0,  0.d0,&
       0.d0/)
    CASE DEFAULT
       WRITE(6,*) 'Warning: cooling factor not defined for impurity ', TRIM(ADJUSTL(impurity_name))
       STOP
    END SELECT
  END SUBROUTINE load_impurity_cooling_coefficients

  SUBROUTINE adimensionalize_impurity_radiation_model()
    REAL*8 :: log_temp_shift, log_energy_rate_scale
    INTEGER :: i

    IF (.NOT. ALLOCATED(phys%alpha_cooling_factor_impurities)) RETURN

    log_temp_shift = LOG(simpar%refval_temperature)
    log_energy_rate_scale = LOG(simpar%refval_density*simpar%refval_time* &
      &simpar%refval_charge/simpar%refval_mass*simpar%refval_time**2/simpar%refval_length**2)

    DO i = 1, phys%n_impurities
       CALL shift_logpoly_1d_17(phys%alpha_cooling_factor_impurities(:, i), log_temp_shift, log_energy_rate_scale)
    END DO
  END SUBROUTINE adimensionalize_impurity_radiation_model

  SUBROUTINE compute_impurity_weighted_cooling(cooling_factors, cooling_total)
    REAL*8, INTENT(IN) :: cooling_factors(:)
    REAL*8, INTENT(OUT) :: cooling_total

    cooling_total = SUM(phys%impurity_concentrations(1:phys%n_impurities)*cooling_factors(1:phys%n_impurities))
  END SUBROUTINE compute_impurity_weighted_cooling

  SUBROUTINE compute_impurity_weighted_dcooling_dU(dcooling_factors_dU, dcooling_total_dU)
    REAL*8, INTENT(IN) :: dcooling_factors_dU(:, :)
    REAL*8, INTENT(OUT) :: dcooling_total_dU(:)
    INTEGER :: i

    dcooling_total_dU = 0.d0
    DO i = 1, phys%n_impurities
       dcooling_total_dU = dcooling_total_dU + phys%impurity_concentrations(i)*dcooling_factors_dU(:, i)
    END DO
  END SUBROUTINE compute_impurity_weighted_dcooling_dU

END MODULE impurity_radiation_model
