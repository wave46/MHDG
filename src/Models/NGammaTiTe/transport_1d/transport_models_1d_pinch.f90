SUBMODULE (transport_models_1d) transport_models_1d_pinch
IMPLICIT NONE

CONTAINS

  MODULE REAL*8 FUNCTION tm1d_smoothstep01(s)
    REAL*8, INTENT(IN) :: s
    REAL*8 :: sc

    sc = MAX(MIN(s, 1.d0), 0.d0)
    tm1d_smoothstep01 = sc*sc*(3.d0 - 2.d0*sc)
  END FUNCTION tm1d_smoothstep01

  MODULE REAL*8 FUNCTION tm1d_axis_ramp(rho, width)
    REAL*8, INTENT(IN) :: rho, width

    IF (rho <= 0.d0) THEN
       tm1d_axis_ramp = 0.d0
    ELSEIF (width <= model_tol .OR. rho >= width) THEN
       tm1d_axis_ramp = 1.d0
    ELSE
       tm1d_axis_ramp = tm1d_smoothstep01(rho/width)
    END IF
  END FUNCTION tm1d_axis_ramp

  MODULE REAL*8 FUNCTION tm1d_edge_cutoff(rho, rho_max, width)
    REAL*8, INTENT(IN) :: rho, rho_max, width
    REAL*8 :: rho_start

    IF (rho >= rho_max) THEN
       tm1d_edge_cutoff = 0.d0
    ELSEIF (width <= model_tol) THEN
       tm1d_edge_cutoff = 1.d0
    ELSE
       rho_start = rho_max - width
       IF (rho <= rho_start) THEN
          tm1d_edge_cutoff = 1.d0
       ELSE
          tm1d_edge_cutoff = 1.d0 - tm1d_smoothstep01((rho - rho_start)/width)
       END IF
    END IF
  END FUNCTION tm1d_edge_cutoff

  MODULE SUBROUTINE tm1d_compute_pinch_profile(this, nuestar_fs)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: nuestar_fs(:)

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    SELECT CASE (this%pinch_model)
    CASE (1)
       CALL tm1d_compute_militello_pinch(this, nuestar_fs, this%vpinch_fs)
    CASE (2)
       CALL tm1d_compute_geometric_pinch(this, this%vpinch_fs)
    CASE (3)
       CALL tm1d_compute_constant_pinch(this, this%vpinch_fs)
    CASE DEFAULT
       this%vpinch_fs = 0.d0
    END SELECT
  END SUBROUTINE tm1d_compute_pinch_profile

  MODULE SUBROUTINE tm1d_compute_militello_pinch(this, nuestar_fs, vpinch_fs)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: nuestar_fs(:)
    REAL*8, INTENT(OUT) :: vpinch_fs(:)
    REAL*8 :: pinch_factor_militello_fs(this%nrho)

    pinch_factor_militello_fs = MIN(1.d0, EXP(1.d0 - nuestar_fs/MAX(this%nu_th, model_tol)))
    vpinch_fs = pinch_factor_militello_fs * this%c_pinch * this%d_part_fs * this%rmin_fs / MAX(this%a_minor, model_tol)**2

    WHERE (this%rmin_fs <= model_tol)
       vpinch_fs = 0.d0
    END WHERE
  END SUBROUTINE tm1d_compute_militello_pinch

  MODULE SUBROUTINE tm1d_compute_geometric_pinch(this, vpinch_fs)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(OUT) :: vpinch_fs(:)

    vpinch_fs = this%c_pinch * this%d_part_fs * this%rmin_fs / MAX(this%a_minor, model_tol)**2

    WHERE (this%rmin_fs <= model_tol)
       vpinch_fs = 0.d0
    END WHERE
  END SUBROUTINE tm1d_compute_geometric_pinch

  MODULE SUBROUTINE tm1d_compute_constant_pinch(this, vpinch_fs)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(OUT) :: vpinch_fs(:)

    vpinch_fs = this%vpinch_const

    WHERE (this%rmin_fs <= model_tol)
       vpinch_fs = 0.d0
    END WHERE
  END SUBROUTINE tm1d_compute_constant_pinch

  MODULE SUBROUTINE tm1d_compute_1D_pinch_matrix(this, b, rho, APinch)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: b(:), rho
    REAL*8, INTENT(OUT) :: APinch(:,:)
    REAL*8 :: vpinch, bnorm(2), bnorm_norm, pinch_weight
    REAL*8 :: chi_i, chi_e, d_part, nu_mom

    APinch = 0.d0
    IF (.NOT. this%is_initialized) RETURN

    pinch_weight = tm1d_pinch_window(this, rho)
    IF (pinch_weight <= model_tol) RETURN

    CALL this%interp_transport(rho, chi_i, chi_e, d_part, nu_mom, vpinch)
    vpinch = pinch_weight*vpinch
    IF (ABS(vpinch) <= model_tol) RETURN

    bnorm = b(1:2)
    bnorm_norm = NORM2(bnorm)
    IF (bnorm_norm <= model_tol) RETURN
    bnorm = bnorm/bnorm_norm

    APinch(1,1) = vpinch*bnorm(2)
    APinch(1,2) = -vpinch*bnorm(1)
  END SUBROUTINE tm1d_compute_1D_pinch_matrix

  MODULE REAL*8 FUNCTION tm1d_pinch_window(this, rho)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho

    tm1d_pinch_window = tm1d_axis_ramp(rho, this%rho_pinch_axis_width) * &
         tm1d_edge_cutoff(rho, this%rho_pinch_model_max, this%rho_pinch_edge_width)
  END FUNCTION tm1d_pinch_window

END SUBMODULE transport_models_1d_pinch
