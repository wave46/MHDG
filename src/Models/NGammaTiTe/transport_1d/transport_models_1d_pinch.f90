SUBMODULE (transport_models_1d) transport_models_1d_pinch
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_compute_pinch_profile(this, work)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(transport_model_workspace_t), INTENT(IN) :: work

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    SELECT CASE (this%pinch_model)
    CASE (1)
       CALL tm1d_compute_militello_pinch(this, work)
    CASE (2)
       CALL tm1d_compute_geometric_pinch(this, work)
    CASE (3)
       CALL tm1d_compute_constant_pinch(this, work)
    CASE DEFAULT
       this%vpinch_fs = 0.d0
    END SELECT
  END SUBROUTINE tm1d_compute_pinch_profile

  MODULE SUBROUTINE tm1d_compute_militello_pinch(this, work)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(transport_model_workspace_t), INTENT(IN) :: work
    REAL*8 :: pinch_factor_militello_fs(this%nrho)

    pinch_factor_militello_fs = MIN(1.d0, EXP(1.d0 - work%nuestar_fs/MAX(this%nu_th, model_tol)))
    this%vpinch_fs = pinch_factor_militello_fs * this%c_pinch * this%d_part_fs * work%rmin_fs / MAX(this%a_minor, model_tol)**2

    WHERE (work%rmin_fs <= model_tol)
       this%vpinch_fs = 0.d0
    END WHERE
  END SUBROUTINE tm1d_compute_militello_pinch

  MODULE SUBROUTINE tm1d_compute_geometric_pinch(this, work)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(transport_model_workspace_t), INTENT(IN) :: work

    this%vpinch_fs = this%c_pinch * this%d_part_fs * work%rmin_fs / MAX(this%a_minor, model_tol)**2

    WHERE (work%rmin_fs <= model_tol)
       this%vpinch_fs = 0.d0
    END WHERE
  END SUBROUTINE tm1d_compute_geometric_pinch

  MODULE SUBROUTINE tm1d_compute_constant_pinch(this, work)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    TYPE(transport_model_workspace_t), INTENT(IN) :: work

    this%vpinch_fs = this%vpinch_const

    WHERE (work%rmin_fs <= model_tol)
       this%vpinch_fs = 0.d0
    END WHERE
  END SUBROUTINE tm1d_compute_constant_pinch

  MODULE REAL*8 FUNCTION tm1d_pinch_window(this, rho)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho

    tm1d_pinch_window = tm1d_axis_ramp(rho, this%rho_pinch_axis_width) * &
         tm1d_edge_cutoff(rho, this%rho_pinch_model_max, this%rho_pinch_edge_width)
  END FUNCTION tm1d_pinch_window

END SUBMODULE transport_models_1d_pinch
