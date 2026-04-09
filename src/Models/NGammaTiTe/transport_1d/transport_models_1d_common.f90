SUBMODULE (transport_models_1d) transport_models_1d_common
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

  MODULE REAL*8 FUNCTION tm1d_blend_weight(this, rho)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho

    tm1d_blend_weight = tm1d_edge_cutoff(rho, this%config%rho_diffusion_model_max, this%config%rho_blend_width)
  END FUNCTION tm1d_blend_weight

  MODULE SUBROUTINE tm1d_interp_profile(this, rho, profile, value)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(IN) :: profile(:)
    REAL*8, INTENT(OUT) :: value
    INTEGER :: ilow, ihigh
    REAL*8 :: alpha
    REAL*8 :: rho_grid(this%nrho)

    value = 0.d0

    IF (this%nrho <= 0) RETURN
    IF (SIZE(profile) /= this%nrho) RETURN
    IF (this%nrho == 1) THEN
      value = profile(1)
      RETURN
    END IF

    IF (.NOT. ALLOCATED(this%rho_grid)) RETURN
    rho_grid = this%rho_grid

    CALL find_cell_and_local_coordinate(this%nrho, rho_grid, MIN(MAX(rho, rho_grid(1)), rho_grid(this%nrho)), ilow, alpha)
    ilow = MIN(MAX(ilow, 1), this%nrho)
    ihigh = MIN(ilow + 1, this%nrho)

    IF (ihigh == ilow) THEN
      value = profile(ilow)
    ELSE
      value = (1.d0 - alpha)*profile(ilow) + alpha*profile(ihigh)
    END IF
  END SUBROUTINE tm1d_interp_profile

END SUBMODULE transport_models_1d_common
