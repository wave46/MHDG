SUBMODULE (transport_models_1d) transport_models_1d_runtime
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_interp_transport(this, rho, chi_i, chi_e, d_part, nu_mom, vpinch)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(OUT) :: chi_i, chi_e, d_part, nu_mom, vpinch

    chi_i = 0.d0
    chi_e = 0.d0
    d_part = 0.d0
    nu_mom = 0.d0
    vpinch = 0.d0

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN
    IF (rho <= model_tol) RETURN

    CALL tm1d_interp_profile(this, rho, this%chi_i_fs, chi_i)
    CALL tm1d_interp_profile(this, rho, this%chi_e_fs, chi_e)
    CALL tm1d_interp_profile(this, rho, this%d_part_fs, d_part)
    CALL tm1d_interp_profile(this, rho, this%nu_mom_fs, nu_mom)
    CALL tm1d_interp_profile(this, rho, this%vpinch_fs, vpinch)
  END SUBROUTINE tm1d_interp_transport

  MODULE SUBROUTINE tm1d_apply_1D_diffusion(this, rho_pol_norm, diff_iso, diff_ani)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho_pol_norm(:)
    REAL*8, INTENT(INOUT) :: diff_iso(:, :, :), diff_ani(:, :, :)
    INTEGER :: g
    REAL*8 :: rho_g, w, chi_i, chi_e, d_part, nu_mom, vpinch

    IF (.NOT. this%is_initialized) RETURN
    IF (SIZE(diff_iso, 3) /= SIZE(rho_pol_norm)) RETURN
    IF (SIZE(diff_ani, 3) /= SIZE(rho_pol_norm)) RETURN

    DO g = 1, SIZE(rho_pol_norm)
       rho_g = MAX(rho_pol_norm(g), 0.d0)
       w = tm1d_blend_weight(this, rho_g)
       IF (w <= 0.d0) CYCLE

       CALL this%interp_transport(rho_g, chi_i, chi_e, d_part, nu_mom, vpinch)

       d_part = MAX(d_part, this%diff_n_min)
       nu_mom = MAX(nu_mom, this%diff_u_min)
       chi_i = MAX(chi_i, this%diff_e_min)
       chi_e = MAX(chi_e, this%diff_ee_min)

       diff_iso(1, 1, g) = diff_iso(1, 1, g) + w*(d_part - diff_iso(1, 1, g))
       diff_iso(2, 2, g) = diff_iso(2, 2, g) + w*(nu_mom - diff_iso(2, 2, g))
       diff_iso(3, 3, g) = diff_iso(3, 3, g) + w*(chi_i - diff_iso(3, 3, g))
       diff_iso(4, 4, g) = diff_iso(4, 4, g) + w*(chi_e - diff_iso(4, 4, g))

       diff_ani(1, 1, g) = diff_ani(1, 1, g) + w*(d_part - diff_ani(1, 1, g))
       diff_ani(2, 2, g) = diff_ani(2, 2, g) + w*(nu_mom - diff_ani(2, 2, g))
       diff_ani(3, 3, g) = diff_ani(3, 3, g) + w*(chi_i - diff_ani(3, 3, g))
       diff_ani(4, 4, g) = diff_ani(4, 4, g) + w*(chi_e - diff_ani(4, 4, g))
    END DO
  END SUBROUTINE tm1d_apply_1D_diffusion

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

END SUBMODULE transport_models_1d_runtime
