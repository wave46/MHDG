SUBMODULE (transport_models_1d) transport_models_1d_runtime
USE magnetic_geometry_state, ONLY: magnetic_geometry_cache
USE magnetic_topology, ONLY: magnetic_region_core, magnetic_region_main_sol, magnetic_region_private_flux
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_interp_transport(this, rho, chi_i, chi_e, d, nu_mom, vpinch)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho
    REAL*8, INTENT(OUT) :: chi_i, chi_e, d, nu_mom, vpinch

    chi_i = 0.d0
    chi_e = 0.d0
    d = 0.d0
    nu_mom = 0.d0
    vpinch = 0.d0

    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN
    IF (rho <= model_tol) RETURN

    CALL tm1d_interp_profile(this, rho, this%chi_i_fs, chi_i)
    CALL tm1d_interp_profile(this, rho, this%chi_e_fs, chi_e)
    CALL tm1d_interp_profile(this, rho, this%d_fs, d)
    CALL tm1d_interp_profile(this, rho, this%nu_mom_fs, nu_mom)
    CALL tm1d_interp_profile(this, rho, this%vpinch_fs, vpinch)
  END SUBROUTINE tm1d_interp_transport

  MODULE SUBROUTINE tm1d_apply_1D_diffusion(this, rho_pol_norm, diff_iso, &
       diff_ani, region)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: rho_pol_norm(:)
    REAL*8, INTENT(INOUT) :: diff_iso(:, :, :), diff_ani(:, :, :)
    INTEGER, INTENT(IN) :: region(:)
    INTEGER :: g, equation
    REAL*8 :: rho_g, w, chi_i, chi_e, d, nu_mom, vpinch
    REAL*8 :: band_start, plateau_start, plateau_end, band_end, coefficient

    IF (.NOT. this%is_initialized) RETURN
    IF (SIZE(diff_iso, 3) /= SIZE(rho_pol_norm)) RETURN
    IF (SIZE(diff_ani, 3) /= SIZE(rho_pol_norm)) RETURN
    IF (SIZE(region) /= SIZE(rho_pol_norm)) RETURN

    DO g = 1, SIZE(rho_pol_norm)
       IF (.NOT. tm1d_region_is_included(&
            this%config%transport_region_policy, region(g))) CYCLE
       rho_g = MAX(rho_pol_norm(g), 0.d0)
       w = tm1d_blend_weight(this, rho_g)
       IF (w <= 0.d0) CYCLE

       CALL this%interp_transport(rho_g, chi_i, chi_e, d, nu_mom, vpinch)

       d = MAX(d, this%config%diff_n_min)
       nu_mom = MAX(nu_mom, this%config%diff_u_min)
       chi_i = MAX(chi_i, this%config%diff_e_min)
       chi_e = MAX(chi_e, this%config%diff_ee_min)

       diff_iso(1, 1, g) = diff_iso(1, 1, g) + w*(d - diff_iso(1, 1, g))
       diff_iso(2, 2, g) = diff_iso(2, 2, g) + w*(nu_mom - diff_iso(2, 2, g))
       diff_iso(3, 3, g) = diff_iso(3, 3, g) + w*(chi_i - diff_iso(3, 3, g))
       diff_iso(4, 4, g) = diff_iso(4, 4, g) + w*(chi_e - diff_iso(4, 4, g))

       diff_ani(1, 1, g) = diff_ani(1, 1, g) + w*(d - diff_ani(1, 1, g))
       diff_ani(2, 2, g) = diff_ani(2, 2, g) + w*(nu_mom - diff_ani(2, 2, g))
       diff_ani(3, 3, g) = diff_ani(3, 3, g) + w*(chi_i - diff_ani(3, 3, g))
       diff_ani(4, 4, g) = diff_ani(4, 4, g) + w*(chi_e - diff_ani(4, 4, g))
    END DO

    ! The band also acts outside the background model's cutoff/topology policy.
    IF (.NOT. ANY(this%config%band_enabled)) RETURN
    IF (.NOT. magnetic_geometry_cache%is_initialized) &
         ERROR STOP 'A transport band requires the magnetic geometry cache'
    CALL this%band_edges(band_start, plateau_start, plateau_end, band_end)
    DO g = 1, SIZE(rho_pol_norm)
       SELECT CASE (region(g))
       CASE (magnetic_region_core, magnetic_region_main_sol, magnetic_region_private_flux)
          w = smoothstep01((rho_pol_norm(g)-band_start)/(plateau_start-band_start)) * &
               smoothstep01((band_end-rho_pol_norm(g))/(band_end-plateau_end))
          IF (w <= 0.d0) CYCLE
          DO equation = 1, 4
             IF (.NOT. this%config%band_enabled(equation)) CYCLE
             coefficient = (1.d0-w)*diff_ani(equation,equation,g) + w*this%config%band_diff(equation)
             ! Preserve the residual isotropic contribution, including at W=1.
             diff_iso(equation,equation,g) = &
                  (diff_iso(equation,equation,g)-diff_ani(equation,equation,g)) + coefficient
             diff_ani(equation,equation,g) = coefficient
          ENDDO
       END SELECT
    ENDDO
  END SUBROUTINE tm1d_apply_1D_diffusion

  MODULE SUBROUTINE tm1d_compute_1D_pinch_matrix(this, b, rho, APinch, &
       region, outward_normal)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: b(:), rho
    REAL*8, INTENT(OUT) :: APinch(:,:)
    INTEGER, INTENT(IN) :: region
    REAL*8, INTENT(IN) :: outward_normal(:)
    REAL*8 :: vpinch, velocity(2), pinch_weight
    REAL*8 :: chi_i, chi_e, d, nu_mom
    INTEGER :: equation

    APinch = 0.d0
    IF (.NOT. this%is_initialized) RETURN
    IF (.NOT. tm1d_region_is_included(&
         this%config%transport_region_policy, region)) RETURN

    pinch_weight = tm1d_pinch_window(this, rho)
    IF (pinch_weight <= model_tol) RETURN

    CALL this%interp_transport(rho, chi_i, chi_e, d, nu_mom, vpinch)
    vpinch = pinch_weight*vpinch
    IF (ABS(vpinch) <= model_tol) RETURN

    CALL tm1d_build_pinch_velocity(vpinch, &
         this%config%transport_region_policy, b, outward_normal, velocity)
    DO equation = 1, SIZE(this%config%pinch_enabled)
       IF (this%config%pinch_enabled(equation)) APinch(equation,1:2) = velocity
    ENDDO
  END SUBROUTINE tm1d_compute_1D_pinch_matrix

END SUBMODULE transport_models_1d_runtime
