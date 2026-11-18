SUBMODULE (transport_models_1d) transport_models_1d_common
USE magnetic_geometry_state, ONLY: magnetic_equilibrium, magnetic_geometry_cache
USE MPI_OMP, ONLY: MPIvar
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_prepare_band(this)
    CLASS(transport_model_1d_t), INTENT(INOUT) :: this
    REAL*8 :: band_start, plateau_start, plateau_end, band_end, ray_radius, r_sep, span
    LOGICAL :: found

    IF (.NOT. ANY(this%config%band_enabled)) RETURN
    IF (this%config%diffusion_coordinate /= 'omp_distance') RETURN
    IF (.NOT. magnetic_equilibrium%is_initialized .OR. .NOT. magnetic_geometry_cache%is_initialized) &
         ERROR STOP 'An OMP transport band requires the magnetic equilibrium and geometry cache'
    IF (this%band_geometry_generation == magnetic_equilibrium%generation) RETURN
    IF (.NOT. ieee_is_finite(phys%lscale) .OR. phys%lscale <= 0.d0) &
         ERROR STOP 'OMP band requires a finite positive equilibrium length scale'
    span = magnetic_equilibrium%psi_lcfs - magnetic_equilibrium%psi_axis
    IF (.NOT. ieee_is_finite(span) .OR. span == 0.d0) ERROR STOP 'OMP band requires a valid flux span'

    CALL magnetic_equilibrium%closed_surface_radius(magnetic_equilibrium%psi_lcfs, 0.d0, ray_radius, found)
    IF (.NOT. found) ERROR STOP 'Could not locate the OMP separatrix for the transport band'
    r_sep = magnetic_equilibrium%r_axis + ray_radius
    CALL tm1d_band_edges(this%config, band_start, plateau_start, plateau_end, band_end)
    CALL omp_distance_to_rho(band_start, r_sep, this%band_start_rho)
    CALL omp_distance_to_rho(plateau_start, r_sep, this%plateau_start_rho)
    CALL omp_distance_to_rho(this%config%band_center, r_sep, this%band_center_rho)
    CALL omp_distance_to_rho(plateau_end, r_sep, this%plateau_end_rho)
    CALL omp_distance_to_rho(band_end, r_sep, this%band_end_rho)
    IF (this%band_start_rho >= this%plateau_start_rho .OR. &
         this%plateau_start_rho >= this%band_center_rho .OR. &
         this%band_center_rho >= this%plateau_end_rho .OR. this%plateau_end_rho >= this%band_end_rho) &
         ERROR STOP 'Mapped OMP band boundaries must be strictly increasing'
    this%band_geometry_generation = magnetic_equilibrium%generation

    IF (MPIvar%glob_id == 0) THEN
       WRITE(6,'(A)') 'Transport band: OMP distance [m] -> rho_pol_norm'
       WRITE(6,'(A,ES12.4,F10.4)') '  band start:    ', band_start, this%band_start_rho
       WRITE(6,'(A,ES12.4,F10.4)') '  plateau start: ', plateau_start, this%plateau_start_rho
       WRITE(6,'(A,ES12.4,F10.4)') '  center:        ', this%config%band_center, this%band_center_rho
       WRITE(6,'(A,ES12.4,F10.4)') '  plateau end:   ', plateau_end, this%plateau_end_rho
       WRITE(6,'(A,ES12.4,F10.4)') '  band end:      ', band_end, this%band_end_rho
    ENDIF
  END SUBROUTINE tm1d_prepare_band

  MODULE SUBROUTINE tm1d_band_rho_edges(this, band_start, plateau_start, plateau_end, band_end, center)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    REAL*8, INTENT(OUT) :: band_start, plateau_start, plateau_end, band_end
    REAL*8, INTENT(OUT), OPTIONAL :: center

    IF (this%config%diffusion_coordinate == 'rho_pol_norm') THEN
       CALL tm1d_band_edges(this%config, band_start, plateau_start, plateau_end, band_end)
       IF (PRESENT(center)) center = this%config%band_center
    ELSE
       IF (.NOT. magnetic_equilibrium%is_initialized .OR. &
            this%band_geometry_generation /= magnetic_equilibrium%generation) &
            ERROR STOP 'OMP transport band must be prepared for the current equilibrium'
       band_start = this%band_start_rho
       plateau_start = this%plateau_start_rho
       plateau_end = this%plateau_end_rho
       band_end = this%band_end_rho
       IF (PRESENT(center)) center = this%band_center_rho
    ENDIF
  END SUBROUTINE tm1d_band_rho_edges

  SUBROUTINE omp_distance_to_rho(distance, r_sep, rho)
    REAL*8, INTENT(IN) :: distance, r_sep
    REAL*8, INTENT(OUT) :: rho
    REAL*8 :: r, psi, psi_n, dpsi_dr, dpsi_dz, span

    r = r_sep + distance/phys%lscale
    IF (.NOT. ieee_is_finite(r) .OR. r < magnetic_equilibrium%r_axis .OR. &
         r > magnetic_equilibrium%r(magnetic_equilibrium%nr)) &
         ERROR STOP 'OMP band boundary lies outside the outboard equilibrium grid'
    CALL magnetic_equilibrium%evaluate_flux(r, magnetic_equilibrium%z_axis, psi, dpsi_dr, dpsi_dz)
    span = magnetic_equilibrium%psi_lcfs - magnetic_equilibrium%psi_axis
    psi_n = (psi-magnetic_equilibrium%psi_axis)/span
    IF (.NOT. ieee_is_finite(psi_n) .OR. psi_n < 0.d0 .OR. &
         dpsi_dr/span < 0.d0) ERROR STOP 'OMP band boundary is outside the increasing normalized-flux branch'
    rho = SQRT(psi_n)
  END SUBROUTINE omp_distance_to_rho

  MODULE REAL*8 FUNCTION tm1d_axis_ramp(rho, width)
    REAL*8, INTENT(IN) :: rho, width

    IF (rho <= 0.d0) THEN
       tm1d_axis_ramp = 0.d0
    ELSEIF (width <= model_tol .OR. rho >= width) THEN
       tm1d_axis_ramp = 1.d0
    ELSE
       tm1d_axis_ramp = smoothstep01(rho/width)
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
          tm1d_edge_cutoff = 1.d0 - smoothstep01((rho - rho_start)/width)
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
