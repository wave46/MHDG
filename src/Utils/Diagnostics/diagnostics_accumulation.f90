SUBMODULE (diagnostics) diagnostics_accumulation
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE diag_reset(this)
    CLASS(diagnostics_type), INTENT(INOUT) :: this

    CALL this%reset_boundary_hdg()
    CALL this%reset_particles_content()
    this%energy = 0.d0
  END SUBROUTINE diag_reset

  MODULE SUBROUTINE diag_reset_boundary_hdg(this)
    CLASS(diagnostics_type), INTENT(INOUT) :: this

    this%boundary_hdg = 0.d0
  END SUBROUTINE diag_reset_boundary_hdg

  MODULE SUBROUTINE diag_reset_particles_content(this)
    CLASS(diagnostics_type), INTENT(INOUT) :: this

    this%particles = 0.d0
    this%content = 0.d0
  END SUBROUTINE diag_reset_particles_content

  MODULE SUBROUTINE diag_add(this, category_id, term_id, value)
    CLASS(diagnostics_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: category_id, term_id
    REAL*8, INTENT(IN) :: value

    SELECT CASE (category_id)
    CASE (diag_category_boundary_hdg)
       IF (term_id >= 1 .AND. term_id <= diag_boundary_term_count) THEN
          this%boundary_hdg(term_id) = this%boundary_hdg(term_id) + value
       ENDIF
    CASE (diag_category_particles)
       IF (term_id >= 1 .AND. term_id <= diag_particle_term_count) THEN
          this%particles(term_id) = this%particles(term_id) + value
       ENDIF
    CASE (diag_category_content)
       IF (term_id >= 1 .AND. term_id <= diag_content_term_count) THEN
          this%content(term_id) = this%content(term_id) + value
       ENDIF
    CASE (diag_category_energy)
       IF (term_id >= 1 .AND. term_id <= diag_energy_term_count) THEN
          this%energy(term_id) = this%energy(term_id) + value
       ENDIF
    END SELECT
  END SUBROUTINE diag_add

  MODULE SUBROUTINE diag_account_boundary_hdg(this, plasma_parallel_flux, plasma_diffusion_flux, plasma_pinch_flux, &
       &neutral_diffusion_flux, neutral_pressure_flux, neutral_convection_flux, neutral_total_flux, tau_numerical_flux, &
       &wall_puff_source, wall_pump_sink, include_wall_sources)
    CLASS(diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_parallel_flux, plasma_diffusion_flux, plasma_pinch_flux
    REAL*8, INTENT(IN) :: neutral_diffusion_flux, neutral_pressure_flux, neutral_convection_flux, neutral_total_flux
    REAL*8, INTENT(IN) :: tau_numerical_flux, wall_puff_source, wall_pump_sink
    LOGICAL, INTENT(IN) :: include_wall_sources

    CALL this%add(diag_category_boundary_hdg, diag_boundary_plasma_parallel_flux, plasma_parallel_flux)
    CALL this%add(diag_category_boundary_hdg, diag_boundary_plasma_diffusion_flux, plasma_diffusion_flux)
    CALL this%add(diag_category_boundary_hdg, diag_boundary_plasma_pinch_flux, plasma_pinch_flux)
    CALL this%add(diag_category_boundary_hdg, diag_boundary_neutral_diffusion_flux, neutral_diffusion_flux)
    CALL this%add(diag_category_boundary_hdg, diag_boundary_neutral_pressure_flux, neutral_pressure_flux)
    CALL this%add(diag_category_boundary_hdg, diag_boundary_neutral_convection_flux, neutral_convection_flux)
    CALL this%add(diag_category_boundary_hdg, diag_boundary_neutral_total_flux, neutral_total_flux)
    CALL this%add(diag_category_boundary_hdg, diag_boundary_tau_numerical_flux, tau_numerical_flux)
    IF (include_wall_sources) THEN
       CALL this%add(diag_category_boundary_hdg, diag_boundary_wall_puff_source, wall_puff_source)
       CALL this%add(diag_category_boundary_hdg, diag_boundary_wall_pump_sink, wall_pump_sink)
    ENDIF
  END SUBROUTINE diag_account_boundary_hdg

  MODULE SUBROUTINE diag_account_volume_particle_reactions(this, ionization_rate, recombination_rate, &
       &charge_exchange_rate)
    CLASS(diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: ionization_rate, recombination_rate, charge_exchange_rate

    CALL this%add(diag_category_particles, diag_particle_plasma_ionization, ionization_rate)
    CALL this%add(diag_category_particles, diag_particle_neutral_ionization, -ionization_rate)
    CALL this%add(diag_category_particles, diag_particle_plasma_recombination, -recombination_rate)
    CALL this%add(diag_category_particles, diag_particle_neutral_recombination, recombination_rate)
    CALL this%add(diag_category_particles, diag_particle_charge_exchange_rate, charge_exchange_rate)
  END SUBROUTINE diag_account_volume_particle_reactions

  MODULE SUBROUTINE diag_account_boundary_particles(this, plasma_boundary_flux, neutral_boundary_flux, &
       &recycling_source, wall_puff_source, wall_pump_sink, include_wall_sources)
    CLASS(diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_boundary_flux, neutral_boundary_flux, recycling_source
    REAL*8, INTENT(IN) :: wall_puff_source, wall_pump_sink
    LOGICAL, INTENT(IN) :: include_wall_sources

    CALL this%add(diag_category_particles, diag_particle_plasma_boundary_flux, plasma_boundary_flux)
    CALL this%add(diag_category_particles, diag_particle_neutral_boundary_flux, neutral_boundary_flux)
    CALL this%add(diag_category_particles, diag_particle_recycling_source, recycling_source)
    IF (include_wall_sources) THEN
       CALL this%add(diag_category_particles, diag_particle_puff_source, wall_puff_source)
       CALL this%add(diag_category_particles, diag_particle_pump_sink, wall_pump_sink)
    ENDIF
  END SUBROUTINE diag_account_boundary_particles

  MODULE SUBROUTINE diag_account_wall_particle_sources(this, puff_source, pump_sink)
    CLASS(diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: puff_source, pump_sink

    CALL this%add(diag_category_particles, diag_particle_puff_source, puff_source)
    CALL this%add(diag_category_particles, diag_particle_pump_sink, pump_sink)
  END SUBROUTINE diag_account_wall_particle_sources

  MODULE SUBROUTINE diag_account_particle_content(this, plasma_particles, neutral_particles)
    CLASS(diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_particles, neutral_particles

    CALL this%add(diag_category_content, diag_content_plasma_particles, plasma_particles)
    CALL this%add(diag_category_content, diag_content_neutral_particles, neutral_particles)
    CALL this%add(diag_category_content, diag_content_total_particles, plasma_particles + neutral_particles)
  END SUBROUTINE diag_account_particle_content

  MODULE FUNCTION diag_get(this, category_id, term_id) RESULT(value)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: category_id, term_id
    REAL*8 :: value

    value = 0.d0
    SELECT CASE (category_id)
    CASE (diag_category_boundary_hdg)
       IF (term_id >= 1 .AND. term_id <= diag_boundary_term_count) value = this%boundary_hdg(term_id)
    CASE (diag_category_particles)
       IF (term_id >= 1 .AND. term_id <= diag_particle_term_count) value = this%particles(term_id)
    CASE (diag_category_content)
       IF (term_id >= 1 .AND. term_id <= diag_content_term_count) value = this%content(term_id)
    CASE (diag_category_energy)
       IF (term_id >= 1 .AND. term_id <= diag_energy_term_count) value = this%energy(term_id)
    END SELECT
  END FUNCTION diag_get

END SUBMODULE diagnostics_accumulation
