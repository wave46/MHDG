SUBMODULE (diagnostics) diagnostics_quantities
  IMPLICIT NONE

CONTAINS

  MODULE FUNCTION diag_boundary_hdg_check(this) RESULT(value)
    CLASS(diagnostics_type), INTENT(IN) :: this
    REAL*8 :: value
    TYPE(diag_boundary_summary_type) :: summary

    summary = diag_boundary_summary(this)
    value = summary%residual
  END FUNCTION diag_boundary_hdg_check

  MODULE FUNCTION diag_boundary_summary(this) RESULT(summary)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_boundary_summary_type) :: summary

    summary%recycled_source = this%boundary_hdg(diag_boundary_plasma_parallel_flux) &
       &+ this%boundary_hdg(diag_boundary_plasma_diffusion_flux) &
       &+ this%boundary_hdg(diag_boundary_plasma_pinch_flux)
    summary%neutral_flux = this%boundary_hdg(diag_boundary_neutral_total_flux)
    summary%wall_flux = this%boundary_hdg(diag_boundary_wall_puff_source) &
       &+ this%boundary_hdg(diag_boundary_wall_pump_sink)
    summary%tau_flux = this%boundary_hdg(diag_boundary_tau_numerical_flux)
    summary%residual = summary%recycled_source - summary%neutral_flux + summary%wall_flux + summary%tau_flux
  END FUNCTION diag_boundary_summary

  MODULE FUNCTION diag_plasma_particle_balance(this) RESULT(value)
    CLASS(diagnostics_type), INTENT(IN) :: this
    REAL*8 :: value

    value = this%particles(diag_particle_plasma_ionization) &
       &+ this%particles(diag_particle_plasma_recombination) &
       &+ this%particles(diag_particle_plasma_boundary_flux)
  END FUNCTION diag_plasma_particle_balance

  MODULE FUNCTION diag_neutral_particle_balance(this) RESULT(value)
    CLASS(diagnostics_type), INTENT(IN) :: this
    REAL*8 :: value

    value = this%particles(diag_particle_neutral_ionization) &
       &+ this%particles(diag_particle_neutral_recombination) &
       &+ this%particles(diag_particle_puff_source) &
       &+ this%particles(diag_particle_pump_sink) &
       &+ this%particles(diag_particle_neutral_boundary_flux)
  END FUNCTION diag_neutral_particle_balance

  MODULE FUNCTION diag_particle_summary(this) RESULT(summary)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_particle_summary_type) :: summary

    summary%plasma_balance = diag_plasma_particle_balance(this)
    summary%neutral_balance = diag_neutral_particle_balance(this)
    summary%total_balance = summary%plasma_balance + summary%neutral_balance
    summary%plasma_content = this%content(diag_content_plasma_particles)
    summary%neutral_content = this%content(diag_content_neutral_particles)
    summary%total_content = this%content(diag_content_total_particles)
  END FUNCTION diag_particle_summary

END SUBMODULE diagnostics_quantities
