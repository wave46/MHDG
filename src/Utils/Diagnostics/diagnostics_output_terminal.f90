SUBMODULE (diagnostics) diagnostics_output_terminal
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE diag_print_boundary_hdg_summary(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_boundary_summary_type) :: summary

    summary = diag_boundary_summary(this)

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,'(A)') '----------------------------------------'
       WRITE(6,'(A)') 'Boundary HDG / BC diagnostics'
       WRITE(6,'(A)') 'neutral closure ['//TRIM(diag_particle_integral_units())//']'
       WRITE(6,'(A,1X,ES11.3)') '  residual           :', summary%residual
       WRITE(6,'(A,1X,ES11.3)') '  neutral flux inward:', summary%neutral_flux
       WRITE(6,'(A,1X,ES11.3)') '  recycled source    :', summary%recycled_source
       WRITE(6,'(A,1X,ES11.3)') '  external sources   :', summary%wall_flux
       WRITE(6,'(A,1X,ES11.3)') '  tau numerical      :', summary%tau_flux
       WRITE(6,'(A)') '  residual = recycled + wall + tau - neutral flux'
       WRITE(6,'(A)') '----------------------------------------'
    ENDIF
  END SUBROUTINE diag_print_boundary_hdg_summary

  MODULE SUBROUTINE diag_print_particle_summary(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_particle_summary_type) :: summary

    summary = diag_particle_summary(this)

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,'(A)') '----------------------------------------'
       WRITE(6,'(A)') 'Physical particle diagnostics'
       WRITE(6,'(A)') 'balance ['//TRIM(diag_particle_integral_units())//']'
       WRITE(6,'(A,1X,ES11.3)') '  plasma :', summary%plasma_balance
       WRITE(6,'(A,1X,ES11.3)') '  neutral:', summary%neutral_balance
       WRITE(6,'(A,1X,ES11.3)') '  total  :', summary%total_balance
       WRITE(6,'(A)') ''
       WRITE(6,'(A)') 'content ['//TRIM(diag_content_units(diag_content_total_particles))//']'
       WRITE(6,'(A,1X,ES11.3)') '  plasma :', summary%plasma_content
       WRITE(6,'(A,1X,ES11.3)') '  neutral:', summary%neutral_content
       WRITE(6,'(A,1X,ES11.3)') '  total  :', summary%total_content
       WRITE(6,'(A)') '----------------------------------------'
    ENDIF
  END SUBROUTINE diag_print_particle_summary

  MODULE SUBROUTINE diag_print_boundary_hdg_detail(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_boundary_summary_type) :: summary

    IF (MPIvar%glob_id .NE. 0) RETURN

    summary = diag_boundary_summary(this)

    WRITE(6,'(A)') '----------------------------------------'
    WRITE(6,'(A)') 'Boundary HDG / BC diagnostics'
    WRITE(6,'(A)') 'neutral closure ['//TRIM(diag_particle_integral_units())//']'
    WRITE(6,'(A,1X,ES11.3)') '  residual           :', summary%residual
    WRITE(6,'(A,1X,ES11.3)') '  neutral flux inward:', summary%neutral_flux
    WRITE(6,'(A,1X,ES11.3)') '  recycled source    :', summary%recycled_source
    WRITE(6,'(A,1X,ES11.3)') '  external sources   :', summary%wall_flux
    WRITE(6,'(A,1X,ES11.3)') '  tau numerical      :', summary%tau_flux
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'components, inward-positive ['//TRIM(diag_particle_integral_units())//']'
    WRITE(6,'(A,1X,ES11.3)') '  recycled plasma:', summary%recycled_source
    WRITE(6,'(A,1X,ES11.3)') '    parallel convection    :', this%boundary_hdg(diag_boundary_plasma_parallel_flux)
    WRITE(6,'(A,1X,ES11.3)') '    perpendicular diffusion:', this%boundary_hdg(diag_boundary_plasma_diffusion_flux)
    WRITE(6,'(A,1X,ES11.3)') '    perpendicular pinch    :', this%boundary_hdg(diag_boundary_plasma_pinch_flux)
    WRITE(6,'(A)') ''
    WRITE(6,'(A,1X,ES11.3)') '  neutral total:', summary%neutral_flux
    WRITE(6,'(A,1X,ES11.3)') '    density diffusion  :', this%boundary_hdg(diag_boundary_neutral_diffusion_flux)
    WRITE(6,'(A,1X,ES11.3)') '    pressure diffusion :', this%boundary_hdg(diag_boundary_neutral_pressure_flux)
    WRITE(6,'(A,1X,ES11.3)') '    parallel convection:', this%boundary_hdg(diag_boundary_neutral_convection_flux)
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '  external sources:'
    WRITE(6,'(A,1X,ES11.3)') '    puff:', this%boundary_hdg(diag_boundary_wall_puff_source)
    WRITE(6,'(A,1X,ES11.3)') '    pump:', this%boundary_hdg(diag_boundary_wall_pump_sink)
    WRITE(6,'(A)') ''
    WRITE(6,'(A,1X,ES11.3)') '  tau numerical:', summary%tau_flux
    WRITE(6,'(A)') '----------------------------------------'
  END SUBROUTINE diag_print_boundary_hdg_detail

  MODULE SUBROUTINE diag_print_particle_detail(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_particle_summary_type) :: summary

    IF (MPIvar%glob_id .NE. 0) RETURN

    summary = diag_particle_summary(this)

    WRITE(6,'(A)') '----------------------------------------'
    WRITE(6,'(A)') 'Physical particle diagnostics'
    WRITE(6,'(A)') 'balance ['//TRIM(diag_particle_integral_units())//']'
    WRITE(6,'(A,1X,ES11.3)') '  plasma :', summary%plasma_balance
    WRITE(6,'(A,1X,ES11.3)') '  neutral:', summary%neutral_balance
    WRITE(6,'(A,1X,ES11.3)') '  total  :', summary%total_balance
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'content ['//TRIM(diag_content_units(diag_content_total_particles))//']'
    WRITE(6,'(A,1X,ES11.3)') '  plasma :', summary%plasma_content
    WRITE(6,'(A,1X,ES11.3)') '  neutral:', summary%neutral_content
    WRITE(6,'(A,1X,ES11.3)') '  total  :', summary%total_content
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') 'components, inward-positive ['//TRIM(diag_particle_integral_units())//']'
    WRITE(6,'(A)') '  plasma:'
    WRITE(6,'(A,1X,ES11.3)') '    ionization   :', this%particles(diag_particle_plasma_ionization)
    WRITE(6,'(A,1X,ES11.3)') '    recombination:', this%particles(diag_particle_plasma_recombination)
    WRITE(6,'(A,1X,ES11.3)') '    boundary flux:', this%particles(diag_particle_plasma_boundary_flux)
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '  neutral:'
    WRITE(6,'(A,1X,ES11.3)') '    ionization   :', this%particles(diag_particle_neutral_ionization)
    WRITE(6,'(A,1X,ES11.3)') '    recombination:', this%particles(diag_particle_neutral_recombination)
    WRITE(6,'(A,1X,ES11.3)') '    boundary flux:', this%particles(diag_particle_neutral_boundary_flux)
    WRITE(6,'(A)') ''
    WRITE(6,'(A)') '  external sources:'
    WRITE(6,'(A,1X,ES11.3)') '    puff:', this%particles(diag_particle_puff_source)
    WRITE(6,'(A,1X,ES11.3)') '    pump:', this%particles(diag_particle_pump_sink)
    WRITE(6,'(A,1X,ES11.3)') '    recycling:', this%particles(diag_particle_recycling_source)
    WRITE(6,'(A,1X,ES11.3)') '      parallel convection    :', this%particles(diag_particle_recycling_parallel_source)
    WRITE(6,'(A,1X,ES11.3)') '      perpendicular diffusion:', this%particles(diag_particle_recycling_diffusion_source)
    WRITE(6,'(A,1X,ES11.3)') '      perpendicular pinch    :', this%particles(diag_particle_recycling_pinch_source)
    WRITE(6,'(A)') '----------------------------------------'
  END SUBROUTINE diag_print_particle_detail

  MODULE FUNCTION diag_particle_integral_units() RESULT(units)
    CHARACTER(LEN=32) :: units

    units = 'particles/s'
  END FUNCTION diag_particle_integral_units

  MODULE FUNCTION diag_content_units(term_id) RESULT(units)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=32) :: units

    SELECT CASE (term_id)
    CASE (diag_content_plasma_particles, diag_content_neutral_particles, diag_content_total_particles)
       units = 'particles'
    CASE (diag_content_total_energy)
       units = 'not_implemented'
    CASE DEFAULT
       units = 'unknown'
    END SELECT
  END FUNCTION diag_content_units

END SUBMODULE diagnostics_output_terminal
