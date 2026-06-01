SUBMODULE (diagnostics) diagnostics_output_terminal
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE diag_print_summary(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_boundary_summary_type) :: boundary
    TYPE(diag_particle_summary_type) :: particles

    boundary = diag_boundary_summary(this)
    particles = diag_particle_summary(this)

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,'(A)') '--- Balance diagnostics summary ---'
       WRITE(6,'(A)') '  particle balances ['//TRIM(diag_particle_integral_units())//']'
       WRITE(6,'(A,1X,ES16.8)') 'boundary HDG / BC check =', boundary%residual
       WRITE(6,'(A,1X,ES16.8)') 'plasma particle balance =', particles%plasma_balance
       WRITE(6,'(A,1X,ES16.8)') 'neutral particle balance =', particles%neutral_balance
       WRITE(6,'(A,1X,ES16.8)') 'total particle balance =', particles%total_balance
       WRITE(6,'(A)') '  particle content ['//TRIM(diag_content_units(diag_content_total_particles))//']'
       WRITE(6,'(A,1X,ES16.8)') 'plasma particle content =', particles%plasma_content
       WRITE(6,'(A,1X,ES16.8)') 'neutral particle content =', particles%neutral_content
       WRITE(6,'(A,1X,ES16.8)') 'total particle content =', particles%total_content
    ENDIF
  END SUBROUTINE diag_print_summary

  MODULE SUBROUTINE diag_print_boundary_hdg_summary(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    TYPE(diag_boundary_summary_type) :: summary

    summary = diag_boundary_summary(this)

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,'(A)') '----------------------------------------'
       WRITE(6,'(A)') 'Boundary HDG / BC diagnostics ['//TRIM(diag_particle_integral_units())//']'
       WRITE(6,'(A,1X,ES11.3)') '  neutral closure residual:', summary%residual
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
       WRITE(6,'(A)') 'Physical particle diagnostics ['//TRIM(diag_particle_integral_units())//']'
       WRITE(6,'(A,1X,ES11.3)') '  plasma balance:', summary%plasma_balance
       WRITE(6,'(A,1X,ES11.3)') '  neutral balance:', summary%neutral_balance
       WRITE(6,'(A,1X,ES11.3)') '  total balance  :', summary%total_balance
       WRITE(6,'(A)') '  particle content ['//TRIM(diag_content_units(diag_content_total_particles))//']'
       WRITE(6,'(A,1X,ES11.3)') '  plasma content :', summary%plasma_content
       WRITE(6,'(A,1X,ES11.3)') '  neutral content:', summary%neutral_content
       WRITE(6,'(A,1X,ES11.3)') '  total content  :', summary%total_content
       WRITE(6,'(A)') '----------------------------------------'
    ENDIF
  END SUBROUTINE diag_print_particle_summary

  MODULE SUBROUTINE diag_print_detail(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER :: i

    IF (MPIvar%glob_id .NE. 0) RETURN

    WRITE(6,'(A)') '--- Boundary HDG / BC components ['//TRIM(diag_particle_integral_units())//'] ---'
    DO i = 1, diag_boundary_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(diag_boundary_label(i))//' =', this%boundary_hdg(i)
    ENDDO
    WRITE(6,'(A)') '--- Physical particle components ['//TRIM(diag_particle_integral_units())//'] ---'
    DO i = 1, diag_particle_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(diag_particle_label(i))//' =', this%particles(i)
    ENDDO
    WRITE(6,'(A)') '--- Global content components ['//TRIM(diag_content_units(diag_content_total_particles))//'] ---'
    DO i = 1, diag_content_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(diag_content_label(i))//' =', diag_content_value(this, i)
    ENDDO
    WRITE(6,'(A)') '--- Reserved energy components ---'
    DO i = 1, diag_energy_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(diag_energy_label(i))//' =', this%energy(i)
    ENDDO
  END SUBROUTINE diag_print_detail

  MODULE SUBROUTINE diag_print_boundary_hdg_detail(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER :: i
    TYPE(diag_boundary_summary_type) :: summary

    IF (MPIvar%glob_id .NE. 0) RETURN

    summary = diag_boundary_summary(this)

    WRITE(6,'(A)') '----------------------------------------'
    WRITE(6,'(A)') 'Boundary HDG / BC diagnostics ['//TRIM(diag_particle_integral_units())//']'
    WRITE(6,'(A,1X,ES11.3)') '  neutral closure residual:', summary%residual
    WRITE(6,'(A,1X,ES11.3)') '  neutral flux into domain:', summary%neutral_flux
    WRITE(6,'(A,1X,ES11.3)') '  recycled source         :', summary%recycled_source
    WRITE(6,'(A,1X,ES11.3)') '  wall sources/sinks      :', summary%wall_flux
    WRITE(6,'(A,1X,ES11.3)') '  tau correction          :', summary%tau_flux
    WRITE(6,'(A)') '  components, inward-positive:'
    DO i = 1, diag_boundary_term_count
       WRITE(6,'(A,1X,ES11.3)') '    '//TRIM(diag_boundary_label(i))//':', this%boundary_hdg(i)
    ENDDO
    WRITE(6,'(A)') '----------------------------------------'
  END SUBROUTINE diag_print_boundary_hdg_detail

  MODULE SUBROUTINE diag_print_particle_detail(this)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER :: i
    TYPE(diag_particle_summary_type) :: summary

    IF (MPIvar%glob_id .NE. 0) RETURN

    summary = diag_particle_summary(this)

    WRITE(6,'(A)') '----------------------------------------'
    WRITE(6,'(A)') 'Physical particle diagnostics ['//TRIM(diag_particle_integral_units())//']'
    WRITE(6,'(A,1X,ES11.3)') '  plasma balance:', summary%plasma_balance
    WRITE(6,'(A,1X,ES11.3)') '  neutral balance:', summary%neutral_balance
    WRITE(6,'(A,1X,ES11.3)') '  total balance  :', summary%total_balance
    WRITE(6,'(A,1X,ES11.3)') '  plasma content :', summary%plasma_content
    WRITE(6,'(A,1X,ES11.3)') '  neutral content:', summary%neutral_content
    WRITE(6,'(A,1X,ES11.3)') '  total content  :', summary%total_content
    WRITE(6,'(A)') '  particle components, inward-positive:'
    DO i = 1, diag_particle_term_count
       WRITE(6,'(A,1X,ES11.3)') '    '//TRIM(diag_particle_label(i))//':', this%particles(i)
    ENDDO
    WRITE(6,'(A)') '  particle content ['//TRIM(diag_content_units(diag_content_total_particles))//']:'
    DO i = 1, diag_content_term_count
       IF (i .EQ. diag_content_total_energy) CYCLE
       WRITE(6,'(A,1X,ES11.3)') '    '//TRIM(diag_content_label(i))//':', diag_content_value(this, i)
    ENDDO
    WRITE(6,'(A)') '----------------------------------------'
  END SUBROUTINE diag_print_particle_detail

  MODULE FUNCTION diag_boundary_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (diag_boundary_plasma_parallel_flux)
       label = 'recycled plasma parallel source'
    CASE (diag_boundary_plasma_diffusion_flux)
       label = 'recycled plasma diffusion source'
    CASE (diag_boundary_plasma_pinch_flux)
       label = 'recycled plasma pinch source'
    CASE (diag_boundary_neutral_diffusion_flux)
       label = 'neutral diffusion boundary flux'
    CASE (diag_boundary_neutral_pressure_flux)
       label = 'neutral pressure boundary flux'
    CASE (diag_boundary_neutral_convection_flux)
       label = 'neutral convection boundary flux'
    CASE (diag_boundary_neutral_total_flux)
       label = 'neutral total boundary flux'
    CASE (diag_boundary_tau_numerical_flux)
       label = 'tau numerical boundary flux'
    CASE (diag_boundary_wall_puff_source)
       label = 'wall puff source'
    CASE (diag_boundary_wall_pump_sink)
       label = 'wall pump sink'
    CASE DEFAULT
       label = 'boundary HDG unknown term'
    END SELECT
  END FUNCTION diag_boundary_label

  MODULE FUNCTION diag_particle_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (diag_particle_plasma_ionization)
       label = 'particle plasma ionization'
    CASE (diag_particle_plasma_recombination)
       label = 'particle plasma recombination'
    CASE (diag_particle_neutral_ionization)
       label = 'particle neutral ionization'
    CASE (diag_particle_neutral_recombination)
       label = 'particle neutral recombination'
    CASE (diag_particle_puff_source)
       label = 'particle puff source'
    CASE (diag_particle_pump_sink)
       label = 'particle pump sink'
    CASE (diag_particle_recycling_source)
       label = 'particle recycling source'
    CASE (diag_particle_plasma_boundary_flux)
       label = 'particle plasma boundary flux'
    CASE (diag_particle_neutral_boundary_flux)
       label = 'particle neutral boundary flux'
    CASE (diag_particle_charge_exchange_rate)
       label = 'particle charge exchange rate'
    CASE DEFAULT
       label = 'particle unknown term'
    END SELECT
  END FUNCTION diag_particle_label

  MODULE FUNCTION diag_content_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (diag_content_plasma_particles)
       label = 'content plasma particles'
    CASE (diag_content_neutral_particles)
       label = 'content neutral particles'
    CASE (diag_content_total_particles)
       label = 'content total particles'
    CASE (diag_content_total_energy)
       label = 'content total energy'
    CASE DEFAULT
       label = 'content unknown term'
    END SELECT
  END FUNCTION diag_content_label

  MODULE FUNCTION diag_energy_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (diag_energy_ion_balance)
       label = 'ion energy balance'
    CASE (diag_energy_electron_balance)
       label = 'electron energy balance'
    CASE (diag_energy_neutral_balance)
       label = 'neutral energy balance'
    CASE (diag_energy_total_balance)
       label = 'total energy balance'
    CASE DEFAULT
       label = 'energy unknown term'
    END SELECT
  END FUNCTION diag_energy_label

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
