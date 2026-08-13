!**********************************************************************
! Central lifecycle and storage for balance diagnostics.
!
! All rates use the physical sign convention "positive into the domain".
! The packed layout and its MPI reduction are owned here so later balance
! categories can be appended without changing the solver driver.
!**********************************************************************
MODULE balance_diagnostics
  USE MPI_OMP
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: balance_mode_invalid = -1
  INTEGER, PARAMETER, PUBLIC :: balance_mode_off = 0
  INTEGER, PARAMETER, PUBLIC :: balance_mode_summary = 1
  INTEGER, PARAMETER, PUBLIC :: balance_mode_detailed = 2

  ! Particle content and conservation terms.
  INTEGER, PARAMETER, PUBLIC :: balance_particle_content_plasma = 1
  INTEGER, PARAMETER, PUBLIC :: balance_particle_content_neutral = 2
  INTEGER, PARAMETER, PUBLIC :: balance_particle_temporal_plasma = 3
  INTEGER, PARAMETER, PUBLIC :: balance_particle_temporal_neutral = 4
  INTEGER, PARAMETER, PUBLIC :: balance_particle_volume_plasma = 5
  INTEGER, PARAMETER, PUBLIC :: balance_particle_volume_neutral = 6
  INTEGER, PARAMETER, PUBLIC :: balance_particle_boundary_physical_plasma = 7
  INTEGER, PARAMETER, PUBLIC :: balance_particle_boundary_physical_neutral = 8
  INTEGER, PARAMETER, PUBLIC :: balance_particle_boundary_tau_plasma = 9
  INTEGER, PARAMETER, PUBLIC :: balance_particle_boundary_tau_neutral = 10

  ! Detailed particle volume components.
  INTEGER, PARAMETER, PUBLIC :: balance_particle_plasma_ionization = 11
  INTEGER, PARAMETER, PUBLIC :: balance_particle_plasma_recombination = 12
  INTEGER, PARAMETER, PUBLIC :: balance_particle_plasma_other_volume_source = 13
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_ionization = 14
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_recombination = 15
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_other_volume_source = 16
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_puff_volume_source = 17
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_pump_volume_source = 18
  INTEGER, PARAMETER, PUBLIC :: balance_particle_charge_exchange = 19

  ! Detailed physical boundary components.
  INTEGER, PARAMETER, PUBLIC :: balance_particle_plasma_parallel_flux = 20
  INTEGER, PARAMETER, PUBLIC :: balance_particle_plasma_diffusion_flux = 21
  INTEGER, PARAMETER, PUBLIC :: balance_particle_plasma_pinch_flux = 22
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_diffusion_flux = 23
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_pressure_flux = 24
  INTEGER, PARAMETER, PUBLIC :: balance_particle_neutral_convection_flux = 25
  INTEGER, PARAMETER :: balance_particle_term_count = 25

  ! Wall/HDG closure terms. Plasma recycling components are recorded as
  ! neutral wall sources; they are not independent global volume sources.
  INTEGER, PARAMETER, PUBLIC :: balance_wall_plasma_diffusion_flux = 1
  INTEGER, PARAMETER, PUBLIC :: balance_wall_plasma_tau_flux = 2
  INTEGER, PARAMETER, PUBLIC :: balance_wall_recycling_parallel_source = 3
  INTEGER, PARAMETER, PUBLIC :: balance_wall_recycling_diffusion_source = 4
  INTEGER, PARAMETER, PUBLIC :: balance_wall_recycling_pinch_source = 5
  INTEGER, PARAMETER, PUBLIC :: balance_wall_puff_source = 6
  INTEGER, PARAMETER, PUBLIC :: balance_wall_pump_source = 7
  INTEGER, PARAMETER, PUBLIC :: balance_wall_neutral_diffusion_flux = 8
  INTEGER, PARAMETER, PUBLIC :: balance_wall_neutral_pressure_flux = 9
  INTEGER, PARAMETER, PUBLIC :: balance_wall_neutral_convection_flux = 10
  INTEGER, PARAMETER, PUBLIC :: balance_wall_neutral_tau_flux = 11
  INTEGER, PARAMETER :: balance_wall_term_count = 11

  INTEGER, PARAMETER, PUBLIC :: balance_particle_offset = 0
  INTEGER, PARAMETER, PUBLIC :: balance_wall_offset = balance_particle_term_count
  INTEGER, PARAMETER, PUBLIC :: balance_diagnostics_value_count = &
       balance_particle_term_count + balance_wall_term_count

  INTEGER, SAVE :: configured_mode = balance_mode_off
  LOGICAL, SAVE :: finalized_values_valid = .FALSE.
  REAL*8, SAVE :: local_values(balance_diagnostics_value_count) = 0.d0
  REAL*8, SAVE :: finalized_values(balance_diagnostics_value_count) = 0.d0

  PUBLIC :: parse_balance_diagnostics_mode
  PUBLIC :: balance_diagnostics_mode_name
  PUBLIC :: balance_diagnostics_configure
  PUBLIC :: balance_diagnostics_get_mode
  PUBLIC :: balance_diagnostics_enabled
  PUBLIC :: balance_diagnostics_detailed
  PUBLIC :: balance_diagnostics_reset
  PUBLIC :: balance_diagnostics_merge
  PUBLIC :: balance_diagnostics_finalize
  PUBLIC :: balance_diagnostics_has_values
  PUBLIC :: balance_diagnostics_copy_values
  PUBLIC :: balance_particle_index
  PUBLIC :: balance_wall_index

CONTAINS

  INTEGER FUNCTION parse_balance_diagnostics_mode(name) RESULT(mode)
    CHARACTER(LEN=*), INTENT(IN) :: name

    SELECT CASE (TRIM(ADJUSTL(name)))
    CASE ('off')
       mode = balance_mode_off
    CASE ('summary')
       mode = balance_mode_summary
    CASE ('detailed')
       mode = balance_mode_detailed
    CASE DEFAULT
       mode = balance_mode_invalid
    END SELECT
  END FUNCTION parse_balance_diagnostics_mode

  CHARACTER(LEN=8) FUNCTION balance_diagnostics_mode_name(mode) RESULT(name)
    INTEGER, INTENT(IN) :: mode

    SELECT CASE (mode)
    CASE (balance_mode_off)
       name = 'off'
    CASE (balance_mode_summary)
       name = 'summary'
    CASE (balance_mode_detailed)
       name = 'detailed'
    CASE DEFAULT
       name = 'invalid'
    END SELECT
  END FUNCTION balance_diagnostics_mode_name

  SUBROUTINE balance_diagnostics_configure(mode)
    INTEGER, INTENT(IN) :: mode

    IF (mode < balance_mode_off .OR. mode > balance_mode_detailed) THEN
       WRITE (6, *) 'Invalid balance diagnostics mode id: ', mode
       STOP
    ENDIF
    configured_mode = mode
    finalized_values_valid = .FALSE.
    IF (configured_mode /= balance_mode_off) THEN
       local_values = 0.d0
       finalized_values = 0.d0
    ENDIF
  END SUBROUTINE balance_diagnostics_configure

  INTEGER FUNCTION balance_diagnostics_get_mode() RESULT(mode)
    mode = configured_mode
  END FUNCTION balance_diagnostics_get_mode

  LOGICAL FUNCTION balance_diagnostics_enabled() RESULT(enabled)
    enabled = configured_mode /= balance_mode_off
  END FUNCTION balance_diagnostics_enabled

  LOGICAL FUNCTION balance_diagnostics_detailed() RESULT(detailed)
    detailed = configured_mode == balance_mode_detailed
  END FUNCTION balance_diagnostics_detailed

  SUBROUTINE balance_diagnostics_reset()
    IF (.NOT. balance_diagnostics_enabled()) RETURN

    local_values = 0.d0
    finalized_values = 0.d0
    finalized_values_valid = .FALSE.
  END SUBROUTINE balance_diagnostics_reset

  SUBROUTINE balance_diagnostics_merge(values)
    REAL*8, INTENT(IN) :: values(balance_diagnostics_value_count)

    ! Call once outside the OpenMP parallel loop. The assembly caller first
    ! combines its thread-private arrays into this rank-local contribution.
    IF (.NOT. balance_diagnostics_enabled()) RETURN
    local_values = local_values + values
    finalized_values_valid = .FALSE.
  END SUBROUTINE balance_diagnostics_merge

  SUBROUTINE balance_diagnostics_finalize()
#ifdef PARALL
    INTEGER :: ierr
#endif

    IF (.NOT. balance_diagnostics_enabled()) RETURN

    finalized_values = local_values
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, finalized_values, &
         balance_diagnostics_value_count, MPI_REAL8, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
#endif
    finalized_values_valid = .TRUE.
  END SUBROUTINE balance_diagnostics_finalize

  LOGICAL FUNCTION balance_diagnostics_has_values() RESULT(has_values)
    has_values = balance_diagnostics_enabled() .AND. finalized_values_valid
  END FUNCTION balance_diagnostics_has_values

  SUBROUTINE balance_diagnostics_copy_values(values)
    REAL*8, INTENT(OUT) :: values(balance_diagnostics_value_count)

    IF (.NOT. balance_diagnostics_has_values()) RETURN
    values = finalized_values
  END SUBROUTINE balance_diagnostics_copy_values

  INTEGER FUNCTION balance_particle_index(term) RESULT(index)
    INTEGER, INTENT(IN) :: term

    index = balance_particle_offset + term
  END FUNCTION balance_particle_index

  INTEGER FUNCTION balance_wall_index(term) RESULT(index)
    INTEGER, INTENT(IN) :: term

    index = balance_wall_offset + term
  END FUNCTION balance_wall_index

END MODULE balance_diagnostics
