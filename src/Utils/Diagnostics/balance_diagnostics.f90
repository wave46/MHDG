MODULE balance_diagnostics
  USE MPI_OMP
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_boundary_hdg = 1
  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_particles = 2
  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_content = 3
  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_energy = 4

  ! All accumulated terms use the same sign convention:
  ! positive is directed into the domain, negative is loss from the domain.
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_plasma_parallel_flux = 1
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_plasma_diffusion_flux = 2
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_plasma_pinch_flux = 3
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_neutral_diffusion_flux = 4
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_neutral_pressure_flux = 5
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_neutral_convection_flux = 6
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_neutral_total_flux = 7
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_tau_numerical_flux = 8
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_wall_puff_source = 9
  INTEGER, PARAMETER, PUBLIC :: balance_diag_boundary_wall_pump_sink = 10
  INTEGER, PARAMETER :: balance_diag_boundary_term_count = 10

  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_plasma_ionization = 1
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_plasma_recombination = 2
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_plasma_boundary_flux = 3
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_neutral_ionization = 4
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_neutral_recombination = 5
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_puff_source = 6
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_pump_sink = 7
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_recycling_source = 8
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_neutral_boundary_flux = 9
  INTEGER, PARAMETER, PUBLIC :: balance_diag_particle_charge_exchange_rate = 10
  INTEGER, PARAMETER :: balance_diag_particle_term_count = 10

  INTEGER, PARAMETER, PUBLIC :: balance_diag_content_plasma_particles = 1
  INTEGER, PARAMETER, PUBLIC :: balance_diag_content_neutral_particles = 2
  INTEGER, PARAMETER, PUBLIC :: balance_diag_content_total_particles = 3
  INTEGER, PARAMETER, PUBLIC :: balance_diag_content_total_energy = 4
  INTEGER, PARAMETER :: balance_diag_content_term_count = 4

  INTEGER, PARAMETER, PUBLIC :: balance_diag_energy_ion_balance = 1
  INTEGER, PARAMETER, PUBLIC :: balance_diag_energy_electron_balance = 2
  INTEGER, PARAMETER, PUBLIC :: balance_diag_energy_neutral_balance = 3
  INTEGER, PARAMETER, PUBLIC :: balance_diag_energy_total_balance = 4
  INTEGER, PARAMETER :: balance_diag_energy_term_count = 4

  TYPE, PUBLIC :: balance_diagnostics_type
     REAL*8 :: boundary_hdg(balance_diag_boundary_term_count) = 0.d0
     REAL*8 :: particles(balance_diag_particle_term_count) = 0.d0
     REAL*8 :: content(balance_diag_content_term_count) = 0.d0
     REAL*8 :: energy(balance_diag_energy_term_count) = 0.d0
   CONTAINS
     PROCEDURE :: reset => balance_diag_reset
     PROCEDURE :: add => balance_diag_add
     PROCEDURE :: get => balance_diag_get
     PROCEDURE :: mpi_reduce => balance_diag_mpi_reduce
     PROCEDURE :: print_summary => balance_diag_print_summary
     PROCEDURE :: print_detail => balance_diag_print_detail
     PROCEDURE :: write_hdf5 => balance_diag_write_hdf5
  END TYPE balance_diagnostics_type

  TYPE(balance_diagnostics_type), PUBLIC :: balance_diag

CONTAINS

  SUBROUTINE balance_diag_reset(this)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this

    this%boundary_hdg = 0.d0
    this%particles = 0.d0
    this%content = 0.d0
    this%energy = 0.d0
  END SUBROUTINE balance_diag_reset

  SUBROUTINE balance_diag_add(this, category_id, term_id, value)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: category_id, term_id
    REAL*8, INTENT(IN) :: value

    SELECT CASE (category_id)
    CASE (balance_diag_category_boundary_hdg)
       IF (term_id >= 1 .AND. term_id <= balance_diag_boundary_term_count) THEN
          this%boundary_hdg(term_id) = this%boundary_hdg(term_id) + value
       ENDIF
    CASE (balance_diag_category_particles)
       IF (term_id >= 1 .AND. term_id <= balance_diag_particle_term_count) THEN
          this%particles(term_id) = this%particles(term_id) + value
       ENDIF
    CASE (balance_diag_category_content)
       IF (term_id >= 1 .AND. term_id <= balance_diag_content_term_count) THEN
          this%content(term_id) = this%content(term_id) + value
       ENDIF
    CASE (balance_diag_category_energy)
       IF (term_id >= 1 .AND. term_id <= balance_diag_energy_term_count) THEN
          this%energy(term_id) = this%energy(term_id) + value
       ENDIF
    END SELECT
  END SUBROUTINE balance_diag_add

  FUNCTION balance_diag_get(this, category_id, term_id) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: category_id, term_id
    REAL*8 :: value

    value = 0.d0
    SELECT CASE (category_id)
    CASE (balance_diag_category_boundary_hdg)
       IF (term_id >= 1 .AND. term_id <= balance_diag_boundary_term_count) value = this%boundary_hdg(term_id)
    CASE (balance_diag_category_particles)
       IF (term_id >= 1 .AND. term_id <= balance_diag_particle_term_count) value = this%particles(term_id)
    CASE (balance_diag_category_content)
       IF (term_id >= 1 .AND. term_id <= balance_diag_content_term_count) value = this%content(term_id)
    CASE (balance_diag_category_energy)
       IF (term_id >= 1 .AND. term_id <= balance_diag_energy_term_count) value = this%energy(term_id)
    END SELECT
  END FUNCTION balance_diag_get

  SUBROUTINE balance_diag_mpi_reduce(this)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
#ifdef PARALL
    INTEGER :: ierr

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%boundary_hdg, balance_diag_boundary_term_count, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%particles, balance_diag_particle_term_count, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%content, balance_diag_content_term_count, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%energy, balance_diag_energy_term_count, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE balance_diag_mpi_reduce

  SUBROUTINE balance_diag_print_summary(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: boundary_residual, plasma_balance, neutral_balance, total_balance

    boundary_residual = SUM(this%boundary_hdg)

    plasma_balance = balance_diag_plasma_particle_balance(this)
    neutral_balance = balance_diag_neutral_particle_balance(this)
    total_balance = plasma_balance + neutral_balance

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,'(A)') '--- Balance diagnostics summary ---'
       WRITE(6,'(A,1X,ES16.8)') 'boundary HDG / BC check =', boundary_residual
       WRITE(6,'(A,1X,ES16.8)') 'plasma particle balance =', plasma_balance
       WRITE(6,'(A,1X,ES16.8)') 'neutral particle balance =', neutral_balance
       WRITE(6,'(A,1X,ES16.8)') 'total particle balance =', total_balance
       WRITE(6,'(A,1X,ES16.8)') 'plasma particle content =', this%content(balance_diag_content_plasma_particles)
       WRITE(6,'(A,1X,ES16.8)') 'neutral particle content =', this%content(balance_diag_content_neutral_particles)
       WRITE(6,'(A,1X,ES16.8)') 'total particle content =', this%content(balance_diag_content_total_particles)
    ENDIF
  END SUBROUTINE balance_diag_print_summary

  SUBROUTINE balance_diag_print_detail(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER :: i

    IF (MPIvar%glob_id .NE. 0) RETURN

    WRITE(6,'(A)') '--- Boundary HDG / BC components ---'
    DO i = 1, balance_diag_boundary_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(balance_diag_boundary_label(i))//' =', this%boundary_hdg(i)
    ENDDO
    WRITE(6,'(A)') '--- Physical particle components ---'
    DO i = 1, balance_diag_particle_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(balance_diag_particle_label(i))//' =', this%particles(i)
    ENDDO
    WRITE(6,'(A)') '--- Global content components ---'
    DO i = 1, balance_diag_content_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(balance_diag_content_label(i))//' =', this%content(i)
    ENDDO
    WRITE(6,'(A)') '--- Reserved energy components ---'
    DO i = 1, balance_diag_energy_term_count
       WRITE(6,'(A,1X,ES16.8)') TRIM(balance_diag_energy_label(i))//' =', this%energy(i)
    ENDDO
  END SUBROUTINE balance_diag_print_detail

  SUBROUTINE balance_diag_write_hdf5(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    ! HDF5 output is wired in a later stage once live terms are migrated.
  END SUBROUTINE balance_diag_write_hdf5

  FUNCTION balance_diag_plasma_particle_balance(this) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: value

    value = this%particles(balance_diag_particle_plasma_ionization) &
       &+ this%particles(balance_diag_particle_plasma_recombination) &
       &+ this%particles(balance_diag_particle_plasma_boundary_flux)
  END FUNCTION balance_diag_plasma_particle_balance

  FUNCTION balance_diag_neutral_particle_balance(this) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: value

    value = this%particles(balance_diag_particle_neutral_ionization) &
       &+ this%particles(balance_diag_particle_neutral_recombination) &
       &+ this%particles(balance_diag_particle_puff_source) &
       &+ this%particles(balance_diag_particle_pump_sink) &
       &+ this%particles(balance_diag_particle_recycling_source) &
       &+ this%particles(balance_diag_particle_neutral_boundary_flux)
  END FUNCTION balance_diag_neutral_particle_balance

  FUNCTION balance_diag_boundary_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (balance_diag_boundary_plasma_parallel_flux)
       label = 'boundary HDG plasma parallel flux'
    CASE (balance_diag_boundary_plasma_diffusion_flux)
       label = 'boundary HDG plasma diffusion flux'
    CASE (balance_diag_boundary_plasma_pinch_flux)
       label = 'boundary HDG plasma pinch flux'
    CASE (balance_diag_boundary_neutral_diffusion_flux)
       label = 'boundary HDG neutral diffusion flux'
    CASE (balance_diag_boundary_neutral_pressure_flux)
       label = 'boundary HDG neutral pressure flux'
    CASE (balance_diag_boundary_neutral_convection_flux)
       label = 'boundary HDG neutral convection flux'
    CASE (balance_diag_boundary_neutral_total_flux)
       label = 'boundary HDG neutral total flux'
    CASE (balance_diag_boundary_tau_numerical_flux)
       label = 'boundary HDG tau numerical flux'
    CASE (balance_diag_boundary_wall_puff_source)
       label = 'boundary HDG wall puff source'
    CASE (balance_diag_boundary_wall_pump_sink)
       label = 'boundary HDG wall pump sink'
    CASE DEFAULT
       label = 'boundary HDG unknown term'
    END SELECT
  END FUNCTION balance_diag_boundary_label

  FUNCTION balance_diag_particle_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (balance_diag_particle_plasma_ionization)
       label = 'particle plasma ionization'
    CASE (balance_diag_particle_plasma_recombination)
       label = 'particle plasma recombination'
    CASE (balance_diag_particle_neutral_ionization)
       label = 'particle neutral ionization'
    CASE (balance_diag_particle_neutral_recombination)
       label = 'particle neutral recombination'
    CASE (balance_diag_particle_puff_source)
       label = 'particle puff source'
    CASE (balance_diag_particle_pump_sink)
       label = 'particle pump sink'
    CASE (balance_diag_particle_recycling_source)
       label = 'particle recycling source'
    CASE (balance_diag_particle_plasma_boundary_flux)
       label = 'particle plasma boundary flux'
    CASE (balance_diag_particle_neutral_boundary_flux)
       label = 'particle neutral boundary flux'
    CASE (balance_diag_particle_charge_exchange_rate)
       label = 'particle charge exchange rate'
    CASE DEFAULT
       label = 'particle unknown term'
    END SELECT
  END FUNCTION balance_diag_particle_label

  FUNCTION balance_diag_content_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (balance_diag_content_plasma_particles)
       label = 'content plasma particles'
    CASE (balance_diag_content_neutral_particles)
       label = 'content neutral particles'
    CASE (balance_diag_content_total_particles)
       label = 'content total particles'
    CASE (balance_diag_content_total_energy)
       label = 'content total energy'
    CASE DEFAULT
       label = 'content unknown term'
    END SELECT
  END FUNCTION balance_diag_content_label

  FUNCTION balance_diag_energy_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (balance_diag_energy_ion_balance)
       label = 'ion energy balance'
    CASE (balance_diag_energy_electron_balance)
       label = 'electron energy balance'
    CASE (balance_diag_energy_neutral_balance)
       label = 'neutral energy balance'
    CASE (balance_diag_energy_total_balance)
       label = 'total energy balance'
    CASE DEFAULT
       label = 'energy unknown term'
    END SELECT
  END FUNCTION balance_diag_energy_label

END MODULE balance_diagnostics
