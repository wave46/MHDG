MODULE balance_diagnostics
  USE HDF5
  USE HDF5_io_module
  USE GLOBALS, ONLY: Mesh, phys
  USE MPI_OMP
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_boundary_hdg = 1
  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_particles = 2
  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_content = 3
  INTEGER, PARAMETER, PUBLIC :: balance_diag_category_energy = 4

  ! All accumulated terms use the same sign convention:
  ! positive is directed into the domain, negative is loss from the domain.
  ! Boundary plasma terms below are recycled neutral source terms in the
  ! neutral-wall closure; physical plasma losses are tracked separately.
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
     PROCEDURE :: reset_boundary_hdg => balance_diag_reset_boundary_hdg
     PROCEDURE :: reset_particles_content => balance_diag_reset_particles_content
     PROCEDURE :: add => balance_diag_add
     PROCEDURE :: account_boundary_hdg => balance_diag_account_boundary_hdg
     PROCEDURE :: account_volume_particle_reactions => balance_diag_account_volume_particle_reactions
     PROCEDURE :: account_boundary_particles => balance_diag_account_boundary_particles
     PROCEDURE :: account_wall_particle_sources => balance_diag_account_wall_particle_sources
     PROCEDURE :: account_particle_content => balance_diag_account_particle_content
     PROCEDURE :: get => balance_diag_get
     PROCEDURE :: mpi_reduce => balance_diag_mpi_reduce
     PROCEDURE :: mpi_reduce_boundary_hdg => balance_diag_mpi_reduce_boundary_hdg
     PROCEDURE :: mpi_reduce_particles_content => balance_diag_mpi_reduce_particles_content
     PROCEDURE :: print_summary => balance_diag_print_summary
     PROCEDURE :: print_boundary_hdg_summary => balance_diag_print_boundary_hdg_summary
     PROCEDURE :: print_particle_summary => balance_diag_print_particle_summary
     PROCEDURE :: print_detail => balance_diag_print_detail
     PROCEDURE :: print_boundary_hdg_detail => balance_diag_print_boundary_hdg_detail
     PROCEDURE :: print_particle_detail => balance_diag_print_particle_detail
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

  SUBROUTINE balance_diag_reset_boundary_hdg(this)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this

    this%boundary_hdg = 0.d0
  END SUBROUTINE balance_diag_reset_boundary_hdg

  SUBROUTINE balance_diag_reset_particles_content(this)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this

    this%particles = 0.d0
    this%content = 0.d0
  END SUBROUTINE balance_diag_reset_particles_content

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

  SUBROUTINE balance_diag_account_boundary_hdg(this, plasma_parallel_flux, plasma_diffusion_flux, plasma_pinch_flux, &
       &neutral_diffusion_flux, neutral_pressure_flux, neutral_convection_flux, neutral_total_flux, tau_numerical_flux, &
       &wall_puff_source, wall_pump_sink, include_wall_sources)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_parallel_flux, plasma_diffusion_flux, plasma_pinch_flux
    REAL*8, INTENT(IN) :: neutral_diffusion_flux, neutral_pressure_flux, neutral_convection_flux, neutral_total_flux
    REAL*8, INTENT(IN) :: tau_numerical_flux, wall_puff_source, wall_pump_sink
    LOGICAL, INTENT(IN) :: include_wall_sources

    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_plasma_parallel_flux, plasma_parallel_flux)
    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_plasma_diffusion_flux, plasma_diffusion_flux)
    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_plasma_pinch_flux, plasma_pinch_flux)
    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_neutral_diffusion_flux, neutral_diffusion_flux)
    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_neutral_pressure_flux, neutral_pressure_flux)
    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_neutral_convection_flux, neutral_convection_flux)
    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_neutral_total_flux, neutral_total_flux)
    CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_tau_numerical_flux, tau_numerical_flux)
    IF (include_wall_sources) THEN
       CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_wall_puff_source, wall_puff_source)
       CALL this%add(balance_diag_category_boundary_hdg, balance_diag_boundary_wall_pump_sink, wall_pump_sink)
    ENDIF
  END SUBROUTINE balance_diag_account_boundary_hdg

  SUBROUTINE balance_diag_account_volume_particle_reactions(this, ionization_rate, recombination_rate, charge_exchange_rate)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: ionization_rate, recombination_rate, charge_exchange_rate

    CALL this%add(balance_diag_category_particles, balance_diag_particle_plasma_ionization, ionization_rate)
    CALL this%add(balance_diag_category_particles, balance_diag_particle_neutral_ionization, -ionization_rate)
    CALL this%add(balance_diag_category_particles, balance_diag_particle_plasma_recombination, -recombination_rate)
    CALL this%add(balance_diag_category_particles, balance_diag_particle_neutral_recombination, recombination_rate)
    CALL this%add(balance_diag_category_particles, balance_diag_particle_charge_exchange_rate, charge_exchange_rate)
  END SUBROUTINE balance_diag_account_volume_particle_reactions

  SUBROUTINE balance_diag_account_boundary_particles(this, plasma_boundary_flux, neutral_boundary_flux, recycling_source, &
       &wall_puff_source, wall_pump_sink, include_wall_sources)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_boundary_flux, neutral_boundary_flux, recycling_source
    REAL*8, INTENT(IN) :: wall_puff_source, wall_pump_sink
    LOGICAL, INTENT(IN) :: include_wall_sources

    CALL this%add(balance_diag_category_particles, balance_diag_particle_plasma_boundary_flux, plasma_boundary_flux)
    CALL this%add(balance_diag_category_particles, balance_diag_particle_neutral_boundary_flux, neutral_boundary_flux)
    CALL this%add(balance_diag_category_particles, balance_diag_particle_recycling_source, recycling_source)
    IF (include_wall_sources) THEN
       CALL this%add(balance_diag_category_particles, balance_diag_particle_puff_source, wall_puff_source)
       CALL this%add(balance_diag_category_particles, balance_diag_particle_pump_sink, wall_pump_sink)
    ENDIF
  END SUBROUTINE balance_diag_account_boundary_particles

  SUBROUTINE balance_diag_account_wall_particle_sources(this, puff_source, pump_sink)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: puff_source, pump_sink

    CALL this%add(balance_diag_category_particles, balance_diag_particle_puff_source, puff_source)
    CALL this%add(balance_diag_category_particles, balance_diag_particle_pump_sink, pump_sink)
  END SUBROUTINE balance_diag_account_wall_particle_sources

  SUBROUTINE balance_diag_account_particle_content(this, plasma_particles, neutral_particles)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_particles, neutral_particles

    CALL this%add(balance_diag_category_content, balance_diag_content_plasma_particles, plasma_particles)
    CALL this%add(balance_diag_category_content, balance_diag_content_neutral_particles, neutral_particles)
    CALL this%add(balance_diag_category_content, balance_diag_content_total_particles, plasma_particles + neutral_particles)
  END SUBROUTINE balance_diag_account_particle_content

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

  SUBROUTINE balance_diag_mpi_reduce_boundary_hdg(this)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
#ifdef PARALL
    INTEGER :: ierr

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%boundary_hdg, balance_diag_boundary_term_count, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE balance_diag_mpi_reduce_boundary_hdg

  SUBROUTINE balance_diag_mpi_reduce_particles_content(this)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
#ifdef PARALL
    INTEGER :: ierr

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%particles, balance_diag_particle_term_count, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%content, balance_diag_content_term_count, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE balance_diag_mpi_reduce_particles_content

  SUBROUTINE balance_diag_print_summary(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: boundary_residual, plasma_balance, neutral_balance, total_balance

    boundary_residual = balance_diag_boundary_hdg_check(this)

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

  SUBROUTINE balance_diag_print_boundary_hdg_summary(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: boundary_residual

    boundary_residual = balance_diag_boundary_hdg_check(this)

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,'(A)') '----------------------------------------'
       WRITE(6,'(A)') 'Boundary HDG / BC diagnostics'
       WRITE(6,'(A,1X,ES11.3)') '  neutral closure residual:', boundary_residual
       WRITE(6,'(A)') '  residual = recycled + wall + tau - neutral flux'
       WRITE(6,'(A)') '----------------------------------------'
    ENDIF
  END SUBROUTINE balance_diag_print_boundary_hdg_summary

  SUBROUTINE balance_diag_print_particle_summary(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: plasma_balance, neutral_balance, total_balance

    plasma_balance = balance_diag_plasma_particle_balance(this)
    neutral_balance = balance_diag_neutral_particle_balance(this)
    total_balance = plasma_balance + neutral_balance

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,'(A)') '----------------------------------------'
       WRITE(6,'(A)') 'Physical particle diagnostics'
       WRITE(6,'(A,1X,ES11.3)') '  plasma balance:', plasma_balance
       WRITE(6,'(A,1X,ES11.3)') '  neutral balance:', neutral_balance
       WRITE(6,'(A,1X,ES11.3)') '  total balance  :', total_balance
       WRITE(6,'(A,1X,ES11.3)') '  plasma content :', this%content(balance_diag_content_plasma_particles)
       WRITE(6,'(A,1X,ES11.3)') '  neutral content:', this%content(balance_diag_content_neutral_particles)
       WRITE(6,'(A,1X,ES11.3)') '  total content  :', this%content(balance_diag_content_total_particles)
       WRITE(6,'(A)') '----------------------------------------'
    ENDIF
  END SUBROUTINE balance_diag_print_particle_summary

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

  SUBROUTINE balance_diag_print_boundary_hdg_detail(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER :: i
    REAL*8 :: boundary_residual, recycled_source, neutral_flux, wall_flux, tau_flux

    IF (MPIvar%glob_id .NE. 0) RETURN

    boundary_residual = balance_diag_boundary_hdg_check(this)
    recycled_source = this%boundary_hdg(balance_diag_boundary_plasma_parallel_flux) &
       &+ this%boundary_hdg(balance_diag_boundary_plasma_diffusion_flux) &
       &+ this%boundary_hdg(balance_diag_boundary_plasma_pinch_flux)
    neutral_flux = this%boundary_hdg(balance_diag_boundary_neutral_total_flux)
    wall_flux = this%boundary_hdg(balance_diag_boundary_wall_puff_source) &
       &+ this%boundary_hdg(balance_diag_boundary_wall_pump_sink)
    tau_flux = this%boundary_hdg(balance_diag_boundary_tau_numerical_flux)

    WRITE(6,'(A)') '----------------------------------------'
    WRITE(6,'(A)') 'Boundary HDG / BC diagnostics'
    WRITE(6,'(A,1X,ES11.3)') '  neutral closure residual:', boundary_residual
    WRITE(6,'(A,1X,ES11.3)') '  neutral flux into domain:', neutral_flux
    WRITE(6,'(A,1X,ES11.3)') '  recycled source         :', recycled_source
    WRITE(6,'(A,1X,ES11.3)') '  wall sources/sinks      :', wall_flux
    WRITE(6,'(A,1X,ES11.3)') '  tau correction          :', tau_flux
    WRITE(6,'(A)') '  components, inward-positive:'
    DO i = 1, balance_diag_boundary_term_count
       WRITE(6,'(A,1X,ES11.3)') '    '//TRIM(balance_diag_boundary_label(i))//':', this%boundary_hdg(i)
    ENDDO
    WRITE(6,'(A)') '----------------------------------------'
  END SUBROUTINE balance_diag_print_boundary_hdg_detail

  SUBROUTINE balance_diag_print_particle_detail(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER :: i
    REAL*8 :: plasma_balance, neutral_balance, total_balance

    IF (MPIvar%glob_id .NE. 0) RETURN

    plasma_balance = balance_diag_plasma_particle_balance(this)
    neutral_balance = balance_diag_neutral_particle_balance(this)
    total_balance = plasma_balance + neutral_balance

    WRITE(6,'(A)') '----------------------------------------'
    WRITE(6,'(A)') 'Physical particle diagnostics'
    WRITE(6,'(A,1X,ES11.3)') '  plasma balance:', plasma_balance
    WRITE(6,'(A,1X,ES11.3)') '  neutral balance:', neutral_balance
    WRITE(6,'(A,1X,ES11.3)') '  total balance  :', total_balance
    WRITE(6,'(A,1X,ES11.3)') '  plasma content :', this%content(balance_diag_content_plasma_particles)
    WRITE(6,'(A,1X,ES11.3)') '  neutral content:', this%content(balance_diag_content_neutral_particles)
    WRITE(6,'(A,1X,ES11.3)') '  total content  :', this%content(balance_diag_content_total_particles)
    WRITE(6,'(A)') '  particle components, inward-positive:'
    DO i = 1, balance_diag_particle_term_count
       WRITE(6,'(A,1X,ES11.3)') '    '//TRIM(balance_diag_particle_label(i))//':', this%particles(i)
    ENDDO
    WRITE(6,'(A)') '  particle content:'
    DO i = 1, balance_diag_content_term_count
       IF (i .EQ. balance_diag_content_total_energy) CYCLE
       WRITE(6,'(A,1X,ES11.3)') '    '//TRIM(balance_diag_content_label(i))//':', this%content(i)
    ENDDO
    WRITE(6,'(A)') '----------------------------------------'
  END SUBROUTINE balance_diag_print_particle_detail

  SUBROUTINE balance_diag_write_hdf5(this, file_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: diagnostics_group_id, boundary_group_id, balance_group_id, content_group_id
    INTEGER(HID_T) :: neutrals_group_id, plasma_group_id, total_group_id, particles_group_id, terms_group_id, summary_group_id
    INTEGER :: ierr

    diagnostics_group_id = -1
    boundary_group_id = -1
    balance_group_id = -1
    content_group_id = -1
    neutrals_group_id = -1
    plasma_group_id = -1
    total_group_id = -1
    particles_group_id = -1
    terms_group_id = -1
    summary_group_id = -1

    IF (MPIvar%glob_id .EQ. 0) THEN
       CALL HDF5_group_create('diagnostics', file_id, diagnostics_group_id, ierr)

       CALL HDF5_group_create('boundary_conditions', diagnostics_group_id, boundary_group_id, ierr)
       CALL HDF5_group_create('neutrals', boundary_group_id, neutrals_group_id, ierr)
       CALL HDF5_group_create('terms', neutrals_group_id, terms_group_id, ierr)
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_plasma_parallel_flux), &
          &'recycled_plasma_parallel_source')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_plasma_diffusion_flux), &
          &'recycled_plasma_diffusion_source')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_plasma_pinch_flux), &
          &'recycled_plasma_pinch_source')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_plasma_parallel_flux) &
          &+ this%boundary_hdg(balance_diag_boundary_plasma_diffusion_flux) &
          &+ this%boundary_hdg(balance_diag_boundary_plasma_pinch_flux), 'recycling_source')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_neutral_diffusion_flux), &
          &'neutral_diffusion_boundary_flux')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_neutral_pressure_flux), &
          &'neutral_pressure_boundary_flux')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_neutral_convection_flux), &
          &'neutral_convection_boundary_flux')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_neutral_total_flux), &
          &'neutral_total_boundary_flux')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_tau_numerical_flux), &
          &'tau_numerical_boundary_flux')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_wall_puff_source), 'wall_puff_source')
       CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(balance_diag_boundary_wall_pump_sink), 'wall_pump_sink')
       CALL HDF5_group_close(terms_group_id, ierr)

       CALL HDF5_group_create('summary', neutrals_group_id, summary_group_id, ierr)
       CALL HDF5_real_saving(summary_group_id, balance_diag_boundary_hdg_check(this), 'neutral_closure_residual')
       CALL HDF5_group_close(summary_group_id, ierr)
       CALL HDF5_group_close(neutrals_group_id, ierr)
       CALL HDF5_group_close(boundary_group_id, ierr)

       CALL HDF5_group_create('balance', diagnostics_group_id, balance_group_id, ierr)

       CALL HDF5_group_create('plasma', balance_group_id, plasma_group_id, ierr)
       CALL HDF5_group_create('particles', plasma_group_id, particles_group_id, ierr)
       CALL HDF5_group_create('terms', particles_group_id, terms_group_id, ierr)
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_plasma_ionization), 'ionization')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_plasma_recombination), 'recombination')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_plasma_boundary_flux), 'boundary_flux')
       CALL HDF5_group_close(terms_group_id, ierr)
       CALL HDF5_group_create('summary', particles_group_id, summary_group_id, ierr)
       CALL HDF5_real_saving(summary_group_id, balance_diag_plasma_particle_balance(this), 'particle_balance')
       CALL HDF5_group_close(summary_group_id, ierr)
       CALL HDF5_group_close(particles_group_id, ierr)
       CALL HDF5_group_close(plasma_group_id, ierr)

       CALL HDF5_group_create('neutrals', balance_group_id, neutrals_group_id, ierr)
       CALL HDF5_group_create('particles', neutrals_group_id, particles_group_id, ierr)
       CALL HDF5_group_create('terms', particles_group_id, terms_group_id, ierr)
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_neutral_ionization), 'ionization')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_neutral_recombination), 'recombination')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_neutral_boundary_flux), 'boundary_flux')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_puff_source), 'puff_source')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_pump_sink), 'pump_sink')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_recycling_source), 'recycling_source')
       CALL HDF5_real_saving(terms_group_id, this%particles(balance_diag_particle_charge_exchange_rate), 'charge_exchange_rate')
       CALL HDF5_group_close(terms_group_id, ierr)
       CALL HDF5_group_create('summary', particles_group_id, summary_group_id, ierr)
       CALL HDF5_real_saving(summary_group_id, balance_diag_neutral_particle_balance(this), 'particle_balance')
       CALL HDF5_group_close(summary_group_id, ierr)
    ENDIF
    CALL balance_diag_write_wall_source_nodal(particles_group_id)
    IF (MPIvar%glob_id .EQ. 0) THEN
       CALL HDF5_group_close(particles_group_id, ierr)
       CALL HDF5_group_close(neutrals_group_id, ierr)

       CALL HDF5_group_create('total', balance_group_id, total_group_id, ierr)
       CALL HDF5_group_create('particles', total_group_id, particles_group_id, ierr)
       CALL HDF5_group_create('summary', particles_group_id, summary_group_id, ierr)
       CALL HDF5_real_saving(summary_group_id, balance_diag_plasma_particle_balance(this) &
          &+ balance_diag_neutral_particle_balance(this), 'particle_balance')
       CALL HDF5_group_close(summary_group_id, ierr)
       CALL HDF5_group_close(particles_group_id, ierr)
       CALL HDF5_group_close(total_group_id, ierr)
       CALL HDF5_group_close(balance_group_id, ierr)

       CALL HDF5_group_create('content', diagnostics_group_id, content_group_id, ierr)
       CALL HDF5_group_create('plasma', content_group_id, plasma_group_id, ierr)
       CALL HDF5_real_saving(plasma_group_id, this%content(balance_diag_content_plasma_particles), 'particles')
       CALL HDF5_group_close(plasma_group_id, ierr)
       CALL HDF5_group_create('neutrals', content_group_id, neutrals_group_id, ierr)
       CALL HDF5_real_saving(neutrals_group_id, this%content(balance_diag_content_neutral_particles), 'particles')
       CALL HDF5_group_close(neutrals_group_id, ierr)
       CALL HDF5_group_create('total', content_group_id, total_group_id, ierr)
       CALL HDF5_real_saving(total_group_id, this%content(balance_diag_content_total_particles), 'particles')
       CALL HDF5_group_close(total_group_id, ierr)
       CALL HDF5_group_close(content_group_id, ierr)

       CALL HDF5_group_close(diagnostics_group_id, ierr)
    ENDIF
  END SUBROUTINE balance_diag_write_hdf5

  SUBROUTINE balance_diag_write_wall_source_nodal(group_id)
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER(HID_T) :: nodal_group_id
    INTEGER :: expected_size, ierr
#ifdef PARALL
    INTEGER :: iel, g, ind_local, ind_global, available_local, available_global
    REAL*8, ALLOCATABLE :: puff_glob(:), pump_glob(:), net_glob(:)
#endif

    expected_size = Mesh%Nelems*Mesh%Nnodesperelem
#ifdef PARALL
    available_local = 0
    IF (ALLOCATED(phys%neutral_wall_source_puff_Nod) .AND. &
       &ALLOCATED(phys%neutral_wall_source_pump_Nod) .AND. &
       &ALLOCATED(phys%neutral_wall_source_net_Nod)) THEN
       IF ((SIZE(phys%neutral_wall_source_puff_Nod) .EQ. expected_size) .AND. &
          &(SIZE(phys%neutral_wall_source_pump_Nod) .EQ. expected_size) .AND. &
          &(SIZE(phys%neutral_wall_source_net_Nod) .EQ. expected_size)) available_local = 1
    ENDIF
    CALL MPI_Allreduce(available_local, available_global, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    IF (available_global .EQ. 0) RETURN
#else
    IF (.NOT. ALLOCATED(phys%neutral_wall_source_puff_Nod)) RETURN
    IF (.NOT. ALLOCATED(phys%neutral_wall_source_pump_Nod)) RETURN
    IF (.NOT. ALLOCATED(phys%neutral_wall_source_net_Nod)) RETURN
    IF (SIZE(phys%neutral_wall_source_puff_Nod) .NE. expected_size) RETURN
    IF (SIZE(phys%neutral_wall_source_pump_Nod) .NE. expected_size) RETURN
    IF (SIZE(phys%neutral_wall_source_net_Nod) .NE. expected_size) RETURN
#endif

#ifdef PARALL
    ALLOCATE(puff_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
    ALLOCATE(pump_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
    ALLOCATE(net_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
    puff_glob = 0.d0
    pump_glob = 0.d0
    net_glob = 0.d0

    DO iel = 1, Mesh%Nelems
       IF (Mesh%ghostElems(iel) .EQ. 0) THEN
          DO g = 1, Mesh%Nnodesperelem
             ind_local = (iel - 1)*Mesh%Nnodesperelem + g
             ind_global = (Mesh%loc2glob_el(iel) - 1)*Mesh%Nnodesperelem + g
             puff_glob(ind_global) = phys%neutral_wall_source_puff_Nod(ind_local)
             pump_glob(ind_global) = phys%neutral_wall_source_pump_Nod(ind_local)
             net_glob(ind_global) = phys%neutral_wall_source_net_Nod(ind_local)
          ENDDO
       ENDIF
    ENDDO

    CALL MPI_Allreduce(MPI_IN_PLACE, puff_glob, SIZE(puff_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, pump_glob, SIZE(pump_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, net_glob, SIZE(net_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    IF (MPIvar%glob_id .EQ. 0) THEN
       CALL HDF5_group_create('nodal_wall_sources', group_id, nodal_group_id, ierr)
       CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_puff_total, 'element_puff_total')
       CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_pump_total, 'element_pump_total')
       CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_puff_total - phys%neutral_wall_source_pump_total, &
          &'element_net_total')
       CALL HDF5_array1D_saving(nodal_group_id, puff_glob, SIZE(puff_glob), 'puff_flux_density')
       CALL HDF5_array1D_saving(nodal_group_id, pump_glob, SIZE(pump_glob), 'pump_flux_density')
       CALL HDF5_array1D_saving(nodal_group_id, net_glob, SIZE(net_glob), 'net_flux_density')
       CALL HDF5_group_close(nodal_group_id, ierr)
    ENDIF

    DEALLOCATE(puff_glob, pump_glob, net_glob)
#else
    CALL HDF5_group_create('nodal_wall_sources', group_id, nodal_group_id, ierr)
    CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_puff_total, 'element_puff_total')
    CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_pump_total, 'element_pump_total')
    CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_puff_total - phys%neutral_wall_source_pump_total, &
       &'element_net_total')
    CALL HDF5_array1D_saving(nodal_group_id, phys%neutral_wall_source_puff_Nod, &
       &SIZE(phys%neutral_wall_source_puff_Nod), 'puff_flux_density')
    CALL HDF5_array1D_saving(nodal_group_id, phys%neutral_wall_source_pump_Nod, &
       &SIZE(phys%neutral_wall_source_pump_Nod), 'pump_flux_density')
    CALL HDF5_array1D_saving(nodal_group_id, phys%neutral_wall_source_net_Nod, &
       &SIZE(phys%neutral_wall_source_net_Nod), 'net_flux_density')
    CALL HDF5_group_close(nodal_group_id, ierr)
#endif
  END SUBROUTINE balance_diag_write_wall_source_nodal

  FUNCTION balance_diag_boundary_hdg_check(this) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: value

    value = this%boundary_hdg(balance_diag_boundary_plasma_parallel_flux) &
       &+ this%boundary_hdg(balance_diag_boundary_plasma_diffusion_flux) &
       &+ this%boundary_hdg(balance_diag_boundary_plasma_pinch_flux) &
       &- this%boundary_hdg(balance_diag_boundary_neutral_total_flux) &
       &+ this%boundary_hdg(balance_diag_boundary_tau_numerical_flux) &
       &+ this%boundary_hdg(balance_diag_boundary_wall_puff_source) &
       &+ this%boundary_hdg(balance_diag_boundary_wall_pump_sink)
  END FUNCTION balance_diag_boundary_hdg_check

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
       &+ this%particles(balance_diag_particle_neutral_boundary_flux)
  END FUNCTION balance_diag_neutral_particle_balance

  FUNCTION balance_diag_boundary_label(term_id) RESULT(label)
    INTEGER, INTENT(IN) :: term_id
    CHARACTER(LEN=80) :: label

    SELECT CASE (term_id)
    CASE (balance_diag_boundary_plasma_parallel_flux)
       label = 'recycled plasma parallel source'
    CASE (balance_diag_boundary_plasma_diffusion_flux)
       label = 'recycled plasma diffusion source'
    CASE (balance_diag_boundary_plasma_pinch_flux)
       label = 'recycled plasma pinch source'
    CASE (balance_diag_boundary_neutral_diffusion_flux)
       label = 'neutral diffusion boundary flux'
    CASE (balance_diag_boundary_neutral_pressure_flux)
       label = 'neutral pressure boundary flux'
    CASE (balance_diag_boundary_neutral_convection_flux)
       label = 'neutral convection boundary flux'
    CASE (balance_diag_boundary_neutral_total_flux)
       label = 'neutral total boundary flux'
    CASE (balance_diag_boundary_tau_numerical_flux)
       label = 'tau numerical boundary flux'
    CASE (balance_diag_boundary_wall_puff_source)
       label = 'wall puff source'
    CASE (balance_diag_boundary_wall_pump_sink)
       label = 'wall pump sink'
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
