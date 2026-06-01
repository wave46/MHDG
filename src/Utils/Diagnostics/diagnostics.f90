MODULE diagnostics
  USE HDF5
  USE MPI_OMP
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: diag_category_boundary_hdg = 1
  INTEGER, PARAMETER, PUBLIC :: diag_category_particles = 2
  INTEGER, PARAMETER, PUBLIC :: diag_category_content = 3
  INTEGER, PARAMETER, PUBLIC :: diag_category_energy = 4

  ! All accumulated terms use the same sign convention:
  ! positive is directed into the domain, negative is loss from the domain.
  ! Boundary plasma terms below are recycled neutral source terms in the
  ! neutral-wall closure; physical plasma losses are tracked separately.
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_plasma_parallel_flux = 1
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_plasma_diffusion_flux = 2
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_plasma_pinch_flux = 3
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_neutral_diffusion_flux = 4
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_neutral_pressure_flux = 5
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_neutral_convection_flux = 6
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_neutral_total_flux = 7
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_tau_numerical_flux = 8
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_wall_puff_source = 9
  INTEGER, PARAMETER, PUBLIC :: diag_boundary_wall_pump_sink = 10
  INTEGER, PARAMETER :: diag_boundary_term_count = 10

  INTEGER, PARAMETER, PUBLIC :: diag_particle_plasma_ionization = 1
  INTEGER, PARAMETER, PUBLIC :: diag_particle_plasma_recombination = 2
  INTEGER, PARAMETER, PUBLIC :: diag_particle_plasma_boundary_flux = 3
  INTEGER, PARAMETER, PUBLIC :: diag_particle_neutral_ionization = 4
  INTEGER, PARAMETER, PUBLIC :: diag_particle_neutral_recombination = 5
  INTEGER, PARAMETER, PUBLIC :: diag_particle_puff_source = 6
  INTEGER, PARAMETER, PUBLIC :: diag_particle_pump_sink = 7
  INTEGER, PARAMETER, PUBLIC :: diag_particle_recycling_source = 8
  INTEGER, PARAMETER, PUBLIC :: diag_particle_neutral_boundary_flux = 9
  INTEGER, PARAMETER, PUBLIC :: diag_particle_charge_exchange_rate = 10
  INTEGER, PARAMETER :: diag_particle_term_count = 10

  INTEGER, PARAMETER, PUBLIC :: diag_content_plasma_particles = 1
  INTEGER, PARAMETER, PUBLIC :: diag_content_neutral_particles = 2
  INTEGER, PARAMETER, PUBLIC :: diag_content_total_particles = 3
  INTEGER, PARAMETER, PUBLIC :: diag_content_total_energy = 4
  INTEGER, PARAMETER :: diag_content_term_count = 4

  INTEGER, PARAMETER, PUBLIC :: diag_energy_ion_balance = 1
  INTEGER, PARAMETER, PUBLIC :: diag_energy_electron_balance = 2
  INTEGER, PARAMETER, PUBLIC :: diag_energy_neutral_balance = 3
  INTEGER, PARAMETER, PUBLIC :: diag_energy_total_balance = 4
  INTEGER, PARAMETER :: diag_energy_term_count = 4

  TYPE :: diag_boundary_summary_type
     REAL*8 :: residual = 0.d0
     REAL*8 :: recycled_source = 0.d0
     REAL*8 :: neutral_flux = 0.d0
     REAL*8 :: wall_flux = 0.d0
     REAL*8 :: tau_flux = 0.d0
  END TYPE diag_boundary_summary_type

  TYPE :: diag_particle_summary_type
     REAL*8 :: plasma_balance = 0.d0
     REAL*8 :: neutral_balance = 0.d0
     REAL*8 :: total_balance = 0.d0
     REAL*8 :: plasma_content = 0.d0
     REAL*8 :: neutral_content = 0.d0
     REAL*8 :: total_content = 0.d0
  END TYPE diag_particle_summary_type

  TYPE, PUBLIC :: diagnostics_type
     REAL*8 :: boundary_hdg(diag_boundary_term_count) = 0.d0
     REAL*8 :: particles(diag_particle_term_count) = 0.d0
     REAL*8 :: content(diag_content_term_count) = 0.d0
     REAL*8 :: energy(diag_energy_term_count) = 0.d0
   CONTAINS
     PROCEDURE :: reset => diag_reset
     PROCEDURE :: reset_boundary_hdg => diag_reset_boundary_hdg
     PROCEDURE :: reset_particles_content => diag_reset_particles_content
     PROCEDURE :: add => diag_add
     PROCEDURE :: account_boundary_hdg => diag_account_boundary_hdg
     PROCEDURE :: account_volume_particle_reactions => diag_account_volume_particle_reactions
     PROCEDURE :: account_boundary_particles => diag_account_boundary_particles
     PROCEDURE :: account_wall_particle_sources => diag_account_wall_particle_sources
     PROCEDURE :: account_particle_content => diag_account_particle_content
     PROCEDURE :: mpi_reduce_boundary_hdg => diag_mpi_reduce_boundary_hdg
     PROCEDURE :: mpi_reduce_particles_content => diag_mpi_reduce_particles_content
     PROCEDURE :: print_boundary_hdg_summary => diag_print_boundary_hdg_summary
     PROCEDURE :: print_particle_summary => diag_print_particle_summary
     PROCEDURE :: print_boundary_hdg_detail => diag_print_boundary_hdg_detail
     PROCEDURE :: print_particle_detail => diag_print_particle_detail
     PROCEDURE :: write_hdf5 => diag_write_hdf5
  END TYPE diagnostics_type

  TYPE(diagnostics_type), PUBLIC :: diag
  LOGICAL :: diag_wall_source_nodal_ready = .FALSE.
  REAL*8, ALLOCATABLE :: diag_wall_source_puff_nodal(:)
  REAL*8, ALLOCATABLE :: diag_wall_source_pump_nodal(:)
  REAL*8, ALLOCATABLE :: diag_wall_source_net_nodal(:)

  INTERFACE
     MODULE SUBROUTINE diag_reset(this)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
     END SUBROUTINE diag_reset

     MODULE SUBROUTINE diag_reset_boundary_hdg(this)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
     END SUBROUTINE diag_reset_boundary_hdg

     MODULE SUBROUTINE diag_reset_particles_content(this)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
     END SUBROUTINE diag_reset_particles_content

     MODULE SUBROUTINE diag_add(this, category_id, term_id, value)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
       INTEGER, INTENT(IN) :: category_id, term_id
       REAL*8, INTENT(IN) :: value
     END SUBROUTINE diag_add

     MODULE SUBROUTINE diag_account_boundary_hdg(this, plasma_parallel_flux, plasma_diffusion_flux, plasma_pinch_flux, &
          &neutral_diffusion_flux, neutral_pressure_flux, neutral_convection_flux, neutral_total_flux, tau_numerical_flux, &
          &wall_puff_source, wall_pump_sink, include_wall_sources)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: plasma_parallel_flux, plasma_diffusion_flux, plasma_pinch_flux
       REAL*8, INTENT(IN) :: neutral_diffusion_flux, neutral_pressure_flux, neutral_convection_flux, neutral_total_flux
       REAL*8, INTENT(IN) :: tau_numerical_flux, wall_puff_source, wall_pump_sink
       LOGICAL, INTENT(IN) :: include_wall_sources
     END SUBROUTINE diag_account_boundary_hdg

     MODULE SUBROUTINE diag_account_volume_particle_reactions(this, ionization_rate, recombination_rate, &
          &charge_exchange_rate)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: ionization_rate, recombination_rate, charge_exchange_rate
     END SUBROUTINE diag_account_volume_particle_reactions

     MODULE SUBROUTINE diag_account_boundary_particles(this, plasma_boundary_flux, neutral_boundary_flux, &
          &recycling_source, wall_puff_source, wall_pump_sink, include_wall_sources)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: plasma_boundary_flux, neutral_boundary_flux, recycling_source
       REAL*8, INTENT(IN) :: wall_puff_source, wall_pump_sink
       LOGICAL, INTENT(IN) :: include_wall_sources
     END SUBROUTINE diag_account_boundary_particles

     MODULE SUBROUTINE diag_account_wall_particle_sources(this, puff_source, pump_sink)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: puff_source, pump_sink
     END SUBROUTINE diag_account_wall_particle_sources

     MODULE SUBROUTINE diag_account_particle_content(this, plasma_particles, neutral_particles)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: plasma_particles, neutral_particles
     END SUBROUTINE diag_account_particle_content

     MODULE SUBROUTINE diag_mpi_reduce_boundary_hdg(this)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
     END SUBROUTINE diag_mpi_reduce_boundary_hdg

     MODULE SUBROUTINE diag_mpi_reduce_particles_content(this)
       CLASS(diagnostics_type), INTENT(INOUT) :: this
     END SUBROUTINE diag_mpi_reduce_particles_content

     MODULE SUBROUTINE diag_print_boundary_hdg_summary(this)
       CLASS(diagnostics_type), INTENT(IN) :: this
     END SUBROUTINE diag_print_boundary_hdg_summary

     MODULE SUBROUTINE diag_print_particle_summary(this)
       CLASS(diagnostics_type), INTENT(IN) :: this
     END SUBROUTINE diag_print_particle_summary

     MODULE SUBROUTINE diag_print_boundary_hdg_detail(this)
       CLASS(diagnostics_type), INTENT(IN) :: this
     END SUBROUTINE diag_print_boundary_hdg_detail

     MODULE SUBROUTINE diag_print_particle_detail(this)
       CLASS(diagnostics_type), INTENT(IN) :: this
     END SUBROUTINE diag_print_particle_detail

     MODULE SUBROUTINE diag_write_hdf5(this, file_id)
       CLASS(diagnostics_type), INTENT(IN) :: this
       INTEGER(HID_T), INTENT(IN) :: file_id
     END SUBROUTINE diag_write_hdf5

     MODULE SUBROUTINE diag_prepare_wall_source_nodal_output()
     END SUBROUTINE diag_prepare_wall_source_nodal_output

     MODULE SUBROUTINE diag_write_wall_source_nodal(group_id)
       INTEGER(HID_T), INTENT(IN) :: group_id
     END SUBROUTINE diag_write_wall_source_nodal

     MODULE FUNCTION diag_boundary_summary(this) RESULT(summary)
       CLASS(diagnostics_type), INTENT(IN) :: this
       TYPE(diag_boundary_summary_type) :: summary
     END FUNCTION diag_boundary_summary

     MODULE FUNCTION diag_plasma_particle_balance(this) RESULT(value)
       CLASS(diagnostics_type), INTENT(IN) :: this
       REAL*8 :: value
     END FUNCTION diag_plasma_particle_balance

     MODULE FUNCTION diag_neutral_particle_balance(this) RESULT(value)
       CLASS(diagnostics_type), INTENT(IN) :: this
       REAL*8 :: value
     END FUNCTION diag_neutral_particle_balance

     MODULE FUNCTION diag_particle_summary(this) RESULT(summary)
       CLASS(diagnostics_type), INTENT(IN) :: this
       TYPE(diag_particle_summary_type) :: summary
     END FUNCTION diag_particle_summary

     MODULE FUNCTION diag_content_value(this, term_id) RESULT(value)
       CLASS(diagnostics_type), INTENT(IN) :: this
       INTEGER, INTENT(IN) :: term_id
       REAL*8 :: value
     END FUNCTION diag_content_value

     MODULE FUNCTION diag_particle_integral_units() RESULT(units)
       CHARACTER(LEN=32) :: units
     END FUNCTION diag_particle_integral_units

     MODULE FUNCTION diag_content_units(term_id) RESULT(units)
       INTEGER, INTENT(IN) :: term_id
       CHARACTER(LEN=32) :: units
     END FUNCTION diag_content_units
  END INTERFACE

END MODULE diagnostics
