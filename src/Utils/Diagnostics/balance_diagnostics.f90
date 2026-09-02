!**********************************************************************
! Public contract and lifecycle for equation-oriented balance diagnostics.
!
! Physical fluxes and HDG stabilization are inward-positive. Primary
! equations follow the conservative model order; total particle and total
! plasma-energy balances are derived after reduction and are never stored.
!**********************************************************************
MODULE balance_diagnostics
  USE HDF5, ONLY: HID_T
  USE MPI_OMP
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: balance_mode_invalid = -1
  INTEGER, PARAMETER, PUBLIC :: balance_mode_off = 0
  INTEGER, PARAMETER, PUBLIC :: balance_mode_summary = 1
  INTEGER, PARAMETER, PUBLIC :: balance_mode_detailed = 2

  ! Stored equation columns mirror U1:U5. Neutral momentum is deliberately
  ! excluded from the independent diagnostic scope.
  INTEGER, PARAMETER :: equation_n = 1
  INTEGER, PARAMETER :: equation_nu = 2
  INTEGER, PARAMETER :: equation_nEi = 3
  INTEGER, PARAMETER :: equation_nEe = 4
  INTEGER, PARAMETER :: equation_nn = 5
  INTEGER, PARAMETER :: balance_equation_count = 5

  ! Each accumulated value is addressed as values(equation,term,section). Terms
  ! describe physics and are reused between the physical and independent-BC
  ! sections; this avoids separate numeric namespaces for the same flux kind.
  INTEGER, PARAMETER :: section_physical = 1
  INTEGER, PARAMETER :: section_bc = 2
  INTEGER, PARAMETER :: balance_section_count = 2

  ! Common aggregate fields.
  INTEGER, PARAMETER :: term_content = 1
  INTEGER, PARAMETER :: term_temporal = 2
  INTEGER, PARAMETER :: term_volume = 3
  INTEGER, PARAMETER :: term_boundary_physical_inward = 4
  INTEGER, PARAMETER :: term_tau_inward = 5

  ! Detailed volume fields. A field is used only for equations to which its
  ! physical meaning applies.
  INTEGER, PARAMETER :: term_ionization = 6
  INTEGER, PARAMETER :: term_recombination = 7
  INTEGER, PARAMETER :: term_prescribed_source = 8
  INTEGER, PARAMETER :: term_puff = 9
  INTEGER, PARAMETER :: term_pump = 10
  INTEGER, PARAMETER :: term_charge_exchange = 11

  ! Detailed exterior physical-flux fields.
  INTEGER, PARAMETER :: term_parallel = 12
  INTEGER, PARAMETER :: term_diffusion = 13
  INTEGER, PARAMETER :: term_pinch = 14
  INTEGER, PARAMETER :: term_pressure = 15
  INTEGER, PARAMETER :: term_convection = 16

  ! Independent boundary-condition fields.
  INTEGER, PARAMETER :: term_recycling_parallel = 17
  INTEGER, PARAMETER :: term_recycling_diffusion = 18
  INTEGER, PARAMETER :: term_recycling_pinch = 19

  ! Detailed non-boundary terms used by momentum and energy equations.
  INTEGER, PARAMETER :: term_parallel_electric_work = 20
  INTEGER, PARAMETER :: term_temperature_exchange = 21
  INTEGER, PARAMETER :: term_pressure_divergence = 22
  INTEGER, PARAMETER :: term_radiation = 23
  INTEGER, PARAMETER :: term_ohmic = 24
  INTEGER, PARAMETER :: balance_term_count = 24

  ! Evaluated physical inputs shared by the particle, momentum, and energy
  ! source mappings. Equation ownership and signs remain private below.
  TYPE, PUBLIC :: balance_atomic_volume_type
     REAL*8 :: ionization_density = 0.d0
     REAL*8 :: recombination_density = 0.d0
     REAL*8 :: ionization_rate_coefficient = 0.d0
     REAL*8 :: recombination_rate_coefficient = 0.d0
     REAL*8 :: charge_exchange_rate_coefficient = 0.d0
  END TYPE balance_atomic_volume_type

  TYPE, PUBLIC :: balance_momentum_volume_type
     REAL*8 :: pressure_divergence = 0.d0
     REAL*8 :: charge_exchange_factor = 0.d0
     REAL*8 :: recombination_factor = 0.d0
     REAL*8 :: neutral_factor = 0.d0
  END TYPE balance_momentum_volume_type

  TYPE, PUBLIC :: balance_energy_volume_type
     REAL*8 :: ion_ionization_factor = 0.d0
     REAL*8 :: ion_recombination_factor = 0.d0
     REAL*8 :: ion_charge_exchange_factor = 0.d0
     REAL*8 :: neutral_factor = 0.d0
     REAL*8 :: electron_ionization_loss_coefficient = 0.d0
     REAL*8 :: electron_recombination_loss_coefficient = 0.d0
     REAL*8 :: recombination_energy = 0.d0
     REAL*8 :: impurity_cooling_factor = 0.d0
     REAL*8 :: ohmic_heating = 0.d0
     REAL*8 :: parallel_electric_transfer = 0.d0
     REAL*8 :: temperature_transfer = 0.d0
  END TYPE balance_energy_volume_type

  ! Read-only equation results shared by terminal and HDF5 presenters.
  TYPE :: physical_balance_type
     REAL*8 :: content = 0.d0
     REAL*8 :: temporal = 0.d0
     REAL*8 :: volume = 0.d0
     REAL*8 :: boundary_inward = 0.d0
     REAL*8 :: physical_imbalance = 0.d0
     REAL*8 :: tau_inward = 0.d0
     REAL*8 :: discrete_residual = 0.d0
  END TYPE physical_balance_type

  TYPE :: particle_bc_balance_type
     REAL*8 :: imposed_source_inward = 0.d0
     REAL*8 :: physical_flux_inward = 0.d0
     REAL*8 :: tau_inward = 0.d0
     REAL*8 :: residual = 0.d0
  END TYPE particle_bc_balance_type

  ! One instance is owned by each OpenMP worker. The equation-oriented matrix
  ! is contiguous and can be merged and reduced without packing.
  TYPE, PUBLIC :: balance_accumulator_type
     PRIVATE
     REAL*8 :: values(balance_equation_count,balance_term_count, &
          balance_section_count) = 0.d0
     LOGICAL :: detailed = .FALSE.
     REAL*8 :: content_scale(balance_equation_count) = 0.d0
     REAL*8 :: rate_scale(balance_equation_count) = 0.d0
   CONTAINS
     PROCEDURE, PUBLIC :: accumulate_volume
     PROCEDURE, PUBLIC :: accumulate_relocated_sources
     PROCEDURE, PUBLIC :: accumulate_particle_face
     PROCEDURE, PUBLIC :: accumulate_particle_tau
     PROCEDURE, PUBLIC :: accumulate_particle_bc
     PROCEDURE, PRIVATE :: reset => accumulator_reset
  END TYPE balance_accumulator_type

  ! Global facade: configuration, thread merge, one MPI reduction, and output
  ! ownership. Accumulation and output implementations live in submodules.
  TYPE, PUBLIC :: balance_diagnostics_type
     PRIVATE
     INTEGER :: mode = balance_mode_off
     LOGICAL :: finalized_values_valid = .FALSE.
     REAL*8 :: content_scale(balance_equation_count) = 0.d0
     REAL*8 :: rate_scale(balance_equation_count) = 0.d0
     REAL*8 :: local_values(balance_equation_count,balance_term_count, &
          balance_section_count) = 0.d0
     REAL*8 :: finalized_values(balance_equation_count,balance_term_count, &
          balance_section_count) = 0.d0
   CONTAINS
     PROCEDURE, PUBLIC :: configure => balance_diagnostics_configure
     PROCEDURE, PUBLIC :: get_mode => balance_diagnostics_get_mode
     PROCEDURE, PUBLIC :: enabled => balance_diagnostics_enabled
     PROCEDURE, PUBLIC :: detailed => balance_diagnostics_detailed
     PROCEDURE, PUBLIC :: begin_assembly => balance_diagnostics_begin_assembly
     PROCEDURE, PUBLIC :: initialize_accumulator => &
          balance_diagnostics_initialize_accumulator
     PROCEDURE, PUBLIC :: merge => balance_diagnostics_merge
     PROCEDURE, PUBLIC :: finalize => balance_diagnostics_finalize
     PROCEDURE, PUBLIC :: has_values => balance_diagnostics_has_values
     PROCEDURE, PUBLIC :: report => balance_diagnostics_report
     PROCEDURE, PUBLIC :: write_hdf5 => balance_diagnostics_write_hdf5
  END TYPE balance_diagnostics_type

  TYPE(balance_diagnostics_type), PUBLIC, SAVE :: balance_diag

  PUBLIC :: parse_balance_diagnostics_mode
  PUBLIC :: balance_diagnostics_mode_name

  INTERFACE
     MODULE SUBROUTINE accumulator_reset(this, detailed, content_scale, &
          &rate_scale)
       CLASS(balance_accumulator_type), INTENT(INOUT) :: this
       LOGICAL, INTENT(IN) :: detailed
       REAL*8, INTENT(IN) :: content_scale(balance_equation_count)
       REAL*8, INTENT(IN) :: rate_scale(balance_equation_count)
     END SUBROUTINE accumulator_reset

     MODULE SUBROUTINE accumulate_volume(this, measure, state, history, &
          &time_coefficients, time_step, steady, prescribed_source, &
          &neutral_state_index, atomic, momentum, energy)
       CLASS(balance_accumulator_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: measure, state(:), history(:,:)
       REAL*8, INTENT(IN) :: time_coefficients(:), time_step
       LOGICAL, INTENT(IN) :: steady
       REAL*8, INTENT(IN) :: prescribed_source(:)
       INTEGER, INTENT(IN) :: neutral_state_index
       TYPE(balance_atomic_volume_type), INTENT(IN) :: atomic
       TYPE(balance_momentum_volume_type), INTENT(IN), OPTIONAL :: momentum
       TYPE(balance_energy_volume_type), INTENT(IN), OPTIONAL :: energy
     END SUBROUTINE accumulate_volume

     MODULE SUBROUTINE accumulate_relocated_sources(this, puff_source, &
          &pump_sink)
       CLASS(balance_accumulator_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: puff_source, pump_sink
     END SUBROUTINE accumulate_relocated_sources

     MODULE SUBROUTINE accumulate_particle_face(this, measure, &
          &plasma_equation, neutral_equation, trace_state, flux_jacobian, &
          &pinch_matrix, gradient, normal, magnetic_direction, diffusion_iso, &
          &diffusion_ani, neutral_perpendicular_diffusion, &
          &neutral_pressure_vector)
       CLASS(balance_accumulator_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: measure
       INTEGER, INTENT(IN) :: plasma_equation, neutral_equation
       REAL*8, INTENT(IN) :: trace_state(:), flux_jacobian(:,:)
       REAL*8, INTENT(IN) :: pinch_matrix(:,:), gradient(:,:), normal(:)
       REAL*8, INTENT(IN) :: magnetic_direction(:), diffusion_iso(:,:)
       REAL*8, INTENT(IN) :: diffusion_ani(:,:)
       LOGICAL, INTENT(IN) :: neutral_perpendicular_diffusion
       REAL*8, INTENT(IN), OPTIONAL :: neutral_pressure_vector(:)
     END SUBROUTINE accumulate_particle_face

     MODULE SUBROUTINE accumulate_particle_tau(this, measure, &
          &plasma_tau_inward, neutral_tau_inward)
       CLASS(balance_accumulator_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: measure, plasma_tau_inward, neutral_tau_inward
     END SUBROUTINE accumulate_particle_tau

     MODULE SUBROUTINE accumulate_particle_bc(this, integration_weight, &
          &density_equation, neutral_equation, trace_state, exterior_state, &
          &gradient, normal, magnetic_direction, magnetic_normal, tau, &
          &diffusion_iso, diffusion_ani, pinch_matrix, flux_jacobian, &
          &recycling_coefficient, puff_source, pump_coefficient, &
          &neutral_perpendicular_diffusion, neutral_pressure_vector, &
          &neutral_momentum_equation)
       CLASS(balance_accumulator_type), INTENT(INOUT) :: this
       REAL*8, INTENT(IN) :: integration_weight
       INTEGER, INTENT(IN) :: density_equation, neutral_equation
       REAL*8, INTENT(IN) :: trace_state(:), exterior_state(:)
       REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
       REAL*8, INTENT(IN) :: magnetic_normal, tau(:,:)
       REAL*8, INTENT(IN) :: diffusion_iso(:,:), diffusion_ani(:,:)
       REAL*8, INTENT(IN) :: pinch_matrix(:,:), flux_jacobian(:,:)
       REAL*8, INTENT(IN) :: recycling_coefficient, puff_source
       REAL*8, INTENT(IN) :: pump_coefficient
       LOGICAL, INTENT(IN) :: neutral_perpendicular_diffusion
       REAL*8, INTENT(IN), OPTIONAL :: neutral_pressure_vector(:)
       INTEGER, INTENT(IN), OPTIONAL :: neutral_momentum_equation
     END SUBROUTINE accumulate_particle_bc

     MODULE SUBROUTINE balance_diagnostics_report(this)
       CLASS(balance_diagnostics_type), INTENT(IN) :: this
     END SUBROUTINE balance_diagnostics_report

     MODULE SUBROUTINE balance_diagnostics_write_hdf5(this, &
          &diagnostics_group_id)
       CLASS(balance_diagnostics_type), INTENT(IN) :: this
       INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
     END SUBROUTINE balance_diagnostics_write_hdf5

     MODULE FUNCTION balance_value(this, equation, term, section) RESULT(value)
       CLASS(balance_diagnostics_type), INTENT(IN) :: this
       INTEGER, INTENT(IN) :: equation, term, section
       REAL*8 :: value
     END FUNCTION balance_value

     MODULE FUNCTION physical_balance(this, equation) RESULT(balance)
       CLASS(balance_diagnostics_type), INTENT(IN) :: this
       INTEGER, INTENT(IN) :: equation
       TYPE(physical_balance_type) :: balance
     END FUNCTION physical_balance

     MODULE FUNCTION add_physical_balances(left, right) RESULT(total)
       TYPE(physical_balance_type), INTENT(IN) :: left, right
       TYPE(physical_balance_type) :: total
     END FUNCTION add_physical_balances

     MODULE FUNCTION neutral_density_bc(this) RESULT(balance)
       CLASS(balance_diagnostics_type), INTENT(IN) :: this
       TYPE(particle_bc_balance_type) :: balance
     END FUNCTION neutral_density_bc

     MODULE FUNCTION perpendicular_diffusive_flux(equation, gradient, &
          &normal, magnetic_direction, magnetic_normal, diffusion_iso, &
          &diffusion_ani) RESULT(flux)
       INTEGER, INTENT(IN) :: equation
       REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
       REAL*8, INTENT(IN) :: magnetic_normal
       REAL*8, INTENT(IN) :: diffusion_iso(:,:), diffusion_ani(:,:)
       REAL*8 :: flux
     END FUNCTION perpendicular_diffusive_flux

     MODULE FUNCTION optional_neutral_pressure_flux(gradient, normal, &
          &magnetic_direction, magnetic_normal, perpendicular_enabled, &
          &pressure_vector) RESULT(flux)
       REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
       REAL*8, INTENT(IN) :: magnetic_normal
       LOGICAL, INTENT(IN) :: perpendicular_enabled
       REAL*8, INTENT(IN), OPTIONAL :: pressure_vector(:)
       REAL*8 :: flux
     END FUNCTION optional_neutral_pressure_flux
  END INTERFACE

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

  SUBROUTINE balance_diagnostics_configure(this, mode)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: mode

    IF (mode < balance_mode_off .OR. mode > balance_mode_detailed) THEN
       WRITE (6, *) 'Invalid balance diagnostics mode id: ', mode
       STOP
    ENDIF
    this%mode = mode
    this%finalized_values_valid = .FALSE.
  END SUBROUTINE balance_diagnostics_configure

  INTEGER FUNCTION balance_diagnostics_get_mode(this) RESULT(mode)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    mode = this%mode
  END FUNCTION balance_diagnostics_get_mode

  LOGICAL FUNCTION balance_diagnostics_enabled(this) RESULT(enabled)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    enabled = this%mode /= balance_mode_off
  END FUNCTION balance_diagnostics_enabled

  LOGICAL FUNCTION balance_diagnostics_detailed(this) RESULT(detailed)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    detailed = this%mode == balance_mode_detailed
  END FUNCTION balance_diagnostics_detailed

  SUBROUTINE balance_diagnostics_begin_assembly(this, reference_density, &
       &reference_length, reference_speed, reference_mass)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: reference_density, reference_length
    REAL*8, INTENT(IN) :: reference_speed, reference_mass
    REAL*8 :: particle_content_scale, particle_rate_scale

    IF (.NOT. this%enabled()) RETURN

    particle_content_scale = &
         2.d0*ACOS(-1.d0)*reference_density*reference_length**3
    particle_rate_scale = &
         particle_content_scale*reference_speed/reference_length
    this%content_scale = 0.d0
    this%rate_scale = 0.d0
    this%content_scale(equation_n) = particle_content_scale
    this%rate_scale(equation_n) = particle_rate_scale
    this%content_scale(equation_nu) = &
         particle_content_scale*reference_mass*reference_speed
    this%rate_scale(equation_nu) = &
         this%content_scale(equation_nu)*reference_speed/reference_length
    this%content_scale(equation_nEi:equation_nEe) = &
         particle_content_scale*reference_mass*reference_speed**2
    this%rate_scale(equation_nEi:equation_nEe) = &
         this%content_scale(equation_nEi:equation_nEe)* &
         reference_speed/reference_length
    this%content_scale(equation_nn) = particle_content_scale
    this%rate_scale(equation_nn) = particle_rate_scale
    this%local_values = 0.d0
    this%finalized_values_valid = .FALSE.
  END SUBROUTINE balance_diagnostics_begin_assembly

  SUBROUTINE balance_diagnostics_initialize_accumulator(this, accumulator)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(balance_accumulator_type), INTENT(OUT) :: accumulator

    CALL accumulator%reset(this%detailed(),this%content_scale,this%rate_scale)
  END SUBROUTINE balance_diagnostics_initialize_accumulator

  SUBROUTINE balance_diagnostics_merge(this, accumulator)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    TYPE(balance_accumulator_type), INTENT(IN) :: accumulator

    this%local_values = this%local_values + accumulator%values
    this%finalized_values_valid = .FALSE.
  END SUBROUTINE balance_diagnostics_merge

  SUBROUTINE balance_diagnostics_finalize(this)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
#ifdef PARALL
    INTEGER :: ierr
#endif

    IF (.NOT. this%enabled()) RETURN

    this%finalized_values = this%local_values
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, this%finalized_values, &
         balance_equation_count*balance_term_count*balance_section_count, &
         MPI_REAL8, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
#endif
    this%finalized_values_valid = .TRUE.
  END SUBROUTINE balance_diagnostics_finalize

  LOGICAL FUNCTION balance_diagnostics_has_values(this) RESULT(has_values)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    has_values = this%enabled() .AND. this%finalized_values_valid
  END FUNCTION balance_diagnostics_has_values

END MODULE balance_diagnostics
