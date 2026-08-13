!**********************************************************************
! Central lifecycle and accumulation for balance diagnostics.
!
! Physical rates use the convention "positive into the domain". Global
! particle diagnostics retain the physical flux and the stabilization term
! from the element numerical flux separately. Wall/BC closure uses its own
! boundary-condition evaluation of those terms.
!**********************************************************************
MODULE balance_diagnostics
  USE HDF5, ONLY: HID_T
  USE HDF5_io_module, ONLY: HDF5_group_create, HDF5_group_close, &
       &HDF5_real_saving, HDF5_string_saving
  USE MPI_OMP
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: balance_mode_invalid = -1
  INTEGER, PARAMETER, PUBLIC :: balance_mode_off = 0
  INTEGER, PARAMETER, PUBLIC :: balance_mode_summary = 1
  INTEGER, PARAMETER, PUBLIC :: balance_mode_detailed = 2

  INTEGER, PARAMETER :: plasma_slot = 0
  INTEGER, PARAMETER :: neutral_slot = 1
  INTEGER, PARAMETER :: species_count = 2

  INTEGER, PARAMETER :: particle_content_start = 1
  INTEGER, PARAMETER :: particle_content_plasma = &
       particle_content_start + plasma_slot
  INTEGER, PARAMETER :: particle_content_neutral = &
       particle_content_start + neutral_slot

  INTEGER, PARAMETER :: particle_temporal_start = &
       particle_content_start + species_count
  INTEGER, PARAMETER :: particle_temporal_plasma = &
       particle_temporal_start + plasma_slot
  INTEGER, PARAMETER :: particle_temporal_neutral = &
       particle_temporal_start + neutral_slot

  INTEGER, PARAMETER :: particle_volume_start = &
       particle_temporal_start + species_count
  INTEGER, PARAMETER :: particle_volume_plasma = &
       particle_volume_start + plasma_slot
  INTEGER, PARAMETER :: particle_volume_neutral = &
       particle_volume_start + neutral_slot

  INTEGER, PARAMETER :: particle_boundary_start = &
       particle_volume_start + species_count
  INTEGER, PARAMETER :: particle_boundary_plasma = &
       particle_boundary_start + plasma_slot
  INTEGER, PARAMETER :: particle_boundary_neutral = &
       particle_boundary_start + neutral_slot

  INTEGER, PARAMETER :: particle_tau_start = &
       particle_boundary_start + species_count
  INTEGER, PARAMETER :: particle_tau_plasma = particle_tau_start + plasma_slot
  INTEGER, PARAMETER :: particle_tau_neutral = &
       particle_tau_start + neutral_slot

  INTEGER, PARAMETER :: particle_volume_components_start = &
       particle_tau_start + species_count
  INTEGER, PARAMETER :: particle_plasma_ionization = &
       particle_volume_components_start
  INTEGER, PARAMETER :: particle_plasma_recombination = &
       particle_plasma_ionization + 1
  INTEGER, PARAMETER :: particle_plasma_other_source = &
       particle_plasma_recombination + 1
  INTEGER, PARAMETER :: particle_neutral_ionization = &
       particle_plasma_other_source + 1
  INTEGER, PARAMETER :: particle_neutral_recombination = &
       particle_neutral_ionization + 1
  INTEGER, PARAMETER :: particle_neutral_other_source = &
       particle_neutral_recombination + 1
  INTEGER, PARAMETER :: particle_neutral_puff_source = &
       particle_neutral_other_source + 1
  INTEGER, PARAMETER :: particle_neutral_pump_source = &
       particle_neutral_puff_source + 1
  INTEGER, PARAMETER :: particle_charge_exchange = &
       particle_neutral_pump_source + 1

  INTEGER, PARAMETER :: particle_boundary_components_start = &
       particle_charge_exchange + 1
  INTEGER, PARAMETER :: particle_plasma_parallel_flux = &
       particle_boundary_components_start
  INTEGER, PARAMETER :: particle_plasma_diffusion_flux = &
       particle_plasma_parallel_flux + 1
  INTEGER, PARAMETER :: particle_plasma_pinch_flux = &
       particle_plasma_diffusion_flux + 1
  INTEGER, PARAMETER :: particle_neutral_diffusion_flux = &
       particle_plasma_pinch_flux + 1
  INTEGER, PARAMETER :: particle_neutral_pressure_flux = &
       particle_neutral_diffusion_flux + 1
  INTEGER, PARAMETER :: particle_neutral_convection_flux = &
       particle_neutral_pressure_flux + 1

  INTEGER, PARAMETER :: wall_plasma_diffusion_flux = &
       particle_neutral_convection_flux + 1
  INTEGER, PARAMETER :: wall_plasma_tau_flux = &
       wall_plasma_diffusion_flux + 1
  INTEGER, PARAMETER :: wall_recycling_parallel_source = &
       wall_plasma_tau_flux + 1
  INTEGER, PARAMETER :: wall_recycling_diffusion_source = &
       wall_recycling_parallel_source + 1
  INTEGER, PARAMETER :: wall_recycling_pinch_source = &
       wall_recycling_diffusion_source + 1
  INTEGER, PARAMETER :: wall_puff_source = &
       wall_recycling_pinch_source + 1
  INTEGER, PARAMETER :: wall_pump_source = wall_puff_source + 1
  INTEGER, PARAMETER :: wall_neutral_diffusion_flux = wall_pump_source + 1
  INTEGER, PARAMETER :: wall_neutral_pressure_flux = &
       wall_neutral_diffusion_flux + 1
  INTEGER, PARAMETER :: wall_neutral_convection_flux = &
       wall_neutral_pressure_flux + 1
  INTEGER, PARAMETER :: wall_neutral_tau_flux = &
       wall_neutral_convection_flux + 1

  INTEGER, PARAMETER :: balance_diagnostics_value_count = &
       wall_neutral_tau_flux

  ! Nondimensional flux densities evaluated from one owned Bohm-boundary
  ! Gauss point. Sources are positive into the neutral population; all
  ! physical and tau fluxes are positive into the computational domain.
  TYPE, PUBLIC :: particle_wall_flux_type
     REAL*8 :: plasma_diffusion_inward = 0.d0
     REAL*8 :: plasma_tau_inward = 0.d0
     REAL*8 :: recycling_parallel_source = 0.d0
     REAL*8 :: recycling_diffusion_source = 0.d0
     REAL*8 :: recycling_pinch_source = 0.d0
     REAL*8 :: puff_source = 0.d0
     REAL*8 :: pump_sink = 0.d0
     REAL*8 :: neutral_diffusion_inward = 0.d0
     REAL*8 :: neutral_pressure_inward = 0.d0
     REAL*8 :: neutral_convection_inward = 0.d0
     REAL*8 :: neutral_tau_inward = 0.d0
  END TYPE particle_wall_flux_type

  ! Final dimensional terms of the two particle wall equations. Keeping the
  ! equation terms together avoids positional scalar interfaces in reporters
  ! and future HDF5 category writers.
  TYPE :: particle_wall_balance_type
     REAL*8 :: plasma_diffusion = 0.d0
     REAL*8 :: plasma_tau = 0.d0
     REAL*8 :: plasma_residual = 0.d0
     REAL*8 :: recycled_plasma = 0.d0
     REAL*8 :: puff_source = 0.d0
     REAL*8 :: pump_sink = 0.d0
     REAL*8 :: neutral_physical_flux = 0.d0
     REAL*8 :: neutral_tau = 0.d0
     REAL*8 :: neutral_residual = 0.d0
  END TYPE particle_wall_balance_type

  ! Thread-local, fixed-size accumulator for every diagnostic category. It
  ! hides term indices, dimensional scaling, signs, and component bookkeeping
  ! from assembly code and is merged once per OpenMP worker.
  TYPE, PUBLIC :: balance_accumulator_type
     PRIVATE
     REAL*8 :: values(balance_diagnostics_value_count) = 0.d0
     LOGICAL :: detailed = .FALSE.
     REAL*8 :: content_scale = 0.d0
     REAL*8 :: rate_scale = 0.d0
   CONTAINS
     PROCEDURE, PUBLIC :: accumulate_particle_volume
     PROCEDURE, PUBLIC :: accumulate_relocated_sources
     PROCEDURE, PUBLIC :: accumulate_particle_face
     PROCEDURE, PUBLIC :: accumulate_particle_tau
     PROCEDURE, PUBLIC :: accumulate_particle_wall
     PROCEDURE, PRIVATE :: reset => accumulator_reset
     PROCEDURE, PRIVATE :: add_content
     PROCEDURE, PRIVATE :: add_temporal
     PROCEDURE, PRIVATE :: add_reactions
     PROCEDURE, PRIVATE :: add_sources
     PROCEDURE, PRIVATE :: add_plasma_flux
     PROCEDURE, PRIVATE :: add_neutral_flux
  END TYPE balance_accumulator_type

  ! Runtime diagnostics service. It owns mode handling, assembly lifecycle,
  ! rank-local and finalized arrays, and the single MPI reduction.
  TYPE, PUBLIC :: balance_diagnostics_type
     PRIVATE
     INTEGER :: mode = balance_mode_off
     LOGICAL :: finalized_values_valid = .FALSE.
     REAL*8 :: content_scale = 0.d0
     REAL*8 :: rate_scale = 0.d0
     REAL*8 :: local_values(balance_diagnostics_value_count) = 0.d0
     REAL*8 :: finalized_values(balance_diagnostics_value_count) = 0.d0
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
       &reference_length, reference_speed)
    CLASS(balance_diagnostics_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: reference_density, reference_length
    REAL*8, INTENT(IN) :: reference_speed

    IF (.NOT. this%enabled()) RETURN

    this%content_scale = 2.d0*ACOS(-1.d0)*reference_density*reference_length**3
    this%rate_scale = this%content_scale*reference_speed/reference_length
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
         balance_diagnostics_value_count, MPI_REAL8, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
#endif
    this%finalized_values_valid = .TRUE.
  END SUBROUTINE balance_diagnostics_finalize

  LOGICAL FUNCTION balance_diagnostics_has_values(this) RESULT(has_values)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    has_values = this%enabled() .AND. this%finalized_values_valid
  END FUNCTION balance_diagnostics_has_values

  SUBROUTINE balance_diagnostics_report(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    IF (.NOT. this%has_values()) RETURN
    IF (MPIvar%glob_id /= 0) RETURN

    SELECT CASE (this%mode)
    CASE (balance_mode_summary)
       CALL print_summary(this)
    CASE (balance_mode_detailed)
       CALL print_detailed(this)
    END SELECT
  END SUBROUTINE balance_diagnostics_report

  SUBROUTINE balance_diagnostics_write_hdf5(this, diagnostics_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: particles_group_id, wall_group_id
    INTEGER :: ierr

    IF (.NOT. this%has_values()) RETURN
    IF (MPIvar%glob_id /= 0) RETURN

    CALL HDF5_group_create('particles',diagnostics_group_id, &
         particles_group_id,ierr)
    CALL write_particle_summary_hdf5(this,particles_group_id)
    IF (this%detailed()) &
         CALL write_particle_detail_hdf5(this,particles_group_id)
    CALL HDF5_group_close(particles_group_id,ierr)

    CALL HDF5_group_create('wall_closure',diagnostics_group_id,wall_group_id, &
         ierr)
    CALL write_wall_hdf5(this,wall_group_id)
    CALL HDF5_group_close(wall_group_id,ierr)
  END SUBROUTINE balance_diagnostics_write_hdf5

  SUBROUTINE write_particle_summary_hdf5(this, particles_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: particles_group_id
    REAL*8 :: content(3), temporal(3), volume(3), boundary(3)
    REAL*8 :: physical_imbalance(3), tau(3), discrete_residual(3)

    CALL particle_aggregates(this,content,temporal,volume,boundary, &
         physical_imbalance,tau,discrete_residual)
    CALL write_species_group(particles_group_id,'content',content,'particles')
    CALL write_species_group(particles_group_id,'conservation', &
         discrete_residual,'particles/s')
    CALL write_species_group(particles_group_id,'physical_imbalance', &
         physical_imbalance,'particles/s')
    CALL write_species_group(particles_group_id,'hdg_tau_inward',tau, &
         'particles/s')
  END SUBROUTINE write_particle_summary_hdf5

  SUBROUTINE write_particle_detail_hdf5(this, particles_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: particles_group_id
    INTEGER(HID_T) :: components_group_id, exchange_group_id
    REAL*8 :: content(3), temporal(3), volume(3), boundary(3)
    REAL*8 :: physical_imbalance(3), tau(3), discrete_residual(3)
    INTEGER :: ierr

    CALL particle_aggregates(this,content,temporal,volume,boundary, &
         physical_imbalance,tau,discrete_residual)
    CALL write_species_group(particles_group_id,'temporal',temporal, &
         'particles/s')
    CALL write_species_group(particles_group_id,'volume',volume, &
         'particles/s')
    CALL write_species_group(particles_group_id,'boundary_physical_inward', &
         boundary,'particles/s')

    CALL HDF5_group_create('components',particles_group_id, &
         components_group_id,ierr)
    CALL write_particle_components_hdf5(this,components_group_id)
    CALL HDF5_group_close(components_group_id,ierr)

    CALL HDF5_group_create('exchange',particles_group_id,exchange_group_id, &
         ierr)
    CALL HDF5_string_saving(exchange_group_id,'particles/s','units')
    CALL HDF5_real_saving(exchange_group_id, &
         this%finalized_values(particle_charge_exchange),'charge_exchange')
    CALL HDF5_group_close(exchange_group_id,ierr)
  END SUBROUTINE write_particle_detail_hdf5

  SUBROUTINE write_species_group(parent_group_id, name, values, units)
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    REAL*8, INTENT(IN) :: values(3)
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create(name,parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    CALL HDF5_real_saving(group_id,values(1),'plasma')
    CALL HDF5_real_saving(group_id,values(2),'neutral')
    CALL HDF5_real_saving(group_id,values(3),'total')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_species_group

  SUBROUTINE write_particle_components_hdf5(this, components_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: components_group_id
    INTEGER(HID_T) :: plasma_group_id, neutral_group_id
    INTEGER(HID_T) :: volume_group_id, boundary_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('plasma',components_group_id,plasma_group_id,ierr)
    CALL HDF5_group_create('volume',plasma_group_id,volume_group_id,ierr)
    CALL HDF5_string_saving(volume_group_id,'particles/s','units')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_plasma_ionization),'ionization')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_plasma_recombination),'recombination')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_plasma_other_source),'other_source')
    CALL HDF5_group_close(volume_group_id,ierr)
    CALL HDF5_group_create('boundary_inward',plasma_group_id, &
         boundary_group_id,ierr)
    CALL HDF5_string_saving(boundary_group_id,'particles/s','units')
    CALL HDF5_real_saving(boundary_group_id, &
         this%finalized_values(particle_plasma_parallel_flux),'parallel')
    CALL HDF5_real_saving(boundary_group_id, &
         this%finalized_values(particle_plasma_diffusion_flux),'diffusion')
    CALL HDF5_real_saving(boundary_group_id, &
         this%finalized_values(particle_plasma_pinch_flux),'pinch')
    CALL HDF5_group_close(boundary_group_id,ierr)
    CALL HDF5_group_close(plasma_group_id,ierr)

    CALL HDF5_group_create('neutral',components_group_id,neutral_group_id,ierr)
    CALL HDF5_group_create('volume',neutral_group_id,volume_group_id,ierr)
    CALL HDF5_string_saving(volume_group_id,'particles/s','units')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_neutral_ionization),'ionization')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_neutral_recombination),'recombination')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_neutral_other_source),'other_source')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_neutral_puff_source),'puff_source')
    CALL HDF5_real_saving(volume_group_id, &
         this%finalized_values(particle_neutral_pump_source),'pump_source')
    CALL HDF5_group_close(volume_group_id,ierr)
    CALL HDF5_group_create('boundary_inward',neutral_group_id, &
         boundary_group_id,ierr)
    CALL HDF5_string_saving(boundary_group_id,'particles/s','units')
    CALL HDF5_real_saving(boundary_group_id, &
         this%finalized_values(particle_neutral_diffusion_flux), &
         'limited_diffusion')
#ifdef NEUTRALP
    CALL HDF5_real_saving(boundary_group_id, &
         this%finalized_values(particle_neutral_pressure_flux), &
         'limited_pressure')
#endif
#ifdef NEUTRALGAMMA
    CALL HDF5_real_saving(boundary_group_id, &
         this%finalized_values(particle_neutral_convection_flux), &
         'neutral_gamma_convection')
#elif defined(NEUTRALCONVECTION)
    CALL HDF5_real_saving(boundary_group_id, &
         this%finalized_values(particle_neutral_convection_flux), &
         'neutral_convection')
#endif
    CALL HDF5_group_close(boundary_group_id,ierr)
    CALL HDF5_group_close(neutral_group_id,ierr)
  END SUBROUTINE write_particle_components_hdf5

  SUBROUTINE write_wall_hdf5(this, wall_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: wall_group_id
    TYPE(particle_wall_balance_type) :: wall
    INTEGER(HID_T) :: plasma_group_id, neutral_group_id
    INTEGER :: ierr

    wall = compute_wall_balance(this)
    CALL HDF5_group_create('plasma_particles',wall_group_id, &
         plasma_group_id,ierr)
    CALL HDF5_string_saving(plasma_group_id,'particles/s','units')
    CALL HDF5_real_saving(plasma_group_id,wall%plasma_diffusion, &
         'diffusion_inward')
    CALL HDF5_real_saving(plasma_group_id,wall%plasma_tau, &
         'stabilization_inward')
    CALL HDF5_real_saving(plasma_group_id,wall%plasma_residual,'residual')
    CALL HDF5_group_close(plasma_group_id,ierr)

    CALL HDF5_group_create('neutral',wall_group_id,neutral_group_id,ierr)
    CALL HDF5_string_saving(neutral_group_id,'particles/s','units')
    CALL HDF5_real_saving(neutral_group_id, &
         wall%recycled_plasma+wall%puff_source-wall%pump_sink,'source_inward')
    CALL HDF5_real_saving(neutral_group_id,wall%neutral_physical_flux, &
         'physical_flux_inward')
    CALL HDF5_real_saving(neutral_group_id,wall%neutral_tau, &
         'stabilization_inward')
    CALL HDF5_real_saving(neutral_group_id,wall%neutral_residual,'residual')
    IF (this%detailed()) THEN
       CALL HDF5_real_saving(neutral_group_id,wall%recycled_plasma, &
            'recycled_plasma_inward')
       CALL HDF5_real_saving(neutral_group_id,wall%puff_source,'puff_source')
       CALL HDF5_real_saving(neutral_group_id,wall%pump_sink,'pump_sink')
       CALL write_wall_components_hdf5(this,neutral_group_id)
    ENDIF
    CALL HDF5_group_close(neutral_group_id,ierr)
  END SUBROUTINE write_wall_hdf5

  SUBROUTINE write_wall_components_hdf5(this, neutral_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: neutral_group_id
    INTEGER(HID_T) :: physical_group_id, recycling_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('physical_flux_components',neutral_group_id, &
         physical_group_id,ierr)
    CALL HDF5_string_saving(physical_group_id,'particles/s','units')
    CALL HDF5_real_saving(physical_group_id, &
         this%finalized_values(wall_neutral_diffusion_flux), &
         'limited_diffusion_inward')
#ifdef NEUTRALP
    CALL HDF5_real_saving(physical_group_id, &
         this%finalized_values(wall_neutral_pressure_flux), &
         'limited_pressure_inward')
#endif
#ifdef NEUTRALGAMMA
    CALL HDF5_real_saving(physical_group_id, &
         this%finalized_values(wall_neutral_convection_flux), &
         'neutral_gamma_inward')
#elif defined(NEUTRALCONVECTION)
    CALL HDF5_real_saving(physical_group_id, &
         this%finalized_values(wall_neutral_convection_flux), &
         'neutral_convection_inward')
#endif
    CALL HDF5_group_close(physical_group_id,ierr)

    CALL HDF5_group_create('recycled_plasma_components',neutral_group_id, &
         recycling_group_id,ierr)
    CALL HDF5_string_saving(recycling_group_id,'particles/s','units')
    CALL HDF5_real_saving(recycling_group_id, &
         this%finalized_values(wall_recycling_parallel_source), &
         'parallel_source')
    CALL HDF5_real_saving(recycling_group_id, &
         this%finalized_values(wall_recycling_diffusion_source), &
         'diffusion_source')
    CALL HDF5_real_saving(recycling_group_id, &
         this%finalized_values(wall_recycling_pinch_source),'pinch_source')
    CALL HDF5_group_close(recycling_group_id,ierr)
  END SUBROUTINE write_wall_components_hdf5

  SUBROUTINE particle_aggregates(this, content, temporal, volume, boundary, &
       &physical_imbalance, tau, discrete_residual)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8, INTENT(OUT) :: content(3), temporal(3), volume(3), boundary(3)
    REAL*8, INTENT(OUT) :: physical_imbalance(3), tau(3)
    REAL*8, INTENT(OUT) :: discrete_residual(3)

    content(1) = this%finalized_values(particle_content_plasma)
    content(2) = this%finalized_values(particle_content_neutral)
    temporal(1) = this%finalized_values(particle_temporal_plasma)
    temporal(2) = this%finalized_values(particle_temporal_neutral)
    volume(1) = this%finalized_values(particle_volume_plasma)
    volume(2) = this%finalized_values(particle_volume_neutral)
    boundary(1) = this%finalized_values(particle_boundary_plasma)
    boundary(2) = this%finalized_values(particle_boundary_neutral)
    tau(1) = this%finalized_values(particle_tau_plasma)
    tau(2) = this%finalized_values(particle_tau_neutral)
    content(3) = content(1) + content(2)
    temporal(3) = temporal(1) + temporal(2)
    volume(3) = volume(1) + volume(2)
    boundary(3) = boundary(1) + boundary(2)
    tau(3) = tau(1) + tau(2)
    physical_imbalance = temporal-volume-boundary
    discrete_residual = physical_imbalance-tau
  END SUBROUTINE particle_aggregates

  FUNCTION compute_wall_balance(this) RESULT(balance)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(particle_wall_balance_type) :: balance

    balance%plasma_diffusion = &
         this%finalized_values(wall_plasma_diffusion_flux)
    balance%plasma_tau = this%finalized_values(wall_plasma_tau_flux)
    balance%plasma_residual = &
         balance%plasma_diffusion+balance%plasma_tau
    balance%recycled_plasma = &
         this%finalized_values(wall_recycling_parallel_source) + &
         this%finalized_values(wall_recycling_diffusion_source) + &
         this%finalized_values(wall_recycling_pinch_source)
    balance%puff_source = this%finalized_values(wall_puff_source)
    balance%pump_sink = -this%finalized_values(wall_pump_source)
    balance%neutral_physical_flux = &
         this%finalized_values(wall_neutral_diffusion_flux) + &
         this%finalized_values(wall_neutral_pressure_flux) + &
         this%finalized_values(wall_neutral_convection_flux)
    balance%neutral_tau = this%finalized_values(wall_neutral_tau_flux)
    balance%neutral_residual = balance%recycled_plasma+ &
         balance%puff_source-balance%pump_sink- &
         balance%neutral_physical_flux-balance%neutral_tau
  END FUNCTION compute_wall_balance

  SUBROUTINE print_summary(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: content(3), temporal(3), volume(3), boundary(3)
    REAL*8 :: physical_imbalance(3), tau(3), discrete_residual(3)
    TYPE(particle_wall_balance_type) :: wall

    CALL particle_aggregates(this,content,temporal,volume,boundary, &
         physical_imbalance,tau,discrete_residual)
    wall = compute_wall_balance(this)

    WRITE(6,'(A)') 'Particle diagnostics'
    WRITE(6,'(A)') '  Content [particles]'
    WRITE(6,'(A,1X,ES12.3)') '    plasma n   ', content(1)
    WRITE(6,'(A,1X,ES12.3)') '    neutral nn ', content(2)
    WRITE(6,'(A,1X,ES12.3)') '    total n+nn ', content(3)
    WRITE(6,'(A)') '  Physical / discrete conservation [particles/s]'
    WRITE(6,'(A)') '                   physical imbalance  HDG tau inward  discrete residual'
    WRITE(6,'(A,3(1X,ES12.3))') '    plasma          ', &
         physical_imbalance(1),tau(1),discrete_residual(1)
    WRITE(6,'(A,3(1X,ES12.3))') '    neutral         ', &
         physical_imbalance(2),tau(2),discrete_residual(2)
    WRITE(6,'(A,3(1X,ES12.3))') '    total           ', &
         physical_imbalance(3),tau(3),discrete_residual(3)
    WRITE(6,'(A)') '  Wall / HDG boundary-condition residuals [particles/s]'
    WRITE(6,'(A,3(1X,ES12.3))') &
         '    plasma: diffusion, stabilization, residual', &
         wall%plasma_diffusion,wall%plasma_tau,wall%plasma_residual
    WRITE(6,'(A,4(1X,ES12.3))') &
         '    neutral: recycled+puff-pump, physical flux, stabilization, residual', &
         wall%recycled_plasma+wall%puff_source-wall%pump_sink, &
         wall%neutral_physical_flux,wall%neutral_tau,wall%neutral_residual
  END SUBROUTINE print_summary

  SUBROUTINE print_detailed(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: content(3), temporal(3), volume(3), boundary(3)
    REAL*8 :: physical_imbalance(3), tau(3), discrete_residual(3)
    TYPE(particle_wall_balance_type) :: wall

    CALL particle_aggregates(this,content,temporal,volume,boundary, &
         physical_imbalance,tau,discrete_residual)
    wall = compute_wall_balance(this)

    WRITE(6,'(A)') 'Particle diagnostics (detailed)'
    WRITE(6,'(A)') '  Content [particles]'
    CALL print_detail_value('plasma n',content(1))
    CALL print_detail_value('neutral nn',content(2))
    CALL print_detail_value('total n+nn',content(3))
    CALL print_particle_detail(this,'plasma',1,temporal(1),volume(1), &
         boundary(1),physical_imbalance(1),tau(1),discrete_residual(1))
    CALL print_particle_detail(this,'neutral',2,temporal(2),volume(2), &
         boundary(2),physical_imbalance(2),tau(2),discrete_residual(2))
    WRITE(6,'(A)') '  total conservation [particles/s]'
    CALL print_detail_value('temporal',temporal(3))
    CALL print_detail_value('volume',volume(3))
    CALL print_detail_value('boundary physical, inward',boundary(3))
    CALL print_detail_value('physical imbalance',physical_imbalance(3))
    CALL print_detail_value('HDG tau inward',tau(3))
    CALL print_detail_value('discrete residual',discrete_residual(3))
    WRITE(6,'(A)') '  exchange information [particles/s]'
    CALL print_detail_value('charge exchange', &
         this%finalized_values(particle_charge_exchange))
    WRITE(6,'(A)') '  Wall / HDG boundary-condition residuals [particles/s]'
    CALL print_wall_detail(this,wall)
  END SUBROUTINE print_detailed

  SUBROUTINE print_particle_detail(this, name, species, temporal, volume, &
       &boundary, physical_imbalance, tau, discrete_residual)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: species
    REAL*8, INTENT(IN) :: temporal, volume, boundary, physical_imbalance
    REAL*8, INTENT(IN) :: tau, discrete_residual

    WRITE(6,'(A)') '  '//TRIM(name)//' conservation [particles/s]'
    CALL print_detail_value('temporal',temporal)
    CALL print_detail_value('volume',volume)
    CALL print_detail_value('boundary physical, inward',boundary)
    CALL print_detail_value('physical imbalance',physical_imbalance)
    CALL print_detail_value('HDG tau inward',tau)
    CALL print_detail_value('discrete residual',discrete_residual)
    WRITE(6,'(A)') '    volume components'
    IF (species == 1) THEN
       CALL print_detail_value('ionization', &
            this%finalized_values(particle_plasma_ionization))
       CALL print_detail_value('recombination', &
            this%finalized_values(particle_plasma_recombination))
       CALL print_detail_value('other source', &
            this%finalized_values(particle_plasma_other_source))
       WRITE(6,'(A)') '    boundary components, inward-positive'
       CALL print_detail_value('parallel', &
            this%finalized_values(particle_plasma_parallel_flux))
       CALL print_detail_value('diffusion', &
            this%finalized_values(particle_plasma_diffusion_flux))
       CALL print_detail_value('pinch', &
            this%finalized_values(particle_plasma_pinch_flux))
    ELSE
       CALL print_detail_value('ionization', &
            this%finalized_values(particle_neutral_ionization))
       CALL print_detail_value('recombination', &
            this%finalized_values(particle_neutral_recombination))
       CALL print_detail_value('other source', &
            this%finalized_values(particle_neutral_other_source))
       CALL print_detail_value('puff source', &
            this%finalized_values(particle_neutral_puff_source))
       CALL print_detail_value('pump source', &
            this%finalized_values(particle_neutral_pump_source))
       WRITE(6,'(A)') '    boundary components, inward-positive'
       CALL print_detail_value('limited diffusion', &
            this%finalized_values(particle_neutral_diffusion_flux))
#ifdef NEUTRALP
       CALL print_detail_value('limited pressure', &
            this%finalized_values(particle_neutral_pressure_flux))
#endif
#ifdef NEUTRALGAMMA
       CALL print_detail_value('NeutralGamma convection', &
            this%finalized_values(particle_neutral_convection_flux))
#elif defined(NEUTRALCONVECTION)
       CALL print_detail_value('neutral advective flux', &
            this%finalized_values(particle_neutral_convection_flux))
#endif
    ENDIF
  END SUBROUTINE print_particle_detail

  SUBROUTINE print_wall_detail(this, wall)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(particle_wall_balance_type), INTENT(IN) :: wall

    WRITE(6,'(A)') '    Plasma density BC: diffusion + stabilization = 0'
    CALL print_detail_value('plasma diffusive flux inward', &
         wall%plasma_diffusion)
    CALL print_detail_value('HDG stabilization flux inward',wall%plasma_tau)
    CALL print_detail_value('residual',wall%plasma_residual)
    WRITE(6,'(A)') '    Neutral density BC:'
    WRITE(6,'(A)') '      recycled + puff - pump - neutral flux - stabilization = 0'
    CALL print_detail_value('neutral physical flux inward', &
         wall%neutral_physical_flux)
    WRITE(6,'(A)') '      physical-flux components, inward-positive'
    CALL print_detail_component('limited neutral diffusion', &
         this%finalized_values(wall_neutral_diffusion_flux))
#ifdef NEUTRALP
    CALL print_detail_component('limited neutral-pressure flux', &
         this%finalized_values(wall_neutral_pressure_flux))
#endif
#ifdef NEUTRALGAMMA
    CALL print_detail_component('NeutralGamma flux', &
         this%finalized_values(wall_neutral_convection_flux))
#elif defined(NEUTRALCONVECTION)
    CALL print_detail_component('neutral advective flux', &
         this%finalized_values(wall_neutral_convection_flux))
#endif
    CALL print_detail_value('recycled plasma flux into neutrals', &
         wall%recycled_plasma)
    WRITE(6,'(A)') '      recycled-flux components'
    CALL print_detail_component('parallel contribution', &
         this%finalized_values(wall_recycling_parallel_source))
    CALL print_detail_component('diffusion contribution', &
         this%finalized_values(wall_recycling_diffusion_source))
    CALL print_detail_component('pinch contribution', &
         this%finalized_values(wall_recycling_pinch_source))
    CALL print_detail_value('puff source',wall%puff_source)
    CALL print_detail_value('pump sink',wall%pump_sink)
    CALL print_detail_value('HDG stabilization flux inward',wall%neutral_tau)
    CALL print_detail_value('residual',wall%neutral_residual)
  END SUBROUTINE print_wall_detail

  SUBROUTINE print_detail_value(label, value)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: value

    WRITE(6,'(A,1X,ES12.3)') '      '//TRIM(label), value
  END SUBROUTINE print_detail_value

  SUBROUTINE print_detail_component(label, value)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: value

    WRITE(6,'(A,1X,ES12.3)') '        '//TRIM(label), value
  END SUBROUTINE print_detail_component

  SUBROUTINE accumulator_reset(this, detailed, content_scale, rate_scale)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: detailed
    REAL*8, INTENT(IN) :: content_scale, rate_scale

    this%values = 0.d0
    this%detailed = detailed
    this%content_scale = content_scale
    this%rate_scale = rate_scale
  END SUBROUTINE accumulator_reset

  SUBROUTINE accumulate_particle_volume(this, measure, plasma_density, &
       &neutral_density, plasma_history, neutral_history, time_coefficients, &
       &time_step, steady, ionization_rate, recombination_rate, &
       &plasma_other_source, neutral_other_source, charge_exchange_rate)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: measure, plasma_density, neutral_density
    REAL*8, INTENT(IN) :: plasma_history(:), neutral_history(:)
    REAL*8, INTENT(IN) :: time_coefficients(:), time_step
    LOGICAL, INTENT(IN) :: steady
    REAL*8, INTENT(IN) :: ionization_rate, recombination_rate
    REAL*8, INTENT(IN) :: plasma_other_source, neutral_other_source
    REAL*8, INTENT(IN) :: charge_exchange_rate
    REAL*8 :: content_weight, rate_weight
    REAL*8 :: plasma_temporal, neutral_temporal

    content_weight = measure*this%content_scale
    rate_weight = measure*this%rate_scale
    CALL this%add_content(plasma_density*content_weight, &
         neutral_density*content_weight)
    IF (.NOT. steady) THEN
       plasma_temporal = discrete_time_derivative(plasma_density, &
            plasma_history,time_coefficients,time_step)*rate_weight
       neutral_temporal = discrete_time_derivative(neutral_density, &
            neutral_history,time_coefficients,time_step)*rate_weight
       CALL this%add_temporal(plasma_temporal,neutral_temporal)
    ENDIF
    CALL this%add_reactions(ionization_rate*rate_weight, &
         recombination_rate*rate_weight,charge_exchange_rate*rate_weight)
    CALL this%add_sources(plasma_other_source*rate_weight, &
         neutral_other_source*rate_weight)
  END SUBROUTINE accumulate_particle_volume

  SUBROUTINE add_content(this, plasma_content, neutral_content)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_content, neutral_content

    this%values(particle_content_plasma) = &
         this%values(particle_content_plasma) + plasma_content
    this%values(particle_content_neutral) = &
         this%values(particle_content_neutral) + neutral_content
  END SUBROUTINE add_content

  SUBROUTINE add_temporal(this, plasma_temporal, neutral_temporal)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_temporal, neutral_temporal

    this%values(particle_temporal_plasma) = &
         this%values(particle_temporal_plasma) + plasma_temporal
    this%values(particle_temporal_neutral) = &
         this%values(particle_temporal_neutral) + neutral_temporal
  END SUBROUTINE add_temporal

  REAL*8 FUNCTION discrete_time_derivative(current, history, coefficients, &
       &time_step) RESULT(derivative)
    REAL*8, INTENT(IN) :: current, history(:), coefficients(:), time_step
    INTEGER :: history_index

    derivative = coefficients(1)*current
    DO history_index = 1, SIZE(history)
       derivative = derivative - coefficients(history_index + 1)* &
            history(history_index)
    ENDDO
    derivative = derivative/time_step
  END FUNCTION discrete_time_derivative

  SUBROUTINE add_reactions(this, ionization_rate, recombination_rate, &
       &charge_exchange_rate)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: ionization_rate, recombination_rate
    REAL*8, INTENT(IN) :: charge_exchange_rate
    REAL*8 :: net_reaction

    net_reaction = ionization_rate-recombination_rate
    this%values(particle_volume_plasma) = &
         this%values(particle_volume_plasma) + net_reaction
    this%values(particle_volume_neutral) = &
         this%values(particle_volume_neutral) - net_reaction
    IF (.NOT. this%detailed) RETURN

    this%values(particle_plasma_ionization) = &
         this%values(particle_plasma_ionization) + ionization_rate
    this%values(particle_plasma_recombination) = &
         this%values(particle_plasma_recombination) - recombination_rate
    this%values(particle_neutral_ionization) = &
         this%values(particle_neutral_ionization) - ionization_rate
    this%values(particle_neutral_recombination) = &
         this%values(particle_neutral_recombination) + recombination_rate
    this%values(particle_charge_exchange) = &
         this%values(particle_charge_exchange) + charge_exchange_rate
  END SUBROUTINE add_reactions

  SUBROUTINE add_sources(this, plasma_source, neutral_source)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_source, neutral_source

    this%values(particle_volume_plasma) = &
         this%values(particle_volume_plasma) + plasma_source
    this%values(particle_volume_neutral) = &
         this%values(particle_volume_neutral) + neutral_source
    IF (.NOT. this%detailed) RETURN

    this%values(particle_plasma_other_source) = &
         this%values(particle_plasma_other_source) + plasma_source
    this%values(particle_neutral_other_source) = &
         this%values(particle_neutral_other_source) + neutral_source
  END SUBROUTINE add_sources

  SUBROUTINE accumulate_relocated_sources(this, puff_source, pump_sink)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: puff_source, pump_sink

    this%values(particle_volume_neutral) = &
         this%values(particle_volume_neutral) + &
         (puff_source-pump_sink)*this%rate_scale
    IF (.NOT. this%detailed) RETURN

    this%values(particle_neutral_puff_source) = &
         this%values(particle_neutral_puff_source) + &
         puff_source*this%rate_scale
    this%values(particle_neutral_pump_source) = &
         this%values(particle_neutral_pump_source) - &
         pump_sink*this%rate_scale
  END SUBROUTINE accumulate_relocated_sources

  SUBROUTINE accumulate_particle_face(this, measure, plasma_equation, &
       &neutral_equation, trace_state, flux_jacobian, pinch_matrix, gradient, &
       &normal, magnetic_direction, diffusion_iso, diffusion_ani, &
       &neutral_perpendicular_diffusion, neutral_pressure_vector)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: measure
    INTEGER, INTENT(IN) :: plasma_equation, neutral_equation
    REAL*8, INTENT(IN) :: trace_state(:), flux_jacobian(:,:), pinch_matrix(:,:)
    REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
    REAL*8, INTENT(IN) :: diffusion_iso(:,:), diffusion_ani(:,:)
    LOGICAL, INTENT(IN) :: neutral_perpendicular_diffusion
    REAL*8, INTENT(IN), OPTIONAL :: neutral_pressure_vector(:)
    REAL*8 :: scale, magnetic_normal
    REAL*8 :: plasma_parallel, plasma_diffusion, plasma_pinch
    REAL*8 :: neutral_diffusion, neutral_pressure, neutral_convection
    REAL*8 :: neutral_pressure_gradient(SIZE(normal))

    scale = measure*this%rate_scale
    magnetic_normal = DOT_PRODUCT(magnetic_direction,normal)
    plasma_parallel = -DOT_PRODUCT(flux_jacobian(plasma_equation,:), &
         trace_state)*magnetic_normal*scale
    plasma_diffusion = (diffusion_iso(plasma_equation,plasma_equation)* &
         DOT_PRODUCT(gradient(:,plasma_equation),normal) - &
         diffusion_ani(plasma_equation,plasma_equation)*magnetic_normal* &
         DOT_PRODUCT(gradient(:,plasma_equation),magnetic_direction))*scale
    plasma_pinch = -trace_state(plasma_equation)* &
         DOT_PRODUCT(pinch_matrix(plasma_equation,:),normal)*scale
    neutral_diffusion = (diffusion_iso(neutral_equation,neutral_equation)* &
         DOT_PRODUCT(gradient(:,neutral_equation),normal) - &
         diffusion_ani(neutral_equation,neutral_equation)*magnetic_normal* &
         DOT_PRODUCT(gradient(:,neutral_equation),magnetic_direction))*scale
    neutral_pressure = 0.d0
    IF (PRESENT(neutral_pressure_vector)) THEN
       neutral_pressure_gradient = MATMUL(gradient,neutral_pressure_vector)
       neutral_pressure = DOT_PRODUCT(neutral_pressure_gradient,normal)*scale
       IF (neutral_perpendicular_diffusion) neutral_pressure = &
            neutral_pressure - magnetic_normal* &
            DOT_PRODUCT(neutral_pressure_gradient,magnetic_direction)*scale
    ENDIF
    neutral_convection = -DOT_PRODUCT(flux_jacobian(neutral_equation,:), &
         trace_state)*magnetic_normal*scale
    CALL this%add_plasma_flux(plasma_parallel,plasma_diffusion,plasma_pinch)
    CALL this%add_neutral_flux(neutral_diffusion,neutral_pressure, &
         neutral_convection)
  END SUBROUTINE accumulate_particle_face

  SUBROUTINE accumulate_particle_tau(this, measure, plasma_tau_inward, &
       &neutral_tau_inward)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: measure, plasma_tau_inward, neutral_tau_inward
    REAL*8 :: scale

    scale = measure*this%rate_scale
    this%values(particle_tau_plasma) = &
         this%values(particle_tau_plasma) + plasma_tau_inward*scale
    this%values(particle_tau_neutral) = &
         this%values(particle_tau_neutral) + neutral_tau_inward*scale
  END SUBROUTINE accumulate_particle_tau

  SUBROUTINE add_plasma_flux(this, parallel_flux, diffusion_flux, pinch_flux)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: parallel_flux, diffusion_flux, pinch_flux

    this%values(particle_boundary_plasma) = &
         this%values(particle_boundary_plasma) + &
         parallel_flux+diffusion_flux+pinch_flux
    IF (.NOT. this%detailed) RETURN

    this%values(particle_plasma_parallel_flux) = &
         this%values(particle_plasma_parallel_flux) + parallel_flux
    this%values(particle_plasma_diffusion_flux) = &
         this%values(particle_plasma_diffusion_flux) + diffusion_flux
    this%values(particle_plasma_pinch_flux) = &
         this%values(particle_plasma_pinch_flux) + pinch_flux
  END SUBROUTINE add_plasma_flux

  SUBROUTINE add_neutral_flux(this, diffusion_flux, pressure_flux, &
       &convection_flux)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: diffusion_flux, pressure_flux, convection_flux

    this%values(particle_boundary_neutral) = &
         this%values(particle_boundary_neutral) + &
         diffusion_flux+pressure_flux+convection_flux
    IF (.NOT. this%detailed) RETURN

    this%values(particle_neutral_diffusion_flux) = &
         this%values(particle_neutral_diffusion_flux) + diffusion_flux
    this%values(particle_neutral_pressure_flux) = &
         this%values(particle_neutral_pressure_flux) + pressure_flux
    this%values(particle_neutral_convection_flux) = &
         this%values(particle_neutral_convection_flux) + convection_flux
  END SUBROUTINE add_neutral_flux

  SUBROUTINE accumulate_particle_wall(this, integration_weight, fluxes)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: integration_weight
    TYPE(particle_wall_flux_type), INTENT(IN) :: fluxes
    REAL*8 :: weight

    weight = integration_weight*this%rate_scale
    this%values(wall_plasma_diffusion_flux) = &
         this%values(wall_plasma_diffusion_flux) + &
         fluxes%plasma_diffusion_inward*weight
    this%values(wall_plasma_tau_flux) = &
         this%values(wall_plasma_tau_flux) + &
         fluxes%plasma_tau_inward*weight
    this%values(wall_recycling_parallel_source) = &
         this%values(wall_recycling_parallel_source) + &
         fluxes%recycling_parallel_source*weight
    this%values(wall_recycling_diffusion_source) = &
         this%values(wall_recycling_diffusion_source) + &
         fluxes%recycling_diffusion_source*weight
    this%values(wall_recycling_pinch_source) = &
         this%values(wall_recycling_pinch_source) + &
         fluxes%recycling_pinch_source*weight
    this%values(wall_puff_source) = this%values(wall_puff_source) + &
         fluxes%puff_source*weight
    this%values(wall_pump_source) = this%values(wall_pump_source) - &
         fluxes%pump_sink*weight
    this%values(wall_neutral_diffusion_flux) = &
         this%values(wall_neutral_diffusion_flux) + &
         fluxes%neutral_diffusion_inward*weight
    this%values(wall_neutral_pressure_flux) = &
         this%values(wall_neutral_pressure_flux) + &
         fluxes%neutral_pressure_inward*weight
    this%values(wall_neutral_convection_flux) = &
         this%values(wall_neutral_convection_flux) + &
         fluxes%neutral_convection_inward*weight
    this%values(wall_neutral_tau_flux) = &
         this%values(wall_neutral_tau_flux) + fluxes%neutral_tau_inward*weight
  END SUBROUTINE accumulate_particle_wall

END MODULE balance_diagnostics
