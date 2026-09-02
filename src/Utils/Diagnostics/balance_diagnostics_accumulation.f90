SUBMODULE (balance_diagnostics) balance_diagnostics_accumulation
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE accumulator_reset(this, detailed, content_scale, &
       &rate_scale)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: detailed
    REAL*8, INTENT(IN) :: content_scale(balance_equation_count)
    REAL*8, INTENT(IN) :: rate_scale(balance_equation_count)

    this%values = 0.d0
    this%detailed = detailed
    this%content_scale = content_scale
    this%rate_scale = rate_scale
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
    LOGICAL :: active(balance_equation_count)
    INTEGER :: state_index(balance_equation_count)
    REAL*8 :: rate_coefficient(balance_equation_count)

    rate_coefficient = measure*this%rate_scale
    CALL select_active_equations(active,state_index,neutral_state_index, &
         &PRESENT(momentum),PRESENT(energy))
    CALL add_conserved_content(this,measure,state,active,state_index)
    CALL add_conserved_temporal(this,steady,state,history, &
         &time_coefficients,time_step,rate_coefficient,active,state_index)
    CALL add_particle_volume_terms(this,atomic,prescribed_source, &
         &neutral_state_index,rate_coefficient)
    IF (PRESENT(momentum)) CALL add_momentum_volume_terms(this,atomic, &
         &momentum,prescribed_source(equation_nu), &
         &rate_coefficient(equation_nu))
    IF (PRESENT(energy)) CALL add_energy_volume_terms(this,atomic,energy, &
         &prescribed_source,rate_coefficient)
  END SUBROUTINE accumulate_volume

  SUBROUTINE select_active_equations(active, state_index, &
       &neutral_state_index, include_momentum, include_energy)
    LOGICAL, INTENT(OUT) :: active(balance_equation_count)
    INTEGER, INTENT(OUT) :: state_index(balance_equation_count)
    INTEGER, INTENT(IN) :: neutral_state_index
    LOGICAL, INTENT(IN) :: include_momentum, include_energy

    active = .FALSE.
    active(equation_n) = .TRUE.
    active(equation_nu) = include_momentum
    active(equation_nEi:equation_nEe) = include_energy
    active(equation_nn) = .TRUE.
    state_index = (/equation_n,equation_nu,equation_nEi, &
         &equation_nEe,neutral_state_index/)
  END SUBROUTINE select_active_equations

  SUBROUTINE add_conserved_content(this, measure, state, active, state_index)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: measure, state(:)
    LOGICAL, INTENT(IN) :: active(balance_equation_count)
    INTEGER, INTENT(IN) :: state_index(balance_equation_count)
    INTEGER :: equation

    DO equation = equation_n,equation_nn
       IF (.NOT. active(equation)) CYCLE
       CALL add_value(this,equation,term_content,section_physical, &
            &state(state_index(equation))*measure* &
            &this%content_scale(equation))
    ENDDO
  END SUBROUTINE add_conserved_content

  SUBROUTINE add_conserved_temporal(this, steady, state, history, &
       &time_coefficients, time_step, rate_coefficient, active, state_index)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: steady
    REAL*8, INTENT(IN) :: state(:), history(:,:), time_coefficients(:)
    REAL*8, INTENT(IN) :: time_step
    REAL*8, INTENT(IN) :: rate_coefficient(balance_equation_count)
    LOGICAL, INTENT(IN) :: active(balance_equation_count)
    INTEGER, INTENT(IN) :: state_index(balance_equation_count)
    INTEGER :: equation

    IF (steady) RETURN
    DO equation = equation_n,equation_nn
       IF (.NOT. active(equation)) CYCLE
       CALL add_value(this,equation,term_temporal,section_physical, &
            &discrete_time_derivative(state(state_index(equation)), &
            &history(state_index(equation),:),time_coefficients,time_step)* &
            &rate_coefficient(equation))
    ENDDO
  END SUBROUTINE add_conserved_temporal

  SUBROUTINE add_particle_volume_terms(this, atomic, prescribed_source, &
       &neutral_state_index, rate_coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    TYPE(balance_atomic_volume_type), INTENT(IN) :: atomic
    REAL*8, INTENT(IN) :: prescribed_source(:)
    INTEGER, INTENT(IN) :: neutral_state_index
    REAL*8, INTENT(IN) :: rate_coefficient(balance_equation_count)
    REAL*8 :: ionization, recombination

    ionization = atomic%ionization_density* &
         &atomic%ionization_rate_coefficient
    recombination = atomic%recombination_density* &
         &atomic%recombination_rate_coefficient
    CALL add_volume_component(this,equation_n,term_ionization,ionization, &
         &rate_coefficient(equation_n))
    CALL add_volume_component(this,equation_n,term_recombination, &
         &-recombination,rate_coefficient(equation_n))
    CALL add_volume_component(this,equation_n,term_prescribed_source, &
         &prescribed_source(equation_n),rate_coefficient(equation_n))
    CALL add_volume_component(this,equation_nn,term_ionization,-ionization, &
         &rate_coefficient(equation_nn))
    CALL add_volume_component(this,equation_nn,term_recombination, &
         &recombination,rate_coefficient(equation_nn))
    CALL add_volume_component(this,equation_nn,term_prescribed_source, &
         &prescribed_source(neutral_state_index),rate_coefficient(equation_nn))
    IF (this%detailed) CALL add_value(this,equation_n,term_charge_exchange, &
         &section_physical,atomic%ionization_density* &
         &atomic%charge_exchange_rate_coefficient*rate_coefficient(equation_n))
  END SUBROUTINE add_particle_volume_terms

  SUBROUTINE add_momentum_volume_terms(this, atomic, momentum, &
       &prescribed_source, rate_coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    TYPE(balance_atomic_volume_type), INTENT(IN) :: atomic
    TYPE(balance_momentum_volume_type), INTENT(IN) :: momentum
    REAL*8, INTENT(IN) :: prescribed_source, rate_coefficient
    REAL*8 :: ionization, charge_exchange

    ionization = momentum%neutral_factor* &
         &atomic%ionization_rate_coefficient
    charge_exchange = (momentum%neutral_factor- &
         &momentum%charge_exchange_factor)* &
         &atomic%charge_exchange_rate_coefficient
    CALL add_volume_component(this,equation_nu,term_ionization,ionization, &
         &rate_coefficient)
    CALL add_volume_component(this,equation_nu,term_recombination, &
         &-momentum%recombination_factor* &
         &atomic%recombination_rate_coefficient, &
         &rate_coefficient)
    CALL add_volume_component(this,equation_nu,term_charge_exchange, &
         &charge_exchange,rate_coefficient)
    CALL add_volume_component(this,equation_nu,term_pressure_divergence, &
         &momentum%pressure_divergence,rate_coefficient)
    CALL add_volume_component(this,equation_nu,term_prescribed_source, &
         &prescribed_source,rate_coefficient)
  END SUBROUTINE add_momentum_volume_terms

  SUBROUTINE add_energy_volume_terms(this, atomic, energy, prescribed_source, &
       &rate_coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    TYPE(balance_atomic_volume_type), INTENT(IN) :: atomic
    TYPE(balance_energy_volume_type), INTENT(IN) :: energy
    REAL*8, INTENT(IN) :: prescribed_source(:)
    REAL*8, INTENT(IN) :: rate_coefficient(balance_equation_count)
    REAL*8 :: ion_ionization, ion_charge_exchange

    ion_ionization = (energy%ion_ionization_factor+ &
         &energy%neutral_factor)*atomic%ionization_rate_coefficient
    ion_charge_exchange = (energy%neutral_factor- &
         &energy%ion_charge_exchange_factor)* &
         &atomic%charge_exchange_rate_coefficient
    CALL add_volume_component(this,equation_nEi,term_ionization, &
         &ion_ionization,rate_coefficient(equation_nEi))
    CALL add_volume_component(this,equation_nEi,term_recombination, &
         &-energy%ion_recombination_factor* &
         &atomic%recombination_rate_coefficient, &
         &rate_coefficient(equation_nEi))
    CALL add_volume_component(this,equation_nEi,term_charge_exchange, &
         &ion_charge_exchange,rate_coefficient(equation_nEi))
    CALL add_volume_component(this,equation_nEi,term_parallel_electric_work, &
         &-energy%parallel_electric_transfer, &
         &rate_coefficient(equation_nEi))
    CALL add_volume_component(this,equation_nEi,term_temperature_exchange, &
         &-energy%temperature_transfer,rate_coefficient(equation_nEi))
    CALL add_volume_component(this,equation_nEi,term_prescribed_source, &
         &prescribed_source(equation_nEi),rate_coefficient(equation_nEi))

    CALL add_volume_component(this,equation_nEe,term_ionization, &
         &-atomic%ionization_density* &
         &energy%electron_ionization_loss_coefficient, &
         &rate_coefficient(equation_nEe))
    CALL add_volume_component(this,equation_nEe,term_recombination, &
         &atomic%recombination_density*(energy%recombination_energy* &
         &atomic%recombination_rate_coefficient- &
         &energy%electron_recombination_loss_coefficient), &
         &rate_coefficient(equation_nEe))
    CALL add_volume_component(this,equation_nEe,term_radiation, &
         &-atomic%recombination_density*energy%impurity_cooling_factor, &
         &rate_coefficient(equation_nEe))
    CALL add_volume_component(this,equation_nEe,term_ohmic, &
         &energy%ohmic_heating, &
         &rate_coefficient(equation_nEe))
    CALL add_volume_component(this,equation_nEe,term_parallel_electric_work, &
         &energy%parallel_electric_transfer,rate_coefficient(equation_nEe))
    CALL add_volume_component(this,equation_nEe,term_temperature_exchange, &
         &energy%temperature_transfer,rate_coefficient(equation_nEe))
    CALL add_volume_component(this,equation_nEe,term_prescribed_source, &
         &prescribed_source(equation_nEe),rate_coefficient(equation_nEe))
  END SUBROUTINE add_energy_volume_terms

  SUBROUTINE add_volume_component(this, equation, term, source, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation, term
    REAL*8, INTENT(IN) :: source, coefficient

    CALL add_value(this,equation,term_volume,section_physical, &
         &source*coefficient)
    IF (this%detailed) CALL add_value(this,equation,term,section_physical, &
         &source*coefficient)
  END SUBROUTINE add_volume_component

  MODULE SUBROUTINE accumulate_relocated_sources(this, puff_source, &
       &pump_sink)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: puff_source, pump_sink
    REAL*8 :: rate_scale

    rate_scale = this%rate_scale(equation_nn)
    CALL add_value(this,equation_nn,term_volume,section_physical, &
         &(puff_source-pump_sink)*rate_scale)
    IF (.NOT. this%detailed) RETURN
    CALL add_value(this,equation_nn,term_puff,section_physical, &
         &puff_source*rate_scale)
    CALL add_value(this,equation_nn,term_pump,section_physical, &
         &-pump_sink*rate_scale)
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
    REAL*8 :: plasma_coefficient, neutral_coefficient, magnetic_normal

    plasma_coefficient = measure*this%rate_scale(equation_n)
    neutral_coefficient = measure*this%rate_scale(equation_nn)
    magnetic_normal = DOT_PRODUCT(magnetic_direction,normal)
    CALL add_plasma_particle_flux(this,plasma_equation,trace_state, &
         &flux_jacobian,pinch_matrix,gradient,normal,magnetic_direction, &
         &magnetic_normal,diffusion_iso,diffusion_ani,plasma_coefficient)
    CALL add_neutral_particle_flux(this,neutral_equation,trace_state, &
         &flux_jacobian,gradient,normal,magnetic_direction,magnetic_normal, &
         &diffusion_iso,diffusion_ani,neutral_perpendicular_diffusion, &
         &neutral_coefficient,neutral_pressure_vector)
  END SUBROUTINE accumulate_particle_face

  SUBROUTINE add_plasma_particle_flux(this, equation, state, flux_jacobian, &
       &pinch_matrix, gradient, normal, magnetic_direction, magnetic_normal, &
       &diffusion_iso, diffusion_ani, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation
    REAL*8, INTENT(IN) :: state(:), flux_jacobian(:,:), pinch_matrix(:,:)
    REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
    REAL*8, INTENT(IN) :: magnetic_normal, diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: coefficient
    REAL*8 :: parallel, diffusion, pinch

    parallel = -DOT_PRODUCT(flux_jacobian(equation,:),state)* &
         &magnetic_normal*coefficient
    diffusion = perpendicular_diffusive_flux(equation,gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)* &
         &coefficient
    pinch = -state(equation)*DOT_PRODUCT(pinch_matrix(equation,:),normal)* &
         &coefficient
    CALL add_value(this,equation_n,term_boundary_physical_inward, &
         &section_physical,parallel+diffusion+pinch)
    IF (.NOT. this%detailed) RETURN
    CALL add_value(this,equation_n,term_parallel,section_physical,parallel)
    CALL add_value(this,equation_n,term_diffusion,section_physical,diffusion)
    CALL add_value(this,equation_n,term_pinch,section_physical,pinch)
  END SUBROUTINE add_plasma_particle_flux

  SUBROUTINE add_neutral_particle_flux(this, equation, state, flux_jacobian, &
       &gradient, normal, magnetic_direction, magnetic_normal, diffusion_iso, &
       &diffusion_ani, perpendicular_enabled, coefficient, pressure_vector)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation
    REAL*8, INTENT(IN) :: state(:), flux_jacobian(:,:), gradient(:,:)
    REAL*8, INTENT(IN) :: normal(:), magnetic_direction(:), magnetic_normal
    REAL*8, INTENT(IN) :: diffusion_iso(:,:), diffusion_ani(:,:), coefficient
    LOGICAL, INTENT(IN) :: perpendicular_enabled
    REAL*8, INTENT(IN), OPTIONAL :: pressure_vector(:)
    REAL*8 :: diffusion, pressure, convection

    diffusion = perpendicular_diffusive_flux(equation,gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)* &
         &coefficient
    pressure = optional_neutral_pressure_flux(gradient,normal, &
         &magnetic_direction,magnetic_normal,perpendicular_enabled, &
         &pressure_vector)*coefficient
    convection = -DOT_PRODUCT(flux_jacobian(equation,:),state)* &
         &magnetic_normal*coefficient
    CALL add_value(this,equation_nn,term_boundary_physical_inward, &
         &section_physical,diffusion+pressure+convection)
    IF (.NOT. this%detailed) RETURN
    CALL add_value(this,equation_nn,term_diffusion,section_physical,diffusion)
    CALL add_value(this,equation_nn,term_pressure,section_physical,pressure)
    CALL add_value(this,equation_nn,term_convection,section_physical,convection)
  END SUBROUTINE add_neutral_particle_flux

  MODULE SUBROUTINE accumulate_particle_tau(this, measure, &
       &plasma_tau_inward, neutral_tau_inward)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: measure, plasma_tau_inward, neutral_tau_inward

    CALL add_value(this,equation_n,term_tau_inward,section_physical, &
         &plasma_tau_inward*measure*this%rate_scale(equation_n))
    CALL add_value(this,equation_nn,term_tau_inward,section_physical, &
         &neutral_tau_inward*measure*this%rate_scale(equation_nn))
  END SUBROUTINE accumulate_particle_tau

  SUBROUTINE add_value(this, equation, term, section, value)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation, term, section
    REAL*8, INTENT(IN) :: value

    this%values(equation,term,section) = &
         this%values(equation,term,section)+value
  END SUBROUTINE add_value

  REAL*8 FUNCTION discrete_time_derivative(current, history, coefficients, &
       &time_step) &
       &RESULT(derivative)
    REAL*8, INTENT(IN) :: current, history(:), coefficients(:), time_step
    INTEGER :: history_index

    derivative = coefficients(1)*current
    DO history_index = 1, SIZE(history)
       derivative = derivative-coefficients(history_index+1)* &
            &history(history_index)
    ENDDO
    derivative = derivative/time_step
  END FUNCTION discrete_time_derivative

END SUBMODULE balance_diagnostics_accumulation
