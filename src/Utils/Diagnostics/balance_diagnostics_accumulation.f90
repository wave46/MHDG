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

  MODULE SUBROUTINE accumulate_particle_volume(this, measure, &
       &plasma_density, neutral_density, plasma_history, neutral_history, &
       &time_coefficients, time_step, steady, ionization_rate, &
       &recombination_rate, &
       &plasma_other_source, neutral_other_source, charge_exchange_rate)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: measure, plasma_density, neutral_density
    REAL*8, INTENT(IN) :: plasma_history(:), neutral_history(:)
    REAL*8, INTENT(IN) :: time_coefficients(:), time_step
    LOGICAL, INTENT(IN) :: steady
    REAL*8, INTENT(IN) :: ionization_rate, recombination_rate
    REAL*8, INTENT(IN) :: plasma_other_source, neutral_other_source
    REAL*8, INTENT(IN) :: charge_exchange_rate
    REAL*8 :: plasma_rate_coefficient, neutral_rate_coefficient

    plasma_rate_coefficient = measure*this%rate_scale(equation_n)
    neutral_rate_coefficient = measure*this%rate_scale(equation_nn)

    CALL add_particle_content(this,measure,plasma_density,neutral_density)
    CALL add_particle_temporal(this,steady,time_coefficients,time_step, &
         &plasma_density, &
         &neutral_density,plasma_history,neutral_history, &
         &plasma_rate_coefficient,neutral_rate_coefficient)
    CALL add_particle_reactions(this,ionization_rate,recombination_rate, &
         &plasma_rate_coefficient,neutral_rate_coefficient)
    CALL add_particle_sources(this,plasma_other_source,neutral_other_source, &
         &plasma_rate_coefficient,neutral_rate_coefficient)
    IF (this%detailed) CALL add_value(this,equation_n, &
         &term_charge_exchange,section_physical, &
         &charge_exchange_rate*plasma_rate_coefficient)
  END SUBROUTINE accumulate_particle_volume

  SUBROUTINE add_particle_content(this, measure, plasma_density, &
       &neutral_density)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: measure, plasma_density, neutral_density

    CALL add_value(this,equation_n,term_content,section_physical, &
         &plasma_density*measure*this%content_scale(equation_n))
    CALL add_value(this,equation_nn,term_content,section_physical, &
         &neutral_density*measure*this%content_scale(equation_nn))
  END SUBROUTINE add_particle_content

  SUBROUTINE add_particle_temporal(this, steady, time_coefficients, &
       &time_step, plasma_density, neutral_density, plasma_history, &
       &neutral_history, &
       &plasma_rate_coefficient, neutral_rate_coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: steady
    REAL*8, INTENT(IN) :: time_coefficients(:), time_step
    REAL*8, INTENT(IN) :: plasma_density, neutral_density
    REAL*8, INTENT(IN) :: plasma_history(:), neutral_history(:)
    REAL*8, INTENT(IN) :: plasma_rate_coefficient, neutral_rate_coefficient

    IF (steady) RETURN
    CALL add_value(this,equation_n,term_temporal,section_physical, &
         &discrete_time_derivative(plasma_density,plasma_history, &
         &time_coefficients,time_step)* &
         &plasma_rate_coefficient)
    CALL add_value(this,equation_nn,term_temporal,section_physical, &
         &discrete_time_derivative(neutral_density,neutral_history, &
         &time_coefficients,time_step)* &
         &neutral_rate_coefficient)
  END SUBROUTINE add_particle_temporal

  SUBROUTINE add_particle_reactions(this, ionization_rate, &
       &recombination_rate, plasma_rate_coefficient, neutral_rate_coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: ionization_rate, recombination_rate
    REAL*8, INTENT(IN) :: plasma_rate_coefficient, neutral_rate_coefficient
    REAL*8 :: net_reaction

    net_reaction = ionization_rate-recombination_rate
    CALL add_value(this,equation_n,term_volume,section_physical, &
         &net_reaction*plasma_rate_coefficient)
    CALL add_value(this,equation_nn,term_volume,section_physical, &
         &-net_reaction*neutral_rate_coefficient)
    IF (.NOT. this%detailed) RETURN
    CALL add_value(this,equation_n,term_ionization,section_physical, &
         &ionization_rate*plasma_rate_coefficient)
    CALL add_value(this,equation_n,term_recombination,section_physical, &
         &-recombination_rate*plasma_rate_coefficient)
    CALL add_value(this,equation_nn,term_ionization,section_physical, &
         &-ionization_rate*neutral_rate_coefficient)
    CALL add_value(this,equation_nn,term_recombination,section_physical, &
         &recombination_rate*neutral_rate_coefficient)
  END SUBROUTINE add_particle_reactions

  SUBROUTINE add_particle_sources(this, plasma_source, neutral_source, &
       &plasma_rate_coefficient, neutral_rate_coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: plasma_source, neutral_source
    REAL*8, INTENT(IN) :: plasma_rate_coefficient, neutral_rate_coefficient

    CALL add_value(this,equation_n,term_volume,section_physical, &
         &plasma_source*plasma_rate_coefficient)
    CALL add_value(this,equation_nn,term_volume,section_physical, &
         &neutral_source*neutral_rate_coefficient)
    IF (.NOT. this%detailed) RETURN
    CALL add_value(this,equation_n,term_prescribed_source,section_physical, &
         &plasma_source*plasma_rate_coefficient)
    CALL add_value(this,equation_nn,term_prescribed_source,section_physical, &
         &neutral_source*neutral_rate_coefficient)
  END SUBROUTINE add_particle_sources

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
