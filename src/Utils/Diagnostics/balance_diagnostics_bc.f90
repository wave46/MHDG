SUBMODULE (balance_diagnostics) balance_diagnostics_bc
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE accumulate_bc(this, integration_weight, &
       &density_equation, neutral_equation, trace_state, exterior_state, &
       &gradient, normal, magnetic_direction, magnetic_normal, tau, &
       &diffusion_iso, diffusion_ani, pinch_matrix, flux_jacobian, &
       &recycling_coefficient, puff_source, pump_coefficient, &
       &ion_sheath_coefficient, electron_sheath_coefficient, &
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
    REAL*8, INTENT(IN) :: ion_sheath_coefficient
    REAL*8, INTENT(IN) :: electron_sheath_coefficient
    LOGICAL, INTENT(IN) :: neutral_perpendicular_diffusion
    REAL*8, INTENT(IN), OPTIONAL :: neutral_pressure_vector(:)
    INTEGER, INTENT(IN), OPTIONAL :: neutral_momentum_equation
    REAL*8 :: density_coefficient, neutral_coefficient
    REAL*8 :: density_diffusion

    density_coefficient = integration_weight*this%rate_scale(equation_n)
    neutral_coefficient = integration_weight*this%rate_scale(equation_nn)
    density_diffusion = perpendicular_diffusive_flux(density_equation, &
         &gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)

    CALL accumulate_density_bc(this,density_equation,trace_state, &
         &exterior_state,tau,density_diffusion,density_coefficient)
    CALL accumulate_plasma_physical_targets(this,density_equation, &
         &trace_state,normal,magnetic_normal,pinch_matrix, &
         &ion_sheath_coefficient,electron_sheath_coefficient, &
         &integration_weight)
    CALL accumulate_recycling_bc(this,density_equation,neutral_equation, &
         &trace_state,exterior_state,normal,magnetic_normal,pinch_matrix, &
         &density_diffusion,recycling_coefficient,puff_source,pump_coefficient, &
         &neutral_coefficient)
    CALL accumulate_neutral_transport_bc(this,neutral_equation,trace_state, &
         &exterior_state,gradient,normal,magnetic_direction,magnetic_normal, &
         &tau,diffusion_iso,diffusion_ani,flux_jacobian, &
         &neutral_perpendicular_diffusion,neutral_coefficient, &
         &neutral_pressure_vector,neutral_momentum_equation)
  END SUBROUTINE accumulate_bc

  SUBROUTINE accumulate_plasma_physical_targets(this, density_equation, &
       &state, normal, magnetic_normal, pinch_matrix, ion_sheath_coefficient, &
       &electron_sheath_coefficient, integration_weight)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: density_equation
    REAL*8, INTENT(IN) :: state(:), normal(:), magnetic_normal
    REAL*8, INTENT(IN) :: pinch_matrix(:,:)
    REAL*8, INTENT(IN) :: ion_sheath_coefficient
    REAL*8, INTENT(IN) :: electron_sheath_coefficient, integration_weight
    REAL*8 :: density_parallel, density_pinch
    REAL*8 :: ion_sheath, electron_sheath, ion_pinch, electron_pinch
    REAL*8 :: density_scale, ion_scale, electron_scale

    density_scale = integration_weight*this%rate_scale(equation_n)
    ion_scale = integration_weight*this%rate_scale(equation_nEi)
    electron_scale = integration_weight*this%rate_scale(equation_nEe)
    density_parallel = -state(equation_nu)*magnetic_normal*density_scale
    density_pinch = -state(density_equation)* &
         &DOT_PRODUCT(pinch_matrix(density_equation,:),normal)*density_scale
    ion_pinch = -state(equation_nEi)* &
         &DOT_PRODUCT(pinch_matrix(equation_nEi,:),normal)*ion_scale
    electron_pinch = -state(equation_nEe)* &
         &DOT_PRODUCT(pinch_matrix(equation_nEe,:),normal)*electron_scale
    ion_sheath = -ion_sheath_flux(state,ion_sheath_coefficient)* &
         &magnetic_normal*ion_scale
    electron_sheath = -electron_sheath_flux(state, &
         &electron_sheath_coefficient)*magnetic_normal*electron_scale

    CALL add_bc_value(this,equation_n,term_boundary_physical_inward, &
         &density_parallel+density_pinch)
    CALL add_bc_value(this,equation_n,term_parallel,density_parallel)
    CALL add_bc_value(this,equation_n,term_pinch,density_pinch)
    CALL add_bc_value(this,equation_nEi,term_boundary_physical_inward, &
         &ion_sheath+ion_pinch)
    CALL add_bc_value(this,equation_nEi,term_sheath,ion_sheath)
    CALL add_bc_value(this,equation_nEi,term_pinch,ion_pinch)
    CALL add_bc_value(this,equation_nEe,term_boundary_physical_inward, &
         &electron_sheath+electron_pinch)
    CALL add_bc_value(this,equation_nEe,term_sheath,electron_sheath)
    CALL add_bc_value(this,equation_nEe,term_pinch,electron_pinch)
  END SUBROUTINE accumulate_plasma_physical_targets

  REAL*8 FUNCTION ion_sheath_flux(state, sheath_coefficient) RESULT(flux)
    REAL*8, INTENT(IN) :: state(:), sheath_coefficient
    REAL*8 :: velocity, internal_energy

    velocity = state(equation_nu)/state(equation_n)
    internal_energy = state(equation_nEi)/state(equation_n)- &
         &0.5d0*velocity**2
    flux = 2.d0/3.d0*sheath_coefficient*state(equation_nu)* &
         &internal_energy+0.5d0*state(equation_n)*velocity**3
  END FUNCTION ion_sheath_flux

  REAL*8 FUNCTION electron_sheath_flux(state, sheath_coefficient) RESULT(flux)
    REAL*8, INTENT(IN) :: state(:), sheath_coefficient

    flux = 2.d0/3.d0*sheath_coefficient*state(equation_nu)* &
         &state(equation_nEe)/state(equation_n)
  END FUNCTION electron_sheath_flux

  SUBROUTINE accumulate_density_bc(this, equation, state, exterior_state, &
       &tau, diffusion, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), tau(:,:)
    REAL*8, INTENT(IN) :: diffusion, coefficient
    REAL*8 :: stabilization

    stabilization = DOT_PRODUCT(tau(equation,:),state-exterior_state)
    CALL add_bc_value(this,equation_n,term_diffusion,diffusion*coefficient)
    CALL add_bc_value(this,equation_n,term_tau_inward, &
         &stabilization*coefficient)
  END SUBROUTINE accumulate_density_bc

  SUBROUTINE accumulate_recycling_bc(this, density_equation, &
       &neutral_equation, state, exterior_state, normal, magnetic_normal, &
       &pinch_matrix, density_diffusion, recycling_coefficient, puff_source, &
       &pump_coefficient, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: density_equation, neutral_equation
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), normal(:)
    REAL*8, INTENT(IN) :: magnetic_normal, pinch_matrix(:,:)
    REAL*8, INTENT(IN) :: density_diffusion
    REAL*8, INTENT(IN) :: recycling_coefficient, puff_source, pump_coefficient
    REAL*8, INTENT(IN) :: coefficient
    REAL*8 :: parallel, diffusion, pinch

    parallel = recycling_coefficient*exterior_state(2)*magnetic_normal
    diffusion = -recycling_coefficient*density_diffusion
    pinch = recycling_coefficient*exterior_state(density_equation)* &
         &DOT_PRODUCT(pinch_matrix(density_equation,:),normal)
    CALL add_bc_value(this,equation_nn,term_recycling_parallel, &
         &parallel*coefficient)
    CALL add_bc_value(this,equation_nn,term_recycling_diffusion, &
         &diffusion*coefficient)
    CALL add_bc_value(this,equation_nn,term_recycling_pinch,pinch*coefficient)
    CALL add_bc_value(this,equation_nn,term_puff,puff_source*coefficient)
    CALL add_bc_value(this,equation_nn,term_pump, &
         &pump_coefficient*state(neutral_equation)*coefficient)
    CALL add_bc_value(this,equation_nn,term_boundary_physical_inward, &
         &(parallel+diffusion+pinch+puff_source- &
         &pump_coefficient*state(neutral_equation))*coefficient)
  END SUBROUTINE accumulate_recycling_bc

  SUBROUTINE accumulate_neutral_transport_bc(this, equation, state, &
       &exterior_state, gradient, normal, magnetic_direction, magnetic_normal, &
       &tau, diffusion_iso, diffusion_ani, flux_jacobian, &
       &perpendicular_enabled, coefficient, pressure_vector, &
       &neutral_momentum_equation)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), gradient(:,:)
    REAL*8, INTENT(IN) :: normal(:), magnetic_direction(:), magnetic_normal
    REAL*8, INTENT(IN) :: tau(:,:), diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: flux_jacobian(:,:), coefficient
    LOGICAL, INTENT(IN) :: perpendicular_enabled
    REAL*8, INTENT(IN), OPTIONAL :: pressure_vector(:)
    INTEGER, INTENT(IN), OPTIONAL :: neutral_momentum_equation
    REAL*8 :: diffusion, pressure, convection, stabilization

    diffusion = perpendicular_diffusive_flux(equation,gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)
    pressure = optional_neutral_pressure_flux(gradient,normal, &
         &magnetic_direction, &
         &magnetic_normal,perpendicular_enabled,pressure_vector)
    convection = -DOT_PRODUCT(flux_jacobian(equation,:),state)*magnetic_normal
    IF (PRESENT(neutral_momentum_equation)) THEN
       IF (neutral_momentum_equation > 0) convection = convection- &
            &state(neutral_momentum_equation)*magnetic_normal
    ENDIF
    stabilization = DOT_PRODUCT(tau(equation,:),state-exterior_state)
    CALL add_bc_value(this,equation_nn,term_diffusion,diffusion*coefficient)
    CALL add_bc_value(this,equation_nn,term_pressure,pressure*coefficient)
    CALL add_bc_value(this,equation_nn,term_convection,convection*coefficient)
    CALL add_bc_value(this,equation_nn,term_tau_inward, &
         &stabilization*coefficient)
  END SUBROUTINE accumulate_neutral_transport_bc

  SUBROUTINE add_bc_value(this, equation, term, value)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation, term
    REAL*8, INTENT(IN) :: value

    this%values(equation,term,section_bc) = &
         this%values(equation,term,section_bc)+value
  END SUBROUTINE add_bc_value

END SUBMODULE balance_diagnostics_bc
