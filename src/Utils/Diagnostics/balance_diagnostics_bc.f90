SUBMODULE (balance_diagnostics) balance_diagnostics_bc
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE accumulate_bc(this, integration_weight, &
       &neutral_equation, trace_state, exterior_state, &
       &gradient, normal, magnetic_direction, magnetic_normal, tau, &
       &diffusion_iso, diffusion_ani, pinch_matrix, &
       &recycling_coefficient, puff_source, pump_coefficient, &
       &plasma, &
       &neutral_perpendicular_diffusion, neutral_pressure_vector, &
       &neutral_momentum_equation)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: integration_weight
    INTEGER, INTENT(IN) :: neutral_equation
    REAL*8, INTENT(IN) :: trace_state(:), exterior_state(:)
    REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
    REAL*8, INTENT(IN) :: magnetic_normal, tau(:,:)
    REAL*8, INTENT(IN) :: diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: pinch_matrix(:,:)
    REAL*8, INTENT(IN) :: recycling_coefficient, puff_source
    REAL*8, INTENT(IN) :: pump_coefficient
    TYPE(balance_plasma_bc_type), INTENT(IN) :: plasma
    LOGICAL, INTENT(IN) :: neutral_perpendicular_diffusion
    REAL*8, INTENT(IN), OPTIONAL :: neutral_pressure_vector(:)
    INTEGER, INTENT(IN), OPTIONAL :: neutral_momentum_equation
    REAL*8 :: density_diffusion

    density_diffusion = perpendicular_diffusive_flux(equation_n, &
         &gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)

    CALL accumulate_plasma_bc(this,trace_state,exterior_state, &
         &gradient,normal,magnetic_direction,magnetic_normal,tau, &
         &diffusion_iso,diffusion_ani,pinch_matrix,density_diffusion,plasma, &
         &integration_weight)
    CALL accumulate_neutral_bc(this,neutral_equation, &
         &trace_state,exterior_state,gradient,normal,magnetic_direction, &
         &magnetic_normal, &
         &tau,diffusion_iso,diffusion_ani,pinch_matrix, &
         &density_diffusion,recycling_coefficient,puff_source,pump_coefficient, &
         &neutral_perpendicular_diffusion,integration_weight, &
         &neutral_pressure_vector,neutral_momentum_equation)
  END SUBROUTINE accumulate_bc

  SUBROUTINE accumulate_plasma_bc(this, state, &
       &exterior_state, gradient, normal, magnetic_direction, magnetic_normal, &
       &tau, diffusion_iso, diffusion_ani, pinch_matrix, density_diffusion, &
       &plasma, integration_weight)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), gradient(:,:)
    REAL*8, INTENT(IN) :: normal(:), magnetic_direction(:), magnetic_normal
    REAL*8, INTENT(IN) :: tau(:,:), diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: pinch_matrix(:,:), density_diffusion
    TYPE(balance_plasma_bc_type), INTENT(IN) :: plasma
    REAL*8, INTENT(IN) :: integration_weight
    REAL*8 :: target(SIZE(state))

    target = exterior_state
    target(equation_nu) = plasma%momentum_target

    CALL accumulate_plasma_physical_fluxes(this,state, &
         &normal,magnetic_normal,pinch_matrix,plasma,integration_weight)
    CALL accumulate_density_closure(this,state,target,tau, &
         &density_diffusion,integration_weight*this%rate_scale(equation_n))
    CALL accumulate_momentum_closure(this,state,target,gradient,normal, &
         &magnetic_direction,magnetic_normal,tau,diffusion_iso,diffusion_ani, &
         &plasma%momentum_split_diffusive_flux, &
         &integration_weight*this%rate_scale(equation_nu))
    CALL accumulate_energy_closure(this,equation_nEi,state,target,gradient, &
         &normal, &
         &magnetic_direction,magnetic_normal,tau,diffusion_iso,diffusion_ani, &
         &plasma%ion_energy_split_diffusive_flux, &
         &plasma%ion_parallel_conductive_flux,plasma%ion_sheath_coefficient, &
         &integration_weight*this%rate_scale(equation_nEi))
    CALL accumulate_energy_closure(this,equation_nEe,state,target,gradient, &
         &normal, &
         &magnetic_direction,magnetic_normal,tau,diffusion_iso,diffusion_ani, &
         &plasma%electron_energy_split_diffusive_flux, &
         &plasma%electron_parallel_conductive_flux, &
         &plasma%electron_sheath_coefficient, &
         &integration_weight*this%rate_scale(equation_nEe))
  END SUBROUTINE accumulate_plasma_bc

  SUBROUTINE accumulate_plasma_physical_fluxes(this, &
       &state, normal, magnetic_normal, pinch_matrix, plasma, &
       &integration_weight)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: state(:), normal(:), magnetic_normal
    REAL*8, INTENT(IN) :: pinch_matrix(:,:)
    TYPE(balance_plasma_bc_type), INTENT(IN) :: plasma
    REAL*8, INTENT(IN) :: integration_weight
    REAL*8 :: density_parallel, density_pinch
    REAL*8 :: ion_sheath, electron_sheath, ion_pinch, electron_pinch
    REAL*8 :: density_scale, ion_scale, electron_scale

    density_scale = integration_weight*this%rate_scale(equation_n)
    ion_scale = integration_weight*this%rate_scale(equation_nEi)
    electron_scale = integration_weight*this%rate_scale(equation_nEe)
    density_parallel = -state(equation_nu)*magnetic_normal*density_scale
    density_pinch = -state(equation_n)* &
         &DOT_PRODUCT(pinch_matrix(equation_n,:),normal)*density_scale
    ion_pinch = -state(equation_nEi)* &
         &DOT_PRODUCT(pinch_matrix(equation_nEi,:),normal)*ion_scale
    electron_pinch = -state(equation_nEe)* &
         &DOT_PRODUCT(pinch_matrix(equation_nEe,:),normal)*electron_scale
    ion_sheath = -ion_sheath_flux(state,plasma%ion_sheath_coefficient)* &
         &magnetic_normal*ion_scale
    electron_sheath = -electron_sheath_flux(state, &
         &plasma%electron_sheath_coefficient)*magnetic_normal*electron_scale

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
  END SUBROUTINE accumulate_plasma_physical_fluxes

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

  SUBROUTINE accumulate_density_closure(this, state, exterior_state, &
       &tau, diffusion, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), tau(:,:)
    REAL*8, INTENT(IN) :: diffusion, coefficient
    REAL*8 :: stabilization

    stabilization = DOT_PRODUCT(tau(equation_n,:),state-exterior_state)
    CALL add_bc_value(this,equation_n,term_diffusion,diffusion*coefficient)
    CALL add_bc_value(this,equation_n,term_tau_stabilization_inward, &
         &stabilization*coefficient)
  END SUBROUTINE accumulate_density_closure

  SUBROUTINE accumulate_momentum_closure(this, state, target, gradient, normal, &
       &magnetic_direction, magnetic_normal, tau, diffusion_iso, &
       &diffusion_ani, split_diffusion, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: state(:), target(:), gradient(:,:), normal(:)
    REAL*8, INTENT(IN) :: magnetic_direction(:), magnetic_normal
    REAL*8, INTENT(IN) :: tau(:,:), diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: split_diffusion, coefficient
    REAL*8 :: diffusion, stabilization

    diffusion = perpendicular_diffusive_flux(equation_nu,gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)
    stabilization = DOT_PRODUCT(tau(equation_nu,:),state-target)
    CALL add_bc_value(this,equation_nu,term_diffusion,diffusion*coefficient)
    CALL add_bc_value(this,equation_nu,term_split_diffusion, &
         &split_diffusion*coefficient)
    CALL add_bc_value(this,equation_nu,term_tau_stabilization_inward, &
         &stabilization*coefficient)
  END SUBROUTINE accumulate_momentum_closure

  SUBROUTINE accumulate_energy_closure(this, equation, state, target, gradient, &
       &normal, magnetic_direction, magnetic_normal, tau, diffusion_iso, &
       &diffusion_ani, split_diffusion, parallel_conduction, &
       &sheath_coefficient, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation
    REAL*8, INTENT(IN) :: state(:), target(:), gradient(:,:), normal(:)
    REAL*8, INTENT(IN) :: magnetic_direction(:), magnetic_normal
    REAL*8, INTENT(IN) :: tau(:,:), diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: split_diffusion, parallel_conduction
    REAL*8, INTENT(IN) :: sheath_coefficient, coefficient
    REAL*8 :: diffusion, sheath_minus_bulk, stabilization

    diffusion = perpendicular_diffusive_flux(equation,gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)
    sheath_minus_bulk = sheath_minus_bulk_flux(state,equation, &
         &sheath_coefficient)*magnetic_normal
    stabilization = DOT_PRODUCT(tau(equation,:),state-target)
    CALL add_bc_value(this,equation,term_diffusion,diffusion*coefficient)
    CALL add_bc_value(this,equation,term_split_diffusion, &
         &split_diffusion*coefficient)
    CALL add_bc_value(this,equation,term_parallel_conduction, &
         &parallel_conduction*magnetic_normal*coefficient)
    CALL add_bc_value(this,equation,term_sheath_minus_bulk, &
         &sheath_minus_bulk*coefficient)
    CALL add_bc_value(this,equation,term_tau_stabilization_inward, &
         &stabilization*coefficient)
  END SUBROUTINE accumulate_energy_closure

  REAL*8 FUNCTION sheath_minus_bulk_flux(state, equation, &
       &sheath_coefficient) RESULT(flux)
    REAL*8, INTENT(IN) :: state(:), sheath_coefficient
    INTEGER, INTENT(IN) :: equation
    REAL*8 :: coefficient

    coefficient = -(5.d0-2.d0*sheath_coefficient)/3.d0
    SELECT CASE (equation)
    CASE (equation_nEi)
       flux = coefficient*state(equation_nu)/state(equation_n)* &
            &(state(equation_nEi)-0.5d0*state(equation_nu)**2/ &
            &state(equation_n))
    CASE (equation_nEe)
       flux = coefficient*state(equation_nu)*state(equation_nEe)/ &
            &state(equation_n)
    CASE DEFAULT
       flux = 0.d0
    END SELECT
  END FUNCTION sheath_minus_bulk_flux

  SUBROUTINE accumulate_neutral_bc(this, neutral_equation, &
       &state, exterior_state, gradient, normal, magnetic_direction, &
       &magnetic_normal, tau, diffusion_iso, diffusion_ani, pinch_matrix, &
       &density_diffusion, recycling_coefficient, puff_source, &
       &pump_coefficient, perpendicular_enabled, integration_weight, &
       &pressure_vector, neutral_momentum_equation)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: neutral_equation
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), gradient(:,:)
    REAL*8, INTENT(IN) :: normal(:), magnetic_direction(:), magnetic_normal
    REAL*8, INTENT(IN) :: tau(:,:), diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: pinch_matrix(:,:)
    REAL*8, INTENT(IN) :: density_diffusion
    REAL*8, INTENT(IN) :: recycling_coefficient, puff_source, pump_coefficient
    LOGICAL, INTENT(IN) :: perpendicular_enabled
    REAL*8, INTENT(IN) :: integration_weight
    REAL*8, INTENT(IN), OPTIONAL :: pressure_vector(:)
    INTEGER, INTENT(IN), OPTIONAL :: neutral_momentum_equation
    REAL*8 :: coefficient

    coefficient = integration_weight*this%rate_scale(equation_nn)
    CALL accumulate_recycling_bc(this,neutral_equation,state, &
         &exterior_state,normal,magnetic_normal,pinch_matrix,density_diffusion, &
         &recycling_coefficient,puff_source,pump_coefficient,coefficient)
    CALL accumulate_neutral_transport_bc(this,neutral_equation,state, &
         &exterior_state,gradient,normal,magnetic_direction,magnetic_normal, &
         &tau,diffusion_iso,diffusion_ani,perpendicular_enabled, &
         &coefficient,pressure_vector,neutral_momentum_equation)
  END SUBROUTINE accumulate_neutral_bc

  SUBROUTINE accumulate_recycling_bc(this, &
       &neutral_equation, state, exterior_state, normal, magnetic_normal, &
       &pinch_matrix, density_diffusion, recycling_coefficient, puff_source, &
       &pump_coefficient, coefficient)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: neutral_equation
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), normal(:)
    REAL*8, INTENT(IN) :: magnetic_normal, pinch_matrix(:,:)
    REAL*8, INTENT(IN) :: density_diffusion
    REAL*8, INTENT(IN) :: recycling_coefficient, puff_source, pump_coefficient
    REAL*8, INTENT(IN) :: coefficient
    REAL*8 :: parallel, diffusion, pinch

    parallel = recycling_coefficient*exterior_state(equation_nu)*magnetic_normal
    diffusion = -recycling_coefficient*density_diffusion
    pinch = recycling_coefficient*exterior_state(equation_n)* &
         &DOT_PRODUCT(pinch_matrix(equation_n,:),normal)
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
       &tau, diffusion_iso, diffusion_ani, &
       &perpendicular_enabled, coefficient, pressure_vector, &
       &neutral_momentum_equation)
    CLASS(balance_accumulator_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: equation
    REAL*8, INTENT(IN) :: state(:), exterior_state(:), gradient(:,:)
    REAL*8, INTENT(IN) :: normal(:), magnetic_direction(:), magnetic_normal
    REAL*8, INTENT(IN) :: tau(:,:), diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8, INTENT(IN) :: coefficient
    LOGICAL, INTENT(IN) :: perpendicular_enabled
    REAL*8, INTENT(IN), OPTIONAL :: pressure_vector(:)
    INTEGER, INTENT(IN), OPTIONAL :: neutral_momentum_equation
    REAL*8 :: diffusion, pressure, convection, stabilization

    diffusion = perpendicular_diffusive_flux(equation,gradient,normal, &
         &magnetic_direction,magnetic_normal,diffusion_iso,diffusion_ani)
    pressure = optional_neutral_pressure_flux(gradient,normal, &
         &magnetic_direction, &
         &magnetic_normal,perpendicular_enabled,pressure_vector)
    ! The default Bohm neutral flux has only the optional NeutralGamma term.
    convection = 0.d0
    IF (PRESENT(neutral_momentum_equation)) THEN
       IF (neutral_momentum_equation > 0) convection = convection- &
            &state(neutral_momentum_equation)*magnetic_normal
    ENDIF
    stabilization = DOT_PRODUCT(tau(equation,:),state-exterior_state)
    CALL add_bc_value(this,equation_nn,term_diffusion,diffusion*coefficient)
    CALL add_bc_value(this,equation_nn,term_pressure,pressure*coefficient)
    CALL add_bc_value(this,equation_nn,term_convection,convection*coefficient)
    CALL add_bc_value(this,equation_nn,term_tau_stabilization_inward, &
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
