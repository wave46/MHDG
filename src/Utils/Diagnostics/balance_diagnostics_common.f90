SUBMODULE (balance_diagnostics) balance_diagnostics_common
  IMPLICIT NONE

CONTAINS

  MODULE FUNCTION balance_value(this, equation, term, section) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation, term, section
    REAL*8 :: value

    value = this%finalized_values(equation,term,section)
  END FUNCTION balance_value

  MODULE FUNCTION physical_balance(this, equation) RESULT(balance)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation
    TYPE(physical_balance_type) :: balance

    balance%content = balance_value(this,equation,term_content,section_physical)
    balance%temporal = &
         balance_value(this,equation,term_temporal,section_physical)
    balance%volume = balance_value(this,equation,term_volume,section_physical)
    balance%boundary_inward = balance_value(this,equation, &
         &term_boundary_physical_inward,section_physical)
    balance%tau_inward = &
         balance_value(this,equation,term_tau_inward,section_physical)
    balance%physical_imbalance = balance%temporal-balance%volume- &
         balance%boundary_inward
    balance%discrete_residual = &
         balance%physical_imbalance-balance%tau_inward
  END FUNCTION physical_balance

  MODULE FUNCTION add_physical_balances(left, right) RESULT(total)
    TYPE(physical_balance_type), INTENT(IN) :: left, right
    TYPE(physical_balance_type) :: total

    total%content = left%content+right%content
    total%temporal = left%temporal+right%temporal
    total%volume = left%volume+right%volume
    total%boundary_inward = left%boundary_inward+right%boundary_inward
    total%physical_imbalance = &
         left%physical_imbalance+right%physical_imbalance
    total%tau_inward = left%tau_inward+right%tau_inward
    total%discrete_residual = &
         left%discrete_residual+right%discrete_residual
  END FUNCTION add_physical_balances

  MODULE FUNCTION neutral_density_bc(this) RESULT(balance)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(particle_bc_balance_type) :: balance

    balance%imposed_source_inward = &
         balance_value(this,equation_nn,term_recycling_parallel,section_bc)+ &
         balance_value(this,equation_nn,term_recycling_diffusion,section_bc)+ &
         balance_value(this,equation_nn,term_recycling_pinch,section_bc)+ &
         balance_value(this,equation_nn,term_puff,section_bc)- &
         balance_value(this,equation_nn,term_pump,section_bc)
    balance%physical_flux_inward = &
         balance_value(this,equation_nn,term_diffusion,section_bc)+ &
         balance_value(this,equation_nn,term_pressure,section_bc)+ &
         balance_value(this,equation_nn,term_convection,section_bc)
    balance%tau_inward = &
         balance_value(this,equation_nn,term_tau_inward,section_bc)
    balance%residual = balance%imposed_source_inward- &
         balance%physical_flux_inward-balance%tau_inward
  END FUNCTION neutral_density_bc

  MODULE FUNCTION perpendicular_diffusive_flux(equation, gradient, normal, &
       &magnetic_direction, magnetic_normal, diffusion_iso, diffusion_ani) &
       &RESULT(flux)
    INTEGER, INTENT(IN) :: equation
    REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
    REAL*8, INTENT(IN) :: magnetic_normal, diffusion_iso(:,:), diffusion_ani(:,:)
    REAL*8 :: flux

    flux = diffusion_iso(equation,equation)* &
         &DOT_PRODUCT(gradient(:,equation),normal)- &
         &diffusion_ani(equation,equation)*magnetic_normal* &
         &DOT_PRODUCT(gradient(:,equation),magnetic_direction)
  END FUNCTION perpendicular_diffusive_flux

  MODULE FUNCTION optional_neutral_pressure_flux(gradient, normal, &
       &magnetic_direction, magnetic_normal, perpendicular_enabled, &
       &pressure_vector) RESULT(flux)
    REAL*8, INTENT(IN) :: gradient(:,:), normal(:), magnetic_direction(:)
    REAL*8, INTENT(IN) :: magnetic_normal
    LOGICAL, INTENT(IN) :: perpendicular_enabled
    REAL*8, INTENT(IN), OPTIONAL :: pressure_vector(:)
    REAL*8 :: flux
    REAL*8 :: pressure_gradient(SIZE(normal))

    flux = 0.d0
    IF (.NOT. PRESENT(pressure_vector)) RETURN
    pressure_gradient = MATMUL(gradient,pressure_vector)
    flux = DOT_PRODUCT(pressure_gradient,normal)
    IF (perpendicular_enabled) flux = flux-magnetic_normal* &
         &DOT_PRODUCT(pressure_gradient,magnetic_direction)
  END FUNCTION optional_neutral_pressure_flux

END SUBMODULE balance_diagnostics_common
