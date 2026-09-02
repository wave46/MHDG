SUBMODULE (balance_diagnostics) balance_diagnostics_reporting
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE balance_diagnostics_report(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    IF (.NOT. this%has_values()) RETURN
    IF (MPIvar%glob_id /= 0) RETURN

    IF (.NOT. this%detailed()) THEN
       CALL print_summary(this)
    ELSE
       CALL print_detailed(this)
    ENDIF
  END SUBROUTINE balance_diagnostics_report

  SUBROUTINE print_summary(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type) :: density, neutral, total
    TYPE(particle_bc_balance_type) :: neutral_bc
    REAL*8 :: density_bc_diffusion, density_bc_tau

    density = physical_balance(this,equation_n)
    neutral = physical_balance(this,equation_nn)
    total = add_physical_balances(density,neutral)
    neutral_bc = neutral_density_bc(this)
    density_bc_diffusion = &
         balance_value(this,equation_n,term_diffusion,section_bc)
    density_bc_tau = &
         balance_value(this,equation_n,term_tau_inward,section_bc)

    WRITE(6,'(A)') 'Balance diagnostics (summary)'
    WRITE(6,'(A)') '  Physical [content, imbalance]'
    CALL print_pair('n',density%content,density%physical_imbalance)
    CALL print_pair('n_n',neutral%content,neutral%physical_imbalance)
    CALL print_pair('total_n',total%content,total%physical_imbalance)
    WRITE(6,'(A)') '  Discrete [HDG tau inward, residual]'
    CALL print_pair('n',density%tau_inward,density%discrete_residual)
    CALL print_pair('n_n',neutral%tau_inward,neutral%discrete_residual)
    CALL print_pair('total_n',total%tau_inward,total%discrete_residual)
    WRITE(6,'(A)') '  Boundary-condition residuals [particles/s]'
    CALL print_value('n',density_bc_diffusion+density_bc_tau)
    CALL print_value('n_n',neutral_bc%residual)
  END SUBROUTINE print_summary

  SUBROUTINE print_detailed(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type) :: density, neutral, total
    TYPE(particle_bc_balance_type) :: neutral_bc
    REAL*8 :: density_bc_diffusion, density_bc_tau

    density = physical_balance(this,equation_n)
    neutral = physical_balance(this,equation_nn)
    total = add_physical_balances(density,neutral)
    neutral_bc = neutral_density_bc(this)
    density_bc_diffusion = &
         balance_value(this,equation_n,term_diffusion,section_bc)
    density_bc_tau = &
         balance_value(this,equation_n,term_tau_inward,section_bc)

    WRITE(6,'(A)') 'Balance diagnostics (detailed)'
    CALL print_physical_detail(this,'n',equation_n,density,'particles', &
         &'particles/s')
    CALL print_physical_detail(this,'n_n',equation_nn,neutral,'particles', &
         &'particles/s')
    CALL print_derived_physical_detail('total_n',total,'particles', &
         &'particles/s')
    WRITE(6,'(A)') '  Independent boundary-condition checks [particles/s]'
    WRITE(6,'(A)') '    n: diffusion + HDG tau = 0'
    CALL print_component('diffusion inward',density_bc_diffusion)
    CALL print_component('HDG tau inward',density_bc_tau)
    CALL print_component('residual',density_bc_diffusion+density_bc_tau)
    CALL print_neutral_bc_detail(this,neutral_bc)
  END SUBROUTINE print_detailed

  SUBROUTINE print_neutral_bc_detail(this, balance)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(particle_bc_balance_type), INTENT(IN) :: balance

    WRITE(6,'(A)') '    n_n: imposed source - physical flux - HDG tau = 0'
    WRITE(6,'(A)') '      Imposed-source components: recycling + puff - pump'
    CALL print_subcomponent('recycling parallel inward', &
         balance_value(this,equation_nn,term_recycling_parallel,section_bc))
    CALL print_subcomponent('recycling diffusion inward', &
         balance_value(this,equation_nn,term_recycling_diffusion,section_bc))
    CALL print_subcomponent('recycling pinch inward', &
         balance_value(this,equation_nn,term_recycling_pinch,section_bc))
    CALL print_subcomponent('puff source', &
         balance_value(this,equation_nn,term_puff,section_bc))
    CALL print_subcomponent('pump sink (subtracted)', &
         balance_value(this,equation_nn,term_pump,section_bc))
    CALL print_component('imposed source inward',balance%imposed_source_inward)
    WRITE(6,'(A)') '      Physical-flux components inward'
    CALL print_subcomponent('limited diffusion', &
         balance_value(this,equation_nn,term_diffusion,section_bc))
#ifdef NEUTRALP
    CALL print_subcomponent('limited pressure', &
         balance_value(this,equation_nn,term_pressure,section_bc))
#endif
#ifdef NEUTRALGAMMA
    CALL print_subcomponent('neutral-gamma convection', &
         balance_value(this,equation_nn,term_convection,section_bc))
#endif
    CALL print_component('physical flux inward',balance%physical_flux_inward)
    CALL print_component('HDG tau inward',balance%tau_inward)
    CALL print_component('residual',balance%residual)
  END SUBROUTINE print_neutral_bc_detail

  SUBROUTINE print_physical_detail(this, name, equation, balance, &
       &content_units, rate_units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name, content_units, rate_units
    INTEGER, INTENT(IN) :: equation
    TYPE(physical_balance_type), INTENT(IN) :: balance

    WRITE(6,'(A)') '  Physical '//TRIM(name)
    CALL print_component('content ['//TRIM(content_units)//']',balance%content)
    CALL print_component('temporal ['//TRIM(rate_units)//']',balance%temporal)
    CALL print_component('volume',balance%volume)
    CALL print_component('boundary physical inward',balance%boundary_inward)
    CALL print_component('physical imbalance',balance%physical_imbalance)
    WRITE(6,'(A)') '    Discrete'
    CALL print_component('HDG tau inward',balance%tau_inward)
    CALL print_component('residual',balance%discrete_residual)
    WRITE(6,'(A)') '    Volume components'
    CALL print_component('ionization', &
         balance_value(this,equation,term_ionization,section_physical))
    CALL print_component('recombination', &
         balance_value(this,equation,term_recombination,section_physical))
    CALL print_component('prescribed source', &
         balance_value(this,equation,term_prescribed_source,section_physical))
    IF (equation == equation_nn) THEN
       CALL print_component('puff source', &
            balance_value(this,equation,term_puff,section_physical))
       CALL print_component('pump source', &
            balance_value(this,equation,term_pump,section_physical))
    ENDIF
  END SUBROUTINE print_physical_detail

  SUBROUTINE print_derived_physical_detail(name, balance, content_units, &
       &rate_units)
    CHARACTER(LEN=*), INTENT(IN) :: name, content_units, rate_units
    TYPE(physical_balance_type), INTENT(IN) :: balance

    WRITE(6,'(A)') '  Physical '//TRIM(name)//' (derived)'
    CALL print_component('content ['//TRIM(content_units)//']',balance%content)
    CALL print_component('temporal ['//TRIM(rate_units)//']',balance%temporal)
    CALL print_component('volume',balance%volume)
    CALL print_component('boundary physical inward',balance%boundary_inward)
    CALL print_component('physical imbalance',balance%physical_imbalance)
    WRITE(6,'(A)') '    Discrete'
    CALL print_component('HDG tau inward',balance%tau_inward)
    CALL print_component('residual',balance%discrete_residual)
  END SUBROUTINE print_derived_physical_detail

  SUBROUTINE print_pair(label, first, second)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: first, second

    WRITE(6,'(A,2(1X,ES12.3))') '    '//TRIM(label),first,second
  END SUBROUTINE print_pair

  SUBROUTINE print_value(label, value)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: value

    WRITE(6,'(A,1X,ES12.3)') '    '//TRIM(label),value
  END SUBROUTINE print_value

  SUBROUTINE print_component(label, value)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: value

    WRITE(6,'(A,1X,ES12.3)') '      '//TRIM(label),value
  END SUBROUTINE print_component

  SUBROUTINE print_subcomponent(label, value)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: value

    WRITE(6,'(A,1X,ES12.3)') '        '//TRIM(label),value
  END SUBROUTINE print_subcomponent

END SUBMODULE balance_diagnostics_reporting
