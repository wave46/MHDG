SUBMODULE (balance_diagnostics) balance_diagnostics_reporting
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE balance_diagnostics_report(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type) :: balances(balance_equation_count)
    TYPE(physical_balance_type) :: total_particles, total_energy

    IF (.NOT. this%has_values()) RETURN
    IF (MPIvar%glob_id /= 0) RETURN

    CALL collect_balances(this,balances,total_particles,total_energy)
    SELECT CASE (this%get_mode())
    CASE (balance_mode_summary)
       CALL print_summary(this,balances,total_particles,total_energy)
    CASE (balance_mode_equations)
       CALL print_equation_tables(this,balances,total_particles,total_energy, &
            &'equations')
    CASE (balance_mode_detailed)
       CALL print_detailed(this,balances,total_particles,total_energy)
    END SELECT
  END SUBROUTINE balance_diagnostics_report

  SUBROUTINE collect_balances(this, balances, total_particles, total_energy)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type), INTENT(OUT) :: &
         &balances(balance_equation_count), total_particles, total_energy
    INTEGER :: equation

    DO equation = 1,balance_equation_count
       balances(equation) = physical_balance(this,equation)
    ENDDO
    total_particles = add_physical_balances(balances(equation_n), &
         &balances(equation_nn))
    total_energy = add_physical_balances(balances(equation_nEi), &
         &balances(equation_nEe))
  END SUBROUTINE collect_balances

  SUBROUTINE print_summary(this, balances, total_particles, total_energy)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type), INTENT(IN) :: &
         &balances(balance_equation_count), total_particles, total_energy

    WRITE(6,'(A)') 'Balance diagnostics (summary)'
    CALL print_inventory(balances,total_particles,total_energy,'Inventory')

    WRITE(6,'(A)') '  External sources'
    WRITE(6,'(A,2(2X,A,1X,ES11.3))') '    particles [particles/s]', &
         &'puff in',external_puff_input(this), &
         &'pump out',external_pump_output(this)
    WRITE(6,'(A,2(2X,A,1X,ES11.3))') '      volume', &
         &'n',prescribed_source(this,equation_n), &
         &'n_n',prescribed_source(this,equation_nn)
    WRITE(6,'(A,2X,A,1X,ES11.3)') '    momentum [N]', &
         &'nu',prescribed_source(this,equation_nu)
    WRITE(6,'(A,2(2X,A,1X,ES11.3))') '    energy [W]', &
         &'nEi',prescribed_source(this,equation_nEi), &
         &'nEe',prescribed_source(this,equation_nEe)

    WRITE(6,'(A)') '  Balances'
    WRITE(6,'(A,2(2X,A,1X,ES11.3))') '    physical imbalance', &
         &'total particles [particles/s]',total_particles%physical_imbalance, &
         &'total plasma energy [W]',total_energy%physical_imbalance
  END SUBROUTINE print_summary

  SUBROUTINE print_inventory(balances, total_particles, total_energy, title)
    TYPE(physical_balance_type), INTENT(IN) :: &
         &balances(balance_equation_count), total_particles, total_energy
    CHARACTER(LEN=*), INTENT(IN) :: title

    WRITE(6,'(A)') '  '//TRIM(title)
    WRITE(6,'(A,3(2X,A,1X,ES11.3))') '    particles [particles]', &
         &'n',balances(equation_n)%content, &
         &'n_n',balances(equation_nn)%content, &
         &'n+n_n',total_particles%content
    WRITE(6,'(A,2X,A,1X,ES11.3)') '    momentum [kg m s^-1]', &
         &'nu',balances(equation_nu)%content
    WRITE(6,'(A,3(2X,A,1X,ES11.3))') '    plasma energy [J]', &
         &'nEi',balances(equation_nEi)%content, &
         &'nEe',balances(equation_nEe)%content, &
         &'nEi+nEe',total_energy%content
  END SUBROUTINE print_inventory

  REAL*8 FUNCTION prescribed_source(this, equation) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation

    value = balance_value(this,equation,term_prescribed_source, &
         &section_physical)
  END FUNCTION prescribed_source

  SUBROUTINE print_equation_tables(this, balances, total_particles, &
       &total_energy, mode_name)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type), INTENT(IN) :: &
         &balances(balance_equation_count), total_particles, total_energy
    CHARACTER(LEN=*), INTENT(IN) :: mode_name

    WRITE(6,'(A)') 'Balance diagnostics ('//TRIM(mode_name)//')'
    WRITE(6,'(A)') '  Particle equations [content: particles; rate: particles/s]'
    CALL print_equation_header()
    CALL print_equation_row(this,'n',equation_n,balances(equation_n))
    CALL print_equation_row(this,'n_n',equation_nn,balances(equation_nn))
    CALL print_derived_row('total_n',total_particles)

    WRITE(6,'(A)') '  Momentum equation [content: kg m s^-1; rate: N]'
    CALL print_equation_header()
    CALL print_equation_row(this,'nu',equation_nu,balances(equation_nu))

    WRITE(6,'(A)') '  Energy equations [content: J; rate: W]'
    CALL print_equation_header()
    CALL print_equation_row(this,'nEi',equation_nEi,balances(equation_nEi))
    CALL print_equation_row(this,'nEe',equation_nEe,balances(equation_nEe))
    CALL print_derived_row('total_E',total_energy)
  END SUBROUTINE print_equation_tables

  SUBROUTINE print_equation_header()
    WRITE(6,'(4X,A10,3(1X,A13),1X,A12)') &
         &'equation','content','physical','discrete','BC'
    WRITE(6,'(4X,A10,3(1X,A13),1X,A12)') &
         &'','','imbalance','residual','residual'
  END SUBROUTINE print_equation_header

  SUBROUTINE print_equation_row(this, name, equation, balance)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: equation
    TYPE(physical_balance_type), INTENT(IN) :: balance
    CHARACTER(LEN=12) :: bc_text
    LOGICAL :: bc_available
    REAL*8 :: bc_residual

    bc_residual = equation_bc_residual(this,equation,bc_available)
    bc_text = ADJUSTR('--')
    IF (bc_available) WRITE(bc_text,'(ES12.4)') bc_residual
    WRITE(6,'(4X,A10,3(1X,ES13.4),1X,A12)') TRIM(name),balance%content, &
         &balance%physical_imbalance,balance%discrete_residual,bc_text
  END SUBROUTINE print_equation_row

  SUBROUTINE print_derived_row(name, balance)
    CHARACTER(LEN=*), INTENT(IN) :: name
    TYPE(physical_balance_type), INTENT(IN) :: balance

    WRITE(6,'(4X,A10,3(1X,ES13.4),1X,A12)') TRIM(name),balance%content, &
         &balance%physical_imbalance,balance%discrete_residual,'--'
  END SUBROUTINE print_derived_row

  SUBROUTINE print_detailed(this, balances, total_particles, total_energy)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type), INTENT(IN) :: balances(balance_equation_count)
    TYPE(physical_balance_type), INTENT(IN) :: total_particles, total_energy

    WRITE(6,'(A)') 'Balance diagnostics (detailed)'
    CALL print_inventory(balances,total_particles,total_energy, &
         &'Conserved contents')
    CALL print_physical_balances(this,balances,total_particles,total_energy)
    CALL print_discrete_equations(this,balances,total_particles,total_energy)
    CALL print_independent_bc_checks(this)
  END SUBROUTINE print_detailed

  SUBROUTINE print_physical_balances(this, balances, total_particles, &
       &total_energy)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type), INTENT(IN) :: &
         &balances(balance_equation_count), total_particles, total_energy
    INTEGER :: equation

    WRITE(6,'(A)') '  Physical balances'
    WRITE(6,'(A)') &
         &'    volume + physical boundary inward - temporal = imbalance'
    DO equation = equation_n,equation_nn
       WRITE(6,'(A)') '    '//TRIM(balance_equation_names(equation))//' ['// &
            &TRIM(balance_rate_units(equation))//']'
       CALL print_physical_aggregate(balances(equation))
       CALL print_volume_components(this,equation)
       CALL print_physical_boundary_components(this,equation)
       IF (equation == equation_n) CALL print_particle_exchange(this)
    ENDDO
    CALL print_physical_total('n+n_n','particles/s',total_particles)
    CALL print_physical_total('nEi+nEe','W',total_energy)
    CALL print_energy_exchange_cancellation(this)
  END SUBROUTINE print_physical_balances

  SUBROUTINE print_physical_aggregate(balance)
    TYPE(physical_balance_type), INTENT(IN) :: balance

    WRITE(6,'(A,4(1X,A13))') '              ', &
         &'temporal','volume','boundary','imbalance'
    WRITE(6,'(A,4(1X,ES13.4))') '              ',balance%temporal, &
         &balance%volume,balance%boundary_inward,balance%physical_imbalance
  END SUBROUTINE print_physical_aggregate

  SUBROUTINE print_physical_total(name, units, balance)
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    TYPE(physical_balance_type), INTENT(IN) :: balance

    WRITE(6,'(A)') '    '//TRIM(name)//' ['//TRIM(units)//'] (derived)'
    CALL print_physical_aggregate(balance)
  END SUBROUTINE print_physical_total

  SUBROUTINE print_energy_exchange_cancellation(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    REAL*8 :: ion, electron

    WRITE(6,'(A)') '    Internal energy exchange [W]'
    WRITE(6,'(A,3(1X,A13))') '              ','ion','electron','sum'
    ion = term_value(this,equation_nEi,term_parallel_electric_work)
    electron = term_value(this,equation_nEe,term_parallel_electric_work)
    WRITE(6,'(A18,3(1X,ES13.4))') 'parallel electric',ion,electron,ion+electron
    ion = term_value(this,equation_nEi,term_temperature_exchange)
    electron = term_value(this,equation_nEe,term_temperature_exchange)
    WRITE(6,'(A18,3(1X,ES13.4))') 'temperature',ion,electron,ion+electron
  END SUBROUTINE print_energy_exchange_cancellation

  SUBROUTINE print_discrete_equations(this, balances, total_particles, &
       &total_energy)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    TYPE(physical_balance_type), INTENT(IN) :: &
         &balances(balance_equation_count), total_particles, total_energy
    INTEGER :: equation

    WRITE(6,'(A)') '  Discrete equations'
    WRITE(6,'(A)') &
         &'    volume + numerical boundary - temporal = residual'
    WRITE(6,'(A)') &
         &'    equation boundary + HDG tau = numerical boundary'
    DO equation = equation_n,equation_nn
       WRITE(6,'(A)') '    '//TRIM(balance_equation_names(equation))//' ['// &
            &TRIM(balance_rate_units(equation))//']'
       CALL print_discrete_aggregate(balances(equation))
       CALL print_equation_boundary_components(this,equation)
    ENDDO
    CALL print_discrete_total('n+n_n','particles/s',total_particles)
    CALL print_discrete_total('nEi+nEe','W',total_energy)
  END SUBROUTINE print_discrete_equations

  SUBROUTINE print_discrete_aggregate(balance)
    TYPE(physical_balance_type), INTENT(IN) :: balance

    WRITE(6,'(A,6(1X,A13))') '              ', &
         &'temporal','volume','equation bnd.','HDG tau','numerical','residual'
    WRITE(6,'(A,6(1X,ES13.4))') '              ',balance%temporal, &
         &balance%volume,balance%equation_boundary_inward, &
         &balance%tau_inward,balance%numerical_boundary_inward, &
         &balance%discrete_residual
  END SUBROUTINE print_discrete_aggregate

  SUBROUTINE print_discrete_total(name, units, balance)
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    TYPE(physical_balance_type), INTENT(IN) :: balance

    WRITE(6,'(A)') '    '//TRIM(name)//' ['//TRIM(units)//'] (derived)'
    CALL print_discrete_aggregate(balance)
  END SUBROUTINE print_discrete_total

  SUBROUTINE print_independent_bc_checks(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    WRITE(6,'(A)') '  Independent boundary-condition checks'
    CALL print_bc_components(this,equation_n)
    CALL print_bc_components(this,equation_nn)
  END SUBROUTINE print_independent_bc_checks

  SUBROUTINE print_volume_components(this, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation

    SELECT CASE (equation)
    CASE (equation_n)
       CALL print_component_pairs('volume', &
            [CHARACTER(LEN=24) :: 'ionization','recombination', &
            &'prescribed source'], &
            (/term_value(this,equation,term_ionization), &
            &term_value(this,equation,term_recombination), &
            &term_value(this,equation,term_prescribed_source)/))
    CASE (equation_nu)
       CALL print_component_pairs('volume', &
            [CHARACTER(LEN=24) :: 'ionization','recombination', &
            &'charge exchange','pressure divergence','prescribed source'], &
            (/term_value(this,equation,term_ionization), &
            &term_value(this,equation,term_recombination), &
            &term_value(this,equation,term_charge_exchange), &
            &term_value(this,equation,term_pressure_divergence), &
            &term_value(this,equation,term_prescribed_source)/))
    CASE (equation_nEi)
       CALL print_component_pairs('volume', &
            [CHARACTER(LEN=24) :: 'ionization','recombination', &
            &'charge exchange','parallel electric work', &
            &'temperature exchange','prescribed source'], &
            (/term_value(this,equation,term_ionization), &
            &term_value(this,equation,term_recombination), &
            &term_value(this,equation,term_charge_exchange), &
            &term_value(this,equation,term_parallel_electric_work), &
            &term_value(this,equation,term_temperature_exchange), &
            &term_value(this,equation,term_prescribed_source)/))
    CASE (equation_nEe)
       CALL print_component_pairs('volume', &
            [CHARACTER(LEN=24) :: 'ionization','recombination','radiation', &
            &'ohmic','parallel electric work','temperature exchange', &
            &'prescribed source'], &
            (/term_value(this,equation,term_ionization), &
            &term_value(this,equation,term_recombination), &
            &term_value(this,equation,term_radiation), &
            &term_value(this,equation,term_ohmic), &
            &term_value(this,equation,term_parallel_electric_work), &
            &term_value(this,equation,term_temperature_exchange), &
            &term_value(this,equation,term_prescribed_source)/))
    CASE (equation_nn)
       CALL print_component_pairs('volume', &
            [CHARACTER(LEN=24) :: 'ionization','recombination', &
            &'prescribed source','puff','pump'], &
            (/term_value(this,equation,term_ionization), &
            &term_value(this,equation,term_recombination), &
            &term_value(this,equation,term_prescribed_source), &
            &term_value(this,equation,term_puff), &
            &term_value(this,equation,term_pump)/))
    END SELECT
  END SUBROUTINE print_volume_components

  SUBROUTINE print_physical_boundary_components(this, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation

    SELECT CASE (equation)
    CASE (equation_n)
       CALL print_component_pairs('wall', &
            [CHARACTER(LEN=24) :: 'parallel convection','pinch'], &
            (/section_value(this,equation,term_parallel,section_bc), &
            &section_value(this,equation,term_pinch,section_bc)/))
    CASE (equation_nu)
       CALL print_transport_boundary_components(this,equation,.FALSE.,'wall')
    CASE (equation_nEi,equation_nEe)
       CALL print_component_pairs('wall', &
            [CHARACTER(LEN=24) :: 'sheath','pinch'], &
            (/section_value(this,equation,term_sheath,section_bc), &
            &section_value(this,equation,term_pinch,section_bc)/))
    CASE (equation_nn)
       CALL print_component_pairs('wall source', &
            [CHARACTER(LEN=24) :: 'recycling parallel', &
            &'recycling diffusion','recycling pinch','puff','pump'], &
            (/section_value(this,equation,term_recycling_parallel,section_bc), &
            &section_value(this,equation,term_recycling_diffusion,section_bc), &
            &section_value(this,equation,term_recycling_pinch,section_bc), &
            &section_value(this,equation,term_puff,section_bc), &
            &-section_value(this,equation,term_pump,section_bc)/))
    END SELECT
  END SUBROUTINE print_physical_boundary_components

  SUBROUTINE print_equation_boundary_components(this, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation

    SELECT CASE (equation)
    CASE (equation_n)
       CALL print_component_pairs('equation boundary', &
            [CHARACTER(LEN=24) :: 'parallel convection','diffusion','pinch'], &
            (/term_value(this,equation,term_parallel), &
            &term_value(this,equation,term_diffusion), &
            &term_value(this,equation,term_pinch)/))
    CASE (equation_nu)
       CALL print_transport_boundary_components(this,equation,.FALSE., &
            &'equation boundary')
    CASE (equation_nEi,equation_nEe)
       CALL print_transport_boundary_components(this,equation,.TRUE., &
            &'equation boundary')
    CASE (equation_nn)
#if defined(NEUTRALP) && defined(NEUTRALGAMMA)
       CALL print_component_pairs('equation boundary', &
            [CHARACTER(LEN=24) :: 'limited diffusion','limited pressure', &
            &'neutral-gamma convection'], &
            (/term_value(this,equation,term_diffusion), &
            &term_value(this,equation,term_pressure), &
            &term_value(this,equation,term_convection)/))
#elif defined(NEUTRALP)
       CALL print_component_pairs('equation boundary', &
            [CHARACTER(LEN=24) :: 'limited diffusion','limited pressure'], &
            (/term_value(this,equation,term_diffusion), &
            &term_value(this,equation,term_pressure)/))
#elif defined(NEUTRALGAMMA)
       CALL print_component_pairs('equation boundary', &
            [CHARACTER(LEN=24) :: 'limited diffusion', &
            &'neutral-gamma convection'], &
            (/term_value(this,equation,term_diffusion), &
            &term_value(this,equation,term_convection)/))
#else
       CALL print_component_pairs('equation boundary', &
            [CHARACTER(LEN=24) :: 'limited diffusion'], &
            (/term_value(this,equation,term_diffusion)/))
#endif
    END SELECT
  END SUBROUTINE print_equation_boundary_components

  SUBROUTINE print_transport_boundary_components(this, equation, energy, title)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation
    LOGICAL, INTENT(IN) :: energy
    CHARACTER(LEN=*), INTENT(IN) :: title

    IF (energy) THEN
       CALL print_component_pairs(title, &
            [CHARACTER(LEN=24) :: 'convection','diffusion', &
            &'parallel conduction','pinch'], &
            (/term_value(this,equation,term_convection), &
            &term_value(this,equation,term_diffusion), &
            &term_value(this,equation,term_parallel_conduction), &
            &term_value(this,equation,term_pinch)/))
    ELSE
       CALL print_component_pairs(title, &
            [CHARACTER(LEN=24) :: 'convection','diffusion','pinch'], &
            (/term_value(this,equation,term_convection), &
            &term_value(this,equation,term_diffusion), &
            &term_value(this,equation,term_pinch)/))
    ENDIF
  END SUBROUTINE print_transport_boundary_components

  SUBROUTINE print_bc_components(this, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation
    TYPE(particle_bc_balance_type) :: neutral

    SELECT CASE (equation)
    CASE (equation_n)
       WRITE(6,'(A)') '    n [particles/s]'
       WRITE(6,'(A)') '      diffusion + HDG tau = residual'
       CALL print_component_pairs('terms', &
            [CHARACTER(LEN=24) :: 'diffusion','HDG tau','residual'], &
            (/section_value(this,equation,term_diffusion,section_bc), &
            &section_value(this,equation,term_tau_inward,section_bc), &
            &section_value(this,equation,term_diffusion,section_bc)+ &
            &section_value(this,equation,term_tau_inward,section_bc)/))
    CASE (equation_nn)
       neutral = neutral_density_bc(this)
       WRITE(6,'(A)') '    n_n [particles/s]'
       WRITE(6,'(A)') &
            &'      imposed source - physical flux - HDG tau = residual'
       CALL print_component_pairs('balance', &
            [CHARACTER(LEN=24) :: 'imposed source','physical flux', &
            &'HDG tau','residual'], &
            (/neutral%imposed_source_inward,neutral%physical_flux_inward, &
            &neutral%tau_inward,neutral%residual/))
       CALL print_component_pairs('imposed source', &
            [CHARACTER(LEN=24) :: 'recycling parallel', &
            &'recycling diffusion','recycling pinch','puff','pump'], &
            (/section_value(this,equation,term_recycling_parallel,section_bc), &
            &section_value(this,equation,term_recycling_diffusion,section_bc), &
            &section_value(this,equation,term_recycling_pinch,section_bc), &
            &section_value(this,equation,term_puff,section_bc), &
            &section_value(this,equation,term_pump,section_bc)/))
       CALL print_neutral_bc_flux_components(this)
    END SELECT
  END SUBROUTINE print_bc_components

  SUBROUTINE print_neutral_bc_flux_components(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

#if defined(NEUTRALP) && defined(NEUTRALGAMMA)
    CALL print_component_pairs('physical flux', &
         [CHARACTER(LEN=24) :: 'limited diffusion','limited pressure', &
         &'neutral-gamma convection'], &
         (/section_value(this,equation_nn,term_diffusion,section_bc), &
         &section_value(this,equation_nn,term_pressure,section_bc), &
         &section_value(this,equation_nn,term_convection,section_bc)/))
#elif defined(NEUTRALP)
    CALL print_component_pairs('physical flux', &
         [CHARACTER(LEN=24) :: 'limited diffusion','limited pressure'], &
         (/section_value(this,equation_nn,term_diffusion,section_bc), &
         &section_value(this,equation_nn,term_pressure,section_bc)/))
#elif defined(NEUTRALGAMMA)
    CALL print_component_pairs('physical flux', &
         [CHARACTER(LEN=24) :: 'limited diffusion', &
         &'neutral-gamma convection'], &
         (/section_value(this,equation_nn,term_diffusion,section_bc), &
         &section_value(this,equation_nn,term_convection,section_bc)/))
#else
    CALL print_component_pairs('physical flux', &
         [CHARACTER(LEN=24) :: 'limited diffusion'], &
         (/section_value(this,equation_nn,term_diffusion,section_bc)/))
#endif
  END SUBROUTINE print_neutral_bc_flux_components

  SUBROUTINE print_particle_exchange(this)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this

    CALL print_component_pairs('exchange', &
         [CHARACTER(LEN=24) :: 'charge exchange rate'], &
         (/term_value(this,equation_n,term_charge_exchange)/))
  END SUBROUTINE print_particle_exchange

  SUBROUTINE print_component_pairs(title, labels, values)
    CHARACTER(LEN=*), INTENT(IN) :: title
    CHARACTER(LEN=*), INTENT(IN) :: labels(:)
    REAL*8, INTENT(IN) :: values(:)
    INTEGER :: item

    WRITE(6,'(A)') '      '//TRIM(title)
    DO item = 1,SIZE(values),2
       IF (item < SIZE(values)) THEN
          WRITE(6,'(8X,2(A24,1X,ES12.4,2X))') labels(item),values(item), &
               &labels(item+1),values(item+1)
       ELSE
          WRITE(6,'(8X,A24,1X,ES12.4)') labels(item),values(item)
       ENDIF
    ENDDO
  END SUBROUTINE print_component_pairs

  REAL*8 FUNCTION term_value(this, equation, term) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation, term

    value = balance_value(this,equation,term,section_physical)
  END FUNCTION term_value

  REAL*8 FUNCTION section_value(this, equation, term, section) RESULT(value)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: equation, term, section

    value = balance_value(this,equation,term,section)
  END FUNCTION section_value

END SUBMODULE balance_diagnostics_reporting
