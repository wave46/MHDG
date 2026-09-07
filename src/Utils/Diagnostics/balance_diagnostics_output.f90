SUBMODULE (balance_diagnostics) balance_diagnostics_output
  USE HDF5_io_module, ONLY: HDF5_group_create, HDF5_group_close, &
       &HDF5_real_saving, HDF5_string_saving
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE balance_diagnostics_write_hdf5(this, &
       &diagnostics_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id

    IF (.NOT. this%has_values()) RETURN
    IF (MPIvar%glob_id /= 0) RETURN

    SELECT CASE (this%get_mode())
    CASE (balance_mode_summary)
       CALL write_summary_hdf5(this,diagnostics_group_id)
    CASE (balance_mode_equations,balance_mode_detailed)
       CALL write_equations_hdf5(this,diagnostics_group_id)
    END SELECT
  END SUBROUTINE balance_diagnostics_write_hdf5

  SUBROUTINE write_summary_hdf5(this, diagnostics_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: summary_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('summary',diagnostics_group_id, &
         &summary_group_id,ierr)
    CALL write_summary_content(this,summary_group_id)
    CALL write_summary_external_sources(this,summary_group_id)
    CALL write_summary_balances(this,summary_group_id)
    CALL HDF5_group_close(summary_group_id,ierr)
  END SUBROUTINE write_summary_hdf5

  SUBROUTINE write_summary_content(this, summary_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: summary_group_id
    INTEGER(HID_T) :: content_group_id, family_group_id
    TYPE(physical_balance_type) :: n, nn, nu, nEi, nEe
    INTEGER :: ierr

    n = physical_balance(this,equation_n)
    nn = physical_balance(this,equation_nn)
    nu = physical_balance(this,equation_nu)
    nEi = physical_balance(this,equation_nEi)
    nEe = physical_balance(this,equation_nEe)
    CALL HDF5_group_create('content',summary_group_id,content_group_id,ierr)

    CALL HDF5_group_create('particles',content_group_id,family_group_id,ierr)
    CALL HDF5_string_saving(family_group_id,'particles','units')
    CALL HDF5_real_saving(family_group_id,n%content,'n')
    CALL HDF5_real_saving(family_group_id,nn%content,'n_n')
    CALL HDF5_real_saving(family_group_id,n%content+nn%content,'total_n')
    CALL HDF5_group_close(family_group_id,ierr)

    CALL HDF5_group_create('momentum',content_group_id,family_group_id,ierr)
    CALL HDF5_string_saving(family_group_id,'kg m s^-1','units')
    CALL HDF5_real_saving(family_group_id,nu%content,'nu')
    CALL HDF5_group_close(family_group_id,ierr)

    CALL HDF5_group_create('plasma_energy',content_group_id, &
         &family_group_id,ierr)
    CALL HDF5_string_saving(family_group_id,'J','units')
    CALL HDF5_real_saving(family_group_id,nEi%content,'nEi')
    CALL HDF5_real_saving(family_group_id,nEe%content,'nEe')
    CALL HDF5_real_saving(family_group_id,nEi%content+nEe%content,'total_E')
    CALL HDF5_group_close(family_group_id,ierr)
    CALL HDF5_group_close(content_group_id,ierr)
  END SUBROUTINE write_summary_content

  SUBROUTINE write_summary_external_sources(this, summary_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: summary_group_id
    INTEGER(HID_T) :: sources_group_id, family_group_id
    TYPE(plasma_heating_summary_type) :: heating
    INTEGER :: ierr

    heating = plasma_heating_summary(this)
    CALL HDF5_group_create('external_sources',summary_group_id, &
         &sources_group_id,ierr)
    CALL HDF5_group_create('particles',sources_group_id,family_group_id,ierr)
    CALL HDF5_string_saving(family_group_id,'particles/s','units')
    CALL HDF5_real_saving(family_group_id,external_puff_input(this),'puff_in')
    CALL HDF5_real_saving(family_group_id,external_pump_output(this),'pump_out')
    CALL write_prescribed_source(this,family_group_id,equation_n,'volume_n')
    CALL write_prescribed_source(this,family_group_id,equation_nn,'volume_n_n')
    CALL HDF5_group_close(family_group_id,ierr)

    CALL HDF5_group_create('momentum',sources_group_id,family_group_id,ierr)
    CALL HDF5_string_saving(family_group_id,'N','units')
    CALL write_prescribed_source(this,family_group_id,equation_nu,'volume_nu')
    CALL HDF5_group_close(family_group_id,ierr)

    CALL HDF5_group_create('energy',sources_group_id,family_group_id,ierr)
    CALL HDF5_string_saving(family_group_id,'W','units')
    CALL HDF5_real_saving(family_group_id,heating%total,'total_heating')
    CALL HDF5_real_saving(family_group_id,heating%ion,'ion_heating')
    CALL HDF5_real_saving(family_group_id,heating%electron,'electron_heating')
    CALL HDF5_real_saving(family_group_id,heating%ohmic,'ohmic_heating')
    CALL HDF5_real_saving(family_group_id,heating%external,'external')
    CALL HDF5_group_close(family_group_id,ierr)
    CALL HDF5_group_close(sources_group_id,ierr)
  END SUBROUTINE write_summary_external_sources

  SUBROUTINE write_prescribed_source(this, group_id, equation, name)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation
    CHARACTER(LEN=*), INTENT(IN) :: name

    CALL HDF5_real_saving(group_id,balance_value(this,equation, &
         &term_prescribed_source,section_physical),name)
  END SUBROUTINE write_prescribed_source

  SUBROUTINE write_summary_balances(this, summary_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: summary_group_id
    INTEGER(HID_T) :: balances_group_id
    TYPE(physical_balance_type) :: total
    INTEGER :: ierr

    CALL HDF5_group_create('balances',summary_group_id, &
         &balances_group_id,ierr)
    total = add_physical_balances(physical_balance(this,equation_n), &
         &physical_balance(this,equation_nn))
    CALL write_summary_balance(balances_group_id,'total_n', &
         &'particles/s',total%physical_imbalance)
    total = add_physical_balances(physical_balance(this,equation_nEi), &
         &physical_balance(this,equation_nEe))
    CALL write_summary_balance(balances_group_id,'total_E','W', &
         &total%physical_imbalance)
    CALL HDF5_group_close(balances_group_id,ierr)
  END SUBROUTINE write_summary_balances

  SUBROUTINE write_summary_balance(parent_group_id, name, units, imbalance)
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    REAL*8, INTENT(IN) :: imbalance
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create(name,parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    CALL HDF5_real_saving(group_id,imbalance,'physical_imbalance')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_summary_balance

  SUBROUTINE write_equations_hdf5(this, diagnostics_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: equations_group_id
    TYPE(physical_balance_type) :: total
    INTEGER :: equation, ierr

    CALL HDF5_group_create('equations',diagnostics_group_id, &
         &equations_group_id,ierr)
    DO equation = equation_n,equation_nn
       CALL write_primary_equation(this,equations_group_id,equation)
    ENDDO
    total = add_physical_balances(physical_balance(this,equation_n), &
         &physical_balance(this,equation_nn))
    CALL write_equation(this,equations_group_id,'total_n',total, &
         &'particles','particles/s')
    total = add_physical_balances(physical_balance(this,equation_nEi), &
         &physical_balance(this,equation_nEe))
    CALL write_equation(this,equations_group_id,'total_E',total,'J','W', &
         &bc_equations=(/equation_nEi,equation_nEe/))
    CALL HDF5_group_close(equations_group_id,ierr)
  END SUBROUTINE write_equations_hdf5

  SUBROUTINE write_primary_equation(this, equations_group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equations_group_id
    INTEGER, INTENT(IN) :: equation

    CALL write_equation(this,equations_group_id, &
         &balance_equation_names(equation),physical_balance(this,equation), &
         &balance_content_units(equation),balance_rate_units(equation),equation)
  END SUBROUTINE write_primary_equation

  SUBROUTINE write_equation(this, parent_group_id, name, balance, &
       &content_units, rate_units, equation, bc_equations)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    CHARACTER(LEN=*), INTENT(IN) :: name, content_units, rate_units
    TYPE(physical_balance_type), INTENT(IN) :: balance
    INTEGER, INTENT(IN), OPTIONAL :: equation
    INTEGER, INTENT(IN), OPTIONAL :: bc_equations(:)
    INTEGER(HID_T) :: equation_group_id
    INTEGER :: ierr

    CALL HDF5_group_create(name,parent_group_id,equation_group_id,ierr)
    CALL HDF5_string_saving(equation_group_id,content_units,'content_units')
    CALL HDF5_string_saving(equation_group_id,rate_units,'rate_units')
    CALL HDF5_real_saving(equation_group_id,balance%content,'content')
    CALL write_physical_view(this,equation_group_id,balance,rate_units,equation)
    CALL write_discrete_view(this,equation_group_id,balance,rate_units,equation)
    IF (PRESENT(equation)) CALL write_bc_view(this,equation_group_id,equation, &
         &rate_units)
    IF (PRESENT(bc_equations)) CALL write_combined_bc_view(this, &
         &equation_group_id,bc_equations,rate_units)
    CALL HDF5_group_close(equation_group_id,ierr)
  END SUBROUTINE write_equation

  SUBROUTINE write_physical_view(this, equation_group_id, balance, units, &
       &equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    TYPE(physical_balance_type), INTENT(IN) :: balance
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER, INTENT(IN), OPTIONAL :: equation
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('physical',equation_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    CALL HDF5_real_saving(group_id,balance%physical_imbalance,'imbalance')
    IF (this%detailed()) THEN
       CALL HDF5_real_saving(group_id,balance%temporal,'temporal')
       CALL HDF5_real_saving(group_id,balance%volume,'volume')
       CALL HDF5_real_saving(group_id,balance%boundary_inward, &
            &'boundary_inward')
       IF (PRESENT(equation)) CALL write_physical_components(this,group_id, &
            &equation,units)
    ENDIF
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_physical_view

  SUBROUTINE write_discrete_view(this, equation_group_id, balance, units, &
       &equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    TYPE(physical_balance_type), INTENT(IN) :: balance
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER, INTENT(IN), OPTIONAL :: equation
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('discrete',equation_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    CALL HDF5_real_saving(group_id,balance%discrete_residual,'residual')
    IF (this%detailed()) THEN
       CALL HDF5_real_saving(group_id,balance%equation_boundary_inward, &
            &'equation_boundary_inward')
       CALL HDF5_real_saving(group_id,balance%tau_stabilization_inward, &
            &'tau_stabilization_inward')
       CALL HDF5_real_saving(group_id,balance%numerical_boundary_inward, &
            &'numerical_boundary_inward')
       IF (PRESENT(equation)) CALL write_equation_boundary_components( &
            &this,group_id,equation,units)
    ENDIF
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_discrete_view

  SUBROUTINE write_physical_components(this, physical_group_id, equation, &
       &units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: physical_group_id
    INTEGER, INTENT(IN) :: equation
    CHARACTER(LEN=*), INTENT(IN) :: units

    CALL write_volume_components(this,physical_group_id,equation,units)
    CALL write_physical_boundary_components(this,physical_group_id,equation, &
         &units)
    IF (equation == equation_n) CALL write_particle_exchange(this, &
         &physical_group_id)
  END SUBROUTINE write_physical_components

  SUBROUTINE write_volume_components(this, parent_group_id, equation, units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER, INTENT(IN) :: equation
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('volume_components',parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    SELECT CASE (equation)
    CASE (equation_n)
       CALL write_term(this,group_id,equation,term_ionization,'ionization')
       CALL write_term(this,group_id,equation,term_recombination,'recombination')
       CALL write_term(this,group_id,equation,term_prescribed_source, &
            &'prescribed_source')
    CASE (equation_nu)
       CALL write_term(this,group_id,equation,term_ionization,'ionization')
       CALL write_term(this,group_id,equation,term_recombination,'recombination')
       CALL write_term(this,group_id,equation,term_charge_exchange, &
            &'charge_exchange')
       CALL write_term(this,group_id,equation,term_pressure_divergence, &
            &'pressure_divergence')
       CALL write_term(this,group_id,equation,term_prescribed_source, &
            &'prescribed_source')
    CASE (equation_nEi)
       CALL write_ion_energy_volume_terms(this,group_id,equation)
    CASE (equation_nEe)
       CALL write_electron_energy_volume_terms(this,group_id,equation)
    CASE (equation_nn)
       CALL write_term(this,group_id,equation,term_ionization,'ionization')
       CALL write_term(this,group_id,equation,term_recombination,'recombination')
       CALL write_term(this,group_id,equation,term_prescribed_source, &
            &'prescribed_source')
       CALL write_term(this,group_id,equation,term_puff,'puff')
       CALL write_term(this,group_id,equation,term_pump,'pump')
    END SELECT
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_volume_components

  SUBROUTINE write_ion_energy_volume_terms(this, group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation

    CALL write_term(this,group_id,equation,term_ionization,'ionization')
    CALL write_term(this,group_id,equation,term_recombination,'recombination')
    CALL write_term(this,group_id,equation,term_charge_exchange,'charge_exchange')
    CALL write_term(this,group_id,equation,term_parallel_electric_work, &
         &'parallel_electric_work')
    CALL write_term(this,group_id,equation,term_temperature_exchange, &
         &'temperature_exchange')
    CALL write_term(this,group_id,equation,term_prescribed_source, &
         &'prescribed_source')
  END SUBROUTINE write_ion_energy_volume_terms

  SUBROUTINE write_electron_energy_volume_terms(this, group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation

    CALL write_term(this,group_id,equation,term_ionization,'ionization')
    CALL write_term(this,group_id,equation,term_recombination,'recombination')
    CALL write_term(this,group_id,equation,term_radiation,'radiation')
    CALL write_term(this,group_id,equation,term_ohmic,'ohmic')
    CALL write_term(this,group_id,equation,term_parallel_electric_work, &
         &'parallel_electric_work')
    CALL write_term(this,group_id,equation,term_temperature_exchange, &
         &'temperature_exchange')
    CALL write_term(this,group_id,equation,term_prescribed_source, &
         &'prescribed_source')
  END SUBROUTINE write_electron_energy_volume_terms

  SUBROUTINE write_physical_boundary_components(this, parent_group_id, &
       &equation, units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER, INTENT(IN) :: equation
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    IF (equation == equation_nu) RETURN
    CALL HDF5_group_create('boundary_components_inward',parent_group_id, &
         &group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    SELECT CASE (equation)
    CASE (equation_n)
       CALL write_section_term(this,group_id,equation,term_parallel, &
            &section_bc,'parallel_convection')
       CALL write_section_term(this,group_id,equation,term_pinch, &
            &section_bc,'pinch')
    CASE (equation_nEi,equation_nEe)
       CALL write_section_term(this,group_id,equation,term_sheath, &
            &section_bc,'sheath')
       CALL write_section_term(this,group_id,equation,term_pinch, &
            &section_bc,'pinch')
    CASE (equation_nn)
       CALL write_neutral_physical_boundary(this,group_id,equation)
    END SELECT
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_physical_boundary_components

  SUBROUTINE write_neutral_physical_boundary(this, group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation

    CALL write_section_term(this,group_id,equation,term_recycling_parallel, &
         &section_bc,'recycling_parallel')
    CALL write_section_term(this,group_id,equation,term_recycling_diffusion, &
         &section_bc,'recycling_diffusion')
    CALL write_section_term(this,group_id,equation,term_recycling_pinch, &
         &section_bc,'recycling_pinch')
    CALL write_section_term(this,group_id,equation,term_puff,section_bc,'puff')
    CALL HDF5_real_saving(group_id,-balance_value(this,equation,term_pump, &
         &section_bc),'pump')
  END SUBROUTINE write_neutral_physical_boundary

  SUBROUTINE write_particle_exchange(this, parent_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('exchange',parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL HDF5_real_saving(group_id,balance_value(this,equation_n, &
         &term_charge_exchange,section_physical),'charge_exchange_rate')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_particle_exchange

  SUBROUTINE write_equation_boundary_components(this, parent_group_id, &
       &equation, units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER, INTENT(IN) :: equation
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('boundary_components_inward',parent_group_id, &
         &group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    SELECT CASE (equation)
    CASE (equation_n)
       CALL write_term(this,group_id,equation,term_parallel, &
            &'parallel_convection')
       CALL write_term(this,group_id,equation,term_diffusion,'diffusion')
       CALL write_term(this,group_id,equation,term_pinch,'pinch')
    CASE (equation_nu)
       CALL write_transport_boundary(this,group_id,equation,.FALSE.)
    CASE (equation_nEi,equation_nEe)
       CALL write_transport_boundary(this,group_id,equation,.TRUE.)
    CASE (equation_nn)
       CALL write_neutral_equation_boundary(this,group_id,equation)
    END SELECT
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_equation_boundary_components

  SUBROUTINE write_transport_boundary(this, group_id, equation, energy)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation
    LOGICAL, INTENT(IN) :: energy

    CALL write_term(this,group_id,equation,term_convection,'convection')
    CALL write_term(this,group_id,equation,term_diffusion,'diffusion')
    IF (energy) CALL write_term(this,group_id,equation, &
         &term_parallel_conduction,'parallel_conduction')
    CALL write_term(this,group_id,equation,term_pinch,'pinch')
  END SUBROUTINE write_transport_boundary

  SUBROUTINE write_neutral_equation_boundary(this, group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation

    CALL write_term(this,group_id,equation,term_diffusion,'limited_diffusion')
#ifdef NEUTRALP
    CALL write_term(this,group_id,equation,term_pressure,'limited_pressure')
#endif
#ifdef NEUTRALGAMMA
    CALL write_term(this,group_id,equation,term_convection, &
         &'neutral_gamma_convection')
#endif
  END SUBROUTINE write_neutral_equation_boundary

  SUBROUTINE write_bc_view(this, equation_group_id, equation, units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER, INTENT(IN) :: equation
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(HID_T) :: group_id
    LOGICAL :: available
    REAL*8 :: residual
    INTEGER :: ierr

    residual = equation_bc_residual(this,equation,available)
    IF (.NOT. available) RETURN
    CALL HDF5_group_create('bc',equation_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    CALL HDF5_real_saving(group_id,residual,'residual')
    IF (this%detailed()) THEN
       IF (equation == equation_n) CALL write_density_bc_detail(this,group_id)
       IF (equation == equation_nu) CALL write_momentum_bc_detail(this,group_id)
       IF (equation == equation_nEi .OR. equation == equation_nEe) &
            &CALL write_energy_bc_detail(this,group_id,equation)
       IF (equation == equation_nn) CALL write_neutral_bc_detail(this,group_id)
    ENDIF
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_bc_view

  SUBROUTINE write_density_bc_detail(this, group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id

    CALL HDF5_real_saving(group_id,balance_value(this,equation_n, &
         &term_diffusion,section_bc),'diffusion_inward')
    CALL HDF5_real_saving(group_id,balance_value(this,equation_n, &
         &term_tau_stabilization_inward,section_bc), &
         &'tau_stabilization_inward')
  END SUBROUTINE write_density_bc_detail

  SUBROUTINE write_momentum_bc_detail(this, group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id

    CALL write_bc_term(this,group_id,equation_nu,term_diffusion, &
         &'perpendicular_diffusion_inward')
    CALL write_bc_term(this,group_id,equation_nu,term_split_diffusion, &
         &'split_diffusion_inward')
    CALL write_bc_term(this,group_id,equation_nu, &
         &term_tau_stabilization_inward,'tau_stabilization_inward')
  END SUBROUTINE write_momentum_bc_detail

  SUBROUTINE write_energy_bc_detail(this, group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation

    CALL write_bc_term(this,group_id,equation,term_diffusion, &
         &'perpendicular_diffusion_inward')
    CALL write_bc_term(this,group_id,equation,term_split_diffusion, &
         &'split_diffusion_inward')
    CALL write_bc_term(this,group_id,equation,term_parallel_conduction, &
         &'parallel_conduction_inward')
    CALL write_bc_term(this,group_id,equation,term_sheath_minus_bulk, &
         &'sheath_minus_bulk_inward')
    CALL write_bc_term(this,group_id,equation, &
         &term_tau_stabilization_inward,'tau_stabilization_inward')
  END SUBROUTINE write_energy_bc_detail

  SUBROUTINE write_bc_term(this, group_id, equation, term, name)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation, term
    CHARACTER(LEN=*), INTENT(IN) :: name

    CALL HDF5_real_saving(group_id,balance_value(this,equation,term, &
         &section_bc),name)
  END SUBROUTINE write_bc_term

  SUBROUTINE write_combined_bc_view(this, equation_group_id, equations, units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER, INTENT(IN) :: equations(:)
    CHARACTER(LEN=*), INTENT(IN) :: units
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr, item
    LOGICAL :: available
    REAL*8 :: residual

    residual = 0.d0
    DO item = 1,SIZE(equations)
       residual = residual+equation_bc_residual(this,equations(item),available)
       IF (.NOT. available) RETURN
    ENDDO
    CALL HDF5_group_create('bc',equation_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,units,'units')
    CALL HDF5_real_saving(group_id,residual,'residual')
    IF (this%detailed()) CALL write_combined_energy_bc_detail(this,group_id, &
         &equations)
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_combined_bc_view

  SUBROUTINE write_combined_energy_bc_detail(this, group_id, equations)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equations(:)

    CALL write_combined_bc_term(this,group_id,equations,term_diffusion, &
         &'perpendicular_diffusion_inward')
    CALL write_combined_bc_term(this,group_id,equations, &
         &term_split_diffusion,'split_diffusion_inward')
    CALL write_combined_bc_term(this,group_id,equations, &
         &term_parallel_conduction,'parallel_conduction_inward')
    CALL write_combined_bc_term(this,group_id,equations, &
         &term_sheath_minus_bulk,'sheath_minus_bulk_inward')
    CALL write_combined_bc_term(this,group_id,equations, &
         &term_tau_stabilization_inward,'tau_stabilization_inward')
  END SUBROUTINE write_combined_energy_bc_detail

  SUBROUTINE write_combined_bc_term(this, group_id, equations, term, name)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equations(:), term
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL*8 :: value
    INTEGER :: item

    value = 0.d0
    DO item = 1,SIZE(equations)
       value = value+balance_value(this,equations(item),term,section_bc)
    ENDDO
    CALL HDF5_real_saving(group_id,value,name)
  END SUBROUTINE write_combined_bc_term

  SUBROUTINE write_neutral_bc_detail(this, group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    TYPE(particle_bc_balance_type) :: balance

    balance = neutral_density_bc(this)
    CALL HDF5_real_saving(group_id,balance%imposed_source_inward, &
         &'imposed_source_inward')
    CALL HDF5_real_saving(group_id,balance%physical_flux_inward, &
         &'physical_flux_inward')
    CALL HDF5_real_saving(group_id,balance%tau_stabilization_inward, &
         &'tau_stabilization_inward')
    CALL write_neutral_bc_sources(this,group_id)
    CALL write_neutral_bc_fluxes(this,group_id)
  END SUBROUTINE write_neutral_bc_detail

  SUBROUTINE write_neutral_bc_sources(this, parent_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('source_components',parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL write_section_term(this,group_id,equation_nn,term_recycling_parallel, &
         &section_bc,'recycling_parallel_inward')
    CALL write_section_term(this,group_id,equation_nn,term_recycling_diffusion, &
         &section_bc,'recycling_diffusion_inward')
    CALL write_section_term(this,group_id,equation_nn,term_recycling_pinch, &
         &section_bc,'recycling_pinch_inward')
    CALL write_section_term(this,group_id,equation_nn,term_puff,section_bc, &
         &'puff')
    CALL write_section_term(this,group_id,equation_nn,term_pump,section_bc, &
         &'pump')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_neutral_bc_sources

  SUBROUTINE write_neutral_bc_fluxes(this, parent_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('physical_flux_components_inward', &
         &parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL write_section_term(this,group_id,equation_nn,term_diffusion, &
         &section_bc,'limited_diffusion')
#ifdef NEUTRALP
    CALL write_section_term(this,group_id,equation_nn,term_pressure, &
         &section_bc,'limited_pressure')
#endif
#ifdef NEUTRALGAMMA
    CALL write_section_term(this,group_id,equation_nn,term_convection, &
         &section_bc,'neutral_gamma_convection')
#endif
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_neutral_bc_fluxes

  SUBROUTINE write_term(this, group_id, equation, term, name)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation, term
    CHARACTER(LEN=*), INTENT(IN) :: name

    CALL HDF5_real_saving(group_id,balance_value(this,equation,term, &
         &section_physical),name)
  END SUBROUTINE write_term

  SUBROUTINE write_section_term(this, group_id, equation, term, section, name)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: equation, term, section
    CHARACTER(LEN=*), INTENT(IN) :: name

    CALL HDF5_real_saving(group_id,balance_value(this,equation,term,section), &
         &name)
  END SUBROUTINE write_section_term

END SUBMODULE balance_diagnostics_output
