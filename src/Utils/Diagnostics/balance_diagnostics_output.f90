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

    CALL write_physical_hdf5(this,diagnostics_group_id)
    CALL write_discrete_hdf5(this,diagnostics_group_id)
    CALL write_bc_hdf5(this,diagnostics_group_id)
  END SUBROUTINE balance_diagnostics_write_hdf5

  SUBROUTINE write_physical_hdf5(this, diagnostics_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: physical_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('physical',diagnostics_group_id, &
         physical_group_id,ierr)
    CALL write_physical_equation_hdf5(this,physical_group_id,'n', &
         equation_n,'particles','particles/s')
    CALL write_physical_equation_hdf5(this,physical_group_id,'n_n', &
         equation_nn,'particles','particles/s')
    CALL write_total_particles_physical_hdf5(this,physical_group_id)
    CALL HDF5_group_close(physical_group_id,ierr)
  END SUBROUTINE write_physical_hdf5

  SUBROUTINE write_physical_equation_hdf5(this, parent_group_id, name, &
       &equation, content_units, rate_units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    CHARACTER(LEN=*), INTENT(IN) :: name, content_units, rate_units
    INTEGER, INTENT(IN) :: equation
    INTEGER(HID_T) :: equation_group_id
    TYPE(physical_balance_type) :: balance
    INTEGER :: ierr

    balance = physical_balance(this,equation)
    CALL HDF5_group_create(name,parent_group_id,equation_group_id,ierr)
    CALL write_physical_fields(equation_group_id,balance,content_units, &
         &rate_units,this%detailed())
    IF (this%detailed()) CALL write_physical_components_hdf5(this, &
         &equation_group_id,equation)
    CALL HDF5_group_close(equation_group_id,ierr)
  END SUBROUTINE write_physical_equation_hdf5

  SUBROUTINE write_total_particles_physical_hdf5(this, parent_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    TYPE(physical_balance_type) :: balance
    INTEGER :: ierr

    balance = add_physical_balances(physical_balance(this,equation_n), &
         physical_balance(this,equation_nn))
    CALL HDF5_group_create('total_n',parent_group_id,group_id,ierr)
    CALL write_physical_fields(group_id,balance,'particles','particles/s', &
         &this%detailed())
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_total_particles_physical_hdf5

  SUBROUTINE write_physical_fields(group_id, balance, content_units, &
       &rate_units, detailed)
    INTEGER(HID_T), INTENT(IN) :: group_id
    TYPE(physical_balance_type), INTENT(IN) :: balance
    CHARACTER(LEN=*), INTENT(IN) :: content_units, rate_units
    LOGICAL, INTENT(IN) :: detailed

    CALL HDF5_string_saving(group_id,content_units,'content_units')
    CALL HDF5_string_saving(group_id,rate_units,'rate_units')
    CALL HDF5_real_saving(group_id,balance%content,'content')
    CALL HDF5_real_saving(group_id,balance%physical_imbalance, &
         &'physical_imbalance')
    IF (.NOT. detailed) RETURN
    CALL HDF5_real_saving(group_id,balance%temporal,'temporal')
    CALL HDF5_real_saving(group_id,balance%volume,'volume')
    CALL HDF5_real_saving(group_id,balance%boundary_inward, &
         &'boundary_physical_inward')
  END SUBROUTINE write_physical_fields

  SUBROUTINE write_physical_components_hdf5(this, equation_group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER, INTENT(IN) :: equation

    CALL write_volume_components_hdf5(this,equation_group_id,equation)
    CALL write_boundary_components_hdf5(this,equation_group_id,equation)
    IF (equation == equation_n) CALL write_particle_exchange_hdf5(this, &
         &equation_group_id,equation)
  END SUBROUTINE write_physical_components_hdf5

  SUBROUTINE write_volume_components_hdf5(this, equation_group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER, INTENT(IN) :: equation
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('volume_components',equation_group_id, &
         group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation,term_ionization,section_physical), &
         'ionization')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation,term_recombination,section_physical), &
         'recombination')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation,term_prescribed_source,section_physical), &
         'prescribed_source')
    IF (equation == equation_nn) THEN
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_puff,section_physical), &
            'puff_source')
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_pump,section_physical), &
            'pump_source')
    ENDIF
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_volume_components_hdf5

  SUBROUTINE write_boundary_components_hdf5(this, equation_group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER, INTENT(IN) :: equation
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('boundary_components_inward',equation_group_id, &
         group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    IF (equation == equation_n) THEN
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_parallel,section_physical), &
            'parallel')
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_diffusion,section_physical), &
            'diffusion')
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_pinch,section_physical),'pinch')
    ELSE
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_diffusion,section_physical), &
            'limited_diffusion')
#ifdef NEUTRALP
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_pressure,section_physical), &
            'limited_pressure')
#endif
#ifdef NEUTRALGAMMA
       CALL HDF5_real_saving(group_id, &
            balance_value(this,equation,term_convection,section_physical), &
            'neutral_gamma_convection')
#endif
    ENDIF
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_boundary_components_hdf5

  SUBROUTINE write_particle_exchange_hdf5(this, equation_group_id, equation)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER, INTENT(IN) :: equation
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('exchange',equation_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation,term_charge_exchange,section_physical), &
         'charge_exchange_rate')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_particle_exchange_hdf5

  SUBROUTINE write_discrete_hdf5(this, diagnostics_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: discrete_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('discrete',diagnostics_group_id, &
         discrete_group_id,ierr)
    CALL write_discrete_equation_hdf5(this,discrete_group_id,'n',equation_n, &
         &'particles/s')
    CALL write_discrete_equation_hdf5(this,discrete_group_id,'n_n', &
         equation_nn,'particles/s')
    CALL write_total_particles_discrete_hdf5(this,discrete_group_id)
    CALL HDF5_group_close(discrete_group_id,ierr)
  END SUBROUTINE write_discrete_hdf5

  SUBROUTINE write_discrete_equation_hdf5(this, parent_group_id, name, &
       &equation, units)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    INTEGER, INTENT(IN) :: equation
    INTEGER(HID_T) :: group_id
    TYPE(physical_balance_type) :: balance
    INTEGER :: ierr

    balance = physical_balance(this,equation)
    CALL HDF5_group_create(name,parent_group_id,group_id,ierr)
    CALL write_discrete_fields(group_id,balance,units)
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_discrete_equation_hdf5

  SUBROUTINE write_total_particles_discrete_hdf5(this, parent_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    TYPE(physical_balance_type) :: balance
    INTEGER :: ierr

    balance = add_physical_balances(physical_balance(this,equation_n), &
         physical_balance(this,equation_nn))
    CALL HDF5_group_create('total_n',parent_group_id,group_id,ierr)
    CALL write_discrete_fields(group_id,balance,'particles/s')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_total_particles_discrete_hdf5

  SUBROUTINE write_discrete_fields(group_id, balance, units)
    INTEGER(HID_T), INTENT(IN) :: group_id
    TYPE(physical_balance_type), INTENT(IN) :: balance
    CHARACTER(LEN=*), INTENT(IN) :: units

    CALL HDF5_string_saving(group_id,units,'units')
    CALL HDF5_real_saving(group_id,balance%tau_inward,'hdg_tau_inward')
    CALL HDF5_real_saving(group_id,balance%discrete_residual,'residual')
  END SUBROUTINE write_discrete_fields

  SUBROUTINE write_bc_hdf5(this, diagnostics_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: bc_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('bc',diagnostics_group_id,bc_group_id,ierr)
    CALL write_density_bc_hdf5(this,bc_group_id)
    CALL write_neutral_density_bc_hdf5(this,bc_group_id)
    CALL HDF5_group_close(bc_group_id,ierr)
  END SUBROUTINE write_bc_hdf5

  SUBROUTINE write_density_bc_hdf5(this, parent_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    REAL*8 :: diffusion, tau
    INTEGER :: ierr

    diffusion = balance_value(this,equation_n,term_diffusion,section_bc)
    tau = balance_value(this,equation_n,term_tau_inward,section_bc)
    CALL HDF5_group_create('n',parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL HDF5_real_saving(group_id,diffusion,'diffusion_inward')
    CALL HDF5_real_saving(group_id,tau,'hdg_tau_inward')
    CALL HDF5_real_saving(group_id,diffusion+tau,'residual')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_density_bc_hdf5

  SUBROUTINE write_neutral_density_bc_hdf5(this, parent_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: group_id
    TYPE(particle_bc_balance_type) :: balance
    INTEGER :: ierr

    balance = neutral_density_bc(this)
    CALL HDF5_group_create('n_n',parent_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL HDF5_real_saving(group_id,balance%imposed_source_inward, &
         'imposed_source_inward')
    CALL HDF5_real_saving(group_id,balance%physical_flux_inward, &
         'physical_flux_inward')
    CALL HDF5_real_saving(group_id,balance%tau_inward,'hdg_tau_inward')
    CALL HDF5_real_saving(group_id,balance%residual,'residual')
    IF (this%detailed()) THEN
       CALL write_neutral_bc_sources_hdf5(this,group_id)
       CALL write_neutral_bc_fluxes_hdf5(this,group_id)
    ENDIF
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_neutral_density_bc_hdf5

  SUBROUTINE write_neutral_bc_sources_hdf5(this, equation_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('source_components',equation_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_recycling_parallel,section_bc), &
         'recycling_parallel_inward')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_recycling_diffusion,section_bc), &
         'recycling_diffusion_inward')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_recycling_pinch,section_bc), &
         'recycling_pinch_inward')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_puff,section_bc), &
         'puff_source')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_pump,section_bc), &
         'pump_sink')
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_neutral_bc_sources_hdf5

  SUBROUTINE write_neutral_bc_fluxes_hdf5(this, equation_group_id)
    CLASS(balance_diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: equation_group_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create('physical_flux_components_inward', &
         equation_group_id,group_id,ierr)
    CALL HDF5_string_saving(group_id,'particles/s','units')
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_diffusion,section_bc), &
         'limited_diffusion')
#ifdef NEUTRALP
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_pressure,section_bc), &
         'limited_pressure')
#endif
#ifdef NEUTRALGAMMA
    CALL HDF5_real_saving(group_id, &
         balance_value(this,equation_nn,term_convection,section_bc), &
         'neutral_gamma_convection')
#endif
    CALL HDF5_group_close(group_id,ierr)
  END SUBROUTINE write_neutral_bc_fluxes_hdf5

END SUBMODULE balance_diagnostics_output
