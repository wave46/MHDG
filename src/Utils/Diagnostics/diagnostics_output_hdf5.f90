SUBMODULE (diagnostics) diagnostics_output_hdf5
  USE HDF5_io_module
  USE GLOBALS, ONLY: phys
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE diag_write_hdf5(this, file_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: file_id

    CALL diag_prepare_wall_source_nodal_output()

    IF (MPIvar%glob_id .EQ. 0) CALL diag_write_root_hdf5(this, file_id)
  END SUBROUTINE diag_write_hdf5

  SUBROUTINE diag_write_root_hdf5(this, file_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: diagnostics_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('diagnostics', file_id, diagnostics_group_id, ierr)
    CALL diag_write_boundary_hdf5(this, diagnostics_group_id)
    CALL diag_write_balance_hdf5(this, diagnostics_group_id)
    CALL diag_write_content_hdf5(this, diagnostics_group_id)
    CALL HDF5_group_close(diagnostics_group_id, ierr)
  END SUBROUTINE diag_write_root_hdf5

  SUBROUTINE diag_write_balance_hdf5(this, diagnostics_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: balance_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('balance', diagnostics_group_id, balance_group_id, ierr)
    CALL diag_write_plasma_particles_hdf5(this, balance_group_id)
    CALL diag_write_neutral_particles_hdf5(this, balance_group_id)
    CALL diag_write_total_particles_hdf5(this, balance_group_id)
    CALL HDF5_group_close(balance_group_id, ierr)
  END SUBROUTINE diag_write_balance_hdf5

  SUBROUTINE diag_write_boundary_hdf5(this, diagnostics_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: boundary_group_id, neutrals_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('boundary_conditions', diagnostics_group_id, boundary_group_id, ierr)
    CALL HDF5_group_create('neutrals', boundary_group_id, neutrals_group_id, ierr)
    CALL diag_write_boundary_terms_hdf5(this, neutrals_group_id)
    CALL diag_write_boundary_summary_hdf5(this, neutrals_group_id)
    CALL HDF5_group_close(neutrals_group_id, ierr)
    CALL HDF5_group_close(boundary_group_id, ierr)
  END SUBROUTINE diag_write_boundary_hdf5

  SUBROUTINE diag_write_boundary_terms_hdf5(this, neutrals_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: neutrals_group_id
    INTEGER(HID_T) :: terms_group_id
    INTEGER :: ierr
    TYPE(diag_boundary_summary_type) :: summary

    summary = diag_boundary_summary(this)

    CALL HDF5_group_create('terms', neutrals_group_id, terms_group_id, ierr)
    CALL HDF5_string_saving(terms_group_id, diag_particle_integral_units(), 'units')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_plasma_parallel_flux), &
       &'recycled_plasma_parallel_source')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_plasma_diffusion_flux), &
       &'recycled_plasma_diffusion_source')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_plasma_pinch_flux), &
       &'recycled_plasma_pinch_source')
    CALL HDF5_real_saving(terms_group_id, summary%recycled_source, 'recycling_source')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_neutral_diffusion_flux), &
       &'neutral_diffusion_boundary_flux')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_neutral_pressure_flux), &
       &'neutral_pressure_boundary_flux')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_neutral_convection_flux), &
       &'neutral_convection_boundary_flux')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_neutral_total_flux), &
       &'neutral_total_boundary_flux')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_tau_numerical_flux), &
       &'tau_numerical_boundary_flux')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_wall_puff_source), 'wall_puff_source')
    CALL HDF5_real_saving(terms_group_id, this%boundary_hdg(diag_boundary_wall_pump_sink), 'wall_pump_sink')
    CALL HDF5_group_close(terms_group_id, ierr)
  END SUBROUTINE diag_write_boundary_terms_hdf5

  SUBROUTINE diag_write_boundary_summary_hdf5(this, neutrals_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: neutrals_group_id
    INTEGER(HID_T) :: summary_group_id
    INTEGER :: ierr
    TYPE(diag_boundary_summary_type) :: summary

    summary = diag_boundary_summary(this)

    CALL HDF5_group_create('summary', neutrals_group_id, summary_group_id, ierr)
    CALL HDF5_string_saving(summary_group_id, diag_particle_integral_units(), 'units')
    CALL HDF5_real_saving(summary_group_id, summary%residual, 'neutral_closure_residual')
    CALL HDF5_group_close(summary_group_id, ierr)
  END SUBROUTINE diag_write_boundary_summary_hdf5

  SUBROUTINE diag_write_plasma_particles_hdf5(this, balance_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: balance_group_id
    INTEGER(HID_T) :: plasma_group_id, particles_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('plasma', balance_group_id, plasma_group_id, ierr)
    CALL HDF5_group_create('particles', plasma_group_id, particles_group_id, ierr)
    CALL diag_write_plasma_particle_terms_hdf5(this, particles_group_id)
    CALL diag_write_plasma_particle_summary_hdf5(this, particles_group_id)
    CALL HDF5_group_close(particles_group_id, ierr)
    CALL HDF5_group_close(plasma_group_id, ierr)
  END SUBROUTINE diag_write_plasma_particles_hdf5

  SUBROUTINE diag_write_plasma_particle_terms_hdf5(this, particles_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: particles_group_id
    INTEGER(HID_T) :: terms_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('terms', particles_group_id, terms_group_id, ierr)
    CALL HDF5_string_saving(terms_group_id, diag_particle_integral_units(), 'units')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_plasma_ionization), 'ionization')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_plasma_recombination), 'recombination')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_plasma_boundary_flux), 'boundary_flux')
    CALL HDF5_group_close(terms_group_id, ierr)
  END SUBROUTINE diag_write_plasma_particle_terms_hdf5

  SUBROUTINE diag_write_plasma_particle_summary_hdf5(this, particles_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: particles_group_id
    INTEGER(HID_T) :: summary_group_id
    INTEGER :: ierr
    TYPE(diag_particle_summary_type) :: summary

    summary = diag_particle_summary(this)

    CALL HDF5_group_create('summary', particles_group_id, summary_group_id, ierr)
    CALL HDF5_string_saving(summary_group_id, diag_particle_integral_units(), 'units')
    CALL HDF5_real_saving(summary_group_id, summary%plasma_balance, 'particle_balance')
    CALL HDF5_group_close(summary_group_id, ierr)
  END SUBROUTINE diag_write_plasma_particle_summary_hdf5

  SUBROUTINE diag_write_neutral_particles_hdf5(this, balance_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: balance_group_id
    INTEGER(HID_T) :: neutrals_group_id, particles_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('neutrals', balance_group_id, neutrals_group_id, ierr)
    CALL HDF5_group_create('particles', neutrals_group_id, particles_group_id, ierr)
    CALL diag_write_neutral_particle_terms_hdf5(this, particles_group_id)
    CALL diag_write_neutral_particle_summary_hdf5(this, particles_group_id)
    CALL diag_write_wall_source_nodal(particles_group_id)
    CALL HDF5_group_close(particles_group_id, ierr)
    CALL HDF5_group_close(neutrals_group_id, ierr)
  END SUBROUTINE diag_write_neutral_particles_hdf5

  SUBROUTINE diag_write_neutral_particle_terms_hdf5(this, particles_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: particles_group_id
    INTEGER(HID_T) :: terms_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('terms', particles_group_id, terms_group_id, ierr)
    CALL HDF5_string_saving(terms_group_id, diag_particle_integral_units(), 'units')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_neutral_ionization), 'ionization')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_neutral_recombination), 'recombination')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_neutral_boundary_flux), 'boundary_flux')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_puff_source), 'puff_source')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_pump_sink), 'pump_sink')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_recycling_source), 'recycling_source')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_recycling_parallel_source), &
       &'recycling_parallel_source')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_recycling_diffusion_source), &
       &'recycling_diffusion_source')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_recycling_pinch_source), &
       &'recycling_pinch_source')
    CALL HDF5_real_saving(terms_group_id, this%particles(diag_particle_charge_exchange_rate), 'charge_exchange_rate')
    CALL HDF5_group_close(terms_group_id, ierr)
  END SUBROUTINE diag_write_neutral_particle_terms_hdf5

  SUBROUTINE diag_write_neutral_particle_summary_hdf5(this, particles_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: particles_group_id
    INTEGER(HID_T) :: summary_group_id
    INTEGER :: ierr
    TYPE(diag_particle_summary_type) :: summary

    summary = diag_particle_summary(this)

    CALL HDF5_group_create('summary', particles_group_id, summary_group_id, ierr)
    CALL HDF5_string_saving(summary_group_id, diag_particle_integral_units(), 'units')
    CALL HDF5_real_saving(summary_group_id, summary%neutral_balance, 'particle_balance')
    CALL HDF5_group_close(summary_group_id, ierr)
  END SUBROUTINE diag_write_neutral_particle_summary_hdf5

  SUBROUTINE diag_write_total_particles_hdf5(this, balance_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: balance_group_id
    INTEGER(HID_T) :: total_group_id, particles_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('total', balance_group_id, total_group_id, ierr)
    CALL HDF5_group_create('particles', total_group_id, particles_group_id, ierr)
    CALL diag_write_total_particle_summary_hdf5(this, particles_group_id)
    CALL HDF5_group_close(particles_group_id, ierr)
    CALL HDF5_group_close(total_group_id, ierr)
  END SUBROUTINE diag_write_total_particles_hdf5

  SUBROUTINE diag_write_total_particle_summary_hdf5(this, particles_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: particles_group_id
    INTEGER(HID_T) :: summary_group_id
    INTEGER :: ierr
    TYPE(diag_particle_summary_type) :: summary

    summary = diag_particle_summary(this)

    CALL HDF5_group_create('summary', particles_group_id, summary_group_id, ierr)
    CALL HDF5_string_saving(summary_group_id, diag_particle_integral_units(), 'units')
    CALL HDF5_real_saving(summary_group_id, summary%total_balance, 'particle_balance')
    CALL HDF5_group_close(summary_group_id, ierr)
  END SUBROUTINE diag_write_total_particle_summary_hdf5

  SUBROUTINE diag_write_content_hdf5(this, diagnostics_group_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: diagnostics_group_id
    INTEGER(HID_T) :: content_group_id
    INTEGER :: ierr

    CALL HDF5_group_create('content', diagnostics_group_id, content_group_id, ierr)
    CALL diag_write_content_entry_hdf5(this, content_group_id, 'plasma', diag_content_plasma_particles)
    CALL diag_write_content_entry_hdf5(this, content_group_id, 'neutrals', diag_content_neutral_particles)
    CALL diag_write_content_entry_hdf5(this, content_group_id, 'total', diag_content_total_particles)
    CALL HDF5_group_close(content_group_id, ierr)
  END SUBROUTINE diag_write_content_hdf5

  SUBROUTINE diag_write_content_entry_hdf5(this, content_group_id, group_name, term_id)
    CLASS(diagnostics_type), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: content_group_id
    CHARACTER(LEN=*), INTENT(IN) :: group_name
    INTEGER, INTENT(IN) :: term_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    CALL HDF5_group_create(group_name, content_group_id, group_id, ierr)
    CALL HDF5_string_saving(group_id, diag_content_units(term_id), 'units')
    CALL HDF5_real_saving(group_id, diag_content_value(this, term_id), 'particles')
    CALL HDF5_group_close(group_id, ierr)
  END SUBROUTINE diag_write_content_entry_hdf5

  MODULE SUBROUTINE diag_write_wall_source_nodal(group_id)
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER(HID_T) :: nodal_group_id
    INTEGER :: ierr

    IF (.NOT. diag_wall_source_nodal_ready) RETURN

    CALL HDF5_group_create('nodal_wall_sources', group_id, nodal_group_id, ierr)
    CALL HDF5_string_saving(nodal_group_id, diag_particle_integral_units(), 'total_units')
    CALL HDF5_string_saving(nodal_group_id, 'particles/(m^3 s)', 'nodal_units')
    CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_puff_total, 'element_puff_total')
    CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_pump_total, 'element_pump_total')
    CALL HDF5_real_saving(nodal_group_id, phys%neutral_wall_source_puff_total - phys%neutral_wall_source_pump_total, &
       &'element_net_total')
    CALL HDF5_array1D_saving(nodal_group_id, diag_wall_source_puff_nodal, &
       &SIZE(diag_wall_source_puff_nodal), 'puff_flux_density')
    CALL HDF5_array1D_saving(nodal_group_id, diag_wall_source_pump_nodal, &
       &SIZE(diag_wall_source_pump_nodal), 'pump_flux_density')
    CALL HDF5_array1D_saving(nodal_group_id, diag_wall_source_net_nodal, &
       &SIZE(diag_wall_source_net_nodal), 'net_flux_density')
    CALL HDF5_group_close(nodal_group_id, ierr)
  END SUBROUTINE diag_write_wall_source_nodal

END SUBMODULE diagnostics_output_hdf5
