SUBMODULE (transport_models_1d) transport_models_1d_io
IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE tm1d_write_hdf5(this, parent_group_id)
    CLASS(transport_model_1d_t), INTENT(IN) :: this
    INTEGER(HID_T), INTENT(IN) :: parent_group_id
    INTEGER(HID_T) :: coeffs_group_id, params_group_id, diagnostics_group_id
    INTEGER :: ierr
    INTEGER :: pinch_equations(4), band_equations(4)
    REAL*8 :: band_rho(5), band_start, plateau_start, plateau_end, band_end, center

    IF (.NOT. switch%transport_1d) RETURN
    IF (.NOT. this%is_initialized) RETURN
    IF (this%nrho <= 0) RETURN

    CALL HDF5_group_create('coefficients', parent_group_id, coeffs_group_id, ierr)
    CALL HDF5_array1D_saving(coeffs_group_id, this%chi_i_fs, SIZE(this%chi_i_fs), 'chi_i_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%chi_e_fs, SIZE(this%chi_e_fs), 'chi_e_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%d_fs, SIZE(this%d_fs), 'd_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%nu_mom_fs, SIZE(this%nu_mom_fs), 'nu_mom_fs')
    CALL HDF5_array1D_saving(coeffs_group_id, this%vpinch_fs, SIZE(this%vpinch_fs), 'vpinch_fs')
    CALL HDF5_group_close(coeffs_group_id, ierr)

    CALL HDF5_group_create('params', parent_group_id, params_group_id, ierr)
    CALL HDF5_real_saving(params_group_id, this%config%transport_1d_relaxation, 'transport_1d_relaxation')
    CALL HDF5_real_saving(params_group_id, this%config%rho_core, 'rho_core')
    CALL HDF5_real_saving(params_group_id, this%config%rho_edge, 'rho_edge')
    CALL HDF5_real_saving(params_group_id, this%config%rho_diffusion_model_max, 'rho_diffusion_model_max')
    CALL HDF5_string_saving(params_group_id, &
         TRIM(tm1d_region_policy_name(this%config%transport_region_policy)), &
         'transport_region_policy')
    CALL HDF5_integer_saving(params_group_id, this%config%transport_region_policy, &
         'transport_region_policy_id')
    CALL HDF5_real_saving(params_group_id, this%config%rho_blend_width, 'rho_blend_width')
    CALL HDF5_real_saving(params_group_id, this%config%diff_n_min, 'diff_n_min')
    CALL HDF5_real_saving(params_group_id, this%config%diff_u_min, 'diff_u_min')
    CALL HDF5_real_saving(params_group_id, this%config%diff_e_min, 'diff_e_min')
    CALL HDF5_real_saving(params_group_id, this%config%diff_ee_min, 'diff_ee_min')
    CALL HDF5_real_saving(params_group_id, this%config%c_bohm_i, 'c_bohm_i')
    CALL HDF5_real_saving(params_group_id, this%config%c_gyrobohm_i, 'c_gyrobohm_i')
    CALL HDF5_real_saving(params_group_id, this%config%c_bohm_e, 'c_bohm_e')
    CALL HDF5_real_saving(params_group_id, this%config%c_gyrobohm_e, 'c_gyrobohm_e')
    CALL HDF5_real_saving(params_group_id, this%config%c_bohm_n, 'c_bohm_n')
    CALL HDF5_real_saving(params_group_id, this%config%c_bohm_n_rho_slope, &
         'c_bohm_n_rho_slope')
    CALL HDF5_real_saving(params_group_id, this%config%prandtl, 'prandtl')
    CALL HDF5_integer_saving(params_group_id, this%config%pinch_model, 'pinch_model')
    pinch_equations = PACK([1, 2, 3, 4], this%config%pinch_enabled, [0, 0, 0, 0])
    CALL HDF5_array1D_saving_int(params_group_id, pinch_equations, SIZE(pinch_equations), 'pinch_equations')
    CALL HDF5_real_saving(params_group_id, this%config%c_pinch, 'c_pinch')
    CALL HDF5_real_saving(params_group_id, this%config%nu_th, 'nu_th')
    CALL HDF5_real_saving(params_group_id, this%config%vpinch_const, 'vpinch_const')
    CALL HDF5_real_saving(params_group_id, this%config%rho_pinch_axis_width, 'rho_pinch_axis_width')
    CALL HDF5_real_saving(params_group_id, this%config%rho_pinch_model_max, 'rho_pinch_model_max')
    CALL HDF5_real_saving(params_group_id, this%config%rho_pinch_edge_width, 'rho_pinch_edge_width')
    IF (ANY(this%config%band_enabled)) THEN
       band_equations = PACK([1, 2, 3, 4], this%config%band_enabled, [0, 0, 0, 0])
       CALL tm1d_band_edges(this%config, band_start, plateau_start, plateau_end, band_end)
       CALL HDF5_array1D_saving_int(params_group_id, band_equations, 4, 'band_equations')
       CALL HDF5_array1D_saving(params_group_id, this%config%band_diff, 4, 'band_diff')
       CALL HDF5_string_saving(params_group_id, TRIM(this%config%diffusion_coordinate), 'band_coordinate')
       CALL HDF5_real_saving(params_group_id, this%config%band_center, 'band_center')
       CALL HDF5_real_saving(params_group_id, this%config%band_width, 'band_width')
       CALL HDF5_real_saving(params_group_id, plateau_start-band_start, 'band_transition_width')
       CALL this%band_edges(band_start, plateau_start, plateau_end, band_end, center)
       band_rho = [band_start,plateau_start,center,plateau_end,band_end]
       CALL HDF5_array1D_saving(params_group_id, band_rho, 5, 'band_rho')
    ENDIF
    CALL HDF5_group_close(params_group_id, ierr)

    CALL HDF5_group_create('diagnostics', parent_group_id, diagnostics_group_id, ierr)
    CALL HDF5_array1D_saving(diagnostics_group_id, this%profile_relative_defect, 5, 'coefficient_relative_defect')
    CALL HDF5_string_saving(diagnostics_group_id, 'D,nu_mom,chi_i,chi_e,vpinch', 'coefficient_order')
    CALL HDF5_group_close(diagnostics_group_id, ierr)
  END SUBROUTINE tm1d_write_hdf5

END SUBMODULE transport_models_1d_io
