PROGRAM test_transport_band
  USE input_reader, ONLY: read_transport_model_input
  USE globals, ONLY: input, transport_model_input, simpar, switch
  USE MPI_OMP, ONLY: MPIvar
  USE transport_models_1d, ONLY: transport_model_1d_t
  USE magnetic_topology, ONLY: magnetic_region_core, magnetic_region_main_sol, &
       magnetic_region_private_flux, magnetic_region_undefined
  USE magnetic_geometry_state, ONLY: magnetic_geometry_cache
  USE HDF5_io_module, ONLY: HID_T, HDF5_create, HDF5_close, &
       HDF5_group_create, HDF5_group_close, HDF5_array2D_saving
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan, ieee_positive_inf
  IMPLICIT NONE

  TYPE(transport_model_1d_t) :: model, background
  REAL*8 :: iso(6,6,9), ani(6,6,9), base_iso(6,6,9), base_ani(6,6,9)
  REAL*8 :: saved(6,6,9), used(9,6), rho(9), center, half_width, transition, target
  REAL*8, PARAMETER :: weights(9) = [0.d0,.15625d0,1.d0,1.d0,1.d0,.84375d0,0.d0,1.d0,0.d0]
  INTEGER :: region(9), equation, g
  INTEGER(HID_T) :: file_id, group_id
  CHARACTER(LEN=32) :: argument

  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  simpar%refval_time = 2.d0
  simpar%refval_length = 4.d0
  switch%transport_1d = .TRUE.
  CALL get_command_argument(1,argument)
  SELECT CASE (TRIM(argument))
  CASE ('negative')
     CALL model%set_config(band_equations=[-1])
  CASE ('neutral')
     CALL model%set_config(band_equations=[5])
  CASE ('overlong')
     CALL model%set_config(band_equations=[1,2,3,4,0])
  CASE ('nan')
     CALL model%set_config(band_equations=[1], band_width=.004d0, &
          band_center=ieee_value(0.d0,ieee_quiet_nan))
  CASE ('infinity')
     CALL model%set_config(band_equations=[1], band_width=.004d0, &
          band_diff_phys=[ieee_value(0.d0,ieee_positive_inf),0.d0,0.d0,0.d0])
  CASE ('no_geometry')
     CALL model%set_config(band_equations=[1], band_width=.004d0)
     CALL model%init(2)
     rho = 1.d0
     region = magnetic_region_core
     iso = 1.d0
     ani = 1.d0
     CALL model%apply_1D_diffusion(rho,iso,ani,region)
  END SELECT
  IF (LEN_TRIM(argument) > 0) ERROR STOP 'Invalid direct configuration accepted'

  input%transport_model_path = 'transport_model.nml'
  CALL read_transport_model_input()
  CALL model%set_config(band_equations=transport_model_input%band_equations, &
       band_diff_phys=transport_model_input%band_diff_phys, &
       band_center=transport_model_input%band_center, band_width=transport_model_input%band_width, &
       band_transition_width=transport_model_input%band_transition_width, &
       transport_region_policy='core_only', rho_diffusion_model_max=.9d0, &
       diff_n_min_phys=80.d0, diff_u_min_phys=80.d0, &
       diff_e_min_phys=80.d0, diff_ee_min_phys=80.d0)
  CALL initialize_profiles(model)
  background = model
  CALL background%set_config(band_equations=[0])
  region = [magnetic_region_core,magnetic_region_core,magnetic_region_core, &
       magnetic_region_core,magnetic_region_main_sol,magnetic_region_main_sol, &
       magnetic_region_private_flux,magnetic_region_private_flux,magnetic_region_undefined]
  CALL sample_band()

  CALL HDF5_create('transport_band.h5',file_id)
  DO equation = 1,6
     used(:,equation) = ani(equation,equation,:)
  ENDDO
  CALL HDF5_array2D_saving(file_id,used,9,6,'used')
  CALL HDF5_group_create('transport_1d',file_id,group_id)
  CALL model%write_hdf5(group_id)
  CALL HDF5_group_close(group_id)
  CALL HDF5_close(file_id)

  saved = ani
  ! Repeated input transfer must normalize from physical inputs, not prior state.
  CALL model%set_config(band_diff_phys=transport_model_input%band_diff_phys)
  CALL model%init(5)
  CALL initialize_profiles(model)
  CALL sample_band()
  IF (ANY(ani /= saved)) ERROR STOP 'Band changed on reinitialization or repeated normalization'

  IF (ANY(model%config%band_enabled)) THEN
     ! Below the BGB cutoff, the band overrides modeled floors and blending too.
     CALL model%set_config(band_center=.75d0, rho_diffusion_model_max=1.d0, rho_blend_width=.5d0)
     CALL background%set_config(rho_diffusion_model_max=1.d0, rho_blend_width=.5d0)
     CALL sample_band()
     ! The same analytic band also supports a prescribed increase and zero.
     CALL model%set_config(band_equations=[1,2], band_diff_phys=[800.d0,0.d0,0.d0,0.d0])
     CALL sample_band()
     ! Automatic shoulders track a new plateau width; explicit widths persist.
     CALL model%set_config(band_width=.008d0)
     CALL sample_band()
  ENDIF

  CALL model%destroy()
  IF (ANY(model%config%band_enabled)) ERROR STOP 'Destroy must reset band selection'
  magnetic_geometry_cache%is_initialized = .FALSE.
  WRITE (*,'(A)') 'transport band checks: PASS'

CONTAINS

  SUBROUTINE initialize_profiles(this)
    TYPE(transport_model_1d_t), INTENT(INOUT) :: this
    INTEGER :: i
    IF (.NOT. this%is_initialized) CALL this%init(3)
    this%rho_grid = [(1.4d0*(i-1)/(this%nrho-1),i=1,this%nrho)]
    this%d_fs = 6.d0
    this%nu_mom_fs = 7.d0
    this%chi_i_fs = 8.d0
    this%chi_e_fs = 9.d0
    this%vpinch_fs = -3.d0
  END SUBROUTINE initialize_profiles

  SUBROUTINE sample_band()
    REAL*8 :: expected
    center = model%config%band_center
    half_width = .5d0*model%config%band_width
    transition = model%config%band_transition_width
    IF (transition == -1.d0) transition = .5d0*half_width
    ! Quarter-shoulder probes distinguish cubic smoothstep from a linear ramp.
    rho = [center-half_width-transition,center-half_width-.75d0*transition, &
         center-half_width,center,center+half_width,center+half_width+.25d0*transition, &
         center+half_width+transition,center,center]
    iso = .125d0
    ani = .125d0
    DO g = 1,9
       DO equation = 1,6
          ani(equation,equation,g) = equation + .125d0*g
          iso(equation,equation,g) = ani(equation,equation,g) + .375d0
       ENDDO
    ENDDO
    base_iso = iso
    base_ani = ani
    ! Disabled band must work without geometry and preserve background exactly.
    magnetic_geometry_cache%is_initialized = .FALSE.
    CALL background%apply_1D_diffusion(rho,base_iso,base_ani,region)
    magnetic_geometry_cache%is_initialized = ANY(model%config%band_enabled)
    CALL model%apply_1D_diffusion(rho,iso,ani,region)
    DO g = 1,9
       DO equation = 1,6
          expected = base_ani(equation,equation,g)
          IF (equation <= 4) THEN
             IF (model%config%band_enabled(equation)) THEN
                target = model%config%band_diff(equation)
                expected = (1.d0-weights(g))*expected + weights(g)*target
             ENDIF
          ENDIF
          CALL check_close(ani(equation,equation,g),expected,'band target/background blend')
       ENDDO
    ENDDO
    CALL check_array_close(iso-ani,base_iso-base_ani,'isotropic contribution changed')
    IF (ANY(iso(5:,:,:) /= base_iso(5:,:,:)) .OR. ANY(ani(5:,:,:) /= base_ani(5:,:,:))) &
         ERROR STOP 'Non-plasma rows changed'
    DO equation = 1,4
       IF (model%config%band_enabled(equation)) CYCLE
       IF (ANY(iso(equation,:,:) /= base_iso(equation,:,:)) .OR. &
            ANY(ani(equation,:,:) /= base_ani(equation,:,:))) ERROR STOP 'Inactive equation changed'
    ENDDO
    IF (ANY(iso(:,:,9) /= base_iso(:,:,9)) .OR. ANY(ani(:,:,9) /= base_ani(:,:,9))) &
         ERROR STOP 'Undefined topology changed'
    IF (ANY(model%d_fs /= 6.d0) .OR. ANY(model%vpinch_fs /= -3.d0)) &
         ERROR STOP 'Raw transport or pinch profile changed'
  END SUBROUTINE sample_band

  SUBROUTINE check_close(actual,expected,label)
    REAL*8, INTENT(IN) :: actual, expected
    CHARACTER(LEN=*), INTENT(IN) :: label
    IF (ABS(actual-expected) > 1.d-10*MAX(1.d0,ABS(expected))) THEN
       WRITE (*,*) 'FAIL: ',label,actual,expected
       ERROR STOP 1
    ENDIF
  END SUBROUTINE check_close

  SUBROUTINE check_array_close(actual,expected,label)
    REAL*8, INTENT(IN) :: actual(:,:,:), expected(:,:,:)
    CHARACTER(LEN=*), INTENT(IN) :: label
    IF (MAXVAL(ABS(actual-expected)) > 1.d-12) THEN
       WRITE (*,*) 'FAIL: ',label
       ERROR STOP 1
    ENDIF
  END SUBROUTINE check_array_close

END PROGRAM test_transport_band
