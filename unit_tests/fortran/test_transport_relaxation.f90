PROGRAM test_transport_relaxation
  USE input_reader, ONLY: read_transport_model_input
  USE globals, ONLY: input,transport_model_input,simpar,switch,phys,sol,time
  USE MPI_OMP, ONLY: MPIvar
  USE transport_models_1d, ONLY: transport_model_1d_t,transport_model_1d
  USE flux_surface_transport_data, ONLY: flux_surface_transport_t
  USE magnetic_geometry_state, ONLY: magnetic_geometry_cache
  USE Main_utils, ONLY: accepted_step_state_t
  USE HDF5_io_module, ONLY: HID_T,HDF5_create,HDF5_close,HDF5_group_create,HDF5_group_close
  IMPLICIT NONE
  TYPE(transport_model_1d_t) :: fresh,explicit_one,relaxed,independent,output_model
  TYPE(flux_surface_transport_t) :: fs
  TYPE(accepted_step_state_t) :: saved
  REAL*8 :: old(5,3),proposal(5,3),expected(5,3),alpha
  INTEGER :: j,ierr
  INTEGER(HID_T) :: file_id,group_id
  CHARACTER(LEN=40) :: argument

  MPIvar%glob_id=0;MPIvar%glob_size=1
  simpar%refval_time=1.d0;simpar%refval_length=.1d0
  simpar%refval_temperature=1000.d0;simpar%refval_density=1.d19
  switch%transport_1d=.TRUE.
  phys%neq=5;phys%npv=11;phys%idx_rhon_eq=5;phys%idx_rhon_pv=11
  phys%Mref=1.d0;phys%a_minor=10.d0;phys%Zeff=1.d0
  CALL get_command_argument(1,argument)
  IF(LEN_TRIM(argument)>0)THEN
     READ(argument,*,IOSTAT=ierr)alpha
     IF(ierr/=0)ERROR STOP 'Invalid test argument'
     CALL relaxed%set_config(transport_1d_relaxation=alpha)
     ERROR STOP 'Invalid direct relaxation accepted'
  ENDIF
  input%transport_model_path='transport_model.nml'
  CALL read_transport_model_input()
  alpha=transport_model_input%transport_1d_relaxation

  CALL fs%init(5,1.d0,.5d0)
  fs%q_fs=2.d0;fs%omega_fs=100.d0;fs%Rmaj_fs=40.d0;fs%shell_weight=1.d0
  CALL setup(fresh,1.d0);CALL setup(explicit_one,1.d0)
  CALL setup(relaxed,.25d0);CALL setup(independent,.5d0)
  CALL fill_profiles(1.d0)
  CALL fresh%update_from_flux_surfaces(fs)
  CALL relaxed%update_from_flux_surfaces(fs)
  old=coefficients(fresh)
  CALL check_equal(coefficients(relaxed),old,'first update must use full proposal')
  IF(ANY(relaxed%profile_relative_defect/=0.d0))ERROR STOP 'Initial defect must be zero'

  ! Calling set_config every NR with the same alpha must preserve history.
  CALL relaxed%set_config(transport_1d_relaxation=.25d0)
  CALL fill_profiles(1.3d0)
  CALL fresh%update_from_flux_surfaces(fs)
  CALL explicit_one%update_from_flux_surfaces(fs)
  proposal=coefficients(fresh)
  CALL check_equal(coefficients(explicit_one),proposal,'default/full-refresh compatibility')
  CALL relaxed%update_from_flux_surfaces(fs)
  expected=.75d0*old+.25d0*proposal
  CALL check_close(coefficients(relaxed),expected,'all five raw profiles must be mixed')
  IF(.NOT.ANY(relaxed%profile_relative_defect>0.d0))ERROR STOP 'Missing coefficient defect diagnostic'
  IF(ANY(relaxed%vpinch_fs>0.d0))ERROR STOP 'Inward pinch sign changed'
  CALL independent%update_from_flux_surfaces(fs)
  CALL check_equal(coefficients(independent),proposal,'different objects must not share history')

  DO j=1,80
     CALL relaxed%update_from_flux_surfaces(fs)
  ENDDO
  CALL check_close(coefficients(relaxed),proposal,'fixed point must equal unrelaxed model')

  CALL relaxed%reset_profile_history()
  CALL fill_profiles(.8d0)
  CALL fresh%update_from_flux_surfaces(fs);CALL relaxed%update_from_flux_surfaces(fs)
  CALL check_equal(coefficients(relaxed),coefficients(fresh),'explicit history reset')
  fs%rho_grid(2)=.55d0
  CALL fill_profiles(1.1d0)
  CALL fresh%update_from_flux_surfaces(fs);CALL relaxed%update_from_flux_surfaces(fs)
  CALL check_equal(coefficients(relaxed),coefficients(fresh),'same-sized radial regrid resets history')
  magnetic_geometry_cache%generation=magnetic_geometry_cache%generation+1
  CALL fill_profiles(1.4d0)
  CALL fresh%update_from_flux_surfaces(fs);CALL relaxed%update_from_flux_surfaces(fs)
  CALL check_equal(coefficients(relaxed),coefficients(fresh),'new geometry generation resets history')
  CALL relaxed%init(3)
  CALL relaxed%update_from_flux_surfaces(fs)
  CALL check_equal(coefficients(relaxed),coefficients(fresh),'reinitialization resets history')
  CALL relaxed%set_config(transport_1d_relaxation=.5d0)
  IF(relaxed%profile_history_ready)ERROR STOP 'Changed alpha must reset history'
  CALL relaxed%update_from_flux_surfaces(fs)
  CALL relaxed%set_config(transport_region_policy='core_only')
  IF(relaxed%profile_history_ready)ERROR STOP 'Changed region policy must reset history'
  CALL relaxed%update_from_flux_surfaces(fs)

  ! Exercise the real predictor/timestep rollback plasma-restore hook.
  transport_model_1d=relaxed
  ALLOCATE(sol%u(5),sol%q(10),sol%u_tilde(5),sol%u0(5,1))
  sol%u=1.d0;sol%q=0.d0;sol%u_tilde=1.d0;sol%u0=1.d0;sol%Nt=0
  time%bdf_start_step=0
  CALL saved%capture()
  CALL fill_profiles(1.6d0)
  CALL transport_model_1d%update_from_flux_surfaces(fs)
  CALL saved%restore_plasma()
  IF(transport_model_1d%profile_history_ready)ERROR STOP 'Rollback retained rejected transport history'
  CALL transport_model_1d%update_from_flux_surfaces(fs);CALL fresh%update_from_flux_surfaces(fs)
  CALL check_equal(coefficients(transport_model_1d),coefficients(fresh),'rollback must rebuild full proposal')

  CALL relaxed%destroy()
  IF(relaxed%profile_history_ready.OR.relaxed%config%transport_1d_relaxation/=1.d0) &
       ERROR STOP 'Destroy must reset history and default alpha'

  CALL setup(output_model,alpha)
  CALL output_model%update_from_flux_surfaces(fs)
  CALL fill_profiles(1.8d0)
  CALL output_model%update_from_flux_surfaces(fs)
  CALL HDF5_create('transport_relaxation.h5',file_id)
  CALL HDF5_group_create('transport_1d',file_id,group_id)
  CALL output_model%write_hdf5(group_id)
  CALL HDF5_group_close(group_id);CALL HDF5_close(file_id)
  WRITE(*,'(A)')'transport relaxation checks: PASS'
CONTAINS
  SUBROUTINE setup(model,value)
    TYPE(transport_model_1d_t),INTENT(INOUT)::model
    REAL*8,INTENT(IN)::value
    CALL model%set_config(transport_1d_relaxation=value,c_pinch=-.5d0, &
         c_bohm_i=.1d0,c_gyrobohm_i=.2d0,c_bohm_e=.05d0,c_gyrobohm_e=.4d0)
  END SUBROUTINE
  SUBROUTINE fill_profiles(factor)
    REAL*8,INTENT(IN)::factor
    REAL*8::rho,n,te,ti
    INTEGER::k
    DO k=1,fs%nrho
       rho=fs%rho_grid(k);n=1.d0+.2d0*rho;te=factor*(5.d0-3.d0*rho);ti=4.d0-rho
       fs%U_fs(:,k)=[n,.1d0*n,n*(1.5d0*ti+.005d0),1.5d0*n*te,.01d0]
       fs%Q_rad_fs(:,k)=0.d0
       fs%Q_rad_fs(1,k)=.2d0
       fs%Q_rad_fs(4,k)=1.5d0*(.2d0*te-3.d0*factor*n)
    ENDDO
    fs%profiles_built=.TRUE.
  END SUBROUTINE
  FUNCTION coefficients(model)RESULT(c)
    TYPE(transport_model_1d_t),INTENT(IN)::model
    REAL*8::c(5,model%nrho)
    c(1,:)=model%d_fs;c(2,:)=model%nu_mom_fs;c(3,:)=model%chi_i_fs
    c(4,:)=model%chi_e_fs;c(5,:)=model%vpinch_fs
  END FUNCTION
  SUBROUTINE check_equal(a,b,message)
    REAL*8,INTENT(IN)::a(:,:),b(:,:)
    CHARACTER(LEN=*),INTENT(IN)::message
    IF(ANY(a/=b))THEN
       WRITE(*,*)'FAIL: ',message,MAXVAL(ABS(a-b));ERROR STOP 1
    ENDIF
  END SUBROUTINE
  SUBROUTINE check_close(a,b,message)
    REAL*8,INTENT(IN)::a(:,:),b(:,:)
    CHARACTER(LEN=*),INTENT(IN)::message
    IF(MAXVAL(ABS(a-b))>1.d-9*MAX(MAXVAL(ABS(b)),1.d-12))THEN
       WRITE(*,*)'FAIL: ',message,MAXVAL(ABS(a-b));ERROR STOP 1
    ENDIF
  END SUBROUTINE
END PROGRAM
