PROGRAM test_transport_pinch
  USE globals, ONLY: input, transport_model_input, simpar, switch
  USE MPI_OMP, ONLY: MPIvar
  USE transport_models_1d, ONLY: transport_model_1d_t
  USE transport_models_1d_derived, ONLY: transport_model_derived_t
  USE magnetic_topology, ONLY: magnetic_region_core, magnetic_region_main_sol
  USE HDF5_io_module, ONLY: HID_T, HDF5_create, HDF5_close, &
       HDF5_group_create, HDF5_group_close, HDF5_array2D_saving
  IMPLICIT NONE

  TYPE(transport_model_1d_t) :: model
  TYPE(transport_model_derived_t) :: work
  REAL*8 :: pinch(6,2), reference(6,2)
  INTEGER(HID_T) :: file_id, group_id
  INTEGER :: i
  CHARACTER(LEN=32) :: argument

  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  simpar%refval_time = 2.d0
  simpar%refval_length = 4.d0
  switch%transport_1d = .TRUE.

  CALL get_command_argument(1,argument)
  SELECT CASE (TRIM(argument))
  CASE ('negative')
     CALL model%set_config(pinch_equations=[-1])
     ERROR STOP 'negative equation accepted'
  CASE ('neutral')
     CALL model%set_config(pinch_equations=[5])
     ERROR STOP 'neutral equation accepted'
  CASE ('overlong')
     CALL model%set_config(pinch_equations=[1,2,3,4,0])
     ERROR STOP 'overlong equation list accepted'
  END SELECT

  input%transport_model_path = 'transport_model.nml'
  CALL read_transport_model_input()
  CALL model%set_config(pinch_equations=transport_model_input%pinch_equations, &
       pinch_model=transport_model_input%pinch_model, &
       vpinch_const_phys=transport_model_input%vpinch_const_phys)
  CALL model%init(3)
  model%rho_grid = [0.d0,0.5d0,1.d0]
  CALL model%compute_pinch_profile(work)
  CALL evaluate_pinch(0.5d0,magnetic_region_core)
  reference = pinch

  CALL HDF5_create('transport_pinch.h5',file_id)
  CALL HDF5_array2D_saving(file_id,pinch,6,2,'pinch')
  CALL HDF5_group_create('transport_1d',file_id,group_id)
  CALL model%write_hdf5(group_id)
  CALL HDF5_group_close(group_id)
  CALL HDF5_close(file_id)

  ! Updating an unrelated setting and reallocating profiles must keep selection.
  CALL model%set_config(c_pinch=0.25d0)
  CALL model%init(5)
  model%rho_grid = [(0.25d0*(i-1),i=1,5)]
  CALL model%compute_pinch_profile(work)
  CALL evaluate_pinch(0.5d0,magnetic_region_core)
  IF (ANY(pinch /= reference)) ERROR STOP 'pinch selection lost on reinitialization'

  CALL evaluate_pinch(0.d0,magnetic_region_core)
  IF (ANY(pinch /= 0.d0)) ERROR STOP 'axis pinch must be zero'
  CALL evaluate_pinch(1.d0,magnetic_region_core)
  IF (ANY(pinch /= 0.d0)) ERROR STOP 'outer pinch must be zero'
  CALL model%set_config(transport_region_policy='core_only')
  CALL evaluate_pinch(0.5d0,magnetic_region_main_sol)
  IF (ANY(pinch /= 0.d0)) ERROR STOP 'excluded-region pinch must be zero'

  CALL model%destroy()
  CALL model%set_config(pinch_model=3,vpinch_const_phys=-8.d0)
  CALL model%init(3)
  model%rho_grid = [0.d0,0.5d0,1.d0]
  CALL model%compute_pinch_profile(work)
  CALL evaluate_pinch(0.5d0,magnetic_region_core)
  IF (ANY(pinch(1,:) /= [-4.d0,0.d0]) .OR. &
       ANY(pinch(2:,:) /= 0.d0)) ERROR STOP 'reset must restore particle-only pinch'
  WRITE (*,'(A)') 'transport pinch checks: PASS'

CONTAINS

  SUBROUTINE evaluate_pinch(rho,region)
    REAL*8, INTENT(IN) :: rho
    INTEGER, INTENT(IN) :: region

    ! Seed every row so inactive rows and suppression must explicitly zero them.
    pinch = 99.d0
    CALL model%compute_1D_pinch_matrix([0.d0,1.d0],rho,pinch, &
         region,[1.d0,0.d0])
  END SUBROUTINE evaluate_pinch

END PROGRAM test_transport_pinch
