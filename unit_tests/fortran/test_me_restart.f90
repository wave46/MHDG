PROGRAM test_me_restart
  USE initialization, ONLY: init_sim
  USE in_out, ONLY: HDF5_load_solution
  USE globals
  USE MPI_OMP, ONLY: MPIvar
  USE interpolation, ONLY: binarysearch, interpolate
  USE external_heating, ONLY: initialize_external_heating, load_external_heating
  IMPLICIT NONE
  INTEGER :: nts
  REAL :: dt, restored_time
  CHARACTER(LEN=1000) :: checkpoint
  REAL :: knots(3), values(2,2), expected_heating

  ! Exact control endpoints and grid corners must leave index+1 in bounds.
  knots = [0.d0, 0.5d0, 0.56d0]
  IF (binarysearch(3,knots,knots(1),1.d-12) /= 1) ERROR STOP 'first interpolation endpoint'
  IF (binarysearch(3,knots,knots(3),1.d-12) /= 2) ERROR STOP 'last interpolation endpoint'
  values = RESHAPE([1.d0,2.d0,3.d0,4.d0],[2,2])
  IF (interpolate(2,knots(1:2),2,knots(1:2),values,0.5d0,0.5d0,1.d-12) /= 4.d0) &
       ERROR STOP 'last grid corner'

  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  CALL read_input()
  CALL init_sim(nts, dt)

  ! The legacy heating reader follows the same initial/endpoint frame selection.
  Mesh%Nnodes = 1
  ALLOCATE(Mesh%X(1,2))
  Mesh%X = 0.5d0/phys%lscale
  input%external_heating_from_grid = .TRUE.
  input%external_heating_path = 'heat_0000.h5'
  IF (switch%ME) input%external_heating_path = 'heat'
  CALL initialize_external_heating()
  CALL load_external_heating()
  IF (ANY(phys%external_heating_ions /= 0.d0) .OR. ANY(phys%external_heating_electrons /= 0.d0)) &
       ERROR STOP 'initial heating frame'
  IF (switch%ME) THEN
     CALL load_external_heating(1)
     expected_heating = 7.d0*simpar%refval_time/simpar%refval_specenergydens/simpar%refval_mass
     IF (ABS(phys%external_heating_ions(1)-expected_heating) > 1.d-12*ABS(expected_heating) .OR. &
         ABS(phys%external_heating_electrons(1)-2.d0*expected_heating) > 1.d-12*ABS(expected_heating)) &
          ERROR STOP 'endpoint heating frame'
  ENDIF
  DEALLOCATE(Mesh%X,phys%external_heating_ions,phys%external_heating_electrons)

  ! Only the dimensions used by the production solution reader are needed.
  Mesh%Nelems = 1
  Mesh%Nfaces = 1
  Mesh%Nnodesperface = 1
  refElPol%Nnodes2D = 1
  refElPol%Nfaces = 1
  refElPol%Nfacenodes = 1
  phys%feedback_integral_error = -11.d0
  phys%feedback_previous_error = -12.d0
  CALL get_command_argument(1, checkpoint)
  CALL HDF5_load_solution(checkpoint)
  restored_time = 0.d0
  IF (sol%Nt > 0) restored_time = sol%time(sol%Nt)
  IF (time%it /= time%ik .OR. time%it /= sol%Nt) STOP 2
  WRITE (*,'(A,2(1X,I0),5(1X,ES24.16))') 'ME_RESTART_STATE', &
       time%it, sol%Nt, time%t, restored_time, phys%puff, &
       phys%feedback_integral_error, phys%feedback_previous_error
END PROGRAM test_me_restart
