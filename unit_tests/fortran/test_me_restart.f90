PROGRAM test_me_restart
  USE initialization, ONLY: init_sim
  USE in_out, ONLY: HDF5_load_solution
  USE globals
  USE MPI_OMP, ONLY: MPIvar
  IMPLICIT NONE
  INTEGER :: nts
  REAL :: dt, restored_time
  CHARACTER(LEN=1000) :: checkpoint

  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  CALL read_input()
  CALL init_sim(nts, dt)

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
