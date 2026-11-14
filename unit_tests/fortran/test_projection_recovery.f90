PROGRAM test_projection_recovery
  USE input_reader, ONLY: read_input
  USE globals, ONLY: adapt, phys, refElPol, utils
  USE MPI_OMP, ONLY: MPIvar
  USE mpi, ONLY: MPI_INIT, MPI_FINALIZE, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_COMM_WORLD
  USE reference_element, ONLY: create_reference_element
  USE adaptivity_projection_module, ONLY: mesh_transfer_t
  IMPLICIT NONE
  TYPE(mesh_transfer_t) :: transfer
  INTEGER, ALLOCATABLE :: connectivity(:,:)
  REAL*8, ALLOCATABLE :: old_coordinates(:,:), new_coordinates(:,:), u_old(:), u_new(:), q_old(:), q_new(:)
  REAL*8, ALLOCATABLE :: expected(:)
  REAL*8 :: gap
  INTEGER :: ierr, degree, nodes, node
  CHARACTER(LEN=64) :: argument

  CALL get_command_argument(1,argument)
  IF (argument == 'input') THEN
     MPIvar%glob_id = 0
     MPIvar%glob_size = 1
     CALL read_input()
     WRITE(*,'(A,1X,ES24.16)') 'PROJECTION_RECOVERY_RATIO', adapt%projection_recovery_distance_ratio
     STOP
  ENDIF
  IF (argument /= 'default') READ(argument,*) adapt%projection_recovery_distance_ratio
  CALL get_command_argument(2,argument)
  READ(argument,*) gap
  CALL get_command_argument(3,argument)
  READ(argument,*) degree
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIvar%glob_id,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,MPIvar%glob_size,ierr)
  CALL create_reference_element(refElPol,2,degree,verbose=0)
  utils%printint = 0
  utils%timing = .FALSE.
  phys%neq = 1
  nodes = refElPol%Nnodes2D
  ALLOCATE(connectivity(1,nodes),old_coordinates(nodes,2),new_coordinates(nodes,2))
  ALLOCATE(u_old(nodes),u_new(nodes),q_old(2*nodes),q_new(2*nodes),expected(nodes))
  connectivity(1,:) = [(node,node=1,nodes)]

  ! An anisotropic triangle: shortest corner edge is 1, longest exceeds 4.
  ! The new triangle has its base gap units outside the old lower edge.
  old_coordinates(:,1) = 0.5d0*(refElPol%coord2d(:,1)+1.d0)
  old_coordinates(:,2) = 2.d0*(refElPol%coord2d(:,2)+1.d0)
  new_coordinates(:,1) = 0.25d0+0.25d0*old_coordinates(:,1)
  new_coordinates(:,2) = 0.25d0*old_coordinates(:,2)-gap
  u_old = 2.d0+old_coordinates(:,1)+3.d0*old_coordinates(:,2)
  expected = 2.d0+new_coordinates(:,1)+3.d0*MAX(0.d0,new_coordinates(:,2))
  q_old(1::2) = 7.d0
  q_old(2::2) = 11.d0
  CALL transfer%build(connectivity,old_coordinates,connectivity,new_coordinates)
  CALL transfer%apply(u_old,u_new,q_old,q_new)
  IF (MAXVAL(ABS(u_new-expected)) > 1.d-7) ERROR STOP 'Incorrect closest-point interpolation'
  IF (MAXVAL(ABS(q_new(1::2)-7.d0)) > 1.d-7 .OR. &
      MAXVAL(ABS(q_new(2::2)-11.d0)) > 1.d-7) ERROR STOP 'Gradient transfer changed'
  CALL transfer%clear()
  IF (MPIvar%glob_id == 0) WRITE(*,'(A)') 'PROJECTION_RECOVERY_PASS'
  CALL MPI_FINALIZE(ierr)
END PROGRAM test_projection_recovery
