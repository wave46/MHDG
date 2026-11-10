PROGRAM test_density_chord
  USE globals
  USE MPI_OMP
  USE density_diagnostics_module
  USE magnetic_geometry_state, ONLY: magnetic_equilibrium
  USE magnetic_field, ONLY: adjust_puff_to_target_density
  USE input_reader, ONLY: read_input
  USE HDF5_io_module
  IMPLICIT NONE
  INTEGER :: ierr, i, j, e, density_len
  INTEGER(HID_T) :: file_id
  REAL*8 :: full, mesh_length, length, diagonal(4), chord(4)
  CHARACTER(LEN=32) :: mode

  CALL get_command_argument(1,mode)
  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  IF (mode == 'input') THEN
     CALL read_input()
     WRITE(*,'(A,1X,L1)') 'DENSITY_INPUT',utils%line_integrated_density_diagnostics
     STOP
  ENDIF
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIvar%glob_id,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,MPIvar%glob_size,ierr)
  phys%lscale = 1.901d-3
  phys%Neq = 5
  simpar%refval_density = 3.d0
  simpar%refval_time = 1.d0
  refElPol%Nnodes2D = 3
  Mesh%Nelems = 6
  Mesh%Nnodes = 8
  ALLOCATE(Mesh%X(8,2),Mesh%T(6,3),sol%u(90))
#ifdef PARALL
  ALLOCATE(Mesh%ghostElems(6))
  DO e = 1,6
     Mesh%ghostElems(e) = MERGE(0,1,MOD(e-1,MPIvar%glob_size) == MPIvar%glob_id)
  ENDDO
#endif
  DO i = 1,4
     Mesh%X(2*i-1,1) = 1.d0
     Mesh%X(2*i,1) = 3.d0
  ENDDO
  Mesh%X(:,2) = [-1.d0,-1.d0,-.5d0,-.5d0,.5d0,.5d0,1.d0,1.d0]
  Mesh%X = Mesh%X*.25d0/phys%lscale
  DO i = 1,3
     Mesh%T(2*i-1,:) = [2*i-1,2*i,2*i+2]
     Mesh%T(2*i,:) = [2*i-1,2*i+2,2*i+1]
  ENDDO
  magnetic_equilibrium%is_initialized = .TRUE.
  ALLOCATE(magnetic_equilibrium%lcfs_r(4),magnetic_equilibrium%lcfs_z(4))
  magnetic_equilibrium%lcfs_r = [1.5d0,2.5d0,2.5d0,1.5d0]
  magnetic_equilibrium%lcfs_z = [-.5d0,-.5d0,.5d0,.5d0]
  magnetic_equilibrium%lcfs_r = magnetic_equilibrium%lcfs_r*.25d0/phys%lscale
  magnetic_equilibrium%lcfs_z = magnetic_equilibrium%lcfs_z*.25d0/phys%lscale
  input%density_chord = [.5d0,-.25d0,.5d0,.25d0]
  CALL set_density(.FALSE.)
  CALL measure_density_chord(input%density_chord,full,mesh_length,length)
  CALL close(full,2.997d0,'uniform full integral, native endpoint convention')
  CALL close(mesh_length,.5d0,'geometric mesh chord length')
  CALL close(length,.25d0,'geometric LCFS chord length')
  CALL compare_native_integral(input%density_chord,full)
  chord = [.5d0,-.4d0,.5d0,.4d0] ! endpoints outside the mesh
  CALL measure_density_chord(chord,full,mesh_length,length)
  CALL close(mesh_length,.5d0,'exclude outside-mesh chord portions')
  CALL close(length,.25d0,'LCFS length independent of outside endpoints')
  CALL compare_native_integral(chord,full)
  diagonal = [.25d0,-.125d0,.75d0,.125d0] ! exactly on an owned/ghost shared edge
  CALL measure_density_chord(diagonal,full,mesh_length,length)
  CALL close(full/mesh_length,5.994d0,'shared-edge MPI ownership')
  CALL close(length,SQRT(1.25d0)*.25d0,'diagonal polygon clipping')
  CALL set_density(.TRUE.)
  CALL measure_density_chord(input%density_chord,full,mesh_length,length)
  CALL close(length,.25d0,'LCFS length must not depend on SOL density')
  CALL compare_native_integral(input%density_chord,full)
  IF (full <= 8.d0) ERROR STOP 'SOL density missing from full integral'
  chord = input%density_chord([3,4,1,2])
  CALL measure_density_chord(chord,full,mesh_length,length)
  CALL close(mesh_length,.5d0,'reversed chord mesh length')
  CALL close(length,.25d0,'reversed chord LCFS length')
  CALL compare_native_integral(chord,full)
  CALL set_density(.FALSE.)
  IF (mode == 'miss') THEN
     input%density_chord = [.1d0,-.25d0,.1d0,.25d0]
     CALL measure_density_chord(input%density_chord,full,mesh_length,length)
     ERROR STOP 'missing LCFS intersection was accepted'
  ENDIF
  IF (mode == 'coverage') THEN
     magnetic_equilibrium%lcfs_z = [-2.d0,-2.d0,2.d0,2.d0]*.25d0/phys%lscale
     input%density_chord = [.5d0,-.5d0,.5d0,.5d0]
     CALL measure_density_chord(input%density_chord,full,mesh_length,length)
     ERROR STOP 'missing core mesh coverage was accepted'
  ENDIF

  switch%target_variable = 0
  time%legacy_input = .FALSE.
  switch%ME = .FALSE.
  switch%testcase = 54
  IF (mode == 'steady') THEN
     switch%target_variable = 1 ! Dormant selector, as in retained steady cases.
     time%legacy_input = .TRUE.
     input%target_density_path = 'absent.h5'
  ENDIF
  IF (mode == 'missing' .OR. mode == 'partial' .OR. mode == 'different' .OR. &
      mode == 'same' .OR. mode == 'reversed') THEN
     switch%target_variable = 1
     input%target_density_path = 'absent.h5'
  ENDIF
  IF (mode == 'partial' .OR. mode == 'different' .OR. mode == 'same' .OR. mode == 'reversed') THEN
     input%target_density_path = 'target.h5'
     chord = input%density_chord
     IF (mode == 'different') chord = [.5d0,-.125d0,.5d0,.125d0]
     IF (mode == 'reversed') chord = chord([3,4,1,2])
     IF (MPIvar%glob_id == 0) THEN
        CALL HDF5_create('target.h5',file_id)
        CALL HDF5_real_saving(file_id,chord(1),'x_lower')
        IF (mode /= 'partial') THEN
           CALL HDF5_real_saving(file_id,chord(3),'x_upper')
           CALL HDF5_real_saving(file_id,chord(2),'y_lower')
           CALL HDF5_real_saving(file_id,chord(4),'y_upper')
           CALL HDF5_array1D_saving(file_id,[0.d0,1.d0],2,'time')
           CALL HDF5_array1D_saving(file_id,[4.d0,4.d0],2,'target_density')
        ENDIF
        CALL HDF5_close(file_id)
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ENDIF
  time%t = .25d0
  utils%line_integrated_density_diagnostics = .FALSE.
  CALL report_density_diagnostics()
  IF (density_diagnostics_available) ERROR STOP 'reporting default is not inactive'
  utils%line_integrated_density_diagnostics = .TRUE.
  IF (mode == 'different' .OR. mode == 'same' .OR. mode == 'reversed') THEN
     time%it = 0
     time%dt = .02d0
     time%t_target = .25d0
     time%t_ME = .25d0
     phys%puff = 1.d0
     phys%feedback_propotional_gain = 0.d0
     phys%feedback_integral_gain = 0.d0
     phys%feedback_derivative_gain = 0.d0
     CALL adjust_puff_to_target_density(input%target_density_path,density_len)
     CALL adjust_puff_to_target_density(input%target_density_path,density_len)
     IF (mode == 'different') THEN
        ! Allow native MPI boundary quadrature differences, but reject the longer reporting chord.
        IF (phys%n_li < 1.d0 .OR. phys%n_li > 2.d0) ERROR STOP 'feedback did not retain its short HDF5 chord'
     ELSE
        IF (phys%n_li < 2.d0 .OR. phys%n_li > 4.d0) ERROR STOP 'feedback did not retain its full HDF5 chord'
     ENDIF
  ENDIF
  phys%n_li = -42.d0
  CALL report_density_diagnostics()
  CALL report_density_diagnostics()
  CALL close(phys%n_li,-42.d0,'reporting must not mutate feedback measurement')
  ! Later solution updates must not change the last assembled snapshot.
  sol%u = sol%u*2
  time%t = .5d0
  input%density_chord = [.1d0,-.25d0,.1d0,.25d0] ! Saving must retain the measured chord too.
  IF (MPIvar%glob_id == 0) THEN
     CALL HDF5_create('density.h5',file_id)
     CALL save_density_diagnostics(file_id)
     CALL HDF5_close(file_id)
  ENDIF
  IF (mode == 'reset') THEN
     CALL reset_density_diagnostics()
     IF (density_diagnostics_available) ERROR STOP 'reset retained an available snapshot'
     IF (MPIvar%glob_id == 0) THEN
        CALL HDF5_create('density_reset.h5',file_id)
        CALL save_density_diagnostics(file_id)
        CALL HDF5_close(file_id)
     ENDIF
     switch%target_variable = 0
     input%density_chord = [.5d0,-.4d0,.5d0,.4d0]
     time%t = .75d0
     CALL report_density_diagnostics()
     IF (MPIvar%glob_id == 0) THEN
        CALL HDF5_create('density_reinitialized.h5',file_id)
        CALL save_density_diagnostics(file_id)
        CALL HDF5_close(file_id)
     ENDIF
  ENDIF
  IF (MPIvar%glob_id == 0) THEN
     WRITE(*,'(A)') 'DENSITY_CHORD_PASS'
  ENDIF
  CALL MPI_FINALIZE(ierr)

CONTAINS
  SUBROUTINE compare_native_integral(chord,actual)
    USE interpolation, ONLY: lineintegration
    REAL*8, INTENT(IN) :: chord(4),actual
#ifndef PARALL
    REAL*8 :: x(1000),y(1000),reference
    INTEGER :: k
    x = [(chord(1)+(chord(3)-chord(1))/1000.d0*(k-1),k=1,1000)]/phys%lscale
    y = [(chord(2)+(chord(4)-chord(2))/1000.d0*(k-1),k=1,1000)]/phys%lscale
    CALL lineintegration(1000,x,y,sol%u(1::5),Mesh%X,Mesh%T,3,reference)
    CALL close(actual,reference*simpar%refval_density,'native feedback quadrature')
#endif
  END SUBROUTINE compare_native_integral

  SUBROUTINE set_density(with_sol)
    LOGICAL, INTENT(IN) :: with_sol
    sol%u = 0.d0
    DO e = 1,6
       DO j = 1,3
          sol%u(((e-1)*3+j-1)*5+1) = 2.d0
          IF (with_sol .AND. (e < 3 .OR. e > 4)) sol%u(((e-1)*3+j-1)*5+1) = 9.d0
#ifdef PARALL
          IF (Mesh%ghostElems(e) /= 0) sol%u(((e-1)*3+j-1)*5+1) = 1.d6
#endif
       ENDDO
    ENDDO
  END SUBROUTINE set_density

  SUBROUTINE close(actual,expected,label)
    REAL*8, INTENT(IN) :: actual,expected
    CHARACTER(LEN=*), INTENT(IN) :: label
    IF (ABS(actual-expected) <= 1.d-10*MAX(1.d0,ABS(expected))) RETURN
    WRITE(*,*) label,actual,expected
    ERROR STOP 1
  END SUBROUTINE close
END PROGRAM test_density_chord
