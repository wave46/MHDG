PROGRAM test_neutral_wall_diagnostics
  USE balance_diagnostics
  USE HDF5_io_module, ONLY: HID_T, HDF5_create, HDF5_close, &
       HDF5_group_create, HDF5_group_close
  USE MPI_OMP, ONLY: MPIvar
  IMPLICIT NONE
  TYPE(balance_diagnostics_type) :: diagnostics
  TYPE(balance_accumulator_type) :: accumulator
  TYPE(balance_plasma_bc_type) :: plasma
  REAL*8 :: state(5), gradient(2,5), normal(2), matrix(5,5), diffusion(5,5), pinch(5,2)
  REAL*8 :: wall_flux, puff, pump
  INTEGER(HID_T) :: file_id, group_id
  INTEGER :: mode, relocated, i
  CHARACTER(LEN=32) :: argument

  CALL get_command_argument(1,argument)
  READ(argument,*) mode
  CALL get_command_argument(2,argument)
  READ(argument,*) wall_flux
  CALL get_command_argument(3,argument)
  READ(argument,*) relocated
  MPIvar%glob_id = 0
  MPIvar%glob_size = 1
  state = [2.d0,5.d0,16.d0,6.d0,3.d0]
  normal = [1.d0,0.d0]
  matrix = 0.d0
  pinch = 0.d0
  diffusion = 0.d0
  diffusion(5,5) = 2.d0
  puff = 7.d0
  pump = 2.d0
  IF (relocated == 1) THEN
     puff = 0.d0
     pump = 0.d0
  ENDIF
  gradient = 0.d0
  gradient(1,5) = (2.5d0+puff-pump*state(5)-wall_flux)/2.d0
  CALL diagnostics%configure(mode)
  CALL diagnostics%begin_assembly(1.d0,1.d0,1.d0/(2.d0*ACOS(-1.d0)),1.d0)
  DO i = 1,2
     CALL diagnostics%initialize_accumulator(accumulator)
     CALL accumulator%accumulate_bc( &
          integration_weight=1.d0, neutral_equation=5, &
          trace_state=state, exterior_state=state, gradient=gradient, &
          normal=normal, magnetic_direction=normal, magnetic_normal=1.d0, &
          tau=matrix, diffusion_iso=diffusion, diffusion_ani=matrix, &
          pinch_matrix=pinch, recycling_coefficient=0.5d0, &
          puff_source=puff, pump_coefficient=pump, &
          neutral_wall_absorption=wall_flux, plasma=plasma, &
          neutral_perpendicular_diffusion=.FALSE.)
     IF (relocated == 1) CALL accumulator%accumulate_relocated_sources(7.d0,6.d0)
     CALL diagnostics%merge(accumulator)
  ENDDO
  CALL diagnostics%finalize()
  CALL diagnostics%report()
  CALL HDF5_create('diagnostics.h5',file_id)
  IF (diagnostics%has_values()) THEN
     CALL HDF5_group_create('diagnostics',file_id,group_id)
     CALL diagnostics%write_hdf5(group_id)
     CALL HDF5_group_close(group_id)
  ENDIF
  CALL HDF5_close(file_id)
END PROGRAM test_neutral_wall_diagnostics
