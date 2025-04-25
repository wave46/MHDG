MODULE external_heating
  USE prec_const
  USE globals
  USE in_out
  USE HDF5_io_module
  USE HDF5
  USE interpolation
CONTAINS

  !*******************************************************
  ! Allocate storing space for the external heating field
  !*******************************************************
  SUBROUTINE initialize_external_heating
    INTEGER :: nnodes
  
#ifdef TOR3D
    nnodes = Mesh%Nnodes*Mesh%Nnodes_toroidal
#else
    nnodes = Mesh%Nnodes
#endif

    ALLOCATE (phys%external_heating(nnodes))
    ALLOCATE (phys%external_heating_distribtution(2))
    phys%external_heating = 0.0
END SUBROUTINE initialize_external_heating

  !*******************************************************
  ! Read the external heating field from a file
  !*******************************************************
  SUBROUTINE load_external_heating
    IF (input%external_heating_from_grid) THEN
      CALL load_external_heating_from_grid
    ELSE
      CALL load_external_heating_nodes
    END IF 
  END SUBROUTINE load_external_heating

  !*******************************************************
  ! Read the external heating field from a grid file
  !*******************************************************
  SUBROUTINE load_external_heating_from_grid
    CHARACTER(LEN=1000)                           :: fname
    INTEGER(HID_T)                                :: file_id
    INTEGER                                       :: IERR, ip, jp, ind, i, k
    CHARACTER(70)        :: nit
    REAL*8, POINTER, DIMENSION(:, :)              :: r2D => NULL(), z2D => NULL(), heating2D => NULL()
    REAL*8, ALLOCATABLE, DIMENSION(:)             :: xvec, yvec
    REAL*8                                        :: heating, x, y

    IF (utils%printint > 0) THEN
      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (6, *) '*************************************************'
         WRITE (6, *) '*       LOADING EXTERNAL HEATING SOURCES        *'
         WRITE (6, *) '*************************************************'
      ENDIF
    END IF

    IF (switch%ME .EQV. .FALSE.) THEN ! if not a moving equilibrium simulation
      fname = input%external_heating_path
    ELSE
      fname = input%external_heating_path
      WRITE(nit, "(i10)") INT(time%it + 1)
      nit = TRIM(ADJUSTL(nit))
      k = INDEX(nit, " ") - 1
      fname = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
    ENDIF

    CALL HDF5_open(fname, file_id, IERR)
    IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6,*) 'Magnetic field loaded from file: ', TRIM(ADJUSTL(fname))
    ENDIF
    CALL HDF5_integer_reading(file_id, ip, 'ip')
    CALL HDF5_integer_reading(file_id, jp, 'jp')

    ALLOCATE (r2D(ip, jp))
    ALLOCATE (z2D(ip, jp))
    ALLOCATE (heating2D(ip, jp))

    CALL HDF5_array2D_reading(file_id, r2D, 'r2D')
    CALL HDF5_array2D_reading(file_id, z2D, 'z2D')
    CALL HDF5_array2D_reading(file_id, heating2D, 'heating2D')
    CALL HDF5_array1D_reading(file_id, phys%external_heating_distribtution, 'additional_heating_distribtution')

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale

    ! Apply energy scale
    heating2D = heating2D*simpar%refval_time/simpar%refval_specenergydens/ &
                simpar%refval_mass

    ! Interpolate the external heating field to the nodes
    ALLOCATE (xvec(jp))
    ALLOCATE (yvec(ip))
    xvec = r2D(1, :)
    yvec = z2D(:, 1)

    DO i = 1, Mesh%Nnodes
      x = Mesh%X(i, 1)
      y = Mesh%X(i, 2)
      heating = interpolate(ip, yvec, jp, xvec, heating2D, y, x, 1e-12)

    ind = i
#ifdef TOR3D
      DO j = 1, Mesh%Nnodes_toroidal
        ind = (j - 1)*Mesh%Nnodes + i
#endif
      phys%external_heating(ind) = heating
#ifdef TOR3D
      END DO
#endif
    END DO


    ! Deallocate arrays
    DEALLOCATE (r2D, z2D, heating2D)
    NULLIFY (r2D, z2D, heating2D)
    DEALLOCATE (xvec, yvec)

    CALL HDF5_close(file_id)
  END SUBROUTINE load_external_heating_from_grid

  !*******************************************************
  ! Read the external heating field from a file
  !*******************************************************
  SUBROUTINE load_external_heating_nodes
    WRITE (*,*) "Not implemented yet"
    STOP
  END SUBROUTINE load_external_heating_nodes

END MODULE external_heating