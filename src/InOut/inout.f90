!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module to load/save files
! in the code
!************************************************************

MODULE in_out
  USE HDF5
  USE HDF5_io_module
  USE globals
  USE printutils
  USE MPI_OMP
  IMPLICIT NONE

CONTAINS

  !********************************
  ! Loads mesh from an hdf5 file
  ! external file
  !********************************
  SUBROUTINE load_mesh_serial(fname)
    USE MPI_OMP
    CHARACTER(LEN=*) :: fname
    CHARACTER(len=1000) :: fname_complete
    CHARACTER(10)  :: str
    REAL*8, PARAMETER::tol = 1e-6
    REAL*8 :: xmin
    INTEGER :: elemType, ndim, Nnodes, Nelems, Nnodesperelem, Nfaces
    INTEGER :: Nextfaces, Nnodesperface, IERR
    INTEGER(HID_T) :: file_id

    fname_complete = TRIM(ADJUSTL(fname))//'.h5'

    IF (utils%printint > 0) THEN
       PRINT *, 'Loading mesh.'
       PRINT *, '        '
    ENDIF

    CALL HDF5_open(fname_complete, file_id, IERR)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error opening mesh file: ", fname_complete
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, elemType, 'elemType', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: elemType"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nfaces, 'Nfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nfaces"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, ndim, 'Ndim', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Ndim"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodes, 'Nnodes', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodes"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nelems, 'Nelems', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nelems"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperelem, 'Nnodesperelem', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperelem"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperface, 'Nnodesperface', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperface"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nextfaces, 'Nextfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nextfaces"
       STOP
    ENDIF

    ALLOCATE (Mesh%T(Nelems, Nnodesperelem))
    ALLOCATE (Mesh%X(Nnodes, ndim))
    ALLOCATE (Mesh%Tb(Nextfaces, Nnodesperface))
    ALLOCATE (Mesh%boundaryFlag(Nextfaces))

    CALL HDF5_array2D_reading_int(file_id, Mesh%T, 'T', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading mesh connectivity T"
       STOP
    ENDIF
    CALL HDF5_array2D_reading_int(file_id, Mesh%Tb, 'Tb', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundary connectivity Tb"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%boundaryFlag, 'boundaryFlag', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundaryFlag"
       STOP
    ENDIF
    CALL HDF5_array2D_reading(file_id, Mesh%X, 'X', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading coordinate matrix X"
       STOP
    ENDIF

    CALL HDF5_close(file_id)

    WRITE (6, *) "Readed mesh file: ", TRIM(ADJUSTL(fname_complete))

    Mesh%Ndim = ndim
    Mesh%Nnodes = Nnodes
    Mesh%Nelems = Nelems
    Mesh%Nnodesperelem = Nnodesperelem
    Mesh%Nnodesperface = Nnodesperface
    Mesh%elemType = elemType
    Mesh%Nextfaces = Nextfaces

    xmin = MINVAL(Mesh%X(:,1))

    ! Apply shift if axisymmetric case
    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) .OR. (switch%axisym .AND. xmin < tol)) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) "*** Applying translation in axisymmetric case!"
       ENDIF
       Mesh%X(:, 1) = Mesh%X(:, 1) + geom%R0
    END IF

    ! Apply length scale
    Mesh%X = Mesh%X/phys%lscale

    Mesh%xmax = MAXVAL(Mesh%X(:, 1))
    Mesh%xmin = MINVAL(Mesh%X(:, 1))
    Mesh%ymax = MAXVAL(Mesh%X(:, 2))
    Mesh%ymin = MINVAL(Mesh%X(:, 2))

    IF (utils%printint > 0) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (elemType == 0) THEN
             WRITE (str, '(A)') 'triangles'
          ELSEIF (elemType == 1) THEN
             WRITE (str, '(A)') 'quads'
          ELSEIF (elemType == 2) THEN
             WRITE (str, '(A)') 'thetra'
          ELSEIF (elemType == 3) THEN
             WRITE (str, '(A)') 'hexa'
          END IF
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*                    MESH                       *'
          WRITE (6, *) '*************************************************'
          WRITE (6, '(A,I18)') ' Number of dimensions:         ', ndim
          WRITE (6, '(A,A34)') ' Element type: ', TRIM(str)
          WRITE (6, '(A,I18)') ' Number of elements:           ', Nelems
          WRITE (6, '(A,I18)') ' Number of nodes:              ', Nnodes
          WRITE (6, '(A,I18)') ' Number of nodes per element:  ', Nnodesperelem
          WRITE (6, '(A,I18)') ' Number of nodes per face:     ', Nnodesperface
          WRITE (6, '(A,I18)') ' Number of exterior faces:     ', Nextfaces
          WRITE (6, *) ' '
          WRITE (6, *) ' '
          IF (utils%printint > 1) THEN
             WRITE (6, *) "Connectivity matrix T:"
             CALL displayMatrixInt(Mesh%T)
             WRITE (6, *) "Boundary connectivity matrix Tb:"
             CALL displayMatrixInt(Mesh%Tb)
          END IF
       ENDIF
    ENDIF

  END SUBROUTINE load_mesh_serial

  SUBROUTINE load_mesh(fname)
    USE MPI_OMP
    CHARACTER(LEN=*) :: fname
    CHARACTER(len=1000) :: fname_complete
    CHARACTER(10)  :: str
    CHARACTER(70)  :: npr, nid
    REAL*8, PARAMETER::tol = 1e-6
    REAL*8 :: xmin
    INTEGER :: elemType, ndim, Nnodes, Nelems, Nnodesperelem
    INTEGER :: Nextfaces, Nnodesperface, IERR
    INTEGER(HID_T) :: file_id
#ifdef PARALL
    INTEGER :: ghfa, ghel, Nel_glob, Nfa_glob, Ndir_glob, Ngho_glob, Nfaces
#endif
#ifdef TOR3D
    INTEGER :: i
#endif
#ifdef TOR3D
#ifdef PARALL
    IF (MPIvar%npol .GT. 1) THEN
       WRITE (nid, *) MPIvar%ipol
       WRITE (npr, *) MPIvar%npol
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
#else
    fname_complete = TRIM(ADJUSTL(fname))//'.h5'
#endif
#else
    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
#endif
    IF (utils%printint > 0) THEN
       PRINT *, 'Loading mesh.'
       PRINT *, '        '
    ENDIF

    CALL HDF5_open(fname_complete, file_id, IERR)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error opening mesh file: ", fname_complete
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, elemType, 'elemType', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: elemType"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, ndim, 'Ndim', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Ndim"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodes, 'Nnodes', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodes"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nelems, 'Nelems', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nelems"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperelem, 'Nnodesperelem', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperelem"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperface, 'Nnodesperface', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperface"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nextfaces, 'Nextfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nextfaces"
       STOP
    ENDIF
#ifdef PARALL
    CALL HDF5_integer_reading(file_id, Nfaces, 'Nfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nfaces"
       STOP
    ENDIF
#endif
    ALLOCATE (Mesh%T(Nelems, Nnodesperelem))
    ALLOCATE (Mesh%X(Nnodes, ndim))
    ALLOCATE (Mesh%Tb(Nextfaces, Nnodesperface))
    ALLOCATE (Mesh%boundaryFlag(Nextfaces))
#ifdef PARALL
    ALLOCATE (Mesh%ghostFaces(Nfaces))
    ALLOCATE (Mesh%loc2glob_fa(Nfaces))
    ALLOCATE (Mesh%loc2glob_el(Nelems))
    ALLOCATE (Mesh%ghostElems(Nelems))
#endif
    CALL HDF5_array2D_reading_int(file_id, Mesh%T, 'T', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading mesh connectivity T"
       STOP
    ENDIF
    CALL HDF5_array2D_reading_int(file_id, Mesh%Tb, 'Tb', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundary connectivity Tb"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%boundaryFlag, 'boundaryFlag', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundaryFlag"
       STOP
    ENDIF
    CALL HDF5_array2D_reading(file_id, Mesh%X, 'X', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading coordinate matrix X"
       STOP
    ENDIF
#ifdef PARALL
    CALL HDF5_array1D_reading_int(file_id, Mesh%loc2glob_fa, 'loc2glob_fa', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading loc2glob_fa"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%loc2glob_el, 'loc2glob_el', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading loc2glob_el"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostFaces, 'ghostFaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostFaces"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostElems, 'ghostElems', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostElems"
       STOP
    ENDIF
    ! Find the number of ghost faces
    ghfa = SUM(Mesh%ghostFaces)
    Mesh%nghostfaces = ghfa

    ! Find the number of ghost elements
    ghel = SUM(Mesh%ghostElems)
    Mesh%nghostElems = ghel

    ALLOCATE (Mesh%ghostflp(ghfa))
    ALLOCATE (Mesh%ghostpro(ghfa))
    ALLOCATE (Mesh%ghostloc(ghfa))
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostflp, 'ghostFlp', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostFlp"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostLoc, 'ghostLoc', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostLoc"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostPro, 'ghostPro', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostPro"
       STOP
    ENDIF

#ifdef TOR3D
    IF (MPIvar%ntor > 1) THEN
       DO i = 1, SIZE(Mesh%ghostPro)
          IF (Mesh%ghostPro(i) .GT. -1) THEN
             Mesh%ghostPro(i) = Mesh%ghostPro(i) + (MPIvar%itor - 1)*MPIvar%npol
          ENDIF
       END DO
    ENDIF
    ALLOCATE (Mesh%ghelspro(ghel))
    ALLOCATE (Mesh%ghelsloc(ghel))
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghelsLoc, 'ghelsLoc', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghelsLoc"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghelsPro, 'ghelsPro', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghelsPro"
       STOP
    ENDIF
    IF (MPIvar%ntor .GT. 1) THEN
       DO i = 1, SIZE(Mesh%ghelspro)
          IF (Mesh%ghelspro(i) .GT. -1) THEN
             Mesh%ghelsPro(i) = Mesh%ghelsPro(i) + (MPIvar%itor - 1)*MPIvar%npol
          END IF
       END DO
    END IF
#endif
#endif
    CALL HDF5_close(file_id)

    !************************************************************************
    !   CONFIRMATION MESSAGE FOR THE USER
    !************************************************************************
#ifdef PARALL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    WRITE (6, *) "Process: ", MPIvar%glob_id, "-- readed mesh file: ", TRIM(ADJUSTL(fname_complete))
#else
    WRITE (6, *) "Readed mesh file: ", TRIM(ADJUSTL(fname_complete))
#endif

#ifdef PARALL
    CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_el), Nel_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_fa), Nfa_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%ndir, Ndir_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%nghostfaces, Ngho_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    Mesh%Nel_glob = Nel_glob
    Mesh%Nfa_glob = Nfa_glob
    Mesh%Ndir_glob = Ndir_glob
    Mesh%Ngho_glob = Ngho_glob
#endif
    Mesh%Ndim = ndim
    Mesh%Nnodes = Nnodes
    Mesh%Nelems = Nelems
    Mesh%Nnodesperelem = Nnodesperelem
    Mesh%Nnodesperface = Nnodesperface
    Mesh%elemType = elemType
    Mesh%Nextfaces = Nextfaces

    xmin = MINVAL(Mesh%X(:,1))
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
    ! Apply shift if axisymmetric case
    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) .OR. (switch%axisym .AND. xmin < tol)) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) "*** Applying translation in axisymmetric case!"
       ENDIF
       Mesh%X(:, 1) = Mesh%X(:, 1) + geom%R0
    END IF

    ! Apply length scale
    Mesh%X = Mesh%X/phys%lscale

    Mesh%xmax = MAXVAL(Mesh%X(:, 1))
    Mesh%xmin = MINVAL(Mesh%X(:, 1))
    Mesh%ymax = MAXVAL(Mesh%X(:, 2))
    Mesh%ymin = MINVAL(Mesh%X(:, 2))

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    IF (utils%printint > 0) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (elemType == 0) THEN
             WRITE (str, '(A)') 'triangles'
          ELSEIF (elemType == 1) THEN
             WRITE (str, '(A)') 'quads'
          ELSEIF (elemType == 2) THEN
             WRITE (str, '(A)') 'thetra'
          ELSEIF (elemType == 3) THEN
             WRITE (str, '(A)') 'hexa'
          END IF
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*                    MESH                       *'
          WRITE (6, *) '*************************************************'
          WRITE (6, '(A,I18)') ' Number of dimensions:         ', ndim
          WRITE (6, '(A,A34)') ' Element type: ', TRIM(str)
          WRITE (6, '(A,I18)') ' Number of elements:           ', Nelems
          WRITE (6, '(A,I18)') ' Number of nodes:              ', Nnodes
          WRITE (6, '(A,I18)') ' Number of nodes per element:  ', Nnodesperelem
          WRITE (6, '(A,I18)') ' Number of nodes per face:     ', Nnodesperface
          WRITE (6, '(A,I18)') ' Number of exterior faces:     ', Nextfaces
          WRITE (6, *) ' '
          WRITE (6, *) ' '
          IF (utils%printint > 1) THEN
             WRITE (6, *) "Connectivity matrix T:"
             CALL displayMatrixInt(Mesh%T)
             WRITE (6, *) "Boundary connectivity matrix Tb:"
             CALL displayMatrixInt(Mesh%Tb)
          END IF
       ENDIF
    ENDIF

  END SUBROUTINE load_mesh

  !**********************************************************************
  ! Save solution in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_solution(fname)
    USE globals
    IMPLICIT NONE

    CHARACTER(LEN=*) :: fname

#ifdef TOR3D
    CHARACTER(70)  :: nip, nit, ngd
#else
    CHARACTER(70)  :: npr, nid
#endif
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

#ifdef TOR3D
    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nip, *) MPIvar%ipol
       WRITE (nit, *) MPIvar%itor
       WRITE (ngd, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_ip'//TRIM(ADJUSTL(nip))//'_it'//TRIM(ADJUSTL(nit))//'_np'//TRIM(ADJUSTL(ngd))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
#else
    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
#endif
    CALL HDF5_create(fname_complete, file_id, ierr)
    CALL HDF5_array1D_saving(file_id, sol%u, SIZE(sol%u), 'u')
    CALL HDF5_array1D_saving(file_id, sol%u_tilde, SIZE(sol%u_tilde), 'u_tilde')
    IF (phys%Neq .GE. 5) THEN
       IF (switch%saveTau) THEN
          CALL HDF5_array1D_saving(file_id, phys%diff_nn_Vol, SIZE(phys%diff_nn_Vol), 'DnnVol')
          CALL HDF5_array1D_saving(file_id, phys%diff_nn_Fac, SIZE(phys%diff_nn_Fac), 'DnnFac')
          CALL HDF5_array1D_saving(file_id, phys%diff_nn_Bou, SIZE(phys%diff_nn_Bou), 'DnnBou')
          CALL HDF5_array2D_saving(file_id, phys%v_nn_Vol, SIZE(phys%v_nn_Vol, 1), SIZE(phys%v_nn_Vol, 2), 'VnnVol')
          CALL HDF5_array2D_saving(file_id, phys%v_nn_Fac, SIZE(phys%v_nn_Fac, 1), SIZE(phys%v_nn_Fac, 2), 'VnnFac')
          CALL HDF5_array2D_saving(file_id, phys%v_nn_Bou, SIZE(phys%v_nn_Bou, 1), SIZE(phys%v_nn_Bou, 2), 'VnnBou')
          CALL HDF5_array2D_saving(file_id, Mesh%Xg, SIZE(Mesh%Xg, 1), SIZE(Mesh%Xg, 2), 'Xg')
          CALL HDF5_array2D_saving(file_id, Mesh%Xgf, SIZE(Mesh%Xgf, 1), SIZE(Mesh%Xgf, 2), 'Xgf')
          CALL HDF5_array2D_saving(file_id, Mesh%Xgb, SIZE(Mesh%Xgb, 1), SIZE(Mesh%Xgb, 2), 'Xgb')
       ENDIF
    ENDIF
    !      call HDF5_array1D_saving(file_id,sol%tres,sol%Nt,'tres')
    !      call HDF5_array1D_saving(file_id,sol%time,sol%Nt,'time')
    IF (switch%steady .OR. switch%psdtime) THEN
       CALL HDF5_integer_saving(file_id,0,'it')
    ELSE
       CALL HDF5_integer_saving(file_id,time%it,'it')
    ENDIF
!!!      call HDF5_integer_saving(file_id,sol%Nt,'Nt')
    CALL HDF5_array1D_saving(file_id, sol%q, SIZE(sol%q), 'q')
    ! Save magnetic field
    CALL HDF5_array2D_saving(file_id, phys%B, SIZE(phys%B, 1), SIZE(phys%B, 2), 'magnetic_field')
    ! Save normalized psi
    CALL HDF5_array1D_saving(file_id, phys%magnetic_psi, SIZE(phys%magnetic_psi), 'magnetic_psi')
    ! Save toroidal current
    IF (switch%ohmicsrc) THEN
       CALL HDF5_array1D_saving(file_id, phys%Jtor, SIZE(phys%Jtor), 'Jtor')
    ENDIF
    ! Save magnetic perturbation and related fields
    IF ((switch%rmp).OR.(switch%ripple)) THEN
       CALL HDF5_array2D_saving(file_id, phys%Bperturb, SIZE(phys%Bperturb, 1), SIZE(phys%Bperturb, 2), 'magnetic_perturbation')
    ENDIF
    IF (switch%rmp) THEN
       CALL HDF5_array3D_saving(file_id, magn%coils_rmp, SIZE(magn%coils_rmp, 1), SIZE(magn%coils_rmp, 2), SIZE(magn%coils_rmp, 3), 'coils_rmp')
    ENDIF
    IF (switch%ripple) THEN
       CALL HDF5_array2D_saving(file_id, magn%coils_ripple, SIZE(magn%coils_ripple, 1), SIZE(magn%coils_ripple, 2), 'coils_ripple')
    ENDIF
    ! Save boundary structure
    CALL HDF5_array2D_saving_int(file_id, Mesh%extfaces, SIZE(Mesh%extfaces, 1), SIZE(Mesh%extfaces, 2), 'extfaces')
    CALL HDF5_array1D_saving_int(file_id, Mesh%boundaryFlag, SIZE(Mesh%boundaryFlag, 1), 'boundaryFlag')
    ! Save simulation parameters
    CALL save_simulation_parameters()
    IF (switch%shockcp .EQ. 3) THEN
       CALL HDF5_array2D_saving(file_id, Mesh%scdiff_nodes, SIZE(Mesh%scdiff_nodes, 1), SIZE(Mesh%scdiff_nodes, 2), 'scdiff_nodes')
    END IF

    IF(switch%saveMeshSol) THEN
       CALL HDF5_integer_saving(file_id,Mesh%Ndim,'Ndim')
       CALL HDF5_integer_saving(file_id,Mesh%Nnodes,'Nnodes')
       IF(ASSOCIATED(Mesh%X_P1)) THEN
          CALL HDF5_integer_saving(file_id,SIZE(Mesh%X_P1,1),'Nnodes_P1')
       ENDIF
       CALL HDF5_integer_saving(file_id,Mesh%Nnodesperelem,'Nnodesperelem')
       CALL HDF5_integer_saving(file_id,Mesh%Nnodesperface,'Nnodesperface')
       CALL HDF5_integer_saving(file_id,Mesh%Nelems,'Nelems')
       CALL HDF5_integer_saving(file_id,Mesh%Nextfaces,'Nextfaces')
       CALL HDF5_integer_saving(file_id,Mesh%Nintfaces,'Nintfaces')
       CALL HDF5_integer_saving(file_id,Mesh%elemType,'elemType')
       CALL HDF5_integer_saving(file_id,Mesh%Ndir,'Ndir')
       CALL HDF5_integer_saving(file_id,Mesh%ukf,'ukf')
       CALL HDF5_array2D_saving_int(file_id,Mesh%T, SIZE(Mesh%T, 1), SIZE(Mesh%T, 2), 'T')
       IF(ASSOCIATED(Mesh%T_gmsh)) THEN
          CALL HDF5_array2D_saving_int(file_id,Mesh%T_gmsh, SIZE(Mesh%T_gmsh, 1), SIZE(Mesh%T_gmsh, 2), 'T_gmsh')
       ENDIF
       IF(ASSOCIATED(Mesh%Tb_gmsh)) THEN
          CALL HDF5_array2D_saving_int(file_id,Mesh%Tb_gmsh, SIZE(Mesh%Tb_gmsh, 1), SIZE(Mesh%Tb_gmsh, 2), 'Tb_gmsh')
       ENDIF
       IF(ASSOCIATED(Mesh%X_P1)) THEN
          CALL HDF5_array2D_saving(file_id,Mesh%X_P1, SIZE(Mesh%X_P1, 1), SIZE(Mesh%X_P1, 2), 'X_P1')
       ENDIF
       IF(ASSOCIATED(Mesh%Tlin)) THEN
          CALL HDF5_array2D_saving_int(file_id,Mesh%Tlin, SIZE(Mesh%Tlin, 1), SIZE(Mesh%Tlin, 2), 'Tlin')
       ENDIF
       CALL HDF5_array2D_saving_int(file_id,Mesh%Tb, SIZE(Mesh%Tb, 1), SIZE(Mesh%Tb, 2), 'Tb')
       !CALL HDF5_array1D_saving_int(file_id,Mesh%boundaryFlag, SIZE(Mesh%boundaryFlag), 'boundaryFlag')
       CALL HDF5_array2D_saving_int(file_id,Mesh%F, SIZE(Mesh%F,1),SIZE(Mesh%F,2), 'F')
       CALL HDF5_array2D_saving_int(file_id,Mesh%N, SIZE(Mesh%N,1),SIZE(Mesh%N,2), 'N')
       IF(ALLOCATED(Mesh%face_info)) THEN
          CALL HDF5_array2D_saving_int(file_id,Mesh%face_info, SIZE(Mesh%face_info,1),SIZE(Mesh%face_info,2), 'face_info')
       ENDIF
       IF(ALLOCATED(Mesh%faces)) THEN
          CALL HDF5_array3D_saving(file_id,REAL(Mesh%faces), SIZE(Mesh%faces,1),SIZE(Mesh%faces,2),SIZE(Mesh%faces,3), 'faces')
       ENDIF
       CALL HDF5_array2D_saving_int(file_id,Mesh%intfaces, SIZE(Mesh%intfaces,1),SIZE(Mesh%intfaces,2), 'intfaces')
       !call HDF5_array2D_saving_int(file_id,int(Mesh%flipface), SIZE(Mesh%flipface,1),SIZE(Mesh%flipface,2), 'flipface')
       !call HDF5_array2D_saving_logical(file_id,Mesh%Fdir, SIZE(Mesh%Fdir,1),SIZE(Mesh%Fdir,2), 'Fdir')
       IF(ALLOCATED(Mesh%periodic_faces)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%periodic_faces, SIZE(Mesh%periodic_faces), 'periodic_faces')
       ENDIF
       IF(ALLOCATED(Mesh%Diric)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%Diric, SIZE(Mesh%Diric), 'Diric')
       ENDIF
       IF(ALLOCATED(Mesh%numberbcs)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%numberbcs, SIZE(Mesh%numberbcs), 'numberbcs')
       ENDIF
       CALL HDF5_array1D_saving(file_id,Mesh%elemSize,SIZE(Mesh%elemSize), 'elemSize')
       CALL HDF5_array2D_saving(file_id,Mesh%X*phys%lscale, SIZE(Mesh%X, 1), SIZE(Mesh%X, 2), 'X')
#ifdef TOR3D
       CALL HDF5_integer_saving(file_id,Mesh%Nnodes_toroidal,'Nnodes_toroidal')
       CALL HDF5_array1D_saving(file_id,Mesh%toroidal,SIZE(Mesh%toroidal), 'toroidal')
#endif
       IF(ALLOCATED(Mesh%flag_elems_rho)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%flag_elems_rho,SIZE(Mesh%flag_elems_rho), 'flag_elems_rho')
       ENDIF
       IF(ALLOCATED(Mesh%flag_elems_sc)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%flag_elems_sc,SIZE(Mesh%flag_elems_sc), 'flag_elems_sc')
       ENDIF
       IF(ALLOCATED(Mesh%minrho_elems)) THEN
          CALL HDF5_array1D_saving(file_id,Mesh%minrho_elems,SIZE(Mesh%minrho_elems), 'minrho_elems')
       ENDIF
       IF(ALLOCATED(Mesh%sour_elems)) THEN
          CALL HDF5_array1D_saving(file_id,Mesh%sour_elems,SIZE(Mesh%sour_elems), 'sour_elems')
       ENDIF
       IF(ALLOCATED(Mesh%diff_elems)) THEN
          CALL HDF5_array1D_saving(file_id,Mesh%diff_elems,SIZE(Mesh%diff_elems), 'diff_elems')
       ENDIF
       IF(ALLOCATED(Mesh%scdiff_nodes)) THEN
          CALL HDF5_array2D_saving(file_id,Mesh%scdiff_nodes,SIZE(Mesh%scdiff_nodes,1),SIZE(Mesh%scdiff_nodes,2), 'scdiff_nodes')
       ENDIF
       CALL HDF5_real_saving(file_id, Mesh%xmax, 'xmax')
       CALL HDF5_real_saving(file_id, Mesh%xmin, 'xmin')
       CALL HDF5_real_saving(file_id, Mesh%ymax, 'ymax')
       CALL HDF5_real_saving(file_id, Mesh%ymin, 'ymin')
       CALL HDF5_real_saving(file_id, Mesh%puff_area, 'puff_area')
       CALL HDF5_real_saving(file_id, Mesh%core_area, 'core_area')
#ifdef PARALL
       CALL HDF5_integer_saving(file_id,Mesh%nghostfaces,'nghostfaces')
       CALL HDF5_integer_saving(file_id,Mesh%nghostelems,'nghostelems')
       CALL HDF5_integer_saving(file_id,Mesh%Nel_glob,'Nel_glob')
       CALL HDF5_integer_saving(file_id,Mesh%Nfa_glob,'Nfa_glob')
       CALL HDF5_integer_saving(file_id,Mesh%Nno_glob,'Nno_glob')
       CALL HDF5_integer_saving(file_id,Mesh%Ndir_glob,'Ndir_glob')
       CALL HDF5_integer_saving(file_id,Mesh%Ngho_glob,'Ngho_glob')

       IF(ASSOCIATED(Mesh%ghostflp)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%ghostflp,SIZE(Mesh%ghostflp), 'ghostflp')
       ENDIF
       IF(ASSOCIATED(Mesh%ghostloc)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%ghostloc,SIZE(Mesh%ghostloc), 'ghostloc')
       ENDIF
       IF(ASSOCIATED(Mesh%ghostpro)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%ghostpro,SIZE(Mesh%ghostpro), 'ghostpro')
       ENDIF
       IF(ASSOCIATED(Mesh%ghelsloc)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%ghelsloc,SIZE(Mesh%ghelsloc), 'ghelsloc')
       ENDIF
       IF(ASSOCIATED(Mesh%ghelspro)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%ghelspro,SIZE(Mesh%ghelspro), 'ghelspro')
       ENDIF
       IF(ALLOCATED(Mesh%fc2sd)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%fc2sd,SIZE(Mesh%fc2sd), 'fc2sd')
       ENDIF
       IF(ALLOCATED(Mesh%pr2sd)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%pr2sd,SIZE(Mesh%pr2sd), 'pr2sd')
       ENDIF
       IF(ALLOCATED(Mesh%fc2rv)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%fc2rv,SIZE(Mesh%fc2rv), 'fc2rv')
       ENDIF
       IF(ALLOCATED(Mesh%pr2rv)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%pr2rv,SIZE(Mesh%pr2rv), 'pr2rv')
       ENDIF
       IF(ALLOCATED(Mesh%el2sd)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%el2sd,SIZE(Mesh%el2sd), 'el2sd')
       ENDIF
       IF(ALLOCATED(Mesh%pe2sd)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%pe2sd,SIZE(Mesh%pe2sd), 'pe2sd')
       ENDIF
       IF(ALLOCATED(Mesh%el2rv)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%el2rv,SIZE(Mesh%el2rv), 'el2rv')
       ENDIF
       IF(ALLOCATED(Mesh%pe2rv)) THEN
          CALL HDF5_array1D_saving_int(file_id,Mesh%pe2rv,SIZE(Mesh%pe2rv), 'pe2rv')
       ENDIF
#endif

    ENDIF

    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    IF (MPIvar%glob_id .EQ. 0) THEN
       PRINT *, 'Output written to file ', TRIM(ADJUSTL(fname_complete))
       PRINT *, '        '
    END IF

  CONTAINS

    !**********************************************************************
    ! Save simulation parameters
    !**********************************************************************
    SUBROUTINE save_simulation_parameters()
      INTEGER(HID_T) :: group_id1, group_id2

      ! Create simulation parameters group
      CALL HDF5_group_create('simulation_parameters', file_id, group_id1, ierr)
      ! Save model definition
      CALL HDF5_string_saving(group_id1, simpar%model, 'model')
      ! Save Ndim and Neq
      CALL HDF5_integer_saving(group_id1, simpar%Ndim, 'Ndim')
      CALL HDF5_integer_saving(group_id1, simpar%Neq, 'Neq')
      ! Create adimensionalization subgroup
      CALL HDF5_group_create('adimensionalization', group_id1, group_id2, ierr)
      ! Save adimensionalization
      CALL HDF5_string_saving(group_id2, simpar%refval_time_dimensions, 'time_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_mass_dimensions, 'mass_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_length_dimensions, 'length_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_temperature_dimensions, 'temperature_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_density_dimensions, 'density_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_neutral_dimensions, 'density_neutral_dimensions')
#ifdef KEQUATION
      CALL HDF5_string_saving(group_id2, simpar%refval_k_dimensions, 'density_k_dimensions')
#endif
      CALL HDF5_string_saving(group_id2, simpar%refval_speed_dimensions, 'speed_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_potential_dimensions, 'potential_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_vorticity_dimensions, 'vorticity_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_magfield_dimensions, 'magfield_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_current_dimensions, 'current_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_diffusion_dimensions, 'diffusion_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_momentum_dimensions, 'momentum_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_specpress_dimensions, 'specific_pressure_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_specenergy_dimensions, 'specific_energy_scale_dimensions')
      CALL HDF5_string_saving(group_id2, simpar%refval_specenergydens_dimensions, 'specific_energy_density_scale_dimensions')
      CALL HDF5_real_saving(group_id2, simpar%refval_length, 'length_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_time, 'time_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_mass, 'mass_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_temperature, 'temperature_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_density, 'density_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_neutral, 'neutral_scale')
#ifdef KEQUATION
      CALL HDF5_real_saving(group_id2, simpar%refval_k, 'k_scale')
#endif
      CALL HDF5_real_saving(group_id2, simpar%refval_speed, 'speed_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_potential, 'potential_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_vorticity, 'vorticity_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_magfield, 'magfield_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_current, 'current_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_diffusion, 'diffusion_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_momentum, 'momentum_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_specpress, 'specific_pressure_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_specenergy, 'specific_energy_scale')
      CALL HDF5_real_saving(group_id2, simpar%refval_specenergydens, 'specific_energy_density_scale')
      CALL HDF5_array1d_saving(group_id2, simpar%physvar_refval, phys%npv, 'reference_values_physical_variables')
      CALL HDF5_array1d_saving(group_id2, simpar%consvar_refval, phys%Neq, 'reference_values_conservative_variables')
      CALL HDF5_group_close(group_id2, ierr)
      ! Close group adimensionalization

      ! Create physics parameters group
      CALL HDF5_group_create('physics', group_id1, group_id2, ierr)
      CALL HDF5_string_array1D_saving(group_id2, phys%phyVarNam, 'physical_variable_names')
      CALL HDF5_string_array1D_saving(group_id2, phys%conVarNam, 'conservative_variable_names')
      CALL HDF5_real_saving(group_id2, phys%a, 'a')
      CALL HDF5_real_saving(group_id2, phys%Mref, 'Mref')
      CALL HDF5_real_saving(group_id2, phys%c1, 'c1')
      CALL HDF5_real_saving(group_id2, phys%c2, 'c2')
      CALL HDF5_real_saving(group_id2, phys%diff_pari, 'diff_pari')
      CALL HDF5_real_saving(group_id2, phys%diff_pare, 'diff_pare')
      CALL HDF5_real_saving(group_id2, phys%etapar, 'eta_parallel')
      CALL HDF5_real_saving(group_id2, phys%diff_n, 'diff_n')
      CALL HDF5_real_saving(group_id2, phys%diff_u, 'diff_u')
      CALL HDF5_real_saving(group_id2, phys%diff_e, 'diff_e')
      CALL HDF5_real_saving(group_id2, phys%diff_ee, 'diff_ee')
      CALL HDF5_real_saving(group_id2, phys%diff_vort, 'diff_vort')
      CALL HDF5_real_saving(group_id2, phys%diff_pot, 'diff_pot')
      CALL HDF5_real_saving(group_id2, phys%diff_nn, 'diff_nn')
      CALL HDF5_real_saving(group_id2, phys%Re, 'recycling')
      CALL HDF5_real_saving(group_id2, phys%Re_pump, 'recycling_pump')
      CALL HDF5_real_saving(group_id2, phys%puff, 'puff')
      IF (switch%ME) THEN
         CALL HDF5_array1d_saving(group_id2, phys%puff_exp, time%nts, 'puff_exp')
      END IF
      CALL HDF5_real_saving(group_id2, phys%tie, 'tau_ie')
      CALL HDF5_real_saving(group_id2, phys%dfcoef, 'dfcoef')
      CALL HDF5_real_saving(group_id2, phys%dexbcoef, 'dexbcoef')
      CALL HDF5_real_saving(group_id2, phys%bohmth, 'bohmth')
      CALL HDF5_real_saving(group_id2, phys%epn, 'epn')
      CALL HDF5_real_saving(group_id2, phys%Gmbohm, 'Gmbohm')
      CALL HDF5_real_saving(group_id2, phys%Gmbohme, 'Gmbohme')
      CALL HDF5_real_saving(group_id2, phys%Potfloat, 'Potfloat')
      CALL HDF5_array1d_saving_int(group_id2, phys%bcflags, 10, 'boundary_flags')
      CALL HDF5_group_close(group_id2, ierr)

      ! Create switches parameters group
      CALL HDF5_group_create('switches', group_id1, group_id2, ierr)
      CALL HDF5_logical_saving(group_id2, switch%driftdia, 'diamagnetic_drift')
      CALL HDF5_logical_saving(group_id2, switch%driftexb, 'ExB_drift')
      CALL HDF5_logical_saving(group_id2, switch%steady, 'steady')
      CALL HDF5_integer_saving(group_id2, switch%testcase, 'testcase')
      CALL HDF5_logical_saving(group_id2, switch%ohmicsrc, 'ohmicsrc')
      CALL HDF5_logical_saving(group_id2, switch%ME, 'ME')
      CALL HDF5_logical_saving(group_id2, switch%rmp, 'RMP')
      CALL HDF5_logical_saving(group_id2, switch%ripple, 'Ripple')
      CALL HDF5_logical_saving(group_id2, switch%psdtime, 'psdtime')
      CALL HDF5_real_saving(group_id2, switch%diffred, 'diffred')
      CALL HDF5_real_saving(group_id2, switch%diffmin, 'diffmin')
      CALL HDF5_integer_saving(group_id2, switch%shockcp, 'shockcp')
      CALL HDF5_integer_saving(group_id2, switch%limrho, 'limrho')
      CALL HDF5_integer_saving(group_id2, switch%difcor, 'difcor')
      CALL HDF5_integer_saving(group_id2, switch%thresh, 'thresh')
      CALL HDF5_logical_saving(group_id2, switch%filter, 'filter')
      CALL HDF5_logical_saving(group_id2, switch%decoup, 'decoup')
      CALL HDF5_logical_saving(group_id2, switch%ckeramp, 'ckeramp')
      CALL HDF5_logical_saving(group_id2, switch%saveNR, 'saveNR')
      CALL HDF5_logical_saving(group_id2, switch%saveTau, 'saveTau')
      CALL HDF5_logical_saving(group_id2, switch%fixdPotLim, 'fixdPotLim')
      CALL HDF5_logical_saving(group_id2, switch%dirivortcore, 'dirivortcore')
      CALL HDF5_logical_saving(group_id2, switch%dirivortlim, 'dirivortlim')
      CALL HDF5_logical_saving(group_id2, switch%convvort, 'convvort')
      CALL HDF5_logical_saving(group_id2, switch%logrho, 'logrho')
      CALL HDF5_group_close(group_id2, ierr)

      ! Create numerics parameters group
      CALL HDF5_group_create('numerics', group_id1, group_id2, ierr)
      CALL HDF5_integer_saving(group_id2, numer%nrp, 'Max_number_of_NR_iterations')
      CALL HDF5_real_saving(group_id2, numer%tnr, 'NR_convergence_criterium')
      CALL HDF5_real_saving(group_id2, numer%ttm, 'Time_convergence_criterium')
      CALL HDF5_real_saving(group_id2, numer%div, 'Divergence_criterium')
      CALL HDF5_array1d_saving(group_id2, numer%tau, 4, 'Stabilization_parameter')
      CALL HDF5_real_saving(group_id2, numer%sc_coe, 'Shock_capturing_parameter')
      CALL HDF5_real_saving(group_id2, numer%sc_sen, 'Shock_capturing_sensibility')
      CALL HDF5_real_saving(group_id2, numer%minrho, 'Value_of_rho_to_start_applying_limiting')
      CALL HDF5_real_saving(group_id2, numer%so_coe, 'Source_coefficient_for_limiting_rho')
      CALL HDF5_real_saving(group_id2, numer%df_coe, 'Diffusion_coefficient_for_limiting_rho')
      CALL HDF5_real_saving(group_id2, numer%dc_coe, 'Diffusion_coefficient_in_corners')
      CALL HDF5_real_saving(group_id2, numer%thr, 'Threshold_to_limit_rho')
      CALL HDF5_real_saving(group_id2, numer%thrpre, 'Threshold_to_limit_pressure')
      CALL HDF5_integer_saving(group_id2, numer%stab, 'Stabilization_type')
      CALL HDF5_real_saving(group_id2, numer%dumpnr, 'dumping_factor_for_Newton_Raphson')
      CALL HDF5_integer_saving(group_id2, numer%ntor, 'Number_of_elements_in_the_toroidal_direction')
      CALL HDF5_integer_saving(group_id2, numer%ptor, 'Polynomial_degree_in_the_toroidal_direction')
      CALL HDF5_real_saving(group_id2, numer%tmax, 'Max_extention_in_the_toroidal_direction')
      CALL HDF5_integer_saving(group_id2, numer%npartor, 'Number_of_MPI_divisions_in_the_toroidal_direction')
      CALL HDF5_real_saving(group_id2, numer%exbdump, 'Dumping_for_ExB_drift')
      CALL HDF5_group_close(group_id2, ierr)

      ! Create time parameters group
      CALL HDF5_group_create('time', group_id1, group_id2, ierr)
      CALL HDF5_real_saving(group_id2, time%dt0, 'Initial_time_step')
      CALL HDF5_real_saving(group_id2, time%dt, 'Current_time_step')
      CALL HDF5_real_saving(group_id2, time%tfi, 'Final_time')
      CALL HDF5_integer_saving(group_id2, time%it, 'Current_time_step_number')
      CALL HDF5_integer_saving(group_id2, time%ik, 'Current_pseudo_time_step_number')
      CALL HDF5_integer_saving(group_id2, time%nts, 'Number_of_time_steps')
      CALL HDF5_integer_saving(group_id2, time%tis, 'Time_integration_scheme')
      CALL HDF5_real_saving(group_id2, time%t, 'Current_time')
      CALL HDF5_group_close(group_id2, ierr)

      ! Create geometry parameters group
      CALL HDF5_group_create('geometry', group_id1, group_id2, ierr)
      CALL HDF5_real_saving(group_id2, geom%R0, 'Major_radius')
      CALL HDF5_real_saving(group_id2, geom%q, 'Safety_factor')
      CALL HDF5_group_close(group_id2, ierr)

      ! Create time parameters group
      CALL HDF5_group_create('magnetic', group_id1, group_id2, ierr)
      CALL HDF5_real_saving(group_id2, magn%amp_rmp, 'Amplitude_RMP')
      CALL HDF5_integer_saving(group_id2, magn%nbCoils_rmp, 'Number_coils_RMP')
      CALL HDF5_integer_saving(group_id2, magn%parite, 'Parity_RMP')
      CALL HDF5_integer_saving(group_id2, magn%nbRow, 'number_rows_RMP')
      CALL HDF5_real_saving(group_id2, magn%amp_ripple, 'Amplitude_Ripple')
      CALL HDF5_integer_saving(group_id2, magn%nbCoils_ripple, 'Number_coils_Ripple')
      CALL HDF5_real_saving(group_id2, magn%triang, 'Triangularity')
      CALL HDF5_real_saving(group_id2, magn%ellip, 'Ellipticity')
      CALL HDF5_group_close(group_id2, ierr)

      CALL HDF5_group_close(group_id1, ierr)


    END SUBROUTINE save_simulation_parameters

  END SUBROUTINE HDF5_save_solution


  SUBROUTINE HDF5_load_mesh_from_solution(fname)
    !*************************************
    !              2D case
    !*************************************

    USE MPI_OMP
    CHARACTER(LEN=*) :: fname
    CHARACTER(len=1000) :: fname_complete
    CHARACTER(10)  :: str
    CHARACTER(70)  :: npr, nid
    REAL*8, PARAMETER::tol = 1e-6
    REAL*8 :: xmin
    INTEGER :: elemType, ndim, Nnodes, Nelems, Nnodesperelem, Nnodes_P1
    INTEGER :: Nextfaces, Nnodesperface, IERR
    INTEGER(HID_T) :: file_id
#ifdef TOR3D
    INTEGER         :: i
#endif
#ifdef PARALL
    INTEGER :: ghfa, ghel, Nel_glob, Nfa_glob, Ndir_glob, Ngho_glob, Nfaces
#endif


#ifdef TOR3D
#ifdef PARALL
    IF (MPIvar%npol .GT. 1) THEN
       WRITE (nid, *) MPIvar%ipol
       WRITE (npr, *) MPIvar%npol
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
#else
    fname_complete = TRIM(ADJUSTL(fname))//'.h5'
#endif
#else
    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
#endif
    IF (utils%printint > 0) THEN
       PRINT *, 'Loading mesh.'
       PRINT *, '        '
    ENDIF

    CALL HDF5_open(fname_complete, file_id, IERR)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error opening mesh from solution file: ", fname_complete
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, elemType, 'elemType', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: elemType"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, ndim, 'Ndim', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Ndim"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodes, 'Nnodes', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodes"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodes_P1, 'Nnodes_P1', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodes_P1"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nelems, 'Nelems', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nelems"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperelem, 'Nnodesperelem', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperelem"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nnodesperface, 'Nnodesperface', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperface"
       STOP
    ENDIF
    CALL HDF5_integer_reading(file_id, Nextfaces, 'Nextfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nextfaces"
       STOP
    ENDIF
#ifdef PARALL
    CALL HDF5_integer_reading(file_id, Nfaces, 'Nfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nfaces"
       STOP
    ENDIF
#endif
    ALLOCATE (Mesh%T(Nelems, Nnodesperelem))
    IF(elemType .EQ. 0) THEN
       ALLOCATE (Mesh%T_gmsh(Nelems, 5 + 3))
       ALLOCATE (Mesh%Tb_gmsh(Nextfaces, 5 + 2))
    ELSE
       ALLOCATE (Mesh%T_gmsh(Nelems, 5 + 4))
       ALLOCATE (Mesh%Tb_gmsh(Nextfaces, 5 + 2))
    ENDIF
    ALLOCATE (Mesh%X(Nnodes, ndim))
    ALLOCATE (Mesh%X_P1(Nnodes_P1,ndim+2))
    ALLOCATE (Mesh%Tb(Nextfaces, Nnodesperface))
    ALLOCATE (Mesh%boundaryFlag(Nextfaces))
#ifdef PARALL
    ALLOCATE (Mesh%ghostFaces(Nfaces))
    ALLOCATE (Mesh%loc2glob_fa(Nfaces))
    ALLOCATE (Mesh%loc2glob_el(Nelems))
    ALLOCATE (Mesh%ghostElems(Nelems))
#endif
    CALL HDF5_array2D_reading_int(file_id, Mesh%T, 'T', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading mesh connectivity T"
       STOP
    ENDIF
    CALL HDF5_array2D_reading_int(file_id, Mesh%Tb, 'Tb', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundary connectivity Tb"
       STOP
    ENDIF
    CALL HDF5_array2D_reading_int(file_id, Mesh%Tb_gmsh, 'Tb_gmsh', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundary connectivity Tb_gmsh"
       STOP
    ENDIF
    CALL HDF5_array2D_reading_int(file_id, Mesh%T_gmsh, 'T_gmsh', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundary connectivity T_gmsh"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%boundaryFlag, 'boundaryFlag', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundaryFlag"
       STOP
    ENDIF
    CALL HDF5_array2D_reading(file_id, Mesh%X, 'X', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading coordinate matrix X"
       STOP
    ENDIF
    CALL HDF5_array2D_reading(file_id, Mesh%X_P1, 'X_P1', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading coordinate matrix X_P1"
       STOP
    ENDIF
#ifdef PARALL
    CALL HDF5_array1D_reading_int(file_id, Mesh%loc2glob_fa, 'loc2glob_fa', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading loc2glob_fa"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%loc2glob_el, 'loc2glob_el', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading loc2glob_el"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostFaces, 'ghostFaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostFaces"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostElems, 'ghostElems', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostElems"
       STOP
    ENDIF
    ! Find the number of ghost faces
    ghfa = SUM(Mesh%ghostFaces)
    Mesh%nghostfaces = ghfa

    ! Find the number of ghost elements
    ghel = SUM(Mesh%ghostElems)
    Mesh%nghostElems = ghel

    ALLOCATE (Mesh%ghostflp(ghfa))
    ALLOCATE (Mesh%ghostpro(ghfa))
    ALLOCATE (Mesh%ghostloc(ghfa))
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostflp, 'ghostFlp', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostFlp"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostLoc, 'ghostLoc', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostLoc"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghostPro, 'ghostPro', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghostPro"
       STOP
    ENDIF

#ifdef TOR3D
    IF (MPIvar%ntor > 1) THEN
       DO i = 1, SIZE(Mesh%ghostPro)
          IF (Mesh%ghostPro(i) .GT. -1) THEN
             Mesh%ghostPro(i) = Mesh%ghostPro(i) + (MPIvar%itor - 1)*MPIvar%npol
          ENDIF
       END DO
    ENDIF
    ALLOCATE (Mesh%ghelspro(ghel))
    ALLOCATE (Mesh%ghelsloc(ghel))
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghelsLoc, 'ghelsLoc', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghelsLoc"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(file_id, Mesh%ghelsPro, 'ghelsPro', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading ghelsPro"
       STOP
    ENDIF
    IF (MPIvar%ntor .GT. 1) THEN
       DO i = 1, SIZE(Mesh%ghelspro)
          IF (Mesh%ghelspro(i) .GT. -1) THEN
             Mesh%ghelsPro(i) = Mesh%ghelsPro(i) + (MPIvar%itor - 1)*MPIvar%npol
          END IF
       END DO
    END IF
#endif
#endif
    CALL HDF5_close(file_id)

    !************************************************************************
    !   CONFIRMATION MESSAGE FOR THE USER
    !************************************************************************
#ifdef PARALL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    WRITE (6, *) "Process: ", MPIvar%glob_id, "-- readed mesh file: ", TRIM(ADJUSTL(fname_complete))
#else
    WRITE (6, *) "Readed mesh file: ", TRIM(ADJUSTL(fname_complete))
#endif

#ifdef PARALL
    CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_el), Nel_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_fa), Nfa_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%ndir, Ndir_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%nghostfaces, Ngho_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    Mesh%Nel_glob = Nel_glob
    Mesh%Nfa_glob = Nfa_glob
    Mesh%Ndir_glob = Ndir_glob
    Mesh%Ngho_glob = Ngho_glob
#endif
    Mesh%Ndim = ndim
    Mesh%Nnodes = Nnodes
    Mesh%Nelems = Nelems
    Mesh%Nnodesperelem = Nnodesperelem
    Mesh%Nnodesperface = Nnodesperface
    Mesh%elemType = elemType
    Mesh%Nextfaces = Nextfaces

    xmin = MINVAL(Mesh%X(:,1))
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
    ! Apply shift if axisymmetric case
    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) .OR. (switch%axisym .AND. xmin < tol)) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) "*** Applying translation in axisymmetric case!"
       ENDIF
       Mesh%X(:, 1) = Mesh%X(:, 1) + geom%R0
    END IF

    ! Apply length scale
    Mesh%X = Mesh%X/phys%lscale

    Mesh%xmax = MAXVAL(Mesh%X(:, 1))
    Mesh%xmin = MINVAL(Mesh%X(:, 1))
    Mesh%ymax = MAXVAL(Mesh%X(:, 2))
    Mesh%ymin = MINVAL(Mesh%X(:, 2))

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    IF (utils%printint > 0) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          IF (elemType == 0) THEN
             WRITE (str, '(A)') 'triangles'
          ELSEIF (elemType == 1) THEN
             WRITE (str, '(A)') 'quads'
          ELSEIF (elemType == 2) THEN
             WRITE (str, '(A)') 'thetra'
          ELSEIF (elemType == 3) THEN
             WRITE (str, '(A)') 'hexa'
          END IF
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*                    MESH                       *'
          WRITE (6, *) '*************************************************'
          WRITE (6, '(A,I18)') ' Number of dimensions:         ', ndim
          WRITE (6, '(A,A34)') ' Element type: ', TRIM(str)
          WRITE (6, '(A,I18)') ' Number of elements:           ', Nelems
          WRITE (6, '(A,I18)') ' Number of nodes:              ', Nnodes
          WRITE (6, '(A,I18)') ' Number of nodes per element:  ', Nnodesperelem
          WRITE (6, '(A,I18)') ' Number of nodes per face:     ', Nnodesperface
          WRITE (6, '(A,I18)') ' Number of exterior faces:     ', Nextfaces
          WRITE (6, *) ' '
          WRITE (6, *) ' '
          IF (utils%printint > 1) THEN
             WRITE (6, *) "Connectivity matrix T:"
             CALL displayMatrixInt(Mesh%T)
             WRITE (6, *) "Boundary connectivity matrix Tb:"
             CALL displayMatrixInt(Mesh%Tb)
          END IF
       ENDIF
    ENDIF

  endsubroutine HDF5_load_mesh_from_solution

  !**********************************************************************
  ! Load solution in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_load_solution(fname)
    USE globals
    USE LinearAlgebra, ONLY: tensorsumint, colint,col
    IMPLICIT NONE

    CHARACTER(LEN=1000) :: fname
    CHARACTER(len=100), POINTER :: mod_ptr
    CHARACTER(len=100), TARGET :: model_string

    CHARACTER(70)  :: npr, nid
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id, group_id
#ifndef TOR3D
    INTEGER(HID_T) ::group_id2
#endif

#ifdef TOR3D
    CHARACTER(70)  :: nip, nit, ngd
    INTEGER :: ntorloc, iel, ifa, iface, Fi, itor, dd, iel3, i, j, N2d, Nfl, Np1d, Np2d
    INTEGER, ALLOCATABLE :: ind_q_add(:), indu2D(:), indq2D(:), indu3D(:), indq3D(:), indufp(:), induf2D(:), induf3D(:), indul(:), indql(:), indutl(:)
    REAL*8, ALLOCATABLE    :: u2D(:), q2D(:), ut2D(:), q2D_add(:, :)
    REAL, POINTER:: u_aux(:), u_tilde_aux(:), q_aux(:)
#else
    REAL*8               :: t
#endif
    INTEGER :: Neq, Ndim, Nel, Np, Nfg, Nf, sizeutilde, sizeu, it
    REAL*8, ALLOCATABLE       :: uaux(:,:),utaux(:,:),qaux(:,:)
    INTEGER              :: logrho_ptr = 0

    Neq = phys%Neq
    mod_ptr => model_string
#ifdef TOR3D
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       ntorloc = numer%ntor/MPIvar%ntor + 1
    ELSE
       ntorloc = numer%ntor
    ENDIF
#else
    ntorloc = numer%ntor
#endif
    Ndim = 3                             ! N. of dimensions
    N2d = Mesh%Nelems                   ! N. of 2D elements
    Nel = N2d*ntorloc                      ! N. of 3D elements
    Np1d = refElTor%Nnodes1D             ! N. of nodes for each toroidal 1d element
    Np2d = refElPol%Nnodes2D             ! N. of nodes for each poloidal 2D element
    Np = Np2d*Np1d                     ! N. of nodes for each 3D element
    Nfl = refElPol%Nnodes1D*Np1d        ! N. of nodes in the lateral faces
    Nfg = Np2d*2 + refElPol%Nfaces*Nfl    ! N. of nodes in all the faces of a 3D element
    Nf = Mesh%Nfaces                   ! N. of faces in the 2D mesh
    sizeu = Neq*Nel*Np                    ! Size of u
#ifdef PARALL
    IF (MPIvar%ntor .GT. 1) THEN
       sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d) + Neq*Np2d*N2d! Size of utilde
    ELSE
       sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
    ENDIF
#else
    sizeutilde = Neq*ntorloc*(Nfl*Nf + Np2d*N2d)! Size of utilde
#endif
#else
    Ndim = 2
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfg = refElPol%Nfacenodes*Nf
    sizeu = Neq*Nel*Np
    sizeutilde = Neq*Mesh%Nfaces*Mesh%Nnodesperface
#endif

    ALLOCATE (sol%u(sizeu))
    ALLOCATE (sol%u_tilde(sizeutilde))
    ALLOCATE (sol%q(sizeu*Ndim))

#ifdef TOR3D
    !*************************************
    !              3D case
    !*************************************

    IF ((fname(1:5) == 'Sol3D') .OR. (fname(1:7) == './Sol3D')) THEN
       ! Initialization with a 3D solution
       WRITE (6, *) "3D initial solution"

       IF (MPIvar%glob_size .GT. 1) THEN
          WRITE (nip, *) MPIvar%ipol
          WRITE (nit, *) MPIvar%itor
          WRITE (ngd, *) MPIvar%glob_size
          fname_complete = TRIM(ADJUSTL(fname))//'_ip'//TRIM(ADJUSTL(nip))//'_it'//TRIM(ADJUSTL(nit))//'_np'//TRIM(ADJUSTL(ngd))//'.h5'
       ELSE
          fname_complete = TRIM(ADJUSTL(fname))//'.h5'
       END IF

       CALL HDF5_open(fname_complete, file_id, IERR)
       CALL HDF5_group_open(file_id, 'simulation_parameters', group_id, ierr)
       CALL HDF5_string_reading(group_id, mod_ptr, 'model')
#ifndef TEMPERATURE
       CALL HDF5_group_open(group_id, 'switches', group_id2, ierr)
       CALL HDF5_integer_reading(group_id2, logrho_ptr, 'logrho')
       CALL HDF5_group_close(group_id2, ierr)
#else
       IF(switch%logrho) THEN
          logrho_ptr = 1
       ELSE
          logrho_ptr = 0
       ENDIF
#endif
       CALL HDF5_group_close(group_id, ierr)
       ! Check if the readed solution corresponds to the right model
       IF (simpar%model .NE. model_string) THEN
          WRITE (6, *) "Wrong model in loaded solution | Loaded model: ", model_string, " | Current model: ", simpar%model
          STOP
       ENDIF

       CALL HDF5_array1D_reading(file_id, sol%u, 'u')
       CALL HDF5_array1D_reading(file_id, sol%u_tilde, 'u_tilde')
       CALL HDF5_integer_reading(file_id,time%it,'it')
       CALL HDF5_array1D_reading(file_id, sol%q, 'q')
       CALL HDF5_close(file_id)
    ELSEIF ((fname(1:5) == 'Sol2D' .OR. fname(1:7) == './Sol2D')) THEN
       ! Initialization with a 2D solution
       WRITE (6, *) "2D initial solution: propagating in the torus..."

       IF (MPIvar%glob_size .GT. 1) THEN
          WRITE (nid, *) MPIvar%ipol
          WRITE (npr, *) MPIvar%npol
          fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
       ELSE
          fname_complete = TRIM(ADJUSTL(fname))//'.h5'
       END IF

       CALL HDF5_open(fname_complete, file_id, IERR)
       CALL HDF5_group_open(file_id, 'simulation_parameters', group_id, ierr)
       CALL HDF5_string_reading(group_id, mod_ptr, 'model')
#ifndef TEMPERATURE
       CALL HDF5_group_open(group_id, 'switches', group_id2, ierr)
       CALL HDF5_integer_reading(group_id2, logrho_ptr, 'logrho')
       CALL HDF5_group_close(group_id2, ierr)
#else
       IF(switch%logrho) THEN
          logrho_ptr = 1
       ELSE
          logrho_ptr = 0
       ENDIF
#endif
       CALL HDF5_group_close(group_id, ierr)
       ! Check if the readed solution corresponds to the right model
       IF (simpar%model .NE. model_string) THEN
          WRITE (6, *) "Wrong model in loaded solution | Loaded model: ", model_string, " | Current model: ", simpar%model
          STOP
       ENDIF

       ALLOCATE (u_aux(Neq*N2d*Np))
       ALLOCATE (u_tilde_aux(Neq*Mesh%Nfaces*Mesh%Nnodesperface))
       ALLOCATE (u2D(Np2D*Neq))
       ALLOCATE (ut2D(refElPol%Nnodes1D*Neq))
       ALLOCATE (indu2D(Np2D*Neq))
       ALLOCATE (indu3D(Np*Neq))
       ALLOCATE (indufp(Np2D*Neq))
       ALLOCATE (induf2D(refElPol%Nnodes1D*Neq))
       ALLOCATE (induf3D(Nfl*Neq))
       ALLOCATE (indul(Np2D*Neq))
       ALLOCATE (indutl(refElPol%Nnodes1D*Neq))
       ALLOCATE (q_aux(Neq*N2d*Np*(Ndim - 1)))
       ALLOCATE (q2D(Np2D*Neq*(Ndim - 1)))
       ALLOCATE (q2D_add(Np2D, Ndim*Neq))
       ALLOCATE (indq2D(Np2D*Neq*(Ndim - 1)))
       ALLOCATE (indq3D(Np*Neq*Ndim))
       ALLOCATE (indql(Np2D*Neq*Ndim))
       ALLOCATE (ind_q_add(Neq*(Ndim - 1)))
       q2D_add = 0.
       ind_q_add = colint(tensorsumint((/(j, j=1, (ndim - 1))/), ndim*(/(i, i=0, (Neq - 1))/)))

       CALL HDF5_open(fname_complete, file_id, IERR)
       CALL HDF5_array1D_reading(file_id, u_aux, 'u')
       CALL HDF5_array1D_reading(file_id, u_tilde_aux, 'u_tilde')
       CALL HDF5_array1D_reading(file_id, q_aux, 'q')
       CALL HDF5_close(file_id)

#ifdef PARALL
       WRITE (6, *) "Process: ", MPIvar%glob_id, "-- readed solution file: ", TRIM(ADJUSTL(fname_complete))
#else
       WRITE (6, *) "Readed solution file: ", TRIM(ADJUSTL(fname_complete))
#endif

       DO iel = 1, N2D
          indu2D = (iel - 1)*Np2d*Neq + (/(i, i=1, Np2d*Neq)/)
          u2D = u_aux(indu2D)
          indq2D = (iel - 1)*Np2d*Neq*(Ndim - 1) + (/(i, i=1, Np2d*Neq*(Ndim - 1))/)
          q2D = q_aux(indq2D)
          q2d_add(:, ind_q_add) = TRANSPOSE(RESHAPE(q2D, [(Ndim - 1)*Neq, Np2d]))

          DO itor = 1, ntorloc
             iel3 = (itor - 1)*N2d+iel
             indu3D = (iel3 - 1)*Np*Neq + (/(i, i=1, Np*Neq)/)
             indq3D = (iel3 - 1)*Np*Neq*Ndim + (/(i, i=1, Np*Neq*Ndim)/)
             dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (iel - 1)*Np2D*Neq
             indufp = dd + (/(i, i=1, Np2D*Neq)/)
             sol%u_tilde(indufp) = u2d
             DO it = 1, Np1d
                indul = (it - 1)*Np2D*Neq + (/(i, i=1, Np2D*Neq)/)
                sol%u(indu3D(indul)) = u2D
                indql = (it - 1)*Np2D*Neq*Ndim + (/(i, i=1, Np2D*Neq*Ndim)/)
                sol%q(indq3D(indql)) = RESHAPE(TRANSPOSE(q2d_add), (/Np2d*Neq*Ndim/))
             END DO
          END DO
       END DO

       DO iface = 1, Mesh%Nintfaces
          Fi = iface
          induf2D = (Fi - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
          ut2d = u_tilde_aux(induf2D)
          DO itor = 1, ntorloc
             dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
             induf3D = dd + (/(i, i=1, Nfl*Neq)/)
             DO it = 1, Np1d
                indutl = (it - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
                sol%u_tilde(induf3D(indutl)) = ut2d
             END DO
          END DO
       END DO

       DO iface = 1, Mesh%Nextfaces
          iel = Mesh%extfaces(iface, 1)
          ifa = Mesh%extfaces(iface, 2)
          IF (Mesh%Fdir(iel, ifa)) CYCLE
          Fi = iface + Mesh%Nintfaces
          induf2D = (Fi - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
          ut2d = u_tilde_aux(induf2D)
          DO itor = 1, ntorloc
             dd = (itor - 1)*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (N2D*Np2D+(Fi - 1)*Nfl)*Neq
             induf3D = dd + (/(i, i=1, Nfl*Neq)/)
             DO it = 1, Np1d
                indutl = (it - 1)*refElPol%Nnodes1D*Neq + (/(i, i=1, refElPol%Nnodes1D*Neq)/)
                sol%u_tilde(induf3D(indutl)) = ut2d
             END DO
          END DO
       END DO

#ifdef PARALL
       ! Add solution on toroidal ghost faces
       IF (MPIvar%ntor .GT. 1) THEN

          DO iel = 1, N2D
             indu2D = (iel - 1)*Np2d*Neq + (/(i, i=1, Np2d*Neq)/)
             u2D = u_aux(indu2D)
             indq2D = (iel - 1)*Np2d*Neq*Ndim + (/(i, i=1, Np2d*Neq*Ndim)/)
             q2D = q_aux(indq2D)
             DO itor = 1, ntorloc
                dd = ntorloc*(N2D*Np2D+(Mesh%Nfaces - Mesh%Ndir)*Nfl)*Neq + (iel - 1)*Np2D*Neq
                indufp = dd + (/(i, i=1, Np2D*Neq)/)
                sol%u_tilde(indufp) = u2d
             END DO
          END DO
       ENDIF
#endif

       WRITE (6, *) "Done!"

       DEALLOCATE (u_aux, u_tilde_aux, u2D, ut2D, indu2D, indu3D, indufp, induf2D, induf3D, indul, indutl)
       DEALLOCATE (q_aux, q2D, q2D_add, indq2D, indq3D, indql, ind_q_add)
    END IF

#else
    !*************************************
    !              2D case
    !*************************************
    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    CALL HDF5_open(fname_complete, file_id, IERR)
    CALL HDF5_group_open(file_id, 'simulation_parameters', group_id, ierr)
    CALL HDF5_string_reading(group_id, mod_ptr, 'model')
#ifndef TEMPERATURE
    CALL HDF5_group_open(group_id, 'switches', group_id2, ierr)
    CALL HDF5_integer_reading(group_id2, logrho_ptr, 'logrho')
    CALL HDF5_group_close(group_id2, ierr)
#else
    IF(switch%logrho) THEN
       logrho_ptr = 1
    ELSE
       logrho_ptr = 0
    ENDIF
#endif
    IF (switch%ME) THEN
       CALL HDF5_group_open(group_id, 'time', group_id2, ierr)
       CALL HDF5_integer_reading(group_id2, it, 'Current_time_step_number')
       CALL HDF5_real_reading(group_id2, t, 'Current_time')
       IF (it .GT. 1) THEN
          time%it = it
          time%ik = it
          time%t = t
          sol%Nt = it
          sol%time(it) = t
       ELSE
       END IF
       CALL HDF5_group_close(group_id2, ierr)
       IF (time%it .NE. 0) THEN
          CALL HDF5_group_open(group_id, 'physics', group_id2, ierr)
          CALL HDF5_array1D_reading(group_id2, phys%puff_exp, 'puff_exp')
          CALL HDF5_group_close(group_id2, ierr)
       END IF
    END IF
    CALL HDF5_group_close(group_id, ierr)


    ! Check if the readed solution corresponds to the right model
    IF (simpar%model .NE. model_string) THEN
       WRITE (6, *) "Wrong model in loaded solution | Loaded model: ", model_string, " | Current model: ", simpar%model
       STOP
    ENDIF
    CALL HDF5_array1D_reading(file_id, sol%u, 'u')
    CALL HDF5_array1D_reading(file_id, sol%u_tilde, 'u_tilde')
    CALL HDF5_array1D_reading(file_id, sol%q, 'q')
    CALL HDF5_close(file_id)
#endif

    IF (switch%logrho .AND. logrho_ptr.EQ.0 ) THEN
       WRITE(6,*) "Readed solution without logrho but switch logrho set to true: "
       WRITE(6,*) "Computing log of density"
       ALLOCATE(uaux(SIZE(sol%u)/phys%neq,phys%neq))
       ALLOCATE(utaux(SIZE(sol%u_tilde)/phys%neq,phys%neq))
       ALLOCATE(qaux(SIZE(sol%q)/phys%neq/Ndim,phys%neq*Ndim))
       uaux = TRANSPOSE(RESHAPE(sol%u,[phys%neq,SIZE(sol%u)/phys%neq]))
       utaux = TRANSPOSE(RESHAPE(sol%u_tilde,[phys%neq,SIZE(sol%u_tilde)/phys%neq]))
       qaux = TRANSPOSE(RESHAPE(sol%q,[phys%neq*Ndim,SIZE(sol%q)/phys%neq/Ndim]))
       qaux(:,1) = qaux(:,1)/uaux(:,1)
       qaux(:,2) = qaux(:,2)/uaux(:,1)
       uaux(:,1) = LOG(uaux(:,1))
       utaux(:,1) = LOG(utaux(:,1))
       sol%u = col(TRANSPOSE(uaux))
       sol%u_tilde = col(TRANSPOSE(utaux))
       sol%q = col(TRANSPOSE(qaux))
       DEALLOCATE(uaux,utaux,qaux)
    ELSEIF (.NOT.switch%logrho .AND. logrho_ptr.EQ.1 ) THEN
       WRITE(6,*) "Readed solution with logrho but switch logrho set to false: "
       WRITE(6,*) "Computing exp of density"
       ALLOCATE(uaux(SIZE(sol%u)/phys%neq,phys%neq))
       ALLOCATE(utaux(SIZE(sol%u_tilde)/phys%neq,phys%neq))
       ALLOCATE(qaux(SIZE(sol%q)/phys%neq/Ndim,phys%neq*Ndim))
       uaux = TRANSPOSE(RESHAPE(sol%u,[phys%neq,SIZE(sol%u)/phys%neq]))
       utaux = TRANSPOSE(RESHAPE(sol%u_tilde,[phys%neq,SIZE(sol%u_tilde)/phys%neq]))
       qaux = TRANSPOSE(RESHAPE(sol%q,[phys%neq*Ndim,SIZE(sol%q)/phys%neq/Ndim]))
       uaux(:,1) = EXP(uaux(:,1))
       utaux(:,1) = EXP(utaux(:,1))
       qaux(:,1) = qaux(:,1)*uaux(:,1)
       qaux(:,2) = qaux(:,2)*uaux(:,1)
       sol%u = col(TRANSPOSE(uaux))
       sol%u_tilde = col(TRANSPOSE(utaux))
       sol%q = col(TRANSPOSE(qaux))
       DEALLOCATE(uaux,utaux,qaux)
    ENDIF


    ! Message to confirm succesful reading of file
    IF (MPIvar%glob_id .EQ. 0) THEN
       PRINT *, 'Solution read from file ', TRIM(ADJUSTL(fname_complete))
       PRINT *, '        '
    END IF

  END SUBROUTINE HDF5_load_solution

  !**********************************************************************
  ! Save HDG matrix (CSR) in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_CSR_matrix(fname)
    USE globals
    IMPLICIT NONE

    CHARACTER(LEN=*) :: fname
    CHARACTER(70)  :: npr, nid
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    CALL HDF5_create(fname_complete, file_id, ierr)
    CALL HDF5_integer_saving(file_id, MatK%n, 'n')
    CALL HDF5_integer_saving(file_id, MatK%nnz, 'nnz')
    CALL HDF5_array1D_saving_int(file_id, MatK%cols, MatK%nnz, 'cols')
    CALL HDF5_array1D_saving_int(file_id, MatK%rowptr, MatK%n + 1, 'rowptr')
    CALL HDF5_array1D_saving_int(file_id, MatK%loc2glob, MatK%n, 'loc2glob')
    CALL HDF5_array1D_saving(file_id, MatK%vals, MatK%nnz, 'vals')
    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    IF (MPIvar%glob_id.EQ.0) THEN
       PRINT*,'Output written to file ', TRIM(ADJUSTL(fname_complete))
       PRINT*,'        '
    END IF

  END SUBROUTINE HDF5_save_CSR_matrix

  !**********************************************************************
  ! Save HDG vector (CSR) in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_CSR_vector(fname)
    USE globals
    IMPLICIT NONE

    CHARACTER(LEN=*) :: fname
    CHARACTER(70)  :: npr, nid
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    CALL HDF5_create(fname_complete, file_id, ierr)
    CALL HDF5_integer_saving(file_id, rhs%n, 'n')
    CALL HDF5_array1D_saving_int(file_id, rhs%loc2glob, rhs%n, 'loc2glob')
    CALL HDF5_array1D_saving(file_id, rhs%vals, rhs%n, 'vals')
    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  END SUBROUTINE HDF5_save_CSR_vector

  !**********************************************************************
  ! Save 3D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_array(Arr, fname)
    USE globals
    IMPLICIT NONE

    REAL, DIMENSION(:, :, :), INTENT(IN) :: Arr
    CHARACTER(LEN=*) :: fname
    CHARACTER(70)  :: npr, nid
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    CALL HDF5_create(fname_complete, file_id, ierr)
    CALL HDF5_array3D_saving(file_id, Arr, SIZE(Arr, 1), SIZE(Arr, 2), SIZE(Arr, 3), 'array')
    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  END SUBROUTINE HDF5_save_array

  !**********************************************************************
  ! Save 2D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_matrix(Mat, fname)
    USE globals
    IMPLICIT NONE

    REAL, DIMENSION(:, :), INTENT(IN) :: Mat
    CHARACTER(LEN=*) :: fname
    CHARACTER(70)  :: npr, nid
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    CALL HDF5_create(fname_complete, file_id, ierr)
    CALL HDF5_array2D_saving(file_id, Mat, SIZE(Mat, 1), SIZE(Mat, 2), 'mat')
    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  END SUBROUTINE HDF5_save_matrix

  !**********************************************************************
  ! Save 1D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_vector(Vec, fname)
    USE globals
    IMPLICIT NONE

    REAL, DIMENSION(:), INTENT(IN) :: Vec
    CHARACTER(LEN=*) :: fname
    CHARACTER(70)  :: npr, nid
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    CALL HDF5_create(fname_complete, file_id, ierr)
    CALL HDF5_array1D_saving(file_id, Vec, SIZE(Vec, 1), 'vec')
    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  END SUBROUTINE HDF5_save_vector

  !**********************************************************************
  ! Save 1D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_vector_int(Vec, fname)
    USE globals
    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN) :: Vec
    CHARACTER(LEN=*) :: fname
    CHARACTER(70)  :: npr, nid
    INTEGER :: ierr
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    CALL HDF5_create(fname_complete, file_id, ierr)
    CALL HDF5_array1D_saving_int(file_id, Vec, SIZE(Vec, 1), 'vec')
    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    !      IF (MPIvar%glob_id.eq.0) THEN
    !                                                   print*,'Output written to file ', trim(adjustl(fname_complete))
    !                                                   print*,'        '
    !      END IF

  END SUBROUTINE HDF5_save_vector_int


  ! Define subroutine copy_file
  SUBROUTINE copy_file(source_file, destination_file)

    USE GMSH_io_module, ONLY: get_unit

    CHARACTER ( len = * ), INTENT(IN) :: source_file
    CHARACTER ( len = * ), INTENT(IN) :: destination_file
    CHARACTER ( len = 255 ) :: temp_source_file
    CHARACTER ( len = 255 ) :: temp_destination_file
    CHARACTER ( len = 255 ) :: buffer
    INTEGER :: unit_in, unit_out
    INTEGER :: ios

    temp_source_file = TRIM(ADJUSTL(source_file))
    temp_destination_file = TRIM(ADJUSTL(destination_file))

    CALL get_unit ( unit_in )
    ! Open the source file for reading
    OPEN(unit=unit_in, file=source_file, status='old', action='read', iostat=ios)
    IF (ios /= 0) THEN
       PRINT *, "Error opening source file:", source_file
       STOP
       RETURN
    END IF

    CALL get_unit ( unit_out )
    ! Open the destination file for writing
    OPEN(unit=unit_out, file=destination_file, status='replace', action='write', iostat=ios)
    IF (ios /= 0) THEN
       PRINT *, "Error opening destination file:", destination_file
       CLOSE(unit_in)
       RETURN
    END IF

    ! Copy data from source to destination
    DO
       READ (unit_in, '(a)', iostat = ios ) buffer
       IF (ios /= 0) EXIT ! Exit loop when end of file is reached
       WRITE(unit_out, '(a)') buffer(1:100)
    END DO

    ! Close the files
    CLOSE(unit_in)
    CLOSE(unit_out)
  END SUBROUTINE copy_file

END MODULE in_out
