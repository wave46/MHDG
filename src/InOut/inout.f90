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
  USE flux_surface_transport_data
  USE transport_models_1d
  USE GLOBALS
  USE MPI_OMP
  USE printutils

  IMPLICIT NONE

CONTAINS

  !********************************
  ! Loads mesh from an hdf5 file
  ! external file
  !********************************
  SUBROUTINE load_mesh_serial_h5(fname)

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
      IF(MPIvar%glob_id .EQ. 0) THEN
         PRINT *, 'Loading mesh.'
         PRINT *, '        '
       ENDIF
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

  ENDSUBROUTINE load_mesh_serial_h5

  SUBROUTINE load_mesh_h5(fname)

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
    INTEGER :: ghfa, ghel, Nel_glob, Nfa_glob, Ndir_glob, Ngho_glob, Nfaces, Nnodes_glob
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
    ALLOCATE (Mesh%loc2glob_nodes(Nnodes))
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
    CALL HDF5_array1D_reading_int(file_id, Mesh%loc2glob_nodes, 'loc2glob_no', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading loc2glob_no"
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
    Mesh%Nextfaces_nogho = COUNT((Mesh%boundaryFlag .NE. 0) .AND. (Mesh%ghostFaces(Mesh%Nintfaces+1:SIZE(Mesh%ghostFaces)) .EQ. 0))
    Mesh%Nintfaces_nogho = COUNT(Mesh%ghostFaces(1:Mesh%Nintfaces) .EQ. 0)
    CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_el), Nel_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_fa), Nfa_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_nodes), Nnodes_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%ndir, Ndir_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%nghostfaces, Ngho_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%Nextfaces_nogho, Mesh%Nextfaces_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(Mesh%Nintfaces_nogho, Mesh%Nintfaces_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    Mesh%Nel_glob = Nel_glob
    Mesh%Nfa_glob = Nfa_glob
    Mesh%Nno_glob = Nnodes_glob
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
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
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
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
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

  ENDSUBROUTINE load_mesh_h5

  !**********************************************************************
  ! Save solution in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_solution(fname)

#ifdef PARALL
    USE communications, ONLY: gather_mesh, gather_solution, gather_additional, gather_magnetic_field, gather_nodal_values
#endif
    IMPLICIT NONE

    CHARACTER(LEN=*)        :: fname
    CHARACTER(len=1000)     :: fname_complete
    INTEGER(HID_T)          :: file_id, group_id1
    INTEGER                 :: ierr

#ifdef TOR3D
    CHARACTER(70)           :: nip, nit, ngd
#endif
#ifdef PARALL
    INTEGER, POINTER        :: T_glob(:,:), Tb_glob(:,:), extfaces_glob(:,:), intfaces_glob(:,:), boundaryFlag_glob(:), periodic_faces_glob(:), F_glob(:,:), N_glob(:,:), face_info_glob(:,:), Tlin_glob(:,:), flag_elems_sc_glob(:)
    REAL*8, POINTER         :: u_tilde_glob(:), u_glob(:), q_glob(:), magnetic_psi_glob(:), magnetic_flux_glob(:), Jtor_glob(:), elemSize_glob(:), scdiff_nodes_glob(:,:)
    REAL*8, POINTER         :: X_glob(:,:), B_glob(:,:), Bperturb_glob(:,:)
    REAL*8, POINTER         :: external_heating_ions_glob(:), external_heating_electrons_glob(:)

    NULLIFY(T_glob, Tb_glob, extfaces_glob, intfaces_glob, boundaryFlag_glob, periodic_faces_glob, F_glob, N_glob, face_info_glob, Tlin_glob, flag_elems_sc_glob)
    NULLIFY(u_tilde_glob, u_glob, q_glob, magnetic_psi_glob, magnetic_flux_glob, Jtor_glob, elemSize_glob, scdiff_nodes_glob)
    NULLIFY(X_glob, B_glob, Bperturb_glob)
    NULLIFY(external_heating_ions_glob, external_heating_electrons_glob)

#endif

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
    fname_complete = TRIM(ADJUSTL(fname))//'.h5'
#endif

#ifndef PARALL
    ! Create file
    CALL HDF5_create(fname_complete, file_id, ierr)

    ! Save simulation parameters
    CALL save_simulation_parameters()

    ! save time iteration number
    IF (switch%steady .OR. switch%psdtime) THEN
       CALL HDF5_integer_saving(file_id,0,'it')
    ELSE
       CALL HDF5_integer_saving(file_id,time%it,'it')
    ENDIF

    ! save solution arrays
    CALL HDF5_group_create('solution', file_id, group_id1, ierr)
    CALL HDF5_array1D_saving(group_id1, sol%u, SIZE(sol%u), 'u')
    CALL HDF5_array1D_saving(group_id1, sol%u_tilde, SIZE(sol%u_tilde), 'u_tilde')
    CALL HDF5_array1D_saving(group_id1, sol%q, SIZE(sol%q), 'q')
    CALL HDF5_group_close(group_id1, ierr)

    CALL save_neutral_flux_limiter_diagnostics()
    CALL save_neutral_reaction_source_diagnostics()
    CALL save_neutral_wall_source_diagnostics()

    IF (switch%transport_1d) THEN
       CALL HDF5_group_create('transport_1d', file_id, group_id1, ierr)
       CALL fs_transport%write_hdf5(group_id1)
       CALL transport_model_1d%write_hdf5(group_id1)
       CALL HDF5_group_close(group_id1, ierr)
    END IF

    ! save magnetic field and Jtor arrays
    CALL HDF5_group_create('magnetic', file_id, group_id1, ierr)
    ! Save magnetic field
    CALL HDF5_array2D_saving(group_id1, phys%B, SIZE(phys%B, 1), SIZE(phys%B, 2), 'magnetic_field')
    ! Save normalized psi
    CALL HDF5_array1D_saving(group_id1, phys%magnetic_psi, SIZE(phys%magnetic_psi), 'magnetic_psi')
    ! Save magnetic flux
    CALL HDF5_array1D_saving(group_id1, phys%magnetic_flux, SIZE(phys%magnetic_flux), 'magnetic_flux')
    ! Save toroidal current
    IF (switch%ohmicsrc) THEN
       CALL HDF5_array1D_saving(group_id1, phys%Jtor, SIZE(phys%Jtor), 'Jtor')
    ENDIF
    ! Save magnetic perturbation and related fields
    IF ((switch%rmp).OR.(switch%ripple)) THEN
       CALL HDF5_array2D_saving(group_id1, phys%Bperturb, SIZE(phys%Bperturb, 1), SIZE(phys%Bperturb, 2), 'magnetic_perturbation')
    ENDIF
    IF (switch%rmp) THEN
       CALL HDF5_array3D_saving(group_id1, magn%coils_rmp, SIZE(magn%coils_rmp, 1), SIZE(magn%coils_rmp, 2), SIZE(magn%coils_rmp, 3), 'coils_rmp')
    ENDIF
    IF (switch%ripple) THEN
       CALL HDF5_array2D_saving(group_id1, magn%coils_ripple, SIZE(magn%coils_ripple, 1), SIZE(magn%coils_ripple, 2), 'coils_ripple')
    ENDIF
    CALL HDF5_group_close(group_id1, ierr)


    
    ! Save mesh related arrays
    CALL HDF5_group_create('mesh', file_id, group_id1, ierr)
    ! Save boundary structure
    CALL HDF5_array2D_saving_int(group_id1, Mesh%extfaces, SIZE(Mesh%extfaces, 1), SIZE(Mesh%extfaces, 2), 'extfaces')
    CALL HDF5_array1D_saving_int(group_id1, Mesh%boundaryFlag, SIZE(Mesh%boundaryFlag, 1), 'boundaryFlag')
    IF (switch%shockcp .EQ. 3) THEN
       CALL HDF5_array2D_saving(group_id1, Mesh%scdiff_nodes, SIZE(Mesh%scdiff_nodes, 1), SIZE(Mesh%scdiff_nodes, 2), 'scdiff_nodes')
    END IF
    CALL HDF5_integer_saving(group_id1,Mesh%Ndim,'Ndim')
    CALL HDF5_integer_saving(group_id1,Mesh%Nnodes,'Nnodes')
    CALL HDF5_integer_saving(group_id1,Mesh%Nnodesperelem,'Nnodesperelem')
    CALL HDF5_integer_saving(group_id1,Mesh%Nnodesperface,'Nnodesperface')
    CALL HDF5_integer_saving(group_id1,Mesh%Nelems,'Nelems')
    CALL HDF5_integer_saving(group_id1,Mesh%Nfaces,'Nfaces')
    CALL HDF5_integer_saving(group_id1,Mesh%Nextfaces,'Nextfaces')
    CALL HDF5_integer_saving(group_id1,Mesh%Nintfaces,'Nintfaces')
    CALL HDF5_integer_saving(group_id1,Mesh%elemType,'elemType')
    CALL HDF5_integer_saving(group_id1,Mesh%Ndir,'Ndir')
    CALL HDF5_integer_saving(group_id1,Mesh%ukf,'ukf')
    CALL HDF5_array2D_saving_int(group_id1,Mesh%T, SIZE(Mesh%T, 1), SIZE(Mesh%T, 2), 'T')
    CALL HDF5_array2D_saving_int(group_id1,Mesh%Tb, SIZE(Mesh%Tb, 1), SIZE(Mesh%Tb, 2), 'Tb')
    CALL HDF5_array2D_saving_int(group_id1,Mesh%Tlin, SIZE(Mesh%Tlin, 1), SIZE(Mesh%Tlin, 2), 'Tlin')
    CALL HDF5_array2D_saving_int(group_id1,Mesh%F, SIZE(Mesh%F,1),SIZE(Mesh%F,2), 'F')
    CALL HDF5_array2D_saving_int(group_id1,Mesh%N, SIZE(Mesh%N,1),SIZE(Mesh%N,2), 'N')
    IF(ALLOCATED(Mesh%faces)) THEN
       CALL HDF5_array3D_saving(group_id1,REAL(Mesh%faces), SIZE(Mesh%faces,1),SIZE(Mesh%faces,2),SIZE(Mesh%faces,3), 'faces')
    ENDIF
    CALL HDF5_array2D_saving_int(group_id1,Mesh%intfaces, SIZE(Mesh%intfaces,1),SIZE(Mesh%intfaces,2), 'intfaces')
    !call HDF5_array2D_saving_int(file_id,int(Mesh%flipface), SIZE(Mesh%flipface,1),SIZE(Mesh%flipface,2), 'flipface')
    !call HDF5_array2D_saving_logical(file_id,Mesh%Fdir, SIZE(Mesh%Fdir,1),SIZE(Mesh%Fdir,2), 'Fdir')
    IF(ALLOCATED(Mesh%periodic_faces)) THEN
       CALL HDF5_array1D_saving_int(group_id1,Mesh%periodic_faces, SIZE(Mesh%periodic_faces), 'periodic_faces')
    ENDIF
    IF(ALLOCATED(Mesh%Diric)) THEN
       CALL HDF5_array1D_saving_int(group_id1,Mesh%Diric, SIZE(Mesh%Diric), 'Diric')
    ENDIF
    IF(ALLOCATED(Mesh%numberbcs)) THEN
       CALL HDF5_array1D_saving_int(group_id1,Mesh%numberbcs, SIZE(Mesh%numberbcs), 'numberbcs')
    ENDIF
    CALL HDF5_array1D_saving(group_id1,Mesh%elemSize,SIZE(Mesh%elemSize), 'elemSize')
    CALL HDF5_array2D_saving(group_id1,Mesh%X*phys%lscale, SIZE(Mesh%X, 1), SIZE(Mesh%X, 2), 'X')
#ifdef TOR3D
    CALL HDF5_integer_saving(group_id1,Mesh%Nnodes_toroidal,'Nnodes_toroidal')
    CALL HDF5_array1D_saving(group_id1,Mesh%toroidal,SIZE(Mesh%toroidal), 'toroidal')
#endif
    IF(ALLOCATED(Mesh%flag_elems_rho)) THEN
       CALL HDF5_array1D_saving_int(group_id1,Mesh%flag_elems_rho,SIZE(Mesh%flag_elems_rho), 'flag_elems_rho')
    ENDIF
    IF(ALLOCATED(Mesh%flag_elems_sc)) THEN
       CALL HDF5_array1D_saving_int(group_id1,Mesh%flag_elems_sc,SIZE(Mesh%flag_elems_sc), 'flag_elems_sc')
    ENDIF
    IF(ALLOCATED(Mesh%minrho_elems)) THEN
       CALL HDF5_array1D_saving(group_id1,Mesh%minrho_elems,SIZE(Mesh%minrho_elems), 'minrho_elems')
    ENDIF
    IF(ALLOCATED(Mesh%sour_elems)) THEN
       CALL HDF5_array1D_saving(group_id1,Mesh%sour_elems,SIZE(Mesh%sour_elems), 'sour_elems')
    ENDIF
    IF(ALLOCATED(Mesh%diff_elems)) THEN
       CALL HDF5_array1D_saving(group_id1,Mesh%diff_elems,SIZE(Mesh%diff_elems), 'diff_elems')
    ENDIF
    IF(ALLOCATED(Mesh%scdiff_nodes)) THEN
       CALL HDF5_array2D_saving(group_id1,Mesh%scdiff_nodes,SIZE(Mesh%scdiff_nodes,1),SIZE(Mesh%scdiff_nodes,2), 'scdiff_nodes')
    ENDIF
    CALL HDF5_real_saving(group_id1, Mesh%xmax, 'xmax')
    CALL HDF5_real_saving(group_id1, Mesh%xmin, 'xmin')
    CALL HDF5_real_saving(group_id1, Mesh%ymax, 'ymax')
    CALL HDF5_real_saving(group_id1, Mesh%ymin, 'ymin')
    CALL HDF5_real_saving(group_id1, Mesh%puff_area, 'puff_area')
    CALL HDF5_real_saving(group_id1, Mesh%pump_area, 'pump_area')
    CALL HDF5_real_saving(group_id1, Mesh%core_area, 'core_area')
    CALL HDF5_group_close(group_id1, ierr)
    

#else

    ! gather solution, mesh and magnetic related arrays
    CALL gather_solution(Mesh_in = Mesh, Nnodesperelem = Mesh%Nnodesperelem, Nnodesperface = Mesh%Nnodesperface, u_tilde_in = sol%u_tilde, u_in = sol%u, q_in = sol%q,  u_tilde_glob = u_tilde_glob, u_glob = u_glob, q_glob = q_glob)

    IF(switch%OhmicSrc) THEN
       IF ((switch%RMP) .OR. (switch%Ripple)) THEN
          CALL gather_magnetic_field(Mesh, B_glob, magnetic_flux_glob, magnetic_psi_glob, Bperturb_glob, Jtor_glob)
       ELSE
          CALL gather_magnetic_field(Mesh_in = Mesh, B_glob = B_glob, magnetic_flux_glob = magnetic_flux_glob, magnetic_psi_glob = magnetic_psi_glob, Jtor_glob = Jtor_glob)
       ENDIF
    ELSE
       IF ((switch%RMP) .OR. (switch%Ripple)) THEN
          CALL gather_magnetic_field(Mesh_in = Mesh, B_glob = B_glob, magnetic_flux_glob = magnetic_flux_glob, magnetic_psi_glob = magnetic_psi_glob, Bperturb_glob = Bperturb_glob)
       ELSE
          CALL gather_magnetic_field(Mesh_in = Mesh, B_glob = B_glob, magnetic_flux_glob = magnetic_flux_glob, magnetic_psi_glob = magnetic_psi_glob)
       ENDIF
    ENDIF

    IF (switch%external_heating) THEN
      CALL gather_nodal_values(Mesh_in = Mesh, value_in=phys%external_heating_ions, value_glob=external_heating_ions_glob,allgather=.TRUE.)
      CALL gather_nodal_values(Mesh_in = Mesh, value_in=phys%external_heating_electrons, value_glob=external_heating_electrons_glob,allgather=.TRUE.)
    ENDIF

    
    IF ((switch%shockcp .NE. 0) .OR. (adapt%shockcp_adapt .NE. 0)) THEN
      CALL gather_mesh(Mesh, T_glob, X_glob, Tb_glob, F_glob, N_glob, intfaces_glob, extfaces_glob, boundaryFlag_glob, Tlin_glob, periodic_faces_glob, elemSize_glob, flag_elems_sc_glob, scdiff_nodes_glob)
    ELSE
      CALL gather_mesh(Mesh, T_glob, X_glob, Tb_glob, F_glob, N_glob, intfaces_glob, extfaces_glob, boundaryFlag_glob, Tlin_glob, periodic_faces_glob, elemSize_glob)
    ENDIF
    ! save to file
    IF (MPIvar%glob_id .EQ. 0) THEN

       CALL HDF5_create(fname_complete, file_id, ierr)

       IF (switch%steady .OR. switch%psdtime) THEN
          CALL HDF5_integer_saving(file_id,0,'it')
       ELSE
          CALL HDF5_integer_saving(file_id,time%it,'it')
       ENDIF

       ! Save simulation parameters
       CALL save_simulation_parameters()

       CALL HDF5_group_create('solution', file_id, group_id1, ierr)
       CALL HDF5_array1D_saving(group_id1, u_tilde_glob, SIZE(u_tilde_glob), 'u_tilde')
       CALL HDF5_array1D_saving(group_id1, u_glob, SIZE(u_glob), 'u')
       CALL HDF5_array1D_saving(group_id1, q_glob, SIZE(q_glob), 'q')
       CALL HDF5_group_close(group_id1)

       IF (switch%transport_1d) THEN
          CALL HDF5_group_create('transport_1d', file_id, group_id1, ierr)
          CALL fs_transport%write_hdf5(group_id1)
          CALL transport_model_1d%write_hdf5(group_id1)
          CALL HDF5_group_close(group_id1)
       END IF

      CALL HDF5_group_create('mesh', file_id, group_id1, ierr)
      CALL HDF5_integer_saving(group_id1,Mesh%Ndim,'Ndim')
      CALL HDF5_integer_saving(group_id1,Mesh%Nno_glob,'Nnodes')
      CALL HDF5_integer_saving(group_id1,Mesh%Nel_glob,'Nelems')
      CALL HDF5_integer_saving(group_id1,Mesh%Nfa_glob,'Nfaces')
      CALL HDF5_integer_saving(group_id1,Mesh%Nextfaces_glob,'Nextfaces')
      CALL HDF5_integer_saving(group_id1,Mesh%Nintfaces_glob,'Nintfaces')
      CALL HDF5_integer_saving(group_id1,Mesh%Ndir_glob,'Ndir')
      CALL HDF5_integer_saving(group_id1,Mesh%Nnodesperelem,'Nnodesperelem')
      CALL HDF5_integer_saving(group_id1,Mesh%Nnodesperface,'Nnodesperface')
      CALL HDF5_integer_saving(group_id1,Mesh%elemType,'elemType')
      ! these are already reduced in preprocess or load mesh
      CALL HDF5_real_saving(group_id1, Mesh%puff_area, 'puff_area')
      CALL HDF5_real_saving(group_id1, Mesh%pump_area, 'pump_area')
      CALL HDF5_real_saving(group_id1, Mesh%core_area, 'core_area')
      CALL HDF5_real_saving(group_id1, Mesh%xmax, 'xmax')
      CALL HDF5_real_saving(group_id1, Mesh%xmin, 'xmin')
      CALL HDF5_real_saving(group_id1, Mesh%ymax, 'ymax')
      CALL HDF5_real_saving(group_id1, Mesh%ymin, 'ymin')
      CALL HDF5_integer_saving(group_id1,Mesh%Nfa_glob,'ukf')

      CALL HDF5_array2D_saving_int(group_id1,T_glob, SIZE(T_glob, 1), SIZE(T_glob, 2), 'T')
      CALL HDF5_array2D_saving_int(group_id1,Tb_glob, SIZE(Tb_glob, 1), SIZE(Tb_glob, 2), 'Tb')
      CALL HDF5_array2D_saving(group_id1,X_glob*phys%lscale, SIZE(X_glob, 1), SIZE(X_glob, 2), 'X')

      !CALL HDF5_array1D_saving_int(file_id,Mesh%boundaryFlag, SIZE(Mesh%boundaryFlag), 'boundaryFlag')
      CALL HDF5_array2D_saving_int(group_id1,F_glob, SIZE(F_glob,1),SIZE(F_glob,2), 'F')
      CALL HDF5_array2D_saving_int(group_id1,N_glob, SIZE(N_glob,1),SIZE(N_glob,2), 'N')

      CALL HDF5_array1D_saving(group_id1,elemSize_glob,SIZE(elemSize_glob), 'elemSize')

      ! Save boundary structure
      CALL HDF5_array2D_saving_int(group_id1, extfaces_glob, SIZE(extfaces_glob, 1), SIZE(extfaces_glob, 2), 'extfaces')
      CALL HDF5_array2D_saving_int(group_id1, intfaces_glob, SIZE(intfaces_glob,1),SIZE(intfaces_glob,2), 'intfaces')
      CALL HDF5_array1D_saving_int(group_id1, boundaryFlag_glob, SIZE(boundaryFlag_glob, 1), 'boundaryFlag')

      IF(ALLOCATED(Mesh%periodic_faces)) THEN
         CALL HDF5_array1D_saving_int(group_id1,periodic_faces_glob, SIZE(periodic_faces_glob), 'periodic_faces')
      ENDIF
      IF(ASSOCIATED(Mesh%Tlin)) THEN
         CALL HDF5_array2D_saving_int(group_id1,Tlin_glob, SIZE(Tlin_glob, 1), SIZE(Tlin_glob, 2), 'Tlin')
      ENDIF
      IF(ALLOCATED(Mesh%flag_elems_sc)) THEN
        CALL HDF5_array1D_saving_int(group_id1,flag_elems_sc_glob,SIZE(flag_elems_sc_glob), 'flag_elems_sc')
      ENDIF
      IF(ALLOCATED(Mesh%scdiff_nodes)) THEN
        CALL HDF5_array2D_saving(group_id1,scdiff_nodes_glob,SIZE(scdiff_nodes_glob,1),SIZE(scdiff_nodes_glob,2), 'scdiff_nodes')
      ENDIF

#ifdef TOR3D
      CALL HDF5_integer_saving(group_id1,Mesh%Nnodes_toroidal,'Nnodes_toroidal')
      CALL HDF5_array1D_saving(group_id1,Mesh%toroidal,SIZE(Mesh%toroidal), 'toroidal')
#endif
      CALL HDF5_group_close(group_id1, ierr)
       

       CALL HDF5_group_create('magnetic', file_id, group_id1, ierr)
       IF (switch%rmp) THEN
          CALL HDF5_array3D_saving(group_id1, magn%coils_rmp, SIZE(magn%coils_rmp, 1), SIZE(magn%coils_rmp, 2), SIZE(magn%coils_rmp, 3), 'coils_rmp')
       ENDIF
       IF (switch%ripple) THEN
          CALL HDF5_array2D_saving(group_id1, magn%coils_ripple, SIZE(magn%coils_ripple, 1), SIZE(magn%coils_ripple, 2), 'coils_ripple')
       ENDIF
       ! Save magnetic field
       CALL HDF5_array2D_saving(group_id1, B_glob, SIZE(B_glob, 1), SIZE(B_glob, 2), 'magnetic_field')
       ! Save magnetic flux
       CALL HDF5_array1D_saving(group_id1, magnetic_flux_glob, SIZE(magnetic_flux_glob), 'magnetic_flux')
       ! Save normalized psi
       CALL HDF5_array1D_saving(group_id1, magnetic_psi_glob, SIZE(magnetic_psi_glob), 'magnetic_psi')
       ! Save toroidal current
       IF (switch%ohmicsrc) THEN
          CALL HDF5_array1D_saving(group_id1, Jtor_glob, SIZE(Jtor_glob), 'Jtor')
       ENDIF
       ! Save magnetic perturbation and related fields
       IF ((switch%rmp).OR.(switch%ripple)) THEN
          CALL HDF5_array2D_saving(group_id1, Bperturb_glob, SIZE(Bperturb_glob, 1), SIZE(Bperturb_glob, 2), 'magnetic_perturbation')
       ENDIF

       IF (switch%shockcp .EQ. 3) THEN
          CALL HDF5_array2D_saving(group_id1, scdiff_nodes_glob, SIZE(scdiff_nodes_glob, 1), SIZE(scdiff_nodes_glob, 2), 'scdiff_nodes')
       END IF
       CALL HDF5_group_close(group_id1, ierr)

    END IF

    CALL save_neutral_flux_limiter_diagnostics()
    CALL save_neutral_reaction_source_diagnostics()
    CALL save_neutral_wall_source_diagnostics()

    IF (ASSOCIATED(T_glob)) THEN
      DEALLOCATE(T_glob, Tb_glob, extfaces_glob, intfaces_glob, boundaryFlag_glob, periodic_faces_glob, F_glob, N_glob, Tlin_glob)
      DEALLOCATE(u_tilde_glob, u_glob, q_glob, magnetic_psi_glob, magnetic_flux_glob, elemSize_glob, X_glob, B_glob)      
      NULLIFY(T_glob, Tb_glob, extfaces_glob, intfaces_glob, boundaryFlag_glob, periodic_faces_glob, F_glob, N_glob, Tlin_glob)
      NULLIFY(u_tilde_glob, u_glob, q_glob, magnetic_psi_glob, magnetic_flux_glob, elemSize_glob, X_glob, B_glob)
    ENDIF
    IF (ASSOCIATED(external_heating_ions_glob)) THEN
      DEALLOCATE(external_heating_ions_glob, external_heating_electrons_glob)
      NULLIFY(external_heating_ions_glob, external_heating_electrons_glob)
    ENDIF

    IF(ASSOCIATED(Jtor_glob)) THEN
      DEALLOCATE(Jtor_glob)
      NULLIFY(Jtor_glob)
    ENDIF
    IF(ASSOCIATED(Bperturb_glob)) THEN
      DEALLOCATE(Bperturb_glob)
      NULLIFY(Bperturb_glob)
    ENDIF
    IF(ASSOCIATED(flag_elems_sc_glob)) THEN
      DEALLOCATE(flag_elems_sc_glob)
      NULLIFY(flag_elems_sc_glob)
    ENDIF
    IF(ASSOCIATED(scdiff_nodes_glob)) THEN
      DEALLOCATE(scdiff_nodes_glob)
      NULLIFY(scdiff_nodes_glob)
    ENDIF
#endif

  IF(MPIvar%glob_id .eq. 0) THEN
    CALL HDF5_close(file_id)
    ! Message to confirm succesful creation and filling of file
    PRINT *, 'Output written to file ', TRIM(ADJUSTL(fname_complete))
    PRINT *, '        '
  ENDIF
  CONTAINS

    SUBROUTINE save_neutral_flux_limiter_diagnostics()
      INTEGER(HID_T) :: group_id
      INTEGER :: expected_size
#ifdef PARALL
      INTEGER :: iel, g, ind_local, ind_global
      REAL*8, ALLOCATABLE :: Dnn_glob(:), phi_glob(:), Deff_glob(:), Gamma_unlim_glob(:)
      REAL*8, ALLOCATABLE :: Gamma_max_glob(:), activation_ratio_glob(:), Gamma_lim_glob(:)
#endif

      IF (.NOT. ALLOCATED(phys%neutral_flux_limiter_phi_Nod)) RETURN
      expected_size = Mesh%Nelems*Mesh%Nnodesperelem
      IF (SIZE(phys%neutral_flux_limiter_phi_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_flux_limiter_Dnn_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_flux_limiter_Deff_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_flux_limiter_Gamma_unlim_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_flux_limiter_Gamma_max_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_flux_limiter_activation_ratio_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_flux_limiter_Gamma_lim_Nod) .NE. expected_size) RETURN

#ifdef PARALL
      ALLOCATE(Dnn_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(phi_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(Deff_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(Gamma_unlim_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(Gamma_max_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(activation_ratio_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(Gamma_lim_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      Dnn_glob = 0.d0
      phi_glob = 0.d0
      Deff_glob = 0.d0
      Gamma_unlim_glob = 0.d0
      Gamma_max_glob = 0.d0
      activation_ratio_glob = 0.d0
      Gamma_lim_glob = 0.d0

      DO iel = 1, Mesh%Nelems
        IF (Mesh%ghostElems(iel) .EQ. 0) THEN
          DO g = 1, Mesh%Nnodesperelem
            ind_local = (iel - 1)*Mesh%Nnodesperelem + g
            ind_global = (Mesh%loc2glob_el(iel) - 1)*Mesh%Nnodesperelem + g
            Dnn_glob(ind_global) = phys%neutral_flux_limiter_Dnn_Nod(ind_local)
            phi_glob(ind_global) = phys%neutral_flux_limiter_phi_Nod(ind_local)
            Deff_glob(ind_global) = phys%neutral_flux_limiter_Deff_Nod(ind_local)
            Gamma_unlim_glob(ind_global) = phys%neutral_flux_limiter_Gamma_unlim_Nod(ind_local)
            Gamma_max_glob(ind_global) = phys%neutral_flux_limiter_Gamma_max_Nod(ind_local)
            activation_ratio_glob(ind_global) = phys%neutral_flux_limiter_activation_ratio_Nod(ind_local)
            Gamma_lim_glob(ind_global) = phys%neutral_flux_limiter_Gamma_lim_Nod(ind_local)
          ENDDO
        ENDIF
      ENDDO

      CALL MPI_Allreduce(MPI_IN_PLACE, Dnn_glob, SIZE(Dnn_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, phi_glob, SIZE(phi_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, Deff_glob, SIZE(Deff_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, Gamma_unlim_glob, SIZE(Gamma_unlim_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, Gamma_max_glob, SIZE(Gamma_max_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, activation_ratio_glob, SIZE(activation_ratio_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, Gamma_lim_glob, SIZE(Gamma_lim_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      IF (MPIvar%glob_id .EQ. 0) THEN
        CALL HDF5_group_create('neutral_flux_limiter_diagnostics', file_id, group_id, ierr)
        CALL HDF5_array1D_saving(group_id, Dnn_glob, SIZE(Dnn_glob), 'Dnn')
        CALL HDF5_array1D_saving(group_id, phi_glob, SIZE(phi_glob), 'phi')
        CALL HDF5_array1D_saving(group_id, Deff_glob, SIZE(Deff_glob), 'D_eff')
        CALL HDF5_array1D_saving(group_id, Gamma_unlim_glob, SIZE(Gamma_unlim_glob), 'Gamma_unlim')
        CALL HDF5_array1D_saving(group_id, Gamma_max_glob, SIZE(Gamma_max_glob), 'Gamma_max')
        CALL HDF5_array1D_saving(group_id, activation_ratio_glob, SIZE(activation_ratio_glob), 'activation_ratio')
        CALL HDF5_array1D_saving(group_id, Gamma_lim_glob, SIZE(Gamma_lim_glob), 'Gamma_lim')
        CALL HDF5_group_close(group_id, ierr)
      ENDIF

      DEALLOCATE(Dnn_glob, phi_glob, Deff_glob, Gamma_unlim_glob)
      DEALLOCATE(Gamma_max_glob, activation_ratio_glob, Gamma_lim_glob)
#else
      CALL HDF5_group_create('neutral_flux_limiter_diagnostics', file_id, group_id, ierr)
      CALL HDF5_array1D_saving(group_id, phys%neutral_flux_limiter_Dnn_Nod, &
        &SIZE(phys%neutral_flux_limiter_Dnn_Nod), 'Dnn')
      CALL HDF5_array1D_saving(group_id, phys%neutral_flux_limiter_phi_Nod, &
        &SIZE(phys%neutral_flux_limiter_phi_Nod), 'phi')
      CALL HDF5_array1D_saving(group_id, phys%neutral_flux_limiter_Deff_Nod, &
        &SIZE(phys%neutral_flux_limiter_Deff_Nod), 'D_eff')
      CALL HDF5_array1D_saving(group_id, phys%neutral_flux_limiter_Gamma_unlim_Nod, &
        &SIZE(phys%neutral_flux_limiter_Gamma_unlim_Nod), 'Gamma_unlim')
      CALL HDF5_array1D_saving(group_id, phys%neutral_flux_limiter_Gamma_max_Nod, &
        &SIZE(phys%neutral_flux_limiter_Gamma_max_Nod), 'Gamma_max')
      CALL HDF5_array1D_saving(group_id, phys%neutral_flux_limiter_activation_ratio_Nod, &
        &SIZE(phys%neutral_flux_limiter_activation_ratio_Nod), 'activation_ratio')
      CALL HDF5_array1D_saving(group_id, phys%neutral_flux_limiter_Gamma_lim_Nod, &
        &SIZE(phys%neutral_flux_limiter_Gamma_lim_Nod), 'Gamma_lim')
      CALL HDF5_group_close(group_id, ierr)
#endif
    ENDSUBROUTINE save_neutral_flux_limiter_diagnostics

    SUBROUTINE save_neutral_reaction_source_diagnostics()
      INTEGER(HID_T) :: group_id

      IF (MPIvar%glob_id .EQ. 0) THEN
        CALL HDF5_group_create('neutral_reaction_sources_diagnostics', file_id, group_id, ierr)
        CALL HDF5_real_saving(group_id, phys%neutral_ionization_total, 'ionization_sink')
        CALL HDF5_real_saving(group_id, phys%neutral_recombination_total, 'recombination_source')
        CALL HDF5_real_saving(group_id, phys%neutral_charge_exchange_total, 'charge_exchange_rate')
        CALL HDF5_group_close(group_id, ierr)
      ENDIF
    ENDSUBROUTINE save_neutral_reaction_source_diagnostics

    SUBROUTINE save_neutral_wall_source_diagnostics()
      INTEGER(HID_T) :: group_id
      INTEGER :: expected_size
#ifdef PARALL
      INTEGER :: iel, g, ind_local, ind_global
      REAL*8, ALLOCATABLE :: puff_glob(:), pump_glob(:), net_glob(:)
#endif

      IF (.NOT. ALLOCATED(phys%neutral_wall_source_puff_Nod)) RETURN
      expected_size = Mesh%Nelems*Mesh%Nnodesperelem
      IF (SIZE(phys%neutral_wall_source_puff_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_wall_source_pump_Nod) .NE. expected_size) RETURN
      IF (SIZE(phys%neutral_wall_source_net_Nod) .NE. expected_size) RETURN

#ifdef PARALL
      ALLOCATE(puff_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(pump_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      ALLOCATE(net_glob(Mesh%Nel_glob*Mesh%Nnodesperelem))
      puff_glob = 0.d0
      pump_glob = 0.d0
      net_glob = 0.d0

      DO iel = 1, Mesh%Nelems
        IF (Mesh%ghostElems(iel) .EQ. 0) THEN
          DO g = 1, Mesh%Nnodesperelem
            ind_local = (iel - 1)*Mesh%Nnodesperelem + g
            ind_global = (Mesh%loc2glob_el(iel) - 1)*Mesh%Nnodesperelem + g
            puff_glob(ind_global) = phys%neutral_wall_source_puff_Nod(ind_local)
            pump_glob(ind_global) = phys%neutral_wall_source_pump_Nod(ind_local)
            net_glob(ind_global) = phys%neutral_wall_source_net_Nod(ind_local)
          ENDDO
        ENDIF
      ENDDO

      CALL MPI_Allreduce(MPI_IN_PLACE, puff_glob, SIZE(puff_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, pump_glob, SIZE(pump_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce(MPI_IN_PLACE, net_glob, SIZE(net_glob), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      IF (MPIvar%glob_id .EQ. 0) THEN
        CALL HDF5_group_create('neutral_wall_sources_diagnostics', file_id, group_id, ierr)
        CALL HDF5_real_saving(group_id, phys%neutral_wall_source_puff_total, 'element_puff_total')
        CALL HDF5_real_saving(group_id, phys%neutral_wall_source_pump_total, 'element_pump_total')
        CALL HDF5_real_saving(group_id, phys%neutral_wall_source_puff_total - phys%neutral_wall_source_pump_total, &
          &'element_net_total')
        CALL HDF5_array1D_saving(group_id, puff_glob, SIZE(puff_glob), 'puff_flux_density')
        CALL HDF5_array1D_saving(group_id, pump_glob, SIZE(pump_glob), 'pump_flux_density')
        CALL HDF5_array1D_saving(group_id, net_glob, SIZE(net_glob), 'net_flux_density')
        CALL HDF5_group_close(group_id, ierr)
      ENDIF

      DEALLOCATE(puff_glob, pump_glob, net_glob)
#else
      CALL HDF5_group_create('neutral_wall_sources_diagnostics', file_id, group_id, ierr)
      CALL HDF5_real_saving(group_id, phys%neutral_wall_source_puff_total, 'element_puff_total')
      CALL HDF5_real_saving(group_id, phys%neutral_wall_source_pump_total, 'element_pump_total')
      CALL HDF5_real_saving(group_id, phys%neutral_wall_source_puff_total - phys%neutral_wall_source_pump_total, &
        &'element_net_total')
      CALL HDF5_array1D_saving(group_id, phys%neutral_wall_source_puff_Nod, &
        &SIZE(phys%neutral_wall_source_puff_Nod), 'puff_flux_density')
      CALL HDF5_array1D_saving(group_id, phys%neutral_wall_source_pump_Nod, &
        &SIZE(phys%neutral_wall_source_pump_Nod), 'pump_flux_density')
      CALL HDF5_array1D_saving(group_id, phys%neutral_wall_source_net_Nod, &
        &SIZE(phys%neutral_wall_source_net_Nod), 'net_flux_density')
      CALL HDF5_group_close(group_id, ierr)
#endif
    ENDSUBROUTINE save_neutral_wall_source_diagnostics

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
      CALL HDF5_real_saving(group_id2, phys%diff_nn_min, 'diff_nn_min')
      CALL HDF5_real_saving(group_id2, phys%Re, 'recycling')
      CALL HDF5_real_saving(group_id2, phys%recycling_neutral_gamma, 'recycling_neutral_gamma')
      CALL HDF5_real_saving(group_id2, phys%impurity_concentration, 'impurity_concentration')
      CALL HDF5_string_saving(group_id2, phys%impurity_name, 'impurity_name')
      CALL HDF5_logical_saving(group_id2, phys%apply_trim, 'apply_trim')
      CALL HDF5_real_saving(group_id2, phys%Zeff, 'Zeff')
      CALL HDF5_real_saving(group_id2, phys%Pohmic, 'ohmic_coeff')
      CALL HDF5_real_saving(group_id2, phys%Re_pump, 'recycling_pump')
      CALL HDF5_real_saving(group_id2, phys%puff, 'puff')
      CALL HDF5_real_saving(group_id2, phys%cryopump_power, 'cryopump_power')
      CALL HDF5_real_saving(group_id2, phys%r_axis, 'r_axis')
      CALL HDF5_real_saving(group_id2, phys%z_axis, 'z_axis')
      CALL HDF5_real_saving(group_id2, phys%a_minor, 'a_minor')
      IF (switch%ME) THEN
         CALL HDF5_array1d_saving(group_id2, phys%puff_exp, input%puff_dimension, 'puff_exp')
         IF (switch%target_variable /= 0) THEN
            CALL HDF5_real_saving(group_id2, phys%feedback_propotional_gain,'feedback_propotional_gain')
            CALL HDF5_real_saving(group_id2, phys%feedback_integral_gain,'feedback_integral_gain')
            CALL HDF5_real_saving(group_id2, phys%feedback_derivative_gain,'feedback_derivative_gain')
            CALL HDF5_real_saving(group_id2, phys%n_li,'n_li')
            CALL HDF5_real_saving(group_id2, phys%feedback_integral_error,'feedback_integral_error')
            CALL HDF5_real_saving(group_id2, phys%feedback_previous_error,'feedback_previous_error')
            IF (switch%target_variable == 3) THEN
               CALL HDF5_real_saving(group_id2, phys%feedback_propotional_gain_xpr,'feedback_propotional_gain_xpr')
               CALL HDF5_real_saving(group_id2, phys%feedback_integral_gain_xpr,'feedback_integral_gain_xpr')
               CALL HDF5_real_saving(group_id2, phys%feedback_derivative_gain_xpr,'feedback_derivative_gain_xpr')
               CALL HDF5_real_saving(group_id2, phys%feedback_integral_error_xpr,'feedback_integral_error_xpr')
               CALL HDF5_real_saving(group_id2, phys%feedback_previous_error_xpr,'feedback_previous_error_xpr')
            ENDIF
         ENDIF
         IF (switch%diff_reverse_Ip) THEN
            CALL HDF5_real_saving(group_id2, phys%I_0, 'I_0')
            CALL HDF5_real_saving(group_id2, phys%I_p, 'I_p')
            CALL HDF5_real_saving(group_id2, phys%ME_diff_n, 'ME_diff_n')
            CALL HDF5_real_saving(group_id2, phys%ME_diff_u, 'ME_diff_u')
            CALL HDF5_real_saving(group_id2, phys%ME_diff_e, 'ME_diff_e')
            CALL HDF5_real_saving(group_id2, phys%ME_diff_ee, 'ME_diff_ee')
         ENDIF
      END IF
      if (switch%external_heating) THEN
#ifdef PARALL
         CALL HDF5_array1D_saving(group_id2, external_heating_ions_glob, SIZE(external_heating_ions_glob), 'external_heating_ions')
         CALL HDF5_array1D_saving(group_id2, external_heating_electrons_glob, SIZE(external_heating_electrons_glob), 'external_heating_electrons')
#else
         CALL HDF5_array1D_saving(group_id2, phys%external_heating_ions, SIZE(phys%external_heating_ions), 'external_heating_ions')
         CALL HDF5_array1D_saving(group_id2, phys%external_heating_electrons, SIZE(phys%external_heating_electrons), 'external_heating_electrons')
#endif
      ENDIF
      CALL HDF5_real_saving(group_id2, phys%tie, 'tau_ie')
      CALL HDF5_real_saving(group_id2, phys%dfcoef, 'dfcoef')
      CALL HDF5_real_saving(group_id2, phys%dexbcoef, 'dexbcoef')
      CALL HDF5_real_saving(group_id2, phys%bohmth, 'bohmth')
      CALL HDF5_real_saving(group_id2, phys%bohm_energy_thresh, 'bohm_energy_thresh')
      CALL HDF5_real_saving(group_id2, phys%epn, 'epn')
      CALL HDF5_real_saving(group_id2, phys%Gmbohm, 'Gmbohm')
      CALL HDF5_real_saving(group_id2, phys%Gmbohme, 'Gmbohme')
      CALL HDF5_real_saving(group_id2, phys%Potfloat, 'Potfloat')
      CALL HDF5_array1d_saving_int(group_id2, phys%bcflags, 10, 'boundary_flags')
      IF (switch%flux_limiter) THEN
         CALL HDF5_real_saving(group_id2, phys%c_fli, 'c_fli')
         CALL HDF5_real_saving(group_id2, phys%c_fle, 'c_fle')
      ENDIF
      CALL HDF5_real_saving(group_id2, phys%T_fluxlim_maxi, 'T_fluxlim_maxi')
      CALL HDF5_real_saving(group_id2, phys%T_fluxlim_maxe, 'T_fluxlim_maxe')
      CALL HDF5_string_saving(group_id2, phys%neutral_flux_limiter_mode, 'neutral_flux_limiter_mode')
      CALL HDF5_real_saving(group_id2, phys%neutral_flux_limiter_gamma, 'neutral_flux_limiter_gamma')
      CALL HDF5_real_saving(group_id2, phys%neutral_flux_limiter_eps, 'neutral_flux_limiter_eps')
      CALL HDF5_real_saving(group_id2, phys%neutral_flux_limiter_fs_fraction, 'neutral_flux_limiter_fs_fraction')
      CALL HDF5_real_saving(group_id2, phys%neutral_flux_limiter_fs_flux_min, 'neutral_flux_limiter_fs_flux_min')
      IF (switch%import_diffusion_1D) THEN 
         CALL HDF5_array1D_saving(group_id2, phys%rho_1D, SIZE(phys%rho_1D), 'rho_1D')
         CALL HDF5_array1D_saving(group_id2, phys%diff_n_1D, SIZE(phys%diff_n_1D), 'diff_n_1D')
         CALL HDF5_array1D_saving(group_id2, phys%diff_u_1D, SIZE(phys%diff_u_1D), 'diff_u_1D')
         CALL HDF5_array1D_saving(group_id2, phys%diff_e_1D, SIZE(phys%diff_e_1D), 'diff_e_1D')
         CALL HDF5_array1D_saving(group_id2, phys%diff_ee_1D, SIZE(phys%diff_ee_1D), 'diff_ee_1D')
      ENDIF
      CALL HDF5_group_close(group_id2, ierr)
         

      ! Create switches parameters group
      CALL HDF5_group_create('switches', group_id1, group_id2, ierr)
      CALL HDF5_logical_saving(group_id2, switch%driftdia, 'diamagnetic_drift')
      CALL HDF5_logical_saving(group_id2, switch%driftexb, 'ExB_drift')
      CALL HDF5_logical_saving(group_id2, switch%steady, 'steady')
      CALL HDF5_integer_saving(group_id2, switch%testcase, 'testcase')
      CALL HDF5_logical_saving(group_id2, switch%ohmicsrc, 'ohmicsrc')
      CALL HDF5_logical_saving(group_id2, switch%ME, 'ME')
      CALL HDF5_integer_saving(group_id2, switch%target_variable, 'target_variable')
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
      CALL HDF5_logical_saving(group_id2, switch%flux_limiter, 'flux_limiter')
      CALL HDF5_logical_saving(group_id2, switch%impurity_radiation, 'impurity_radiation')
      CALL HDF5_logical_saving(group_id2, switch%external_heating, 'external_heating')
      CALL HDF5_logical_saving(group_id2, switch%import_diffusion_1D, 'import_diffusion_1D')
      CALL HDF5_logical_saving(group_id2, switch%neutral_wall_sources_in_elements, 'neutral_wall_sources_in_elements')
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
      CALL HDF5_real_saving(group_id2, numer%neutralp_lambda, 'NeutralP_lambda')
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


    ENDSUBROUTINE save_simulation_parameters

  ENDSUBROUTINE HDF5_save_solution

  SUBROUTINE HDF5_load_mesh_from_solution(fname)
    !*************************************
    !              2D case
    !*************************************


    CHARACTER(LEN=*) :: fname
    CHARACTER(len=1000) :: fname_complete
    CHARACTER(10)  :: str
    REAL*8, PARAMETER::tol = 1e-6
    REAL*8 :: xmin
    INTEGER :: elemType, ndim, Nnodes, Nelems, Nnodesperelem
    INTEGER :: Nextfaces, Nnodesperface, Nfaces, IERR
    INTEGER(HID_T) :: file_id, group_id
#ifdef TOR3D
    INTEGER         :: i
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
    fname_complete = TRIM(ADJUSTL(fname))//'.h5'
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

    CALL HDF5_group_open(file_id, 'mesh', group_id, ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error opening group 'mesh'"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, elemType, 'elemType', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: elemType"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, ndim, 'Ndim', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Ndim"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, Nnodes, 'Nnodes', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodes"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, Nelems, 'Nelems', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nelems"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, Nnodesperelem, 'Nnodesperelem', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperelem"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, Nnodesperface, 'Nnodesperface', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nnodesperface"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, Nextfaces, 'Nextfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nextfaces"
       STOP
    ENDIF
    CALL HDF5_integer_reading(group_id, Nfaces, 'Nfaces', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading integer: Nfaces"
       STOP
    ENDIF

    ALLOCATE (Mesh%T(Nelems, Nnodesperelem))
    ALLOCATE (Mesh%X(Nnodes, ndim))

    ALLOCATE (Mesh%Tb(Nextfaces, Nnodesperface))
    ALLOCATE (Mesh%boundaryFlag(Nextfaces))

    CALL HDF5_array2D_reading_int(group_id, Mesh%T, 'T', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading mesh connectivity T"
       STOP
    ENDIF
    CALL HDF5_array2D_reading_int(group_id, Mesh%Tb, 'Tb', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundary connectivity Tb"
       STOP
    ENDIF
    CALL HDF5_array1D_reading_int(group_id, Mesh%boundaryFlag, 'boundaryFlag', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading boundaryFlag"
       STOP
    ENDIF
    CALL HDF5_array2D_reading(group_id, Mesh%X, 'X', ierr)
    IF (IERR .NE. 0) THEN
       WRITE (6, *) "Error reading coordinate matrix X"
       STOP
    ENDIF

    CALL HDF5_group_close(group_id, ierr)

    CALL HDF5_close(file_id)

    !************************************************************************
    !   CONFIRMATION MESSAGE FOR THE USER
    !************************************************************************
    IF(MPIvar%glob_id .EQ. 0) THEN
      WRITE (6, *) "Mesh read from solution file: ", TRIM(ADJUSTL(fname_complete))
    ENDIF

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

  endsubroutine HDF5_load_mesh_from_solution

  !**********************************************************************
  ! Load solution in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_load_solution(fname)

    USE LinearAlgebra, ONLY: tensorsumint, colint,col
    IMPLICIT NONE

    CHARACTER(LEN=1000) :: fname
    CHARACTER(len=100), POINTER :: mod_ptr
    CHARACTER(len=100), TARGET :: model_string

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

    fname_complete = TRIM(ADJUSTL(fname))//'.h5'


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
       IF (switch%target_variable /= 0) THEN           
         CALL HDF5_group_open(group_id, 'physics', group_id2, ierr)
         IF (switch%target_variable == 1) THEN
            CALL HDF5_real_reading(group_id2, phys%puff, 'puff')
         ELSEIF (switch%target_variable == 2) THEN
            CALL HDF5_real_reading(group_id2, phys%Re, 'recycling')
         ELSEIF (switch%target_variable == 3) THEN
            CALL HDF5_real_reading(group_id2, phys%puff, 'puff')
            CALL HDF5_real_reading(group_id2, phys%impurity_concentration, 'impurity_concentration')
         ENDIF

         IF (time%it .GT. 1) THEN
            CALL HDF5_real_reading(group_id2, phys%feedback_integral_error, 'feedback_integral_error')
            CALL HDF5_real_reading(group_id2, phys%feedback_previous_error, 'feedback_previous_error')
            IF (switch%target_variable == 3) THEN
               CALL HDF5_real_reading(group_id2, phys%feedback_integral_error_xpr, 'feedback_integral_error_xpr')
               CALL HDF5_real_reading(group_id2, phys%feedback_previous_error_xpr, 'feedback_previous_error_xpr')
            ENDIF
         ENDIF
         CALL HDF5_group_close(group_id2, ierr)
       ENDIF
    END IF
    CALL HDF5_group_close(group_id, ierr)


    ! Check if the readed solution corresponds to the right model
    IF (simpar%model .NE. model_string) THEN
       WRITE (6, *) "Wrong model in loaded solution | Loaded model: ", model_string, " | Current model: ", simpar%model
       STOP
    ENDIF
    CALL HDF5_group_open(file_id, 'solution', group_id, ierr)
    CALL HDF5_array1D_reading(group_id, sol%u, 'u')
    CALL HDF5_array1D_reading(group_id, sol%u_tilde, 'u_tilde')
    CALL HDF5_array1D_reading(group_id, sol%q, 'q')
    CALL HDF5_group_close(group_id, ierr)

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
       PRINT *, 'Solution read from file: ', TRIM(ADJUSTL(fname_complete))
       PRINT *, '        '
    END IF

  ENDSUBROUTINE HDF5_load_solution

  !**********************************************************************
  ! Save HDG matrix (CSR) in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_CSR_matrix(fname)

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

  ENDSUBROUTINE HDF5_save_CSR_matrix

  !**********************************************************************
  ! Save HDG vector (CSR) in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_CSR_vector(fname)

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

  ENDSUBROUTINE HDF5_save_CSR_vector

  !**********************************************************************
  ! Save 3D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_array(Arr, fname)

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

  ENDSUBROUTINE HDF5_save_array

  !**********************************************************************
  ! Save 2D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_matrix(Mat, fname)

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

  ENDSUBROUTINE HDF5_save_matrix

  !**********************************************************************
  ! Save 1D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_vector(Vec, fname)

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

  ENDSUBROUTINE HDF5_save_vector

  !**********************************************************************
  ! Save 1D array in HDF5 file format
  !**********************************************************************
  SUBROUTINE HDF5_save_vector_int(Vec, fname)

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

  ENDSUBROUTINE HDF5_save_vector_int


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
  ENDSUBROUTINE copy_file

END MODULE in_out
