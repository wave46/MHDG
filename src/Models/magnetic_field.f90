!*****************************************
! project: MHDG
! file: magnetic_field.f90
! date: 06/04/2020
! Definition of the magnetic field
! The magnetic field is defined at the mesh
! nodes
!*****************************************

MODULE magnetic_field
  USE prec_const
  USE globals
  USE in_out
  USE HDF5_io_module
  USE HDF5
  USE interpolation
CONTAINS

  !**********************************************
  ! Allocate storing space for the magnetic field
  !**********************************************
  SUBROUTINE initialize_magnetic_field
    INTEGER  :: nnodes
#ifdef TOR3D
    nnodes = Mesh%Nnodes*Mesh%Nnodes_toroidal
#else
    nnodes = Mesh%Nnodes
#endif
    ! Allocate storing space in phys
    ALLOCATE (phys%B(nnodes, 3))
    ALLOCATE (phys%magnetic_flux(nnodes))
    ALLOCATE (phys%magnetic_psi(nnodes))
#ifdef KEQUATION
    ALLOCATE (phys%omega(nnodes))
    ALLOCATE (phys%q_cyl(nnodes))
#endif
    IF (switch%ohmicsrc) THEN
       ALLOCATE (phys%Jtor(nnodes))
    END IF
    IF ((switch%RMP).OR.(switch%Ripple)) THEN
       ALLOCATE (phys%Bperturb(nnodes, 3))
    END IF
  END SUBROUTINE initialize_magnetic_field

  SUBROUTINE load_magnetic_field_Jtor
      CALL load_magnetic_field()
      IF (switch%ohmicsrc) CALL loadJtorMap()
  ENDSUBROUTINE load_magnetic_field_Jtor

  !**********************************************
  ! Load magnetic field in phys
  !**********************************************
  SUBROUTINE load_magnetic_field

    phys%B = 0.
    phys%magnetic_flux = 0.
    phys%magnetic_psi = 0.
#ifdef KEQUATION
    phys%omega = 0.
    phys%q_cyl = 0.
#endif

    SELECT CASE (switch%testcase)
    CASE (1:49)
       ! Analytic definition of the magnetic field
       CALL load_magnetic_field_analytical

    CASE (50:59)
       ! Magnetic field loaded from file in a cartesian grid
       ! Interpolation is needed
       IF (input%field_from_grid) THEN
          CALL load_magnetic_field_grid
       ELSE
          CALL load_magnetic_field_nodes
       ENDIF

    CASE (60:69)

       ! Analytic definition of the magnetic field
       CALL load_magnetic_field_analytical

    CASE (70:79)
       ! Magnetic field loaded from file in the mesh nodes
       CALL load_magnetic_field_nodes

    CASE (80:89)
       ! Magnetic field loaded from file in a cartesian grid
       ! Interpolation is needed
       CALL load_magnetic_field_grid


    END SELECT
#ifdef TOR3D
    ! Magnetic perturbations
    ! RMP part testcase 60-69.
    ! We need to be here to have refElTor%Nodes1d, refElTor%coord1d and numer%ntor (for tdiv)
    IF ((switch%RMP).OR.(switch%Ripple)) THEN
       phys%Bperturb = 0.
       CALL addMagneticPerturbation()
    END IF
#endif
    ! Adimensionalization of the magnetic field
    !phys%B = phys%B/phys%B0
  END SUBROUTINE load_magnetic_field

  !***********************************************************************
  ! Magnetic field defined analytically
  !***********************************************************************
  SUBROUTINE load_magnetic_field_analytical
    REAL*8   :: x(Mesh%Nnodes), y(Mesh%Nnodes)  ! Coordinates in the plane
#ifdef TOR3D
    REAL*8      :: t(Mesh%Nnodes_toroidal), tt       ! Toroidal coordinates
    INTEGER*4   :: j, N1D
#endif
    REAL*8                  :: xc, yc, R0, q, r
    REAL*8                  :: xmax, xmin, ymax, ymin, xm, ym
    REAL*8                  :: xx, yy, B0, xr, yr
    INTEGER*4               :: i, ind, N2D

    x = Mesh%X(:, 1)
    y = Mesh%X(:, 2)
    N2d = SIZE(X, 1)
#ifdef TOR3D
    t = Mesh%toroidal
    N1d = SIZE(t, 1)
#endif
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    xm = 0.5*(xmax + xmin)
    ym = 0.5*(ymax + ymin)

    xc = 0.
    yc = 0.
    DO i = 1, N2d
       xx = x(i)
       yy = y(i)
       ind = i
#ifdef TOR3D
       DO j = 1, N1d
          tt = t(j)
          ind = (j - 1)*N2d+i
#endif
          SELECT CASE (switch%testcase)
          CASE (1)
             IF (switch%axisym) THEN
                WRITE (6, *) "This is NOT an axisymmetric test case!"
                STOP
             END IF
             ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
             phys%B(ind, 1) = (yy - yc)
             phys%B(ind, 2) = (-xx + xc)
             phys%B(ind, 3) = 1.

          CASE (2)
             IF (.NOT. switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
             phys%B(ind, 1) = (yy - ym)/xx
             phys%B(ind, 2) = (-xx + xm)/xx
             phys%B(ind, 3) = 1.

          CASE (5)
             IF (switch%axisym) THEN
                WRITE (6, *) "This is NOT an axisymmetric test case!"
                STOP
             END IF
             !
             phys%B(ind, 1) = 0.
             phys%B(ind, 2) = 0.
             phys%B(ind, 3) = 1.
          CASE (6:7)
             IF (.NOT.switch%axisym) THEN
                WRITE (6, *) "This is an axisymmetric test case!"
                STOP
             END IF
             ! Axysimmetric case
             phys%B(ind, 1) = 0.
             phys%B(ind, 2) = 0.
             phys%B(ind, 3) = 1.+xx
          CASE (50:59)
             WRITE (6, *) "Error in defineMagneticField: you should not be here!"
             STOP
          CASE (60:68)

             ! Circular case with limiter
             R0 = geom%R0
             q = geom%q
             B0 = 2!*0.1522
             xr = xx*phys%lscale
             yr = yy*phys%lscale

             r = SQRT((xr - R0)**2 + yr**2)
             phys%B(ind, 1) = -B0*yr/(xr*q*SQRT(1 - (r/R0)**2))
             phys%B(ind, 2) = B0*(xr - R0)/(xr*q*SQRT(1 - (r/R0)**2))
             phys%B(ind, 3) = B0*R0/xr
          CASE DEFAULT
             WRITE (6, *) "Error! Test case not valid"
             STOP
          END SELECT
#ifdef TOR3D
       END DO
#endif
    END DO

  END SUBROUTINE load_magnetic_field_analytical

  !***********************************************************************
  ! Magnetic field loaded by a hdf5 file in a cartesian grid
  !***********************************************************************
  SUBROUTINE load_magnetic_field_grid()



    USE reference_element

    INTEGER                           :: i, ierr, ip, jp, ind, k
#ifdef TOR3D
    INTEGER                           :: j
#endif
    INTEGER(HID_T)                    :: file_id
    REAL*8, POINTER, DIMENSION(:, :)  :: r2D, z2D, flux2D, Br2D, Bz2D, Bphi2D
!!! variables for computation derivatives of the flux
    INTEGER                           :: iel, inode
    REAL*8                            :: shapeFunctions(refElpol%Nnodes2D,refElpol%Nnodes2D,3)
    REAL*8                            :: Xel(refElpol%Nnodes2D,2)        !only for 2D so far
    REAL*8                            :: J11(refElpol%Nnodes2D),J12(refElpol%Nnodes2D)
    REAL*8                            :: J21(refElpol%Nnodes2D),J22(refElpol%Nnodes2D)
    REAL*8                            :: iJ11(refElpol%Nnodes2D),iJ12(refElpol%Nnodes2D)
    REAL*8                            :: iJ21(refElpol%Nnodes2D),iJ22(refElpol%Nnodes2D)
    REAL*8                            :: Nxn(refElpol%Nnodes2D),Nyn(refElpol%Nnodes2D)
    REAL*8                            :: detJ(refElpol%Nnodes2D)
    REAL*8                            :: coord2D_fixed(refElpol%Nnodes2D,2)      ! applying some shift to the third node of thriangle to avoid infinite derivative
!!! end of variables for computation derivatives of the flux
    REAL*8, ALLOCATABLE, DIMENSION(:) :: xvec, yvec
    REAL*8                            :: x, y
    REAL*8                            :: Br, Bz, Bt, flux, psiSep, dt_ME,t_ME
    CHARACTER(LEN=1000) :: fname
    CHARACTER(50)  :: nit
    INTEGER                            :: min_ind(2)
#ifdef KEQUATION
    REAL*8                            :: q_cyl, omega,a    
#endif


    IF (utils%printint > 0) THEN
       IF(MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*           LOADING MAGNETIC FIELD              *'
          WRITE (6, *) '*************************************************'
       ENDIF
    END IF

    ! Read file
    IF (switch%testcase>=50 .AND. switch%testcase<60) THEN
       ! WEST case
       ! Dimensions of the file storing the magnetic field for West
       ip = input%field_dimensions(1)
       jp = input%field_dimensions(2)
       !ip = 541
       !jp = 391
       IF (switch%ME .EQV. .FALSE.)  THEN !if not a moving equilibrium simulation
          fname = input%field_path
       ELSE
          fname = input%field_path
          WRITE(nit, "(i10)") INT(time%it + 1)
          nit = TRIM(ADJUSTL(nit))
          k = INDEX(nit, " ") -1
          fname = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
       ENDIF
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(6,*) 'Magnetic field loaded from file: ', TRIM(ADJUSTL(fname))
       ENDIF
       CALL HDF5_open(fname, file_id, IERR)
    ELSEIF (switch%testcase>=80 .AND. switch%testcase<90) THEN
       ! ITER case
       IF (switch%ME .EQV. .FALSE.) THEN !if not a moving equilibrium simulation
          !ip = 513
          !jp = 257
          !fname = 'ITER_2008_MagField.h5'
          fname = 'ITER_135011_0000.h5'
       ELSE
          fname = 'B_field_exp/ITER_135011'
          WRITE(nit, "(i10)") INT(time%it + 1)
          nit = TRIM(ADJUSTL(nit))
          k = INDEX(nit, " ") -1
          fname = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
       ENDIF
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(6,*) 'Magnetic field loaded from file: ', TRIM(ADJUSTL(fname))
       ENDIF
       CALL HDF5_open(fname, file_id, IERR)
       CALL HDF5_integer_reading(file_id, ip, 'ip')
       CALL HDF5_integer_reading(file_id, jp, 'jp')
    ENDIF

    ALLOCATE (r2D(ip, jp))
    ALLOCATE (z2D(ip, jp))
    ALLOCATE (flux2D(ip, jp))
    ALLOCATE (Br2D(ip, jp))
    ALLOCATE (Bz2D(ip, jp))
    ALLOCATE (Bphi2D(ip, jp))
    CALL HDF5_array2D_reading(file_id, r2D, 'r2D')
    CALL HDF5_array2D_reading(file_id, z2D, 'z2D')
    CALL HDF5_real_reading(file_id, psiSep, 'psiSep')
    CALL HDF5_array2D_reading(file_id, flux2D, 'flux2D')
    CALL HDF5_array2D_reading(file_id, Br2D, 'Br2D')
    CALL HDF5_array2D_reading(file_id, Bz2D, 'Bz2D')
    CALL HDF5_array2D_reading(file_id, Bphi2D, 'Bphi2D')
    IF (switch%ME ) THEN
       CALL HDF5_real_reading(file_id, dt_ME, 'dt')
       CALL HDF5_real_reading(file_id, t_ME, 'time')
    ENDIF
    CALL HDF5_close(file_id)

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale

    !finding axis
    min_ind = MINLOC(flux2D)
    phys%r_axis = r2D(min_ind(1),min_ind(2))
    phys%z_axis = z2D(min_ind(1),min_ind(2))

    ! Min and Max flux for inizialization
    !phys%Flux2Dmin = minval(flux2D)
    !phys%Flux2Dmax = maxval(flux2D)

    ! Interpolate
    ALLOCATE (xvec(jp))
    ALLOCATE (yvec(ip))
    xvec = r2D(1, :)
    yvec = z2D(:, 1)
    DO i = 1, Mesh%Nnodes
       x = Mesh%X(i, 1)
       y = Mesh%X(i, 2)
       Br = interpolate(ip, yvec, jp, xvec, Br2D, y, x, 1e-12)
       Bz = interpolate(ip, yvec, jp, xvec, Bz2D, y, x, 1e-12)
       Bt = interpolate(ip, yvec, jp, xvec, Bphi2D, y, x, 1e-12)
       flux = interpolate(ip, yvec, jp, xvec, flux2D, y, x, 1e-12)
#ifdef KEQUATION
       omega = simpar%refval_charge/simpar%refval_mass*SQRT(Br**2+Bz**2+Bt**2)*simpar%refval_time
       a = SQRT((x-phys%r_axis)**2+(y-phys%z_axis)**2)
       q_cyl = ABS(Bt)*a/SQRT(Br**2+Bz**2)/x
       IF(q_cyl>1.e4) q_cyl = 1.e4
       IF(q_cyl<1.) q_cyl = 1.

#endif
       ind = i
#ifdef TOR3D
       DO j = 1, Mesh%Nnodes_toroidal
          ind = (j - 1)*Mesh%Nnodes + i
#endif
          phys%B(ind, 1) = Br
          phys%B(ind, 2) = Bz
          phys%B(ind, 3) = Bt
          phys%magnetic_flux(ind) = flux
#ifdef KEQUATION
          phys%omega(ind) = omega
          phys%q_cyl(ind) = q_cyl
#endif
#ifdef TOR3D
       END DO
#endif
    END DO
    ! Field from fluxes (ONLY 2D, ONLY triangles checked)
    ! gives nan at third point of the triangle, because its eta coordinate equal to straight 1.0

    IF (input%compute_from_flux) THEN
       coord2D_fixed =  refElpol%coord2d
       coord2D_fixed(3,2) = coord2D_fixed(3,2)-1.e-10 !! dirty trick, need to solve it later
       CALL compute_shape_functions_at_points(refElpol,coord2D_fixed,shapeFunctions)
       DO iel = 1, Mesh%Nelems
          ! taking coordinates for given element
          Xel = Mesh%X(Mesh%T(iel,:),:)
          !Jacobian computations
          J11 = MATMUL(shapeFunctions(:,:,2),Xel(:,1))                           ! nnodes x 1
          J12 = MATMUL(shapeFunctions(:,:,2),Xel(:,2))                           ! nnodes x 1
          J21 = MATMUL(shapeFunctions(:,:,3),Xel(:,1))                          ! nnodes x 1
          J22 = MATMUL(shapeFunctions(:,:,3),Xel(:,2))                          ! nnodes x 1
          detJ = J11*J22 - J21*J12                    ! determinant of the Jacobian
          iJ11 = J22/detJ
          iJ12 = -J12/detJ
          iJ21 = -J21/detJ
          iJ22 = J11/detJ
          DO inode = 1, Mesh%Nnodesperelem
             ! x and y derivatives of the shape functions
             Nxn = iJ11(inode)*shapeFunctions(inode,:,2) + iJ12(inode)*shapeFunctions(inode,:,3)
             Nyn = iJ21(inode)*shapeFunctions(inode,:,2) + iJ22(inode)*shapeFunctions(inode,:,3)
             ! Remember about 2pi
             Br = -1.*dot_PRODUCT(Nyn,phys%magnetic_flux(Mesh%T(iel,:)))/Xel(inode,1)/simpar%refval_length**2
             Bz = dot_PRODUCT(Nxn,phys%magnetic_flux(Mesh%T(iel,:)))/Xel(inode,1)/simpar%refval_length**2
             phys%B(Mesh%T(iel,inode),1) = Br
             phys%B(Mesh%T(iel,inode),2) = Bz
          ENDDO
       ENDDO
       IF (input%divide_by_2pi) THEN
          phys%B(:,1) = phys%B(:,1)/2./PI
          phys%B(:,2) = phys%B(:,2)/2./PI
       ENDIF
    ENDIF


    ! Min and Max flux for inizialization
    phys%Flux2Dmin = MINVAL(phys%magnetic_flux)
    phys%Flux2Dmax = MAXVAL(phys%magnetic_flux)

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%Flux2Dmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%Flux2Dmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

    ! Magnetic flux normalized to separatrix: PSI
    phys%magnetic_psi = (phys%magnetic_flux - phys%Flux2Dmin)/(psiSep - phys%Flux2Dmin)

    IF (switch%ME) THEN
       time%dt_ME = dt_ME
       time%t_ME = t_ME
    ENDIF
    ! Free memory
    DEALLOCATE (Br2D, Bz2D, Bphi2D, xvec, yvec)
    DEALLOCATE (r2D, z2D, flux2D)
    NULLIFY (r2D, z2D, flux2D, Br2D, Bz2D, Bphi2D)

  END SUBROUTINE load_magnetic_field_grid

  !***********************************************************************
  ! Magnetic field loaded by a hdf5 file in the nodes !TODO modify for 3D
  !***********************************************************************
  SUBROUTINE load_magnetic_field_nodes()


    USE MPI_OMP, only: MPIvar
    INTEGER        ::  ierr, k
#ifdef KEQUATION
    INTEGER        ::  i
#endif
    CHARACTER(LEN=1000) :: fname = 'Evolving_equilibrium'
    CHARACTER(50)  :: npr, nid, nit
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id
    REAL*8, POINTER, DIMENSION(:) :: Br, Bz, Bt, flux
    REAL*8            :: psiSep,dt_ME,t_ME
    INTEGER  :: nnodes
    INTEGER                            :: min_ind(1)
#ifdef PARALL
    REAL*8                              :: minflux_in(2), minflux_out(2)
    INTEGER                           :: my_rank
#endif
#ifdef TOR3D
    nnodes = Mesh%Nnodes*Mesh%Nnodes_toroidal
#else
    nnodes = Mesh%Nnodes
#endif

    ALLOCATE (flux(nnodes))
    ALLOCATE (Br(nnodes))
    ALLOCATE (Bz(nnodes))
    ALLOCATE (Bt(nnodes))

    IF (switch%ME .EQV. .FALSE.)  THEN !if not a moving equilibrium simulation
       fname = TRIM(ADJUSTL(input%field_path))
    ELSE
       fname = TRIM(ADJUSTL(input%field_path))
       WRITE(nit, "(i10)") INT(time%it + 1)
       nit = TRIM(ADJUSTL(nit))
       k = INDEX(nit, " ") -1
       fname = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))
    ENDIF


    IF (MPIvar%glob_size .GT. 1) THEN
       WRITE (nid, *) MPIvar%glob_id + 1
       WRITE (npr, *) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'.h5'
    END IF
    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) 'Magnetic field loaded from file: ', TRIM(ADJUSTL(fname_complete))
    ENDIF
    ! Read file
    CALL HDF5_open(fname_complete, file_id, IERR)
    CALL HDF5_array1D_reading(file_id, Br, 'Br')
    CALL HDF5_array1D_reading(file_id, Bz, 'Bz')
    CALL HDF5_array1D_reading(file_id, Bt, 'Bt')
    CALL HDF5_real_reading(file_id, psiSep, 'psiSep')
    CALL HDF5_array1D_reading(file_id, flux, 'flux')
    IF (switch%ME ) THEN
       CALL HDF5_real_reading(file_id, dt_ME, 'dt')
       CALL HDF5_real_reading(file_id, t_ME, 'time')
    ENDIF
    CALL HDF5_close(file_id)

    phys%B(:, 1) = Br
    phys%B(:, 2) = Bz
    phys%B(:, 3) = Bt
    phys%magnetic_flux = flux

    ! Min and Max flux for inizialization
    phys%Flux2Dmin = MINVAL(phys%magnetic_flux)
    phys%Flux2Dmax = MAXVAL(phys%magnetic_flux)

    !finding magnetic axis
    !finding axis
    min_ind = MINLOC(flux)

    phys%r_axis = Mesh%X(min_ind(1),1)
    phys%z_axis = Mesh%X(min_ind(1),2)

#ifdef PARALL

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

    minflux_in(1) = phys%Flux2Dmin
    minflux_in(2) = my_rank
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    CALL MPI_ALLREDUCE(minflux_in,minflux_out,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(phys%r_axis,1,MPI_REAL8,INT(minflux_out(2)),MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(phys%z_axis,1,MPI_REAL8,INT(minflux_out(2)),MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

#endif

#ifdef PARALL
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%Flux2Dmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%Flux2Dmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

    ! Magnetic flux normalized to separatrix: PSI
    phys%magnetic_psi = (phys%magnetic_flux - phys%Flux2Dmin)/(psiSep - phys%Flux2Dmin)

#ifdef KEQUATION
    DO i = 1, Mesh%Nnodes
       phys%omega(i) = simpar%refval_charge/simpar%refval_mass*SQRT(Br(i)**2+Bz(i)**2+Bt(i)**2)*simpar%refval_time

       phys%q_cyl(i) = ABS(Bt(i))*SQRT((Mesh%X(i,1)-phys%r_axis)**2+(Mesh%X(i,2)-phys%z_axis)**2)/SQRT(Br(i)**2+Bz(i)**2)/Mesh%X(i,1)
       phys%q_cyl(i) = MAX(phys%q_cyl(i),1.)
       phys%q_cyl(i) = MIN(phys%q_cyl(i),1e4)
    ENDDO
    WRITE(6,*) 'r_axis', phys%r_axis*simpar%refval_length
    WRITE(6,*) 'z_axis', phys%z_axis*simpar%refval_length

#endif
    IF (switch%ME) THEN
       time%dt_ME = dt_ME
       time%t_ME = t_ME
    ENDIF

    DEALLOCATE (Br, Bz, Bt, flux)
  END SUBROUTINE load_magnetic_field_nodes

#ifdef TOR3D
  !********************************
  ! A subroutine adding RMP or Ripple to an equilibrium not normalized B field
  ! Input:
  !        x, y, t: coordinates x, y and toroidal one
  !        b: (R,Z,Phi) in: computation of B at equilibrium (NOT normalized). out: normalized b with perturbation
  !********************************
  SUBROUTINE addMagneticPerturbation()
    REAL*8                                   :: xc,yc
    REAL*8                                   :: xmax,xmin,ymax,ymin,xm,ym
    REAL*8                                   :: xx,yy,tt,BB
    INTEGER*4                                :: i,j,ind,N2D,N1D
    REAL*8,  DIMENSION(Mesh%Nnodes)          :: x, y
    REAL*8,  DIMENSION(Mesh%Nnodes_toroidal) :: t
    REAL*8, ALLOCATABLE                      :: brmp(:,:), bripple(:,:)
    REAL*8, DIMENSION(2,2)                   :: coilCoord
    INTEGER                                  :: elDiscr, rowNb


    x = Mesh%X(:,1)
    y = Mesh%X(:,2)
    t = Mesh%toroidal
    N2d = SIZE(x,1)
    N1d = SIZE(t,1)
    xmax = Mesh%xmax
    xmin = Mesh%xmin
    ymax = Mesh%ymax
    ymin = Mesh%ymin
    !xmax = 2314.57127827459!Mesh%xmax
    !xmin = 1262.49342451341!Mesh%xmin
    !ymax = 526.038926880589!Mesh%ymax
    !ymin = -526.038926880589!Mesh%ymin
    xm = 0.5*(xmax+xmin)
    ym = 0.5*(ymax+ymin)

    xc = -0.5
    yc = -0.5

    ! RMP part testcase 61-69
    ALLOCATE(brmp(N2d*N1d,3))
    ALLOCATE(bripple(N2d*N1d,3))
    brmp = 0.
    bripple = 0.
    IF (switch%RMP) THEN ! RMP
       elDiscr = 16
       ALLOCATE(magn%coils_rmp(magn%nbCoils_rmp*4*elDiscr,6,magn%nbRow)) !4 because square, 6: 2 positions (start/stop) and 3 coordinates
       IF (magn%nbRow.EQ.2) THEN
          ! Upper row
          rowNb = 1
          coilCoord(1,1) = 0.95*xmax ! R coordinate top, upper row
          coilCoord(2,1) = ym + 3./4.*(ymax - ym) ! Z coordinate top, upper row
          coilCoord(1,2) = 1.05*xmax ! R coordinate bottom, upper row
          coilCoord(2,2) = ym + 1./4.*(ymax - ym) ! Z coordinate bottom, upper row
          CALL calcRMPField(brmp, coilCoord, rowNb, elDiscr)
          ! Lower row
          rowNb = 2
          coilCoord(1,1) = 1.05*xmax ! R coordinate top, lower row
          coilCoord(2,1) = ym - 1./4.*(ymax - ym) ! Z coordinate top, lower row
          coilCoord(1,2) = 0.95*xmax ! R coordinate bottom, lower row
          coilCoord(2,2) = ym - 3./4.*(ymax - ym) ! Z coordinate bottom, lower row
          CALL calcRMPField(brmp, coilCoord, rowNb, elDiscr)
       ELSE
          WRITE(6, *) 'TODO: not implemented yet'
          STOP
       ENDIF
    ENDIF

    IF (switch%Ripple) THEN ! Ripple
       elDiscr = 32
       ALLOCATE(magn%coils_ripple(magn%nbCoils_ripple*elDiscr,6))
       CALL calcRippleField(bripple, elDiscr)
    ENDIF

    DO i=1,N2d
       DO j=1,N1d
          xx = x(i)
          yy = y(i)
          tt = t(j)
          ind = (j-1)*N2d+i

          phys%B(ind,1) = phys%B(ind,1) + brmp(ind,1) + bripple(ind,1)
          phys%B(ind,2) = phys%B(ind,2) + brmp(ind,2) + bripple(ind,2)
          phys%B(ind,3) = phys%B(ind,3) + brmp(ind,3) + bripple(ind,3)

          BB = SQRT(phys%B(ind,1)**2 + phys%B(ind,2)**2 + phys%B(ind,3)**2)
          !phys%B(ind,1) = phys%B(ind,1)/BB
          !phys%B(ind,2) = phys%B(ind,2)/BB
          !phys%B(ind,3) = phys%B(ind,3)/BB
          phys%bperturb(ind,1) = (brmp(ind,1) + bripple(ind,1))!/BB
          phys%bperturb(ind,2) = (brmp(ind,2) + bripple(ind,2))!/BB
          phys%bperturb(ind,3) = (brmp(ind,3) + bripple(ind,3))!/BB
       END DO
    END DO
    DEALLOCATE(brmp,bripple)
  END SUBROUTINE addMagneticPerturbation


  !********************************
  ! A subroutine drawing the (rectangular) coils for RMP and computing B RMP field
  ! Input:
  !        brmp: (R,Z,Phi) computation of B (not normalized) generated by coils.
  !        coilCoord: (2,2) with in each colum the (R,Z) coordinates of top (bottom) of the coils in a row
  !           ex: (1,1): R top, (2,1): Z top, (1,2): R bottom, (2,2): Z bottom
  !        rowNb: current number of row
  !********************************
  SUBROUTINE calcRMPField(brmp, coilCoord, rowNb, elDiscr)
    REAL*8,  DIMENSION(:,:), INTENT(inout)       :: brmp
    REAL*8,  DIMENSION(2,2), INTENT(in)          :: coilCoord
    INTEGER, INTENT(in)                          :: rowNb,elDiscr

    REAL*8,  DIMENSION(Mesh%Nnodes)              :: x, y
    REAL*8,  DIMENSION(Mesh%Nnodes_toroidal)     :: t
    REAL*8                                       :: tmax, aLcoil, phiCoil
    REAL*8                                       :: spaceBetwCoil, dlxx, dlyy, dlzz
    REAL*8                                       :: x0, y0, z0, xxc, yyc, zzc
    REAL*8                                       :: rrx, rry, rrz, rr2, rr3
    REAL*8                                       :: Bx, By, Bz, dBx, dBy, dBz
    INTEGER                                      :: i, j, k, l, N2d, N1d, i2d, i1d, ind, par
    REAL*8, DIMENSION(magn%nbCoils_RMP, 4)       :: xx, yy, zz

    x = Mesh%X(:,1)
    y = Mesh%X(:,2)
    t = Mesh%toroidal
    ! Size of the domain
    tmax = 2*pi !numer%tmax
    N2d = SIZE(x,1)
    N1d = SIZE(t,1)

    aLcoil = (coilCoord(1,1) + coilCoord(1,2))/2.0*magn%torElongCoils_rmp
    spaceBetwCoil = (tmax - magn%nbCoils_rmp*magn%torElongCoils_rmp)/magn%nbCoils_rmp
    IF (spaceBetwCoil.LE.-0.01) THEN
       WRITE(6,*) "Error in RMP coils: widths of coils in a row to much for the chosen toroidal expansion"
       WRITE(6,*) "Negative space between coils: ", spaceBetwCoil
       STOP
    ENDIF
    ! Be careful to the direction of the coil (4 corners from bottom right then trigo)
    ! Creation of the n_coils_row by rotation (in row/ toroidal direction)
    DO i=1,magn%nbCoils_RMP
       phiCoil = tmax*(i-1)/magn%nbCoils_rmp + spaceBetwCoil

       xx(i,1) = coilCoord(1,2)*COS(phiCoil) - (-aLcoil/2.0)*SIN(phiCoil)
       yy(i,1) = coilCoord(1,2)*SIN(phiCoil) + (-alcoil/2.0)*COS(phiCoil)
       zz(i,1) = coilCoord(2,2)

       xx(i,2) = coilCoord(1,2)*COS(phiCoil) - (aLcoil/2.0)*SIN(phiCoil)
       yy(i,2) = coilCoord(1,2)*SIN(phiCoil) + (alcoil/2.0)*COS(phiCoil)
       zz(i,2) = coilCoord(2,2)

       xx(i,3) = coilCoord(1,1)*COS(phiCoil) - (aLcoil/2.0)*SIN(phiCoil)
       yy(i,3) = coilCoord(1,1)*SIN(phiCoil) + (alcoil/2.0)*COS(phiCoil)
       zz(i,3) = coilCoord(2,1)

       xx(i,4) = coilCoord(1,1)*COS(phiCoil) - (-aLcoil/2.0)*SIN(phiCoil)
       yy(i,4) = coilCoord(1,1)*SIN(phiCoil) + (-alcoil/2.0)*COS(phiCoil)
       zz(i,4) = coilCoord(2,1)
    ENDDO

    !Saving the coils coordinates for drawing
    !Loop on coils
    ind = 0
    DO i=1,magn%nbCoils_rmp
       !Loop on the 4 parts of a coil: 1->2, 2->3, 3->4, 4->1
       DO j=1,4
          IF (j.LE.3) THEN
             k = j + 1
          ELSE
             k = 1
          ENDIF
          dlxx = (xx(i,k)-xx(i,j))/elDiscr
          dlyy = (yy(i,k)-yy(i,j))/elDiscr
          dlzz = (zz(i,k)-zz(i,j))/elDiscr
          ! Loop on elements of coils for writing coils coordinates
          DO l=1,elDiscr
             ind = ind + 1
             magn%coils_rmp(ind,1,rowNb) = (l-1)*dlxx + xx(i,j)
             magn%coils_rmp(ind,3,rowNb) = (l-1)*dlyy + yy(i,j)
             magn%coils_rmp(ind,5,rowNb) = (l-1)*dlzz + zz(i,j)

             magn%coils_rmp(ind,2,rowNb) = l*dlxx + xx(i,j)
             magn%coils_rmp(ind,4,rowNb) = l*dlyy + yy(i,j)
             magn%coils_rmp(ind,6,rowNb) = l*dlzz + zz(i,j)
          ENDDO
       ENDDO
    ENDDO

    ! RMP hard coded with Biot and Savard law (see phd E. Nardon and ERGOS)
    DO i2d=1,N2d
       DO i1d=1,N1d
          ! In cartesian coordinates (for Biot and Savard)
          Bx = 0.0
          By = 0.0
          Bz = 0.0

          x0 = SQRT(x(i2d)**2 + y(i2d)**2)*COS(t(i1d))
          y0 = SQRT(x(i2d)**2 + y(i2d)**2)*SIN(t(i1d))
          z0 = y(i2d)
          ind = (i1d-1)*N2d+i2d

          !Loop on coils
          DO i=1,magn%nbCoils_rmp
             IF ((magn%parite.EQ.1).OR.(magn%parite.EQ.-1)) THEN
                par = magn%parite*(-1)**i
             ELSE
                par = (-1)**i
             ENDIF
             !Loop on the 4 parts of a coil: 1->2, 2->3, 3->4, 4->1
             DO j=1,4
                IF (j.LE.3) THEN
                   k = j + 1
                ELSE
                   k = 1
                ENDIF

                dlxx = (xx(i,k)-xx(i,j))/elDiscr
                dlyy = (yy(i,k)-yy(i,j))/elDiscr
                dlzz = (zz(i,k)-zz(i,j))/elDiscr
                ! Loop on elements of coils for writing coils coordinates
                ! rr is the r vector in Biot and Savard
                DO l=1,elDiscr
                   xxc = (l-1)*dlxx + xx(i,j)
                   yyc = (l-1)*dlyy + yy(i,j)
                   zzc = (l-1)*dlzz + zz(i,j)

                   rrx = x0 - xxc
                   rry = y0 - yyc
                   rrz = z0 - zzc

                   rr2 = rrx**2 + rry**2 + rrz**2
                   rr3 = rr2**(3./2)

                   dBx = par*(dlyy*rrz-dlzz*rry)/rr3/phys%lscale
                   dBy = par*(dlzz*rrx-dlxx*rrz)/rr3/phys%lscale
                   dBz = par*(dlxx*rry-dlyy*rrx)/rr3/phys%lscale

                   Bx = Bx + dBx
                   By = By + dBy
                   Bz = Bz + dBz
                ENDDO ! end elements of one coil
             ENDDO ! end 4 parts of one coil
          ENDDO ! end loop on coils
          ! Br, Bz, Bt
          brmp(ind,1) = brmp(ind,1) + magn%amp_rmp*(-Bx*SIN(t(i1d)) + By*COS(t(i1d)))
          brmp(ind,2) = brmp(ind,2) + magn%amp_rmp*Bz
          brmp(ind,3) = brmp(ind,3) + magn%amp_rmp*(Bx*COS(t(i1d)) + By*SIN(t(i1d)))
       ENDDO
    ENDDO
  END SUBROUTINE calcRMPField

  !********************************
  ! A subroutine drawing the toroidal coils for ripple and computing B ripple field
  ! Input:
  !        bripple: (R,Z,Phi) computation of B (not normalized) generated by coils minus averaged value of B.
  !        amp: amplitude of the magnetic field defined as the ratio of the coil current and the plasma current (in TK3X)
  !        triang: triangularity
  !        ellip: ellipticity
  !********************************
  SUBROUTINE calcRippleField(bripple, elDiscr)
    REAL*8,  DIMENSION(:,:), INTENT(inout)       :: bripple
    INTEGER, INTENT(in)                          :: elDiscr

    REAL*8,  DIMENSION(Mesh%Nnodes)              :: x, y
    REAL*8,  DIMENSION(Mesh%Nnodes_toroidal)     :: t
    REAL*8,  DIMENSION(Mesh%Nnodes,3)            :: bripple_av
    REAL*8                                       :: minRadius, tmax, thetaCoil, phiCoil, csteR
    REAL*8                                       :: dlxx, dlyy, dlzz
    REAL*8                                       :: x0, y0, z0
    REAL*8                                       :: rrx, rry, rrz, rr2, rr3
    REAL*8                                       :: Bx, By, Bz, dBx, dBy, dBz
    INTEGER                                      :: i, j, k, N2d, N1d, i2d, i1d, ind, par
    REAL*8, DIMENSION(:,:), ALLOCATABLE          :: xx, yy, zz

    x = Mesh%X(:,1)
    y = Mesh%X(:,2)
    t = Mesh%toroidal
    ! Size of the domain
    tmax = 2*pi !numer%tmax
    N2d = SIZE(x,1)
    N1d = SIZE(t,1)
    csteR = 1.4
    minRadius = (Mesh%xmax - Mesh%xmin)/2.0

    ALLOCATE(xx(magn%nbCoils_ripple, elDiscr))
    ALLOCATE(yy(magn%nbCoils_ripple, elDiscr))
    ALLOCATE(zz(magn%nbCoils_ripple, elDiscr))

    ! Creation of the N_coils by rotation (in row/ toroidal direction)
    ! Need to create coils a0round full torus to avoid B-field inconsistency
    DO i=1,magn%nbCoils_ripple
       ! Shift from 0 for first phi to avoid non-axisymmetry
       phiCoil = tmax*(i-1)/magn%nbCoils_ripple + tmax/(2*magn%nbCoils_ripple)
       ! Theta discretization on Ndiscr points and minor radius of 2*a for toroidal coils
       DO j = 1,elDiscr
          thetaCoil = 2*PI*j/elDiscr
          xx(i, j) = (geom%R0/phys%lscale + csteR*minRadius*COS(thetaCoil + magn%triang*SIN(thetaCoil)))*COS(phiCoil)
          yy(i, j) = (geom%R0/phys%lscale + csteR*minRadius*COS(thetaCoil + magn%triang*SIN(thetaCoil)))*SIN(phiCoil)
          zz(i, j) = magn%ellip*csteR*minRadius*SIN(thetaCoil)
       ENDDO
    ENDDO

    !Saving the coils coordinates for drawing
    !Loop on coils
    ind = 0
    DO i=1,magn%nbCoils_ripple
       DO j=1,elDiscr
          IF (j.LT.elDiscr) THEN
             k = j + 1
          ELSE
             k = 1
          ENDIF
          ind = ind + 1
          magn%coils_ripple(ind,1) = xx(i,j)
          magn%coils_ripple(ind,3) = yy(i,j)
          magn%coils_ripple(ind,5) = zz(i,j)

          magn%coils_ripple(ind,2) = xx(i,k)
          magn%coils_ripple(ind,4) = yy(i,k)
          magn%coils_ripple(ind,6) = zz(i,k)
       ENDDO
    ENDDO

    ! Ripple hard coded with Biot and Savard law (see phd E. Nardon and ERGOS)
    DO i2d=1,N2d
       DO i1d=1,N1d
          ! In cartesian coordinates (for Biot and Savard)
          Bx = 0.0
          By = 0.0
          Bz = 0.0

          x0 = SQRT(x(i2d)**2 + y(i2d)**2)*COS(t(i1d))
          y0 = SQRT(x(i2d)**2 + y(i2d)**2)*SIN(t(i1d))
          z0 = y(i2d)
          ind = (i1d-1)*N2d+i2d

          !Loop on coils
          DO i=1,magn%nbCoils_ripple
             par = 1 !(-1)**i
             ! Loop on elements of coils for writing coils coordinates
             ! rr is the r vector in Biot and Savard
             DO j=1,elDiscr
                IF (j.LT.elDiscr) THEN
                   k = j + 1
                ELSE
                   k = 1
                ENDIF

                rrx = x0 - xx(i,j)
                rry = y0 - yy(i,j)
                rrz = z0 - zz(i,j)

                dlxx = xx(i,k)-xx(i,j)
                dlyy = yy(i,k)-yy(i,j)
                dlzz = zz(i,k)-zz(i,j)

                rr2 = rrx**2 + rry**2 + rrz**2
                rr3 = rr2**(3./2)

                dBx = par*(dlyy*rrz-dlzz*rry)/rr3/phys%lscale
                dBy = par*(dlzz*rrx-dlxx*rrz)/rr3/phys%lscale
                dBz = par*(dlxx*rry-dlyy*rrx)/rr3/phys%lscale

                Bx = Bx + dBx
                By = By + dBy
                Bz = Bz + dBz
             ENDDO ! end elements of one coil
          ENDDO ! end loop on coils
          ! Br, Bz, Bt
          bripple(ind,1) = bripple(ind,1) + magn%amp_ripple*(-Bx*SIN(t(i1d)) + By*COS(t(i1d)))
          bripple(ind,2) = bripple(ind,2) + magn%amp_ripple*Bz
          bripple(ind,3) = bripple(ind,3) + magn%amp_ripple*(Bx*COS(t(i1d)) + By*SIN(t(i1d)))
       ENDDO
    ENDDO
    ! Compute the ripple average
    bripple_av = 0.
    DO i2d=1,N2d
       DO i1d=1,N1d
          ind = (i1d-1)*N2d+i2d
          bripple_av(i2d,1) = bripple_av(i2d,1) + bripple(ind,1)
          bripple_av(i2d,2) = bripple_av(i2d,2) + bripple(ind,2)
          bripple_av(i2d,3) = bripple_av(i2d,3) + bripple(ind,3)
       ENDDO
    ENDDO
    bripple_av = bripple_av/SIZE(t,1)
    ! Substract the mean along phi to the magnetic field to only keep the ripple
    ! minloc can be necessary to find the correct index in the global matrix bripple_av to substract to the local matrix bripple
    ! When bripple and phys%bripple are the same size, loc must be equal to ind (which unfortunately is NOT the case for a few
    ! index).
    DO i2d=1,N2d
       !loc = minloc(abs(mesh%X(:,1) - x(i2d)) + abs(mesh%X(:,2) - y(i2d)),1)
       ! Fortran 2008 only !!!
       !if (any(mesh%T(1:size(T,1)/2,:).eq.loc)) then
       !   loc = minloc(abs(mesh%X(:,1) - x(i2d)) + abs(mesh%X(:,2) - y(i2d)),1, BACK=.FALSE.)
       !else
       !   loc = minloc(abs(mesh%X(:,1) - x(i2d)) + abs(mesh%X(:,2) - y(i2d)),1, BACK=.TRUE.)
       !endif
       DO i1d=1,N1d
          ind = (i1d-1)*N2d+i2d
          !bripple(ind,1) = bripple(ind,1) - phys%bripple(loc, 1)
          !bripple(ind,2) = bripple(ind,2) - phys%bripple(loc, 2)
          !bripple(ind,3) = bripple(ind,3) - phys%bripple(loc, 3)
          bripple(ind,1) = bripple(ind,1) - bripple_av(i2d, 1)
          bripple(ind,2) = bripple(ind,2) - bripple_av(i2d, 2)
          bripple(ind,3) = bripple(ind,3) - bripple_av(i2d, 3)
       ENDDO
    ENDDO
    DEALLOCATE(xx,yy,zz)
  END SUBROUTINE calcRippleField
#endif

  ! Below are routines from Manuel MHDG v2.1. Copy as it without any check: TODO adapt it to global magnetic field

  SUBROUTINE loadJtorMap()




    INTEGER        :: i,ierr,ip,jp,ind, k
#ifdef TOR3D
    INTEGER        :: j
#endif
    INTEGER(HID_T) :: file_id

    CHARACTER(LEN=1000)    :: fname
    CHARACTER(70)        :: nit

    REAL*8,POINTER,DIMENSION(:,:) :: r2D,z2D,Jtor
    REAL*8,ALLOCATABLE,DIMENSION(:)   :: xvec,yvec
    REAL*8                            :: dt_ME,t_ME
    REAL*8                            :: x,y
    REAL*8,PARAMETER                  :: tol = 1.e-12


    IF (utils%printint > 0) THEN
       IF(MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*          LOADING TOROIDAL CURRENT             *'
          WRITE (6, *) '*************************************************'
       ENDIF
    END IF

    ! Allocate storing space in phys
    phys%Jtor = 0.

    ! Read file
    IF (switch%testcase>=50 .AND. switch%testcase<60) THEN
       ! WEST case
       ! Dimensions of the file storing the magnetic field for West
       ip =  input%jtor_dimensions(1)
       jp =  input%jtor_dimensions(2)
       !ip = 541
       !jp = 391
       IF(switch%ME .EQV. .FALSE.) THEN !if not a moving equilibrium simulation
          fname = input%jtor_path
       ELSE
          fname = input%jtor_path
          WRITE(nit, "(i10)") INT(time%it + 1)
          nit = TRIM(ADJUSTL(nit))
          k = INDEX(nit, " ") -1
          fname = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
       ENDIF
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(6,*) 'Toroidal current loaded from file: ', TRIM(ADJUSTL(fname))
       ENDIF
       CALL HDF5_open(fname, file_id, IERR)
    ELSEIF (switch%testcase>=80 .AND. switch%testcase<90) THEN
       ! ITER case
       IF(switch%ME .EQV. .FALSE.) THEN !if not a moving equilibrium simulation
          !fname = 'ITER_2008_MagField.h5'
          fname = 'ITER_135011_Jtor_0000.h5'
       ELSE
          fname = 'B_field_exp/ITER_135011_Jtor'
          WRITE(nit, "(i10)") INT(time%it + 1)
          nit = TRIM(ADJUSTL(nit))
          k = INDEX(nit, " ") -1
          fname = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
       ENDIF
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(6,*) 'Toroidal current loaded from file: ', TRIM(ADJUSTL(fname))
       ENDIF
       CALL HDF5_open(fname, file_id, IERR)
       CALL HDF5_integer_reading(file_id, ip, 'ip')
       CALL HDF5_integer_reading(file_id, jp, 'jp')
    ENDIF

    ALLOCATE(r2D(ip,jp))
    ALLOCATE(z2D(ip,jp))
    ALLOCATE(Jtor(ip,jp))

    CALL HDF5_array2D_reading(file_id,r2D,'r2D')
    CALL HDF5_array2D_reading(file_id,z2D,'z2D')
    CALL HDF5_array2D_reading(file_id,Jtor,'Jtor')
    IF (switch%ME) THEN
       CALL HDF5_real_reading(file_id, dt_ME, 'dt')
       CALL HDF5_real_reading(file_id, t_ME, 'time')
    ENDIF
    CALL HDF5_close(file_id)

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale

    ! Interpolate
    ALLOCATE(xvec(jp))
    ALLOCATE(yvec(ip))
    xvec = r2D(1,:)
    yvec = z2D(:,1)
    DO i = 1,Mesh%Nnodes
       x = Mesh%X(i,1)
       y = Mesh%X(i,2)
       ind = i
#ifdef TOR3D
       DO j = 1, Mesh%Nnodes_toroidal
          ind = (j - 1)*Mesh%Nnodes + i
#endif
          phys%Jtor(ind) = interpolate(ip, yvec,jp, xvec,Jtor, y,x, 1e-12)
#ifdef TOR3D
       END DO
#endif
    END DO

    !Compute Ip
    CALL computeIplasma()

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,*) 'I_p =  ', phys%I_p, '[MA]'
    ENDIF

    ! check that time is the same
    IF (switch%ME) THEN
       IF ((ABS(dt_ME - time%dt_ME) > TOL) .OR. (ABS(t_ME - time%t_ME) > TOL)) THEN
          WRITE(6,*) 'Time in current and in equilibrium files are different'
          STOP
       ENDIF
    ENDIF

    ! Free memory
    DEALLOCATE(r2D,z2D,Jtor,xvec,yvec)
    NULLIFY(r2D,z2D,Jtor)

  END SUBROUTINE loadJtorMap

  SUBROUTINE loadMagneticFieldFromExperimentalData()



    USE MPI_OMP, only: MPIvar
    INTEGER        :: ierr,k,ip,jp,i,ind
    REAL*8,POINTER,DIMENSION(:,:) :: r2D,z2D,flux2D,Br2D,Bz2D,Bphi2D
    REAL*8,ALLOCATABLE,DIMENSION(:)   :: xvec,yvec,Bmod
    REAL*8                            :: x,y, Br, Bz, Bt, flux
#ifdef TOR3D
    INTEGER        :: j

#endif
    CHARACTER(LEN=25) :: fname = 'WEST_54487'
    CHARACTER(10)  :: npr,nid,nit
    CHARACTER(len=1000) :: fname_complete
    INTEGER(HID_T) :: file_id

    ! Dimensions of the file storing the magnetic field for West
    ip = 457
    jp = 457
    ALLOCATE(r2D(ip,jp))
    ALLOCATE(z2D(ip,jp))
    ALLOCATE(flux2D(ip,jp))
    ALLOCATE(Br2D(ip,jp))
    ALLOCATE(Bz2D(ip,jp))
    ALLOCATE(Bphi2D(ip,jp))
    ALLOCATE(Bmod(Mesh%Nnodes))

    ! File name
    WRITE(nit, "(i10)") time%it+1
    nit = TRIM(ADJUSTL(nit))
    k = INDEX(nit, " ") -1

    IF (MPIvar%glob_size.GT.1) THEN
       WRITE(nid,*) MPIvar%glob_id+1
       WRITE(npr,*) MPIvar%glob_size
       fname_complete = TRIM(ADJUSTL(fname))//'_'//TRIM(ADJUSTL(nid))//'_'//TRIM(ADJUSTL(npr))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
    ELSE
       fname_complete = TRIM(ADJUSTL(fname))//'_'//REPEAT("0", 4 - k)//TRIM(ADJUSTL(nit))//'.h5'
    END IF

    WRITE(6,*) 'Magnetic field loaded from file: ', TRIM(ADJUSTL(fname_complete))

    ! Read file
    CALL HDF5_open(fname_complete,file_id,IERR)
    CALL HDF5_array2D_reading(file_id,r2D,'r2D')
    CALL HDF5_array2D_reading(file_id,z2D,'z2D')
    CALL HDF5_array2D_reading(file_id,flux2D,'flux2D')
    CALL HDF5_array2D_reading(file_id,Br2D,'Br2D')
    CALL HDF5_array2D_reading(file_id,Bz2D,'Bz2D')
    CALL HDF5_array2D_reading(file_id,Bphi2D,'Bphi2D')
    CALL HDF5_close(file_id)

    ! Apply length scale
    r2D = r2D/phys%lscale
    z2D = z2D/phys%lscale



    ! Interpolate
    ALLOCATE (xvec(jp))
    ALLOCATE (yvec(ip))
    xvec = r2D(1, :)
    yvec = z2D(:, 1)
    DO i = 1, Mesh%Nnodes
       x = Mesh%X(i,1)
       y = Mesh%X(i,2)
       Br = interpolate(ip, yvec, jp, xvec, Br2D, y, x, 1e-12)
       Bz = interpolate(ip, yvec, jp, xvec, Bz2D, y, x, 1e-12)
       Bt = interpolate(ip, yvec, jp, xvec, Bphi2D, y, x, 1e-12)
       flux = interpolate(ip, yvec, jp, xvec, flux2D, y, x, 1e-12)
       ind = i
#ifdef TOR3D
       DO j = 1, Mesh%Nnodes_toroidal
          ind = (j - 1)*Mesh%Nnodes + i
#endif
          phys%B(ind, 1) = Br
          phys%B(ind, 2) = Bz
          phys%B(ind, 3) = Bt
          phys%magnetic_flux(ind) = flux
#ifdef TOR3D
       END DO
#endif
    END DO

    ! Free memory
    DEALLOCATE(Br2D,Bz2D,Bphi2D,xvec,yvec)
    DEALLOCATE(r2D,z2D,flux2D)

  END SUBROUTINE loadMagneticFieldFromExperimentalData

  SUBROUTINE computeIplasma()
    REAL*8  	:: Xel(refElPol%Nnodes2D,2),xyg(refElPol%NGauss2D,2),dvolu
    REAL*8		:: Jtorel(refElPol%Nnodes2D),Jtorg(refElPol%NGauss2D)
    REAL*8		:: J11(refElPol%NGauss2D),J12(refElPol%NGauss2D)
    REAL*8		:: J21(refElPol%NGauss2D),J22(refElPol%NGauss2D)
    REAL*8		:: detJ(refElPol%NGauss2D)
    INTEGER	  :: iel,g
#ifdef PARALL
    INTEGER   :: ierr
#endif

    phys%I_p = 0.

    DO iel = 1, Mesh%Nelems
       ! Coordinates of the nodes of the element
       Xel = Mesh%X(Mesh%T(iel,:),:)

       ! Toroidal current of the nodes of the element
       Jtorel = phys%Jtor(Mesh%T(iel,:))

       ! Gauss points position
       xyg = MATMUL(refElPol%N2D,Xel)

       ! Toroidal current at Gauss points
       Jtorg = MATMUL(refElPol%N2D,Jtorel)

       ! Jacobian
       J11 = MATMUL(refElPol%Nxi2D,Xel(:,1))                             ! ng x 1
       J12 = MATMUL(refElPol%Nxi2D,Xel(:,2))                             ! ng x 1
       J21 = MATMUL(refElPol%Neta2D,Xel(:,1))                          ! ng x 1
       J22 = MATMUL(refElPol%Neta2D,Xel(:,2))                          ! ng x 1
       detJ = J11*J22 - J21*J12                 			           ! determinant of the Jacobian

#ifdef PARALL
       IF (Mesh%ghostElems(iel) .EQ. 0) THEN
#endif

          ! Loop in 2D Gauss points
          DO g = 1, refElPol%NGauss2D
             ! Integration weight
             dvolu = refElPol%gauss_weights2D(g)*detJ(g)

             ! Compute I plasma
             phys%I_p = phys%I_p + Jtorg(g)*dvolu*phys%lscale**2

          END DO

#ifdef PARALL
       ENDIF
#endif

    END DO

#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, phys%I_p, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    ! Toroidal current in [MA]
    phys%I_p = phys%I_p/1.e6

  END SUBROUTINE computeIplasma

  SUBROUTINE initialize_puff()

    IF (utils%printint > 0) THEN
       IF(MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) '*************************************************'
          WRITE (6, *) '*             INITIALIZING PUFF                 *'
          WRITE (6, *) '*************************************************'
       ENDIF
    END IF

#ifdef NEUTRAL
    IF (switch%ME .EQV. .FALSE.) THEN
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE(6,*) 'Puff is analytical'
       ENDIF
    ELSE
       CALL SetParticleSource()
    END IF
#endif
  ENDSUBROUTINE initialize_puff


SUBROUTINE SetParticleSource()

   CHARACTER(LEN=1000) :: fname, fname_density, fname_impurity, fname_zeff
   INTEGER(HID_T)    :: file_id
   INTEGER           :: qp, Nn2D
   REAL*8            :: lower, upper, nli, n_Gw, n_la, a = 2.
   REAL*8, POINTER, DIMENSION(:) :: puff_time, target_density_time, target_density_exp
   INTEGER           :: puff_len, density_len, impurity_concentration_len, zeff_len

   NULLIFY(puff_time, target_density_time, target_density_exp)

   fname = input%puff_path
   puff_len = input%puff_dimension
   fname_density = input%target_density_path
   density_len = input%target_density_dimension
   fname_impurity = input%impurity_concentration_path
   impurity_concentration_len = input%impurity_concentration_dimension

   ! Allocate storing space in phys (puff for WEST, 403 entries)
   IF (switch%testcase .GE. 50 .AND. switch%testcase .LE. 59) THEN
      IF (switch%target_variable .EQ. 0) THEN
         CALL load_puff_from_file(fname, puff_len)
         IF (switch%impurity_radiation)THEN
            CALL load_impurity_concentration(fname_impurity, impurity_concentration_len)
         ENDIF
      ELSEIF (switch%target_variable .EQ. 1) THEN
         CALL adjust_puff_to_target_density(fname_density, density_len)
         IF (switch%impurity_radiation)THEN
            CALL load_impurity_concentration(fname_impurity, impurity_concentration_len)
         ENDIF
      ELSEIF (switch%target_variable .EQ. 2) THEN
         !We first load the puff from file
         CALL load_puff_from_file(fname, puff_len)
         !Then we adjust wall recycling
         CALL adjust_recycling_to_target_density(fname_density, density_len)
         IF (switch%impurity_radiation)THEN
            CALL load_impurity_concentration(fname_impurity, impurity_concentration_len)
         ENDIF
      END IF
   END IF

   ! ITER puff: linear increase up to nli = 4.00E+19
   IF (switch%testcase .GE. 80 .AND. switch%testcase .LE. 89) THEN
       ! Puff feedback: check if central line integrated has reached the target value (4e19 for ITER ohmic)
       CALL compute_line_integrated_density(MINVAL(Mesh%X(:, 1)), MAXVAL(Mesh%X(:, 1)), 0.5/phys%lscale, 0.5/phys%lscale, nli)
       CALL adjust_ITER_puff(nli)
   END IF

END SUBROUTINE SetParticleSource

SUBROUTINE load_puff_from_file(fname, puff_len)
   CHARACTER(LEN=1000), INTENT(IN) :: fname
   INTEGER, INTENT(IN) :: puff_len
   INTEGER(HID_T) :: file_id
   REAL*8, POINTER, DIMENSION(:) :: puff_time
   INTEGER :: puff_time_idx
   INTEGER :: ierr

   ALLOCATE(puff_time(puff_len))
   ALLOCATE(phys%puff_exp(puff_len))

   ! Read file
   CALL HDF5_open(fname, file_id, ierr)
   CALL HDF5_array1D_reading(file_id, phys%puff_exp, 'puff')
   CALL HDF5_array1D_reading(file_id, puff_time, 'time')
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'Puff loaded from file: ', TRIM(ADJUSTL(fname))
   END IF
   CALL HDF5_close(file_id)

   ! Linear interpolation of puff
   puff_time_idx = binarySearch(puff_len, puff_time, time%t_ME, 1e-12)
   phys%puff = phys%puff_exp(puff_time_idx)*(puff_time(puff_time_idx+1)-time%t_ME)/(puff_time(puff_time_idx+1)-puff_time(puff_time_idx)) + &
               phys%puff_exp(puff_time_idx+1)*(time%t_ME-puff_time(puff_time_idx))/(puff_time(puff_time_idx+1)-puff_time(puff_time_idx))
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'puff =  ', phys%puff
   END IF
   DEALLOCATE(puff_time)
   NULLIFY(puff_time)
END SUBROUTINE load_puff_from_file

SUBROUTINE adjust_puff_to_target_density(fname_density, density_len)
   CHARACTER(LEN=1000), INTENT(IN) :: fname_density
   INTEGER, INTENT(IN) :: density_len
   REAL*8 :: x_lower, x_upper, y_lower, y_upper
   INTEGER(HID_T) :: file_id
   REAL*8, POINTER, DIMENSION(:) :: target_density_time, target_density_exp
   INTEGER :: target_density_idx
   REAL*8 :: target_density, nli
   INTEGER :: ierr

   ALLOCATE(target_density_time(density_len))
   ALLOCATE(target_density_exp(density_len))

   ! Read file
   CALL HDF5_open(fname_density, file_id, ierr)
   CALL HDF5_array1D_reading(file_id, target_density_exp, 'target_density')
   CALL HDF5_array1D_reading(file_id, target_density_time, 'time')
   CALL HDF5_real_reading(file_id, x_lower, 'x_lower')
   CALL HDF5_real_reading(file_id, x_upper, 'x_upper')
   CALL HDF5_real_reading(file_id, y_lower, 'y_lower')
   CALL HDF5_real_reading(file_id, y_upper, 'y_upper')
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'Target density loaded from file: ', TRIM(ADJUSTL(fname_density))
   END IF
   CALL HDF5_close(file_id)

   ! Linear interpolation of target density
   target_density_idx = binarySearch(density_len, target_density_time, time%t_ME, 1e-12)
   target_density = target_density_exp(target_density_idx)*(target_density_time(target_density_idx+1)-time%t_ME)/(target_density_time(target_density_idx+1)-target_density_time(target_density_idx)) + &
                    target_density_exp(target_density_idx+1)*(time%t_ME-target_density_time(target_density_idx))/(target_density_time(target_density_idx+1)-target_density_time(target_density_idx))
   target_density = target_density / 2.

   ! Compute line integrated density
   CALL compute_line_integrated_density(x_lower, x_upper, y_lower, y_upper, nli)

   ! Adjust puff using feedback
   CALL adjust_puff_feedback(target_density, nli)

   DEALLOCATE(target_density_time, target_density_exp)
   NULLIFY(target_density_time, target_density_exp)
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'puff =  ', phys%puff
   END IF
END SUBROUTINE adjust_puff_to_target_density

SUBROUTINE adjust_recycling_to_target_density(fname_density, density_len)
   CHARACTER(LEN=1000), INTENT(IN) :: fname_density
   INTEGER, INTENT(IN) :: density_len
   REAL*8 :: x_lower, x_upper, y_lower, y_upper
   INTEGER(HID_T) :: file_id
   REAL*8, POINTER, DIMENSION(:) :: target_density_time, target_density_exp
   INTEGER :: target_density_idx
   REAL*8 :: target_density, nli
   INTEGER :: ierr

   ALLOCATE(target_density_time(density_len))
   ALLOCATE(target_density_exp(density_len))

   ! Read file
   CALL HDF5_open(fname_density, file_id, ierr)
   CALL HDF5_array1D_reading(file_id, target_density_exp, 'target_density')
   CALL HDF5_array1D_reading(file_id, target_density_time, 'time')
   CALL HDF5_real_reading(file_id, x_lower, 'x_lower')
   CALL HDF5_real_reading(file_id, x_upper, 'x_upper')
   CALL HDF5_real_reading(file_id, y_lower, 'y_lower')
   CALL HDF5_real_reading(file_id, y_upper, 'y_upper')
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'Target density loaded from file: ', TRIM(ADJUSTL(fname_density))
   END IF
   CALL HDF5_close(file_id)

   ! Linear interpolation of target density
   target_density_idx = binarySearch(density_len, target_density_time, time%t_ME, 1e-12)
   target_density = target_density_exp(target_density_idx)*(target_density_time(target_density_idx+1)-time%t_ME)/(target_density_time(target_density_idx+1)-target_density_time(target_density_idx)) + &
                    target_density_exp(target_density_idx+1)*(time%t_ME-target_density_time(target_density_idx))/(target_density_time(target_density_idx+1)-target_density_time(target_density_idx))
   target_density = target_density / 2.

   ! Compute line integrated density
   CALL compute_line_integrated_density(x_lower, x_upper, y_lower, y_upper, nli)
   ! Adjust recycling using feedback
   CALL adjust_recycling_feedback(target_density, nli)

   DEALLOCATE(target_density_time, target_density_exp)
   NULLIFY(target_density_time, target_density_exp)
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'recycling =  ', phys%Re
   END IF

END SUBROUTINE adjust_recycling_to_target_density

SUBROUTINE load_impurity_concentration(fname_impurity, impurity_concentration_len)
   CHARACTER(LEN=1000), INTENT(IN) :: fname_impurity
   INTEGER, INTENT(IN) :: impurity_concentration_len
   INTEGER(HID_T) :: file_id
   REAL*8, POINTER, DIMENSION(:) :: impurity_concentration_time, impurity_concentration_exp
   INTEGER :: impurity_concentration_idx
   REAL*8 :: impurity_concentration
   INTEGER :: ierr

   ALLOCATE(impurity_concentration_time(impurity_concentration_len))
   ALLOCATE(impurity_concentration_exp(impurity_concentration_len))

   ! Read file
   CALL HDF5_open(fname_impurity, file_id, ierr)
   CALL HDF5_array1D_reading(file_id, impurity_concentration_exp, 'impurity_concentration')
   CALL HDF5_array1D_reading(file_id, impurity_concentration_time, 'time')
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'Impurity concentration loaded from file: ', TRIM(ADJUSTL(fname_impurity))
   END IF
   CALL HDF5_close(file_id)

   ! Linear interpolation of impurity concentration
   impurity_concentration_idx = binarySearch(impurity_concentration_len, impurity_concentration_time, time%t_ME, 1e-12)
   impurity_concentration = impurity_concentration_exp(impurity_concentration_idx)*(impurity_concentration_time(impurity_concentration_idx+1)-time%t_ME)/(impurity_concentration_time(impurity_concentration_idx+1)-impurity_concentration_time(impurity_concentration_idx)) + &
                            impurity_concentration_exp(impurity_concentration_idx+1)*(time%t_ME-impurity_concentration_time(impurity_concentration_idx))/(impurity_concentration_time(impurity_concentration_idx+1)-impurity_concentration_time(impurity_concentration_idx))

   phys%impurity_concentration = impurity_concentration
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6, *) 'impurity_concentration =  ', phys%impurity_concentration
   END IF
   DEALLOCATE(impurity_concentration_time, impurity_concentration_exp)
   NULLIFY(impurity_concentration_time, impurity_concentration_exp)
END SUBROUTINE load_impurity_concentration

SUBROUTINE compute_line_integrated_density(x_lower, x_upper, y_lower, y_upper, nli)
   REAL*8, INTENT(IN) :: x_lower, x_upper, y_lower, y_upper
   REAL*8, INTENT(OUT) :: nli
   INTEGER :: qp, Nn2D
#ifdef PARALL
   INTEGER :: ierr
#endif
   REAL*8 :: linex(1000), liney(1000), n_i(Mesh%Nelems*refElPol%Nnodes2D)
   REAL*8 :: X(Mesh%Nnodes,2), u(Mesh%Nelems*refElPol%Nnodes2D,phys%Neq)
   INTEGER :: T(Mesh%Nelems,refElPol%Nnodes2D)

   qp = SIZE(linex)
   Nn2D = refElPol%Nnodes2D
   X = mesh%X
   T = mesh%T
   nli = 0.
   linex = (/(x_lower + (x_upper - x_lower)/1000.*(i-1), i=1, 1000)/)/phys%lscale
   liney = (/(y_lower + (y_upper - y_lower)/1000.*(i-1), i=1, 1000)/)/phys%lscale
   u = TRANSPOSE(RESHAPE(sol%u,[phys%Neq,SIZE(sol%u)/phys%Neq]))
   n_i = u(:,1)

#ifndef PARALL
   CALL lineintegration(qp, linex, liney, n_i, X, T, Nn2D, nli)
#else
   CALL lineintegration(qp, linex, liney, n_i, X, T, Nn2D,Mesh%ghostElems, nli)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, nli, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   nli = nli*simpar%refval_density
   phys%n_li = nli
END SUBROUTINE compute_line_integrated_density

SUBROUTINE adjust_puff_feedback(target_density, nli)
   REAL*8, INTENT(IN) :: target_density, nli
   REAL*8 :: control_signal, anti_windup_gain

   IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,*) 'n_li = ', nli, ' [m^-2]'
       WRITE(6,*) 'n_litarget = ', target_density, '[m^-2]'
   END IF   

   ! Initialize integral error and previous error on the first timestep
   IF (time%it .EQ. 0) THEN
       phys%feedback_integral_error = 0.0
       phys%feedback_previous_error = target_density - nli
   END IF

   ! Calculate the control signal
   control_signal = phys%puff + phys%feedback_propotional_gain * (target_density - nli) + &
                            phys%feedback_integral_gain * phys%feedback_integral_error + &
                            phys%feedback_derivative_gain * (target_density - nli - phys%feedback_previous_error)/time%dt_ME
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6,*) 'Proportional impact: ', phys%feedback_propotional_gain * (target_density - nli)
      WRITE(6,*) 'Integral impact: ', phys%feedback_integral_gain * phys%feedback_integral_error
      WRITE(6,*) 'Derivative impact: ', phys%feedback_derivative_gain * (target_density - nli - phys%feedback_previous_error)/time%dt_ME
   END IF
   ! Saturate the control signal
   phys%puff = MAX(control_signal, 0.0)

   ! Back-calculate the integral error to prevent windup
   anti_windup_gain = 0.1  ! Tunable parameter
   phys%feedback_integral_error = phys%feedback_integral_error + &
                                                 anti_windup_gain * (phys%puff - control_signal)

   ! Update the integral error only if the output is not saturated
   IF (phys%puff > 0.0) THEN
       phys%feedback_integral_error = phys%feedback_integral_error + (target_density - nli) * time%dt_ME
   END IF

   ! Update the previous error
   phys%feedback_previous_error = target_density - nli
END SUBROUTINE adjust_puff_feedback

SUBROUTINE adjust_recycling_feedback(target_density, nli)
   REAL*8, INTENT(IN) :: target_density, nli
   REAL*8 :: control_signal, anti_windup_gain

   IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE(6,*) 'n_li = ', nli, ' [m^-2]'
       WRITE(6,*) 'n_litarget = ', target_density, '[m^-2]'
   END IF   

   ! Initialize integral error and previous error on the first timestep
   IF (time%it .EQ. 0) THEN
       phys%feedback_integral_error = 0.0
       phys%feedback_previous_error = target_density - nli
   END IF
   ! Print the impact of each part of the feedback on the control signal
   IF (MPIvar%glob_id .EQ. 0) THEN
      WRITE(6,*) 'Proportional impact: ', phys%feedback_propotional_gain * (target_density - nli)
      WRITE(6,*) 'Integral impact: ', phys%feedback_integral_gain * phys%feedback_integral_error
      WRITE(6,*) 'Derivative impact: ', phys%feedback_derivative_gain * (target_density - nli - phys%feedback_previous_error)/time%dt_ME
   END IF
   ! Calculate the control signal
   control_signal = phys%Re + phys%feedback_propotional_gain * (target_density - nli) + &
                            phys%feedback_integral_gain * phys%feedback_integral_error + &
                            phys%feedback_derivative_gain * (target_density - nli - phys%feedback_previous_error)/time%dt_ME
   
   ! Saturate the control signal
   phys%Re = MAX(control_signal, 0.0)

   ! Back-calculate the integral error to prevent windup
   anti_windup_gain = 0.1  ! Tunable parameter
   phys%feedback_integral_error = phys%feedback_integral_error + &
                                                 anti_windup_gain * (phys%Re - control_signal)

   ! Update the integral error only if the output is not saturated
   IF (phys%Re > 0.0) THEN
       phys%feedback_integral_error = phys%feedback_integral_error + (target_density - nli) * time%dt_ME
   END IF

   ! Update the previous error
   phys%feedback_previous_error = target_density - nli
END SUBROUTINE adjust_recycling_feedback

SUBROUTINE adjust_ITER_puff(nli)
   REAL*8, INTENT(IN) :: nli
   REAL*8 :: n_Gw, n_la, a = 2.

   n_la = nli/(8.3659 - 4.04)*simpar%refval_density
   n_Gw = phys%I_p/(pi*a**2)*10.*simpar%refval_density

   IF (time%it .EQ. 0) THEN
       phys%puff_exp(time%it+1) = MAX(10.*(phys%puff_slope*n_Gw - n_la), 0.)
       phys%puff = phys%puff_exp(time%it+1)
   ELSE
       phys%puff_exp(time%it+1) = MAX(phys%puff_exp(time%it) + 50.*(2 - SIGN(1.,phys%puff_slope*n_Gw - n_la))*(phys%puff_slope*n_Gw - n_la), 0.)
       phys%puff = phys%puff_exp(time%it+1)
   END IF

   IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, '(" * Puff = ", E10.3, 27X, " *")')  phys%puff
   END IF
END SUBROUTINE adjust_ITER_puff

END MODULE Magnetic_field
