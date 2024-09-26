!*****************************************
! project: MHDG
! file: HDG_LimitingTechniques.f90
! date: 29/03/2017
! Limiting techniques
!*****************************************
MODULE HDG_LimitingTechniques

  USE globals
  USE analytical
  USE LinearAlgebra
  USE printUtils
  USE MPI_OMP

  IMPLICIT NONE

CONTAINS

  !****************************************************************
  ! Option for limiting the min value of rho:
  ! 0 - No limiting
  ! 1 - Use a constant source in flagged elements
  ! 2 - Use a constant diffusion in flagged elements
  ! 3 - Use a constant source & diffusion in flagged elements
  !****************************************************************
  SUBROUTINE initializeLimiting()
    INTEGER           :: Ne

    IF (switch%limrho .EQ. 0) RETURN

    Ne = Mesh%Nelems
    ALLOCATE (Mesh%flag_elems_rho(Ne))
    ALLOCATE (Mesh%minrho_elems(Ne))
    Mesh%flag_elems_rho = 0
    Mesh%minrho_elems = 1e8

    SELECT CASE (switch%limrho)
    CASE (1)
      ALLOCATE (Mesh%sour_elems(Ne))
      Mesh%sour_elems = 0.
    CASE (2)
      ALLOCATE (Mesh%diff_elems(Ne))
      Mesh%diff_elems = 0.
    CASE (3)
      ALLOCATE (Mesh%sour_elems(Ne))
      ALLOCATE (Mesh%diff_elems(Ne))
      Mesh%sour_elems = 0.
      Mesh%diff_elems = 0.
    CASE DEFAULT
      WRITE (6, *) "Error: limiting case not valid"
      STOP
    END SELECT
  END SUBROUTINE initializeLimiting

  !****************************************************************
  ! Option for shock capturing:
  ! 0 - No shock capturing
  ! 1 - Constant coefficient in each element
  ! 2 - Linear interpolation of the coefficient
  ! 3 - Same as 2 but also save the artificial diffusion in each node
  !****************************************************************
  SUBROUTINE initializeShockCapturing()
    INTEGER           :: Ne, shock

    IF ((switch%shockcp .EQ. 0) .AND. (adapt%shockcp_adapt .EQ. 0)) RETURN

    IF (utils%printint > 0) THEN
      WRITE (6, *) '*************************************************'
      WRITE (6, *) '*  Initializing shock capturing                 *'
      WRITE (6, *) '*************************************************'
    END IF

    Ne = Mesh%Nelems
    IF (ALLOCATED(Mesh%flag_elems_sc)) THEN
       DEALLOCATE(Mesh%flag_elems_sc)
    ENDIF
    ALLOCATE (Mesh%flag_elems_sc(Ne))

    IF (ALLOCATED(Mesh%scdiff_nodes)) THEN
       DEALLOCATE(Mesh%scdiff_nodes)
    ENDIF

    ALLOCATE (Mesh%scdiff_nodes(Mesh%Nelems, Mesh%Nnodesperelem))
    Mesh%flag_elems_sc = 0

    shock = 0

    IF(switch%shockcp .EQ. 0) THEN
       shock = adapt%shockcp_adapt
    ENDIF

    IF(adapt%shockcp_adapt .EQ. 0) THEN
       shock = switch%shockcp
    ENDIF

    IF(ALLOCATED(refElPol%Nlin)) DEALLOCATE(refElPol%Nlin)

    SELECT CASE (shock)
    CASE (1)

    CASE (2)
      ALLOCATE (refElPol%Nlin(Mesh%Nnodesperelem, refElPol%Nvertices))
       IF (refElPol%Nvertices .EQ. 3) THEN
        refElPol%Nlin(:, 1) = -0.5*(refElPol%coord2d(:, 1) + refElPol%coord2d(:, 2))
        refElPol%Nlin(:, 2) = 0.5*(1 + refElPol%coord2d(:, 1))
        refElPol%Nlin(:, 3) = 0.5*(1 + refElPol%coord2d(:, 2))
       ELSEIF (refElPol%Nvertices .EQ. 4) THEN
        refElPol%Nlin(:, 1) = 0.25*((1 - refElPol%coord2d(:, 1))*(1 - refElPol%coord2d(:, 2)))
        refElPol%Nlin(:, 2) = 0.25*((1 + refElPol%coord2d(:, 1))*(1 - refElPol%coord2d(:, 2)))
        refElPol%Nlin(:, 3) = 0.25*((1 + refElPol%coord2d(:, 1))*(1 + refElPol%coord2d(:, 2)))
        refElPol%Nlin(:, 4) = 0.25*((1 - refElPol%coord2d(:, 1))*(1 + refElPol%coord2d(:, 2)))
      ELSE
        WRITE (6, *) "Not coded yet"
        STOP
      END IF

    CASE (3)
      ALLOCATE (refElPol%Nlin(Mesh%Nnodesperelem, refElPol%Nvertices))
      Mesh%scdiff_nodes = 0.
       IF (refElPol%Nvertices .EQ. 3) THEN
        refElPol%Nlin(:, 1) = -0.5*(refElPol%coord2d(:, 1) + refElPol%coord2d(:, 2))
        refElPol%Nlin(:, 2) = 0.5*(1 + refElPol%coord2d(:, 1))
        refElPol%Nlin(:, 3) = 0.5*(1 + refElPol%coord2d(:, 2))
       ELSEIF (refElPol%Nvertices .EQ. 4) THEN
        refElPol%Nlin(:, 1) = 0.25*((1 - refElPol%coord2d(:, 1))*(1 - refElPol%coord2d(:, 2)))
        refElPol%Nlin(:, 2) = 0.25*((1 + refElPol%coord2d(:, 1))*(1 - refElPol%coord2d(:, 2)))
        refElPol%Nlin(:, 3) = 0.25*((1 + refElPol%coord2d(:, 1))*(1 + refElPol%coord2d(:, 2)))
        refElPol%Nlin(:, 4) = 0.25*((1 - refElPol%coord2d(:, 1))*(1 + refElPol%coord2d(:, 2)))
      ELSE
        WRITE (6, *) "Not coded yet"
        STOP
      END IF

    CASE DEFAULT
      WRITE (6, *) "Error: shock capturing case not valid"
      STOP
    END SELECT
  END SUBROUTINE initializeShockCapturing

  !****************************************************************************************
  !
  !                        LIMITING OF RHO MINIMUM
  !
  !****************************************************************************************
  !**********************************************************
  ! This routine try to limit the too small values of rho,
  ! either introducing a source of density, either adding
  ! diffusion, or both
  !**********************************************************
  SUBROUTINE HDG_LimitingRho()
    IMPLICIT NONE
    INTEGER*4             :: iel, i, ct, Ndim, Neq, Nel, Np, Nfp, Nf
    INTEGER*4             :: ind(Mesh%Nnodesperelem), ind_rho(Mesh%Nnodesperelem)
    REAL*8                :: Xe(Mesh%Nnodesperelem, Mesh%Ndim)
    REAL*8                :: rhoe(Mesh%Nnodesperelem)
    REAL*8                :: rhog(refElPol%NGauss2D)
    REAL*8                :: minrhoe, minrhog, minrho_element, minrho_glob
    REAL*8                :: source, diff, maxdiff, maxsource
    REAL*8, ALLOCATABLE    :: Pel(:, :), Qel(:, :), fe(:)

    IF (switch%limrho .EQ. 0) RETURN

    IF (utils%printint > 0) THEN
      WRITE (6, *) '*************************************************'
      WRITE (6, *) '*           Limiting rho                        *'
      WRITE (6, *) '*************************************************'
    END IF

    Ndim = Mesh%ndim
    Neq = phys%Neq
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfp = refElPol%Nfacenodes

    !****************************************************
    ! Determine the number of elements that need limiting
    !****************************************************
    ct = 0
    DO iel = 1, Nel
      ind = (iel - 1)*Np*Neq + (/(i, i=1, Neq*Np, Neq)/)
      rhoe = sol%u(ind)
       rhog = MATMUL(refElPol%N2D, rhoe)
       minrhog = MINVAL(rhog)
       minrhoe = MINVAL(rhoe)
       minrho_element = MIN(minrhog, minrhoe)
       IF ((minrho_element .GT. numer%minrho) .AND. (Mesh%flag_elems_rho(iel) .EQ. 0)) CYCLE
      ct = ct + 1
      Mesh%flag_elems_rho(iel) = ct
    END DO

    ! If no element needs limiting, exit this subroutine
    IF (ct .EQ. 0) RETURN

    IF (switch%limrho .EQ. 1) THEN
      ALLOCATE (fe(Np))
      IF (ALLOCATED(elMat%S_lrho)) THEN
        DEALLOCATE (elMat%S_lrho)
      END IF
      ALLOCATE (elMat%S_lrho(Neq*Np, ct))
      elMat%S_lrho = 0.
    ELSE IF (switch%limrho .EQ. 2) THEN
      IF (ALLOCATED(elMat%P_lrho)) THEN
        DEALLOCATE (elMat%P_lrho)
      END IF
      ALLOCATE (elMat%P_lrho(Neq*Np, Ndim*Neq*Np, ct))
      ALLOCATE (Pel(Neq*Np, Ndim*Neq*Np))
      ALLOCATE (Qel(Neq*Np, Ndim*Neq*Np))
      elMat%P_lrho = 0.
    ELSE
      ALLOCATE (fe(Np))
      IF (ALLOCATED(elMat%S_lrho)) THEN
        DEALLOCATE (elMat%S_lrho)
      END IF
      ALLOCATE (elMat%S_lrho(Neq*Np, ct))
      IF (ALLOCATED(elMat%P_lrho)) THEN
        DEALLOCATE (elMat%P_lrho)
      END IF
      ALLOCATE (elMat%P_lrho(Neq*Np, Ndim*Neq*Np, ct))
      ALLOCATE (Pel(Neq*Np, Ndim*Neq*Np))
      ALLOCATE (Qel(Neq*Np, Ndim*Neq*Np))
      elMat%S_lrho = 0.
      elMat%P_lrho = 0.
    END IF

    ind_rho = (/(i, i=1, Neq*Np, Neq)/)
    minrho_glob = 1.e9
    maxsource = -1.e9
    maxdiff = -1.e9
    ct = 1
    !*****************
    ! Loop in elements
    !*****************
    DO iel = 1, Nel
       IF (Mesh%flag_elems_rho(iel) .EQ. 0) CYCLE
      ind = (iel - 1)*Np*Neq + (/(i, i=1, Neq*Np, Neq)/)
      rhoe = sol%u(ind)
       rhog = MATMUL(refElPol%N2D, rhoe)
       minrhog = MINVAL(rhog)
       minrhoe = MINVAL(rhoe)
       minrho_element = MIN(minrhog, minrhoe)
       minrho_glob = MIN(minrho_glob, minrho_element)

      ! Coordinates of the nodes of the element
      Xe = Mesh%X(Mesh%T(iel, :), :)

      !*****************
      ! CASE ADD SOURCE
      !*****************
       IF (switch%limrho .EQ. 1) THEN

        ! Determine the source to add in this element
          IF (minrho_element .GT. Mesh%minrho_elems(iel)) THEN
          ! Here I use the surce that I used in the previous time step
          source = Mesh%sour_elems(iel)
        ELSE
             IF (minrho_element .LT. 0.) THEN
            source = numer%so_coe*numer%minrho
          ELSE
                source = numer%so_coe*numer%minrho*LOG10(numer%minrho/minrho_element)
          END IF
          Mesh%minrho_elems(iel) = minrho_element
        END IF
        Mesh%sour_elems(iel) = source
          maxsource = MAX(source, maxsource)
        fe = 0.d0

        ! Compute the matrices for the element
        CALL elemental_matrices_source(Xe, fe, source)

        elMat%S_lrho(ind_rho, ct) = fe
        !********************
        ! CASE ADD DIFFUSION
        !********************
       ELSE IF (switch%limrho .EQ. 2) THEN

        ! Determine the diffusion to add in this element
          IF (minrho_element .GT. Mesh%minrho_elems(iel) .AND. Mesh%flag_elems_rho(iel) .NE. 0) THEN
          ! Here I use the diffusion that I used in the previous time step
          diff = Mesh%diff_elems(iel)
        ELSE
             IF (minrho_element .LT. 0.) THEN
            diff = numer%df_coe*phys%diff_n
          ELSE
                diff = numer%df_coe*phys%diff_n*LOG10(numer%minrho/minrho_element)
          END IF
          Mesh%minrho_elems(iel) = minrho_element
        END IF
        Mesh%diff_elems(iel) = diff
          maxdiff = MAX(diff, maxdiff)

        Pel = 0.d0; Qel = 0.d0

        ! Compute the matrices for the element
        CALL elemental_matrices_diff(Xe, Pel, Qel, diff)

        ! Store the matrices for all the elements
        elMat%P_lrho(:, :, ct) = Pel - Qel

        !******************************
        ! CASE ADD SOURCE AND DIFFUSION
        !******************************
      ELSE
        ! Determine the diffusion to add in this element
          IF (minrho_element .GT. Mesh%minrho_elems(iel) .AND. Mesh%flag_elems_rho(iel) .NE. 0) THEN
          ! Here I use the diffusion and source that I used in the previous time step
          diff = Mesh%diff_elems(iel)
          source = Mesh%sour_elems(iel)
        ELSE
             IF (minrho_element .LT. 0.) THEN
            diff = numer%df_coe*phys%diff_n
            source = numer%so_coe*numer%minrho
          ELSE
                diff = numer%df_coe*phys%diff_n*LOG10(numer%minrho/minrho_element)
                source = numer%so_coe*numer%minrho*LOG10(numer%minrho/minrho_element)
          END IF
          Mesh%minrho_elems(iel) = minrho_element
        END IF
        Mesh%diff_elems(iel) = diff
          maxdiff = MAX(diff, maxdiff)
        Mesh%sour_elems(iel) = source
          maxsource = MAX(source, maxsource)

        fe = 0.d0; Pel = 0.d0; Qel = 0.d0

        ! Compute the matrices for the element
        CALL elemental_matrices_diff(Xe, Pel, Qel, diff)
        CALL elemental_matrices_source(Xe, fe, source)

        ! Store the matrices for all the elements
        elMat%S_lrho(ind_rho, ct) = fe
        elMat%P_lrho(:, :, ct) = Pel - Qel
      END IF
      ct = ct + 1
    END DO

    IF (switch%limrho .EQ. 1) THEN
      DEALLOCATE (fe)
    ELSE IF (switch%limrho .EQ. 2) THEN
      DEALLOCATE (Pel, Qel)
    ELSE
      DEALLOCATE (fe, Pel, Qel)
    END IF

    WRITE (6, *) "*********************************************************"
    WRITE (6, '(A,I0,A)') "Limiting min value of rho in ", ct - 1, " elements:"
    WRITE (6, '(A,E11.4)') "Min rho detected: ", minrho_glob
    IF (switch%limrho .EQ. 1) THEN
      WRITE (6, '(A,E11.4)') "Max source applied: ", maxsource
    ELSE IF (switch%limrho .EQ. 2) THEN
      WRITE (6, '(A,E11.4)') "Max diffusion applied: ", maxdiff
    ELSE
      WRITE (6, '(A,E11.4)') "Max source applied: ", maxsource
      WRITE (6, '(A,E11.4)') "Max diffusion applied: ", maxdiff
    END IF
    WRITE (6, *) "*********************************************************"

  CONTAINS

    !*****************************************
    ! Elemental computation: add source
    !*****************************************
    SUBROUTINE elemental_matrices_source(Xe, fe, source)
      REAL*8, INTENT(IN)     :: Xe(1:Mesh%Nnodesperelem, 1:Mesh%Ndim)
      REAL*8, INTENT(INOUT)  :: fe(:)
      REAL*8, INTENT(IN)     :: source
      INTEGER*4             :: g, NGauss
      REAL*8                :: detJ, dvolu
      REAL*8                :: xy(1:refElPol%Ngauss2d, 1:ndim)
      REAL*8                :: Jacob(1:ndim, 1:ndim)

      !***********************************
      !    Volume computation
      !***********************************

      ! Gauss points position
      xy = MATMUL(refElPol%N2D, Xe)

      ! Loop in 2D Gauss points
      Ngauss = refElPol%NGauss2D
      DO g = 1, NGauss

        ! Jacobian
        Jacob = 0.d0
         Jacob(1, 1) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 1))
         Jacob(1, 2) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 2))
         Jacob(2, 1) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 1))
         Jacob(2, 2) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 2))
        detJ = Jacob(1, 1)*Jacob(2, 2) - Jacob(1, 2)*Jacob(2, 1)

        ! Integration weight
        dvolu = refElPol%gauss_weights2D(g)*detJ
        IF (switch%axisym) THEN
          dvolu = dvolu*xy(g, 1)
        END IF
        ! Contribution of the current integration point to the elemental matrix
        fe = fe + source*refElPol%N2D(g, :)*dvolu
      END DO

    END SUBROUTINE elemental_matrices_source

    !*****************************************
    ! Elemental computation: add diffusion
    !*****************************************
    SUBROUTINE elemental_matrices_diff(Xe, Pel, Qel, diff)
      REAL*8, INTENT(IN)     :: Xe(1:Mesh%Nnodesperelem, 1:Mesh%Ndim)
      REAL*8, INTENT(INOUT)  :: Pel(:, :), Qel(:, :)
      REAL*8, INTENT(IN)     :: diff
      INTEGER*4             :: g, NGauss, ifa
      REAL*8                :: detJ, dvolu, dline, xyDerNorm_g
      REAL*8                :: xy(1:refElPol%Ngauss2d, 1:ndim)
      REAL*8                :: xyf(1:refElPol%Ngauss1d, 1:ndim)
      REAL*8                :: xyDer(1:refElPol%Ngauss1d, 1:ndim)
      REAL*8                :: Jacob(1:ndim, 1:ndim)
      INTEGER*4             :: ind_fe(neq*Nfp)
      INTEGER*4             :: ind_fG(neq*ndim*Nfp)
      REAL*8, DIMENSION(Np)  :: Nxg, Nyg, Nx_ax
      REAL*8                :: invJ(1:ndim, 1:ndim), t_g(1:ndim), n_g(1:ndim)
      REAL*8, ALLOCATABLE, DIMENSION(:, :) :: Px, Py, Q_loc, Qx, Qy

      !***********************************
      !    Volume computation
      !***********************************
      ALLOCATE (Px(1:Np, 1:Np))
      ALLOCATE (Py(1:Np, 1:Np))
      Px = 0.; Py = 0.

      ! Gauss points position
      xy = MATMUL(refElPol%N2D, Xe)

      ! Loop in 2D Gauss points
      Ngauss = refElPol%NGauss2D
      DO g = 1, NGauss

        ! Jacobian
        Jacob = 0.d0
         Jacob(1, 1) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 1))
         Jacob(1, 2) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 2))
         Jacob(2, 1) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 1))
         Jacob(2, 2) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 2))
        detJ = Jacob(1, 1)*Jacob(2, 2) - Jacob(1, 2)*Jacob(2, 1)
        IF (detJ < 1e-12) THEN
            error STOP "Negative jacobian"
         END IF

        ! x and y derivatives of the shape functions
         CALL invert_matrix(Jacob, invJ)
        Nxg = invJ(1, 1)*refElPol%Nxi2D(g, :) + invJ(1, 2)*refElPol%Neta2D(g, :)
        Nyg = invJ(2, 1)*refElPol%Nxi2D(g, :) + invJ(2, 2)*refElPol%Neta2D(g, :)
        IF (switch%axisym) THEN
          Nx_ax = Nxg + 1./xy(g, 1)*refElPol%N2D(g, :)
        ELSE
          Nx_ax = Nxg
        END IF
        ! Integration weight
        dvolu = refElPol%gauss_weights2D(g)*detJ
        IF (switch%axisym) THEN
          dvolu = dvolu*xy(g, 1)
        END IF
        ! Contribution of the current integration point to the elemental matrix
        Px = Px + diff*TensorProduct(Nxg, refElPol%N2D(g, :))*dvolu
        Py = Py + diff*TensorProduct(Nyg, refElPol%N2D(g, :))*dvolu
      END DO

      CALL expand_matrix_Bt(Px, Py, Pel)

      !***********************************
      ! Faces computations
      !***********************************
      NGauss = refElPol%Ngauss1D
      ALLOCATE (Q_loc(Neq*Ndim*Nfp, Ndim*Nfp))
      ALLOCATE (Qx(Nfp, Nfp))
      ALLOCATE (Qy(Nfp, Nfp))
      Q_loc = 0.

      ! Loop in the faces of the element
      DO ifa = 1, Nf

         ind_fG = RESHAPE(tensorSumInt((/(i, i=1, neq*ndim)/), neq*ndim*(refElPol%face_nodes(ifa, :) - 1)), (/neq*ndim*Nfp/))
         ind_fe = RESHAPE(tensorSumInt((/(i, i=1, neq)/), neq*(refElPol%face_nodes(ifa, :) - 1)), (/neq*Nfp/))
         xyf = MATMUL(refElPol%N1D, Xe(refElPol%face_nodes(ifa, :), :))
         xyDer = MATMUL(refElPol%Nxi1D, Xe(refElPol%face_nodes(ifa, :), :))

        Qx = 0.; Qy = 0.
        ! Loop in 1D Gauss points
        DO g = 1, NGauss

          ! Calculate the integration weight
            xyDerNorm_g = NORM2(xyDer(g, :))
          dline = refElPol%gauss_weights1D(g)*xyDerNorm_g

          IF (switch%axisym) THEN
            dline = dline*xyf(g, 1)
          END IF
          ! Unit normal to the boundary
          t_g = xyDer(g, :)/xyDerNorm_g
          n_g = [t_g(2), -t_g(1)]
          ! Contribution of the current integration point to the elemental matrix
          Qx = Qx + diff*tensorProduct(refElPol%N1D(g, :), refElPol%N1D(g, :))*n_g(1)*dline
          Qy = Qy + diff*tensorProduct(refElPol%N1D(g, :), refElPol%N1D(g, :))*n_g(2)*dline
        END DO

        ! Elemental assembly
        CALL expand_matrix_B(Qx, Qy, Q_loc)
         Qel(ind_fe, ind_fG) = Qel(ind_fe, ind_fG) + TRANSPOSE(Q_loc)

      END DO
      DEALLOCATE (Px, Py, Q_loc, Qx, Qy)

    END SUBROUTINE elemental_matrices_diff

  END SUBROUTINE HDG_LimitingRho

  !****************************************************************************************
  !
  !                        ADD DIFFUSION IN CORNERS
  !
  !****************************************************************************************
  SUBROUTINE HDG_AddDiffusionCorner()

    INTEGER*4             :: iel, i, j, Ndim, Neq, Nel, Np, Nf, Nfp
    INTEGER*4             :: perm(phys%Neq*Mesh%Nnodesperface), ind_loc(refElPol%Nfaces, phys%Neq*Mesh%Nnodesperface)
    REAL*8                :: Xe(1:Mesh%Nnodesperelem, 1:Mesh%Ndim)
    REAL*8                :: diff(1:Mesh%Nnodesperelem)
    REAL*8, ALLOCATABLE    :: Pel(:, :), Qel(:, :), Lfel(:, :)
    REAL*8                :: tol = 1e-3

    IF (switch%difcor .EQ. 0) RETURN

    IF (MPIvar%glob_id .EQ. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) '*************************************************'
        WRITE (6, *) '*           Adding diffusion in corners         *'
        WRITE (6, *) '*************************************************'
      END IF
    END IF

    Ndim = Mesh%ndim
    Neq = phys%Neq
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfp = refElPol%Nfacenodes

    ALLOCATE (Pel(Neq*Np, Ndim*Neq*Np))
    ALLOCATE (Qel(Neq*Np, Ndim*Neq*Np))
    ALLOCATE (Lfel(Neq*Nf*Nfp, Neq*Ndim*Np))

    perm = 0
    CALL set_permutations(Neq*Nfp, Neq, perm)

    ind_loc = 0
    DO i = 1, Nf
      DO j = 1, Neq*Nfp
        ind_loc(i, j) = Neq*Nfp*(i - 1) + j
      END DO
    END DO

    !*****************
    ! Loop in elements
    !*****************
    DO iel = 1, Nel
      Pel = 0.
      Qel = 0.
      Lfel = 0.

      ! Coordinates of the nodes of the element
      Xe = Mesh%X(Mesh%T(iel, :), :)

      ! Compute the local diffusion in this point
      diff = computeDiffusion(Xe)

       IF (MAXVAL(diff) .LT. tol) CYCLE

      ! Compute the matrices for the element
      CALL elemental_matrices(Xe, Pel, Qel, Lfel, diff)

      ! Flip the faces for the elements for which the face is already been counted
      DO j = 1, Nf
        IF (Mesh%flipface(iel, j)) THEN
          Lfel(ind_loc(j, :), :) = Lfel(ind_loc(j, perm), :)
          END IF
      END DO

      ! Store the matrices for all the elements
      elMat%P(:, :, iel) = elMat%P(:, :, iel) + Pel - Qel
      elMat%Lf(:, :, iel) = elMat%Lf(:, :, iel) + Lfel
    END DO

    DEALLOCATE (Pel, Qel, Lfel)

    IF (MPIvar%glob_id .EQ. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Done!"
      END IF
    END IF

  CONTAINS

    !*********************************************
    ! Elemental matrices: add diffusion in a point
    !*********************************************
    SUBROUTINE elemental_matrices(Xe, Pel, Qel, Lfel, diff)
      REAL*8, INTENT(IN)     :: Xe(Np, ndim)
      REAL*8, INTENT(INOUT)  :: Pel(:, :), Qel(:, :), Lfel(:, :)
      REAL*8, INTENT(IN)     :: diff(1:Mesh%Nnodesperelem)
      REAL*8                :: diffg(1:refElPol%Ngauss2d), diffgf(refElPol%Ngauss1d)
      INTEGER*4             :: g, NGauss, ifa
      REAL*8                :: detJ, dvolu, dline, xyDerNorm_g, isext
      REAL*8                :: xy(1:refElPol%Ngauss2d, 1:ndim)
      REAL*8                :: xyf(1:refElPol%Ngauss1d, 1:ndim)
      REAL*8                :: xyDer(1:refElPol%Ngauss1d, 1:ndim)
      REAL*8                :: Jacob(1:ndim, 1:ndim)
      INTEGER*4             :: ind_fe(neq*Nfp), ind_ff(neq*Nfp)
      INTEGER*4             :: ind_fG(neq*ndim*Nfp)
      REAL*8, DIMENSION(Np)  :: Nxg, Nyg, Nx_ax
      REAL*8                :: invJ(1:ndim, 1:ndim), t_g(1:ndim), n_g(1:ndim)
      REAL*8, ALLOCATABLE, DIMENSION(:, :)    :: Px, Py, Q_loc, Qx, Qy, Lf_t

      !***********************************
      !    Volume computation
      !***********************************
      ALLOCATE (Px(1:Np, 1:Np))
      ALLOCATE (Py(1:Np, 1:Np))
      ALLOCATE (Lf_t(Neq*Ndim*Np, Nf*Neq*Nfp))
      Px = 0.; Py = 0.; Lf_t = 0.

      ! Gauss points position
      xy = MATMUL(refElPol%N2D, Xe)

      ! Diffusion at Gauss points
      diffg = MATMUL(refElPol%N2D, diff)

      ! Loop in 2D Gauss points
      Ngauss = refElPol%NGauss2D
      DO g = 1, NGauss

        ! Jacobian
        Jacob = 0.d0
         Jacob(1, 1) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 1))
         Jacob(1, 2) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 2))
         Jacob(2, 1) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 1))
         Jacob(2, 2) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 2))
        detJ = Jacob(1, 1)*Jacob(2, 2) - Jacob(1, 2)*Jacob(2, 1)
        IF (detJ < tol) THEN
            error STOP "Negative jacobian"
         END IF

        ! x and y derivatives of the shape functions
         CALL invert_matrix(Jacob, invJ)
        Nxg = invJ(1, 1)*refElPol%Nxi2D(g, :) + invJ(1, 2)*refElPol%Neta2D(g, :)
        Nyg = invJ(2, 1)*refElPol%Nxi2D(g, :) + invJ(2, 2)*refElPol%Neta2D(g, :)
        IF (switch%axisym) THEN
          Nx_ax = Nxg + 1./xy(g, 1)*refElPol%N2D(g, :)
        ELSE
          Nx_ax = Nxg
        END IF
        ! Integration weight
        dvolu = refElPol%gauss_weights2D(g)*detJ
        IF (switch%axisym) THEN
          dvolu = dvolu*xy(g, 1)
        END IF
        ! Contribution of the current integration point to the elemental matrix
        Px = Px + diffg(g)*TensorProduct(Nxg, refElPol%N2D(g, :))*dvolu
        Py = Py + diffg(g)*TensorProduct(Nyg, refElPol%N2D(g, :))*dvolu
      END DO

      CALL expand_matrix_Bt(Px, Py, Pel)

      !***********************************
      ! Faces computations
      !***********************************
      NGauss = refElPol%Ngauss1D
      ALLOCATE (Q_loc(Neq*Ndim*Nfp, Neq*Nfp))
      ALLOCATE (Qx(Nfp, Nfp))
      ALLOCATE (Qy(Nfp, Nfp))
      Q_loc = 0.

      ! Loop in the faces of the element
      DO ifa = 1, Nf
        isext = 0.
         IF ((Mesh%F(iel, ifa) .GT. Mesh%Nintfaces)) THEN
            IF (Mesh%boundaryFlag(Mesh%F(iel, ifa) - Mesh%Nintfaces) .NE. 0) THEN
            isext = 1.
          END IF
        ENDIF

         ind_fG = RESHAPE(tensorSumInt((/(i, i=1, neq*ndim)/), neq*ndim*(refElPol%face_nodes(ifa, :) - 1)), (/neq*ndim*Nfp/))
         ind_fe = RESHAPE(tensorSumInt((/(i, i=1, neq)/), neq*(refElPol%face_nodes(ifa, :) - 1)), (/neq*Nfp/))
        ind_ff = (ifa - 1)*neq*Nfp + (/(i, i=1, neq*Nfp)/)
         xyf = MATMUL(refElPol%N1D, Xe(refElPol%face_nodes(ifa, :), :))
         xyDer = MATMUL(refElPol%Nxi1D, Xe(refElPol%face_nodes(ifa, :), :))

        ! Diffusion at Gauss points
         diffgf = MATMUL(refElPol%N1D, diff(refElPol%face_nodes(ifa, :)))

        Qx = 0.; Qy = 0.
        ! Loop in 1D Gauss points
        DO g = 1, NGauss

          ! Calculate the integration weight
            xyDerNorm_g = NORM2(xyDer(g, :))
          dline = refElPol%gauss_weights1D(g)*xyDerNorm_g

          IF (switch%axisym) THEN
            dline = dline*xyf(g, 1)
          END IF
          ! Unit normal to the boundary
          t_g = xyDer(g, :)/xyDerNorm_g
          n_g = [t_g(2), -t_g(1)]
          ! Contribution of the current integration point to the elemental matrix
          Qx = Qx + diffgf(g)*tensorProduct(refElPol%N1D(g, :), refElPol%N1D(g, :))*n_g(1)*dline
          Qy = Qy + diffgf(g)*tensorProduct(refElPol%N1D(g, :), refElPol%N1D(g, :))*n_g(2)*dline
        END DO

        ! Elemental assembly
        CALL expand_matrix_B(Qx, Qy, Q_loc)
         Qel(ind_fe, ind_fG) = Qel(ind_fe, ind_fG) + TRANSPOSE(Q_loc)
        Lf_t(ind_fG, ind_ff) = Lf_t(ind_fG, ind_ff) + Q_loc*(1 - isext)

      END DO

      Lfel = TRANSPOSE(Lf_t)

      DEALLOCATE (Px, Py, Q_loc, Qx, Qy, Lf_t)

    END SUBROUTINE elemental_matrices

    !*******************************************
    ! Compute local diffusion in points
    !*******************************************
    FUNCTION computeDiffusion(Xe) RESULT(res)
      REAL*8, INTENT(IN)  :: Xe(Mesh%Nnodesperelem, Mesh%Ndim)
      REAL*8             :: res(Mesh%Nnodesperelem), A(3, 3)
      REAL*8, ALLOCATABLE :: xycorn(:, :), b(:, :), barcor(:, :), barcormax(:), barcormin(:)
      REAL*8, ALLOCATABLE :: rad(:)
      REAL*8             :: h
      REAL*8, PARAMETER   :: tol = 1.e-6
      INTEGER            :: i
      INTEGER            :: chose

      res = 0.

      chose = 2

      SELECT CASE (switch%difcor)
      CASE (1)
        ! Circular case with infinitely small limiter
        ALLOCATE (xycorn(2, 2))
        xycorn(1, 1) = geom%R0
        xycorn(2, 1) = -0.75
        xycorn(1, 2) = 1e8
        xycorn(2, 2) = 1e8
      CASE (2)
        ! Circular case with infinitely small limiter
        ALLOCATE (xycorn(2, 2))
        xycorn(1, 1) = geom%R0
        xycorn(2, 1) = -0.287
        xycorn(1, 2) = 1e8
        xycorn(2, 2) = 1e8
      CASE (3)
        ! West
        ALLOCATE (xycorn(2, 2))
        xycorn(1, 1) = 2.7977
        xycorn(2, 1) = -0.5128
        xycorn(1, 2) = 1e8
        xycorn(2, 2) = 1e8
      CASE DEFAULT
        WRITE (6, *) "Case not valid"
        STOP
      END SELECT

      SELECT CASE (chose)

      CASE (1)
        !!**********************************************************
        !! Constant diffusion in the elements containing the corner
        !!**********************************************************
        ! Find if the corner belongs to the element
         ALLOCATE (b(3, SIZE(xycorn, 2)))
         ALLOCATE (barcor(3, SIZE(xycorn, 2)))
         ALLOCATE (barcormax(SIZE(xycorn, 2)))
         ALLOCATE (barcormin(SIZE(xycorn, 2)))
        A(1, 1:3) = Xe(1:3, 1)*phys%lscale
        A(2, 1:3) = Xe(1:3, 2)*phys%lscale
        A(3, 1:3) = 1.
        b(1, :) = xycorn(1, :)
        b(2, :) = xycorn(2, :)
        b(3, :) = 1.

        CALL solve_linear_system(A, b, barcor)
         barcormax = MAXVAL(barcor, 1)
         barcormin = MINVAL(barcor, 1)

         DO i = 1, SIZE(barcor, 2)
            IF ((barcormax(i) .GT. (1.+tol)) .OR. (barcormin(i) .LT. (0.-tol))) CYCLE
          res = numer%dc_coe*phys%diff_n
        END DO
        DEALLOCATE (b, barcor, barcormax, barcormin)

      CASE (2)
        !!**********************************************************
        !! Gaussian around the corner
        !!**********************************************************
        ALLOCATE (rad(Mesh%Nnodesperelem))
        h = 10e-3
         DO i = 1, SIZE(xycorn, 2)
            rad = SQRT((Xe(:, 1)*phys%lscale - xycorn(1, i))**2 + (Xe(:, 2)*phys%lscale - xycorn(2, i))**2)
            res = res + numer%dc_coe*phys%diff_n*EXP(-(2*rad/h)**2)
        END DO
        DEALLOCATE (rad)

      END SELECT

      DEALLOCATE (xycorn)

    END FUNCTION computeDiffusion

    !*****************************************
    ! Set permutations for flipping faces
    !****************************************
    SUBROUTINE set_permutations(n, m, perm)
      INTEGER, INTENT(IN)  :: n, m
      INTEGER, INTENT(OUT) :: perm(:)
      INTEGER              :: i
      INTEGER              :: temp(m, n/m), templr(m, n/m)

      IF (MOD(n, m) .NE. 0) THEN
        WRITE (6, *) 'Error! n must be a multiple of m'
        STOP
      END IF

      templr = 0
      temp = RESHAPE((/(i, i=1, n)/), (/m, n/m/))
      DO i = 1, n/m
        templr(:, i) = temp(:, n/m - i + 1)
      END DO
      perm = RESHAPE(templr, (/n/))
    END SUBROUTINE set_permutations

  END SUBROUTINE HDG_AddDiffusionCorner

  !****************************************************************************************
  !
  !                        SHOCK CAPTURING
  !
  !****************************************************************************************
  SUBROUTINE HDG_ShockCapturing()
    USE reference_element, ONLY: vandermonde_2d, vandermonde_qua
    INTEGER*4             :: iel, i, Ndim, Neq, Nel, Np, ct, Nf, Nfp, ctgl, inod
    REAL*8                :: diff(Mesh%Nelems), maxdiff
    REAL*8                :: tol = 1e-12
    INTEGER               :: els(SIZE(Mesh%N, 2))
    REAL*8                :: auxones(Mesh%Nnodesperelem), diffnod(Mesh%Nnodesperelem), diffver(refElPol%Nvertices)
    LOGICAL               :: check
    REAL*8                :: Vand(refElPol%Nnodes2D, refElPol%Nnodes2D), invVand(refElPol%Nnodes2D, refElPol%Nnodes2D)
#ifdef PARALL
    INTEGER               :: ierr
#endif
    IF (switch%shockcp .EQ. 0) RETURN

    IF (MPIvar%glob_id .EQ. 0) THEN
      IF (utils%printint > 1) THEN
        WRITE (6, *) '*************************************************'
        WRITE (6, *) '*           Shock capturing                     *'
        WRITE (6, *) '*************************************************'
      END IF
    END IF

    Ndim = Mesh%ndim
    Neq = phys%Neq
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nf = refElPol%Nfaces
    Nfp = refElPol%Nfacenodes

    maxdiff = 0.

    ! Reset flags for shock capturing
    Mesh%flag_elems_sc = 0

    !******* Find shock capturing coefficient in each element
    ! Vandermonde matrix
    IF (refElPol%elemType == 0) THEN
      ! Triangles
      CALL vandermonde_2d(Vand, refElPol)
    ELSEIF (refElPol%elemType == 1) THEN
      ! Quadrilaterals
      CALL vandermonde_qua(Vand, refElPol)
    ELSE
      WRITE (6, *) "Vandermonde matrix for this element type not coded yet"
      STOP
    END IF
    ! Invert Vandermonde matrix
    CALL invert_matrix(Vand, invVand)
    CALL findCoeffShockCaptur(1.e-1, diff, invVand)
    ! *****
    ! First loop in elements to allocate
    ct = 0
    IF (switch%shockcp .EQ. 1) THEN
      DO iel = 1, Nel
          IF (diff(iel) .LT. tol) CYCLE
        ct = ct + 1
      END DO
    ELSEIF (switch%shockcp .GE. 2) THEN
      DO iel = 1, Nel
          check = .FALSE.
        DO inod = 1, refElPol%Nvertices
          els = Mesh%N(Mesh%Tlin(iel, inod), :)
             DO i = 1, SIZE(els)
                IF (els(i) .EQ. 0) CYCLE
                IF (diff(els(i)) .GE. tol) THEN
                   check = .TRUE.
              EXIT
            END IF
          END DO
        END DO
        IF (check) ct = ct + 1
      END DO
    ELSE
      WRITE (6, *) "Something is wrong"
      STOP
    END IF

    ! Initialize shock capturing diffusion storing
    Mesh%scdiff_nodes = 0.

    IF (ct .GT. 0) THEN

      ! Allocate
      Mesh%flag_elems_sc = 0

      auxones = 1.

      !*****************
      ! Loop in elements
      !*****************
      ct = 1
      DO iel = 1, Nel
          IF (switch%shockcp .EQ. 1) THEN
             IF (diff(iel) .LT. tol) THEN
            CYCLE
          ELSE
            diffnod = diff(iel)*auxones
          END IF
          ELSE IF (switch%shockcp .GE. 2) THEN
          DO inod = 1, refElPol%Nvertices
            els = Mesh%N(Mesh%Tlin(iel, inod), :)
            diffver(inod) = 0.
                DO i = 1, SIZE(els)
                   IF (els(i) .EQ. 0) CYCLE
                   diffver(inod) = MAX(diffver(inod), diff(els(i)))
            END DO
          END DO
             IF (MAXVAL(diffver) .LT. tol) THEN
            CYCLE
          ELSE
                diffnod = MATMUL(refElPol%Nlin, diffver)
          END IF
        END IF

          maxdiff = MAX(maxdiff, MAXVAL(diffnod))

        ! Store the shock capturing diffusion
        Mesh%scdiff_nodes(iel, :) = diffnod
        Mesh%flag_elems_sc(iel) = ct
        ct = ct + 1
      END DO

      ctgl = ct - 1
    ELSE
      ctgl = 0
    END IF

#ifdef PARALL
    CALL mpi_allreduce(MPI_IN_PLACE, ctgl, 1, mpi_integer, mpi_sum, MPI_COMM_WORLD, ierr)
    CALL mpi_allreduce(MPI_IN_PLACE, maxdiff, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
    IF (MPIvar%glob_id .EQ. 0 .AND. ctgl .GT. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) "Added shock capturing diffusion in ", ctgl, " elements. Max diffusion: ", maxdiff
        WRITE (6, *) "Done!"
      END IF
    END IF

  CONTAINS
    !*****************************************
    ! Elemental matrices computation: 2D case
    !*****************************************
    !                                        SUBROUTINE elemental_matrices(Xe,Pel,Qel,Lfel,diff)
    !       real*8,intent(IN)     :: Xe(Mesh%Nnodesperelem,Mesh%Ndim)
    !       real*8,intent(INOUT)  :: Pel(:,:),Qel(:,:),Lfel(:,:)
    !       real*8,intent(IN)     :: diff
    !                                                        integer*4             :: g,NGauss,ifa
    !                                                        real*8                :: detJ,dvolu,dline,xyDerNorm_g,isext
    !                                                        real*8                :: xy(refElPol%Ngauss2d,Mesh%Ndim)
    !                                                        real*8                :: xyf(refElPol%Ngauss1d,Mesh%Ndim)
    !                                                        real*8                :: xyDer(refElPol%Ngauss1d,Mesh%Ndim)
    !                                                        real*8                :: Jacob(Mesh%Ndim,Mesh%Ndim)
    !                                                        integer*4             :: ind_fe(phys%neq*refElPol%Nfacenodes), ind_ff(neq*Nfp)
    !       integer*4             :: ind_fG(phys%neq*Mesh%Ndim*refElPol%Nfacenodes)
    !       real*8                :: invJ(Mesh%Ndim,Mesh%Ndim),t_g(Mesh%Ndim),n_g(Mesh%Ndim)
    !       real*8,dimension(Np)  :: Nxg,Nyg,Nx_ax
    !       real*8,allocatable,dimension(:,:)    :: Px,Py,Q_loc,Qx,Qy,Lf_t

    !       !***********************************
    !       !    Volume computation
    !       !***********************************
    !       ALLOCATE(Px(Np,Np))
    !       ALLOCATE(Py(Np,Np))
    !       ALLOCATE(Lf_t(Neq*Ndim*Np,Nf*Ndim*Nfp))
    !       Px = 0.; Py = 0.; Lf_t = 0.

    !                                                        ! Gauss points position
    !                                                        xy = matmul(refElPol%N2D,Xe)
    !
    !                                                        ! Loop in 2D Gauss points
    !       Ngauss = refElPol%NGauss2D
    !                                                        DO g = 1,NGauss
    !
    !                                                                                ! Jacobian
    !                                                                                Jacob = 0.d0
    !                                                                                Jacob(1,1) = dot_product(refElPol%Nxi2D(g,:),Xe(:,1))
    !                                                                                Jacob(1,2) = dot_product(refElPol%Nxi2D(g,:),Xe(:,2))
    !                                                                                Jacob(2,1) = dot_product(refElPol%Neta2D(g,:),Xe(:,1))
    !                                                                                Jacob(2,2) = dot_product(refElPol%Neta2D(g,:),Xe(:,2))
    !                                                                                detJ = Jacob(1,1)*Jacob(2,2) - Jacob(1,2)*Jacob(2,1)
    !                                                                                IF (detJ < tol) THEN
    !                                                                                   error stop "Negative jacobian"
    !                                                                                END if

    !                                                                                ! x and y derivatives of the shape functions
    !                                                                                call invert_matrix(Jacob,invJ)
    !                                                           Nxg = invJ(1,1)*refElPol%Nxi2D(g,:) + invJ(1,2)*refElPol%Neta2D(g,:)
    !                                                           Nyg = invJ(2,1)*refElPol%Nxi2D(g,:) + invJ(2,2)*refElPol%Neta2D(g,:)
    !                                                           IF (switch%axisym) THEN
    !                                                              Nx_ax = Nxg+1./xy(g,1)*refElPol%N2D(g,:)
    !                                                           ELSE
    !                                                              Nx_ax = Nxg
    !                                                           END IF
    !                                                                                ! Integration weight
    !                                                                                dvolu = refElPol%gauss_weights2D(g)*detJ
    !          IF (switch%axisym) THEN
    !             dvolu = dvolu*xy(g,1)
    !          END IF
    !                                                                                ! Contribution of the current integration point to the elemental matrix
    !          Px = Px + diff*TensorProduct(Nxg,refElPol%N2D(g,:))*dvolu
    !          Py = Py + diff*TensorProduct(Nyg,refElPol%N2D(g,:))*dvolu
    !                                                        END DO

    !       CALL expand_matrix_Bt(Px,Py,Pel)

    !       !***********************************
    !                                          ! Faces computations
    !       !***********************************
    !                                                        NGauss = refElPol%Ngauss1D
    !       ALLOCATE(Q_loc(Neq*Ndim*Nfp,Ndim*Nfp))
    !       ALLOCATE(Qx(Nfp,Nfp))
    !       ALLOCATE(Qy(Nfp,Nfp))
    !       Q_loc = 0.
    !
    !       ! Loop in the faces of the element
    !                                                        DO ifa = 1,Nf
    !          isext = 0.
    !          IF ( (Mesh%F(iel,ifa).gt.Mesh%Nintfaces)  ) THEN
    !             IF (Mesh%boundaryFlag(Mesh%F(iel,ifa)-Mesh%Nintfaces).ne.0) THEN
    !                isext = 1.
    !             END IF
    !          ENDIF
    !          ind_fG = reshape(tensorSumInt( (/ (i, i = 1, neq*ndim) /),neq*ndim*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*ndim*Nfp/))
    !          ind_fe = reshape(tensorSumInt( (/ (i, i = 1, neq) /),neq*(refElPol%face_nodes(ifa,:)-1) ) , (/neq*Nfp/))
    !          ind_ff = (ifa-1)*neq*Nfp + (/ (i, i = 1, neq*Nfp) /)
    !                                                                                xyf    = matmul(refElPol%N1D,Xe(refElPol%face_nodes(ifa,:),:))
    !          xyDer  = matmul(refElPol%Nxi1D,Xe(refElPol%face_nodes(ifa,:),:))
    !
    !          Qx = 0.; Qy = 0.
    !                                                                                ! Loop in 1D Gauss points
    !                                                                                DO g = 1,NGauss
    !
    !                                                                                   ! Calculate the integration weight
    !                                                                                   xyDerNorm_g = norm2(xyDer(g,:))
    !                                                                                   dline = refElPol%gauss_weights1D(g)*xyDerNorm_g

    !                                                       IF (switch%axisym) THEN
    !                                                          dline = dline*xyf(g,1)
    !                                                       END IF
    !                                                                                   ! Unit normal to the boundary
    !                                                                                   t_g = xyDer(g,:)/xyDerNorm_g
    !                                                                                   n_g = [t_g(2), -t_g(1)]
    !                                                                                   ! Contribution of the current integration point to the elemental matrix
    !             Qx = Qx + diff*tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(1)*dline
    !             Qy = Qy + diff*tensorProduct(refElPol%N1D(g,:),refElPol%N1D(g,:))*n_g(2)*dline
    !                                                                                END DO

    !                                                                                ! Elemental assembly
    !                                                                                CALL expand_matrix_B(Qx,Qy,Q_loc)
    !          Qel(ind_fe,ind_fG) = Qel(ind_fe,ind_fG) + transpose(Q_loc)
    !          Lf_t(ind_fG,ind_ff) = Lf_t(ind_fG,ind_ff) + Q_loc*(1-isext)

    !                                                        END DO
    !
    !                                                        Lfel = transpose(Lf_t)

    !      DEALLOCATE(Px,Py,Q_loc,Qx,Qy,Lf_t)

    !                                        END SUBROUTINE elemental_matrices

    SUBROUTINE elemental_matrices(Xe, Pel, Qel, Lfel, diff)
      REAL*8, INTENT(IN)     :: Xe(Np, ndim)
      REAL*8, INTENT(INOUT)  :: Pel(:, :), Qel(:, :), Lfel(:, :)
      REAL*8, INTENT(IN)     :: diff(1:Mesh%Nnodesperelem)
      REAL*8                :: diffg(1:refElPol%Ngauss2d), diffgf(refElPol%Ngauss1d)
      INTEGER*4             :: g, NGauss, ifa
      REAL*8                :: detJ, dvolu, dline, xyDerNorm_g, isext
      REAL*8                :: xy(1:refElPol%Ngauss2d, 1:ndim)
      REAL*8                :: xyf(1:refElPol%Ngauss1d, 1:ndim)
      REAL*8                :: xyDer(1:refElPol%Ngauss1d, 1:ndim)
      REAL*8                :: Jacob(1:ndim, 1:ndim)
      INTEGER*4             :: ind_fe(neq*Nfp), ind_ff(neq*Nfp)
      INTEGER*4             :: ind_fG(neq*ndim*Nfp)
      REAL*8, DIMENSION(Np)  :: Nxg, Nyg, Nx_ax
      REAL*8                :: invJ(1:ndim, 1:ndim), t_g(1:ndim), n_g(1:ndim)
      REAL*8, ALLOCATABLE, DIMENSION(:, :)    :: Px, Py, Q_loc, Qx, Qy, Lf_t

      !***********************************
      !    Volume computation
      !***********************************
      ALLOCATE (Px(1:Np, 1:Np))
      ALLOCATE (Py(1:Np, 1:Np))
      ALLOCATE (Lf_t(Neq*Ndim*Np, Nf*Neq*Nfp))
      Px = 0.; Py = 0.; Lf_t = 0.

      ! Gauss points position
      xy = MATMUL(refElPol%N2D, Xe)

      ! Diffusion at Gauss points
      diffg = MATMUL(refElPol%N2D, diff)

      ! Loop in 2D Gauss points
      Ngauss = refElPol%NGauss2D
      DO g = 1, NGauss

        ! Jacobian
        Jacob = 0.d0
         Jacob(1, 1) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 1))
         Jacob(1, 2) = dot_PRODUCT(refElPol%Nxi2D(g, :), Xe(:, 2))
         Jacob(2, 1) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 1))
         Jacob(2, 2) = dot_PRODUCT(refElPol%Neta2D(g, :), Xe(:, 2))
        detJ = Jacob(1, 1)*Jacob(2, 2) - Jacob(1, 2)*Jacob(2, 1)
        IF (detJ < tol) THEN
            error STOP "Negative jacobian"
         END IF

        ! x and y derivatives of the shape functions
         CALL invert_matrix(Jacob, invJ)
        Nxg = invJ(1, 1)*refElPol%Nxi2D(g, :) + invJ(1, 2)*refElPol%Neta2D(g, :)
        Nyg = invJ(2, 1)*refElPol%Nxi2D(g, :) + invJ(2, 2)*refElPol%Neta2D(g, :)
        IF (switch%axisym) THEN
          Nx_ax = Nxg + 1./xy(g, 1)*refElPol%N2D(g, :)
        ELSE
          Nx_ax = Nxg
        END IF
        ! Integration weight
        dvolu = refElPol%gauss_weights2D(g)*detJ
        IF (switch%axisym) THEN
          dvolu = dvolu*xy(g, 1)
        END IF
        ! Contribution of the current integration point to the elemental matrix
        Px = Px + diffg(g)*TensorProduct(Nxg, refElPol%N2D(g, :))*dvolu
        Py = Py + diffg(g)*TensorProduct(Nyg, refElPol%N2D(g, :))*dvolu
      END DO

      CALL expand_matrix_Bt(Px, Py, Pel)

      !***********************************
      ! Faces computations
      !***********************************
      NGauss = refElPol%Ngauss1D
      ALLOCATE (Q_loc(Neq*Ndim*Nfp, Neq*Nfp))
      ALLOCATE (Qx(Nfp, Nfp))
      ALLOCATE (Qy(Nfp, Nfp))
      Q_loc = 0.

      ! Loop in the faces of the element
      DO ifa = 1, Nf
        isext = 0.
         IF ((Mesh%F(iel, ifa) .GT. Mesh%Nintfaces)) THEN
            IF (Mesh%boundaryFlag(Mesh%F(iel, ifa) - Mesh%Nintfaces) .NE. 0) THEN
            isext = 1.
          END IF
        ENDIF

         ind_fG = RESHAPE(tensorSumInt((/(i, i=1, neq*ndim)/), neq*ndim*(refElPol%face_nodes(ifa, :) - 1)), (/neq*ndim*Nfp/))
         ind_fe = RESHAPE(tensorSumInt((/(i, i=1, neq)/), neq*(refElPol%face_nodes(ifa, :) - 1)), (/neq*Nfp/))
        ind_ff = (ifa - 1)*neq*Nfp + (/(i, i=1, neq*Nfp)/)
         xyf = MATMUL(refElPol%N1D, Xe(refElPol%face_nodes(ifa, :), :))
         xyDer = MATMUL(refElPol%Nxi1D, Xe(refElPol%face_nodes(ifa, :), :))

        ! Diffusion at Gauss points
         diffgf = MATMUL(refElPol%N1D, diff(refElPol%face_nodes(ifa, :)))

        Qx = 0.; Qy = 0.
        ! Loop in 1D Gauss points
        DO g = 1, NGauss

          ! Calculate the integration weight
            xyDerNorm_g = NORM2(xyDer(g, :))
          dline = refElPol%gauss_weights1D(g)*xyDerNorm_g

          IF (switch%axisym) THEN
            dline = dline*xyf(g, 1)
          END IF
          ! Unit normal to the boundary
          t_g = xyDer(g, :)/xyDerNorm_g
          n_g = [t_g(2), -t_g(1)]
          ! Contribution of the current integration point to the elemental matrix
          Qx = Qx + diffgf(g)*tensorProduct(refElPol%N1D(g, :), refElPol%N1D(g, :))*n_g(1)*dline
          Qy = Qy + diffgf(g)*tensorProduct(refElPol%N1D(g, :), refElPol%N1D(g, :))*n_g(2)*dline
        END DO

        ! Elemental assembly
        CALL expand_matrix_B(Qx, Qy, Q_loc)
         Qel(ind_fe, ind_fG) = Qel(ind_fe, ind_fG) + TRANSPOSE(Q_loc)
        Lf_t(ind_fG, ind_ff) = Lf_t(ind_fG, ind_ff) + Q_loc*(1 - isext)

      END DO

      Lfel = TRANSPOSE(Lf_t)
      DEALLOCATE (Px, Py, Q_loc, Qx, Qy, Lf_t)

    END SUBROUTINE elemental_matrices

    !*****************************************
    ! Find shock capturing coefficient in each
    ! element (Persson-Peraire)
    !*****************************************
    SUBROUTINE findCoeffShockCaptur(thresh, eps, invV)
      REAL*8, INTENT(IN)      :: thresh
      REAL*8, INTENT(OUT)     :: eps(:)
      REAL*8, INTENT(IN)      :: invV(refElPol%Nnodes2D, refElPol%Nnodes2D)
      INTEGER*4              :: Ndim, Neq, Nel, Np, Npm1, i
      REAL*8                 :: se(Mesh%Nelems), s0, eps0(Mesh%Nelems)
      REAL*8, ALLOCATABLE     :: up(:, :), udet(:)
      REAL*8, ALLOCATABLE     :: um(:, :), umho(:, :)
      REAL*8, PARAMETER        :: tol = 1e-12

      Ndim = Mesh%ndim
      Neq = phys%Neq
      Nel = Mesh%Nelems
      Np = refElPol%Nnodes2D

      ALLOCATE (up(Mesh%Nelems*refElPol%Nnodes2D, phys%npv))
      ALLOCATE (udet(Mesh%Nelems*refElPol%Nnodes2D))
      ALLOCATE (um(refElPol%Nnodes2D, Mesh%Nelems))
      ALLOCATE (umho(refElPol%Nnodes2D, Mesh%Nelems))
      CALL cons2phys(TRANSPOSE(RESHAPE(sol%u, (/neq, Nel*Np/))), up)
      !        u = transpose(reshape(sol%u,(/neq,Nel*Np/)))
      ! use rho for shock detection
      ! ud = u(:,1);
      ! use gamma for shock detection
      ! ud = u(:,2);
#ifndef TEMPERATURE
      udet = up(:, 1)
#else
      udet = up(:, 9)
#endif

      ! Convert solution into modal expansion
      um = MATMUL(invV, RESHAPE(udet, (/Np, Nel/)))

      ! Solution with only the ho mode
      Npm1 = refElPol%Ndeg*(refElPol%Ndeg + 1)/2
      umho = 0.
      umho(Npm1 + 1:Np, :) = um(Npm1 + 1:Np, :)

      ! Shock detector
      se = LOG10(tol + SUM(umho**2, 1)/(SUM(um**2, 1) + tol))

      ! coefficients
      s0 = LOG10(1./refElPol%Ndeg**4)
      eps0 = numer%sc_coe*Mesh%elemSize/refElPol%Ndeg
      eps = 0.
      ! Loop in elements
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP PRIVATE(i)
      !$OMP DO SCHEDULE(STATIC)
      DO i = 1, Nel
         IF (SUM(um(:, i)**2) .LT. thresh) THEN
          se(i) = -100
        END IF
         IF (se(i) .GE. s0 - numer%sc_sen .AND. se(i) .LE. s0 + numer%sc_sen) THEN
            eps(i) = eps0(i)/2.*(1 + SIN(pi/(2*numer%sc_sen)*(se(i) - s0)))
         ELSE IF (se(i) .GT. s0 - numer%sc_sen) THEN
          eps(i) = eps0(i)
        END IF
      END DO
      !$OMP END DO
      !$OMP END PARALLEL
      DEALLOCATE (up, udet, um, umho)
    END SUBROUTINE findCoeffShockCaptur

    !*****************************************
    ! Set permutations for flipping faces
    !****************************************
    SUBROUTINE set_permutations(n, m, perm)
      INTEGER, INTENT(IN)  :: n, m
      INTEGER, INTENT(OUT) :: perm(:)
      INTEGER              :: i
      INTEGER              :: temp(m, n/m), templr(m, n/m)

      IF (MOD(n, m) .NE. 0) THEN
        WRITE (6, *) 'Error! n must be a multiple of m'
        STOP
      END IF

      templr = 0
      temp = RESHAPE((/(i, i=1, n)/), (/m, n/m/))
      DO i = 1, n/m
        templr(:, i) = temp(:, n/m - i + 1)
      END DO
      perm = RESHAPE(templr, (/n/))
    END SUBROUTINE set_permutations

  END SUBROUTINE HDG_ShockCapturing

  !                                                                !*******************************************
  !                                                                ! Expand matrix B
  !                                                                !*******************************************
  !                                                                SUBROUTINE expand_matrix_B(Bx,By,B)
  !                                                                                real*8,intent(in)    :: Bx(:,:),By(:,:)
  !                                                                                real*8,intent(out)   :: B(1:4*size(Bx,1),1:2*size(Bx,2))
  !                                                                                integer*4            :: i,j,n,m
  !
  !                                                                                n = size(Bx,1)
  !                                                                                m = size(Bx,2)
  !                                                                                B = 0.d0
  !                                                                                DO j= 1,m
  !                                                                                                        DO i = 1,n
  !                                                                                                                                B((i-1)*4+1,(j-1)*2+1) = Bx(i,j)
  !                                                                                                                                B((i-1)*4+2,(j-1)*2+1) = By(i,j)
  !                                                                                                                                B((i-1)*4+3,(j-1)*2+2) = Bx(i,j)
  !                                                                                                                                B((i-1)*4+4,(j-1)*2+2) = By(i,j)
  !                                                                                                        END DO
  !                                                                                END DO
  !                                                                END SUBROUTINE expand_matrix_B
  !
  !
  !

  !                                !*******************************************
  !                                ! Expand matrix B transpose
  !                                !*******************************************
  !                                SUBROUTINE expand_matrix_Bt(Bx,By,B)
  !                                                real*8,intent(in)    :: Bx(:,:),By(:,:)
  !                                                real*8,intent(out)   :: B(1:2*size(Bx,1),1:4*size(Bx,2))
  !                                                integer*4            :: i,j,n,m
  !
  !                                                n = size(Bx,1)
  !                                                m = size(Bx,2)
  !                                                B = 0.d0
  !                                                DO j = 1,m
  !                                                                        DO i = 1,n
  !                                                                                                B((i-1)*2+1,(j-1)*4+1) = Bx(i,j)
  !                                                                                                B((i-1)*2+1,(j-1)*4+2) = By(i,j)
  !                                                                                                B((i-1)*2+2,(j-1)*4+3) = Bx(i,j)
  !                                                                                                B((i-1)*2+2,(j-1)*4+4) = By(i,j)
  !                                                                        END DO
  !                                                END DO
  !                                END SUBROUTINE expand_matrix_Bt

  !*******************************************
  ! Expand matrix B
  !*******************************************
  SUBROUTINE expand_matrix_B(Bx, By, B)
    REAL*8, INTENT(in)    :: Bx(:, :), By(:, :)
    REAL*8, INTENT(out)   :: B(phys%Neq*Mesh%ndim*SIZE(Bx, 1), phys%Neq*SIZE(Bx, 2))
    INTEGER*4            :: i, j, k, n, m

    n = SIZE(Bx, 1)
    m = SIZE(Bx, 2)
    B = 0.d0
    DO j = 1, m
      DO i = 1, n
        DO k = 1, phys%Neq
          B((i - 1)*phys%Neq*Mesh%ndim + 1 + (k - 1)*Mesh%ndim, (j - 1)*phys%Neq + k) = Bx(i, j)
          B((i - 1)*phys%Neq*Mesh%ndim + 2 + (k - 1)*Mesh%ndim, (j - 1)*phys%Neq + k) = By(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE expand_matrix_B

  !*******************************************
  ! Expand matrix B transpose
  !*******************************************
  SUBROUTINE expand_matrix_Bt(Bx, By, B)
    REAL*8, INTENT(in)    :: Bx(:, :), By(:, :)
    REAL*8, INTENT(out)   :: B(phys%Neq*SIZE(Bx, 1), phys%Neq*Mesh%ndim*SIZE(Bx, 2))
    INTEGER*4            :: i, j, k, n, m

    n = SIZE(Bx, 1)
    m = SIZE(Bx, 2)
    B = 0.d0
    DO j = 1, m
      DO i = 1, n
        DO k = 1, phys%Neq
          B((i - 1)*phys%Neq + k, (j - 1)*phys%Neq*Mesh%ndim + 1 + (k - 1)*Mesh%ndim) = Bx(i, j)
          B((i - 1)*phys%Neq + k, (j - 1)*phys%Neq*Mesh%ndim + 2 + (k - 1)*Mesh%ndim) = By(i, j)
        END DO
      END DO
    END DO
  END SUBROUTINE expand_matrix_Bt

  !****************************************************************************************
  !
  !                        APPLY THRESHOLD
  !
  !****************************************************************************************
  SUBROUTINE HDG_applyThreshold(mkelms)
    USE physics, ONLY: cons2phys
    USE reference_element, ONLY: vandermonde_2d, vandermonde_qua
    LOGICAL, INTENT(inout) :: mkelms(:)
    INTEGER*4             :: Neq, Nel, Np, Nfp, iel, ifa, i, nu, nu_tilde, Fa, Nelth
    INTEGER*4             :: ind(refElPol%Nnodes2D), faceNodes(Mesh%Nnodesperface), ind_uf(Mesh%Nnodesperface)
    REAL*8, ALLOCATABLE    :: u(:, :), u_tilde(:, :)
    REAL*8                :: ug(refElPol%NGauss2D), uf(refElPol%Nnodes2D), uav

#ifdef TEMPERATURE
    REAL*8 :: uav_ti, uav_te
#endif

    INTEGER :: limtq, neq_phy
    REAL :: tval

    IF (switch%thresh .EQ. 0) RETURN
    !  IF (MPIvar%glob_id.eq.0) THEN
    !     WRITE(6,*) "Applying threshold"
    !  END IF

    neq = phys%Neq
    neq_phy = phys%Npv
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Nfp = Mesh%Nnodesperface
    nu = SIZE(sol%u)/neq
    nu_tilde = SIZE(sol%u_tilde)/neq

    !****************************************************
    ! Limiting technique:
    ! 1 - based on conservative variables
    ! 2 - based on physical variables
    ! 3 - like filtering
    !****************************************************

    limtq = switch%thresh !3

    IF (limtq == 1) THEN
      !***************************************************
      !         Using limiting on conservative values
      !
      !***************************************************

      !*************************************************
      ! Checking the trace solution for too small values
      !*************************************************
       DO i = 1, SIZE(sol%u_tilde), phys%Neq
          IF (sol%u_tilde(i) .LT. numer%thr) THEN
          sol%u_tilde(i) = numer%thr
        END IF
#ifdef TEMPERATURE
          IF (sol%u_tilde(i + 2) .LT. numer%thr) THEN
          sol%u_tilde(i + 2) = numer%thr
        END IF
          IF (sol%u_tilde(i + 3) .LT. numer%thr) THEN
          sol%u_tilde(i + 3) = numer%thr
        END IF
#endif
      END DO
      !***************************************************
      ! Checking the element solution for too small values
      !***************************************************
       DO i = 1, SIZE(sol%u), phys%Neq
          IF (sol%u(i) .LT. numer%thr) THEN
          sol%u(i) = numer%thr
        END IF
#ifdef TEMPERATURE
          IF (sol%u(i + 2) .LT. numer%thr) THEN
          sol%u(i + 2) = numer%thr
        END IF
          IF (sol%u(i + 3) .LT. numer%thr) THEN
          sol%u(i + 3) = numer%thr
        END IF
#endif
      END DO

    ELSE IF (limtq == 2) THEN
      !***************************************************
      !         Using limiting on physical values
      !
      !***************************************************

      ALLOCATE (u(nu, neq))
      ALLOCATE (u_tilde(nu_tilde, neq))
       u = TRANSPOSE(RESHAPE(sol%u, (/neq, nu/)))
       u_tilde = TRANSPOSE(RESHAPE(sol%u_tilde, (/neq, nu_tilde/)))

      !*************************************************
      ! Checking the trace solution for too small values
      !*************************************************
      DO i = 1, nu_tilde

          IF (u_tilde(i, 1) .LT. numer%thr) THEN
          u_tilde(i, 1) = numer%thr
          !           u_tilde(i,3)= 1.5*phys%Mref*numer%thrpre*u_tilde(i,1) + 0.5*u_tilde(i,2)**2/u_tilde(i,1)
          !           u_tilde(i,4)= 1.5*phys%Mref*numer%thrpre*u_tilde(i,1)
          !           u_tilde(i,2)=0.
        END IF
#ifdef TEMPERATURE
        ! Check if ions temperature is below tolerance in the N-Gamma-Ti-Te model
        tval = 2./(3.*phys%Mref)*(u_tilde(i, 3) - 0.5*u_tilde(i, 2)**2/u_tilde(i, 1))/u_tilde(i, 1)
          IF (tval .LT. numer%thrpre) THEN
          u_tilde(i, 3) = 1.5*phys%Mref*numer%thrpre*u_tilde(i, 1) + 0.5*u_tilde(i, 2)**2/u_tilde(i, 1)
        END IF
        ! Check if electrons temperature is below tolerance in the N-Gamma-Ti-Te model
        tval = 2./(3.*phys%Mref)*u_tilde(i, 4)/u_tilde(i, 1)
          IF (tval .LT. numer%thrpre) THEN
          u_tilde(i, 4) = 1.5*phys%Mref*numer%thrpre*u_tilde(i, 1)
        END IF
#endif
      END DO

      !***************************************************
      ! Checking the element solution for too small values
      !***************************************************
      DO i = 1, nu
          IF (u(i, 1) .LT. numer%thr) THEN
          u(i, 1) = numer%thr
          !           u(i,3)= 1.5*phys%Mref*numer%thrpre*u(i,1) + 0.5*u(i,2)**2/u(i,1)
          !           u(i,4)= 1.5*phys%Mref*numer%thrpre*u(i,1)
          !           u(i,2)=0.
        END IF
#ifdef TEMPERATURE
        ! Check if ions temperature is below tolerance in the N-Gamma-Ti-Te model
        tval = 2./(3.*phys%Mref)*(u(i, 3) - 0.5*u(i, 2)**2/u(i, 1))/u(i, 1)
          IF (tval .LT. numer%thrpre) THEN
          u(i, 3) = 1.5*phys%Mref*numer%thrpre*u(i, 1) + 0.5*u(i, 2)**2/u(i, 1)

        END IF
        ! Check if electrons temperature is below tolerance in the N-Gamma-Ti-Te model
        tval = 2./(3.*phys%Mref)*u(i, 4)/u(i, 1)
          IF (tval .LT. numer%thrpre) THEN
          u(i, 4) = 1.5*phys%Mref*numer%thrpre*u(i, 1)
        END IF
#endif
      END DO

       sol%u = RESHAPE(TRANSPOSE(u), (/neq*nu/))
       sol%u_tilde = RESHAPE(TRANSPOSE(u_tilde), (/neq*nu_tilde/))
      DEALLOCATE (u, u_tilde)

    ELSE IF (limtq == 3) THEN
      !***************************************************
      !         Using filtering style limiting
      !
      !***************************************************

      ALLOCATE (u(nu, neq))
      ALLOCATE (u_tilde(nu_tilde, neq))
       u = TRANSPOSE(RESHAPE(sol%u, (/neq, nu/)))
       u_tilde = TRANSPOSE(RESHAPE(sol%u_tilde, (/neq, nu_tilde/)))

      Nelth = 0
      !*****************
      ! Loop in elements
      !*****************
      DO iel = 1, Nel
        ind = (iel - 1)*Np + (/(i, i=1, Np)/)
          ug = MATMUL(refElPol%N2D, u(ind, 1))

          IF (MINVAL(u(ind, 1)) .GT. numer%thr .AND. MINVAL(ug) .GT. numer%thr .AND. .NOT. mkelms(iel)) CYCLE
        !                                     IF (minval(u(ind,1)).gt.numer%thr .and. minval(ug).gt.numer%thr ) CYCLE
          mkelms(iel) = .TRUE.

        Nelth = Nelth + 1
        ! I set the solution to the max between the
        ! local average and the threshold

          uav = SUM(u(ind, 1))/Np
#ifdef TEMPERATURE
          uav_ti = SUM(u(ind, 3))/Np
          uav_te = SUM(u(ind, 4))/Np
#endif
        uf = 0.
          u(ind, 1) = MAX(numer%thr, uav)
#ifdef TEMPERATURE
          u(ind, 2) = SUM(u(ind, 2))/Np
#else
        u(ind, 2) = 0.
#endif

#ifdef TEMPERATURE
          u(ind, 3) = MAX(1.5*phys%Mref*numer%thrpre*numer%thr, uav_ti)
          u(ind, 4) = MAX(1.5*phys%Mref*numer%thrpre*numer%thr, uav_te)
#endif
        ! Set the face values to the new values
        DO ifa = 1, refElPol%Nfaces
          faceNodes = refElPol%Face_nodes(ifa, :)
          Fa = Mesh%F(iel, ifa)
          ind_uf = (Fa - 1)*Nfp + (/(i, i=1, Nfp)/)
          u_tilde(ind_uf, :) = u(ind(faceNodes), :)
        END DO
      END DO

      ! Store filtered solution and deallcate
       sol%u = RESHAPE(TRANSPOSE(u), (/neq*nu/))
       sol%u_tilde = RESHAPE(TRANSPOSE(u_tilde), (/neq*nu_tilde/))
      DEALLOCATE (u, u_tilde)

       IF (MPIvar%glob_id .EQ. 0) THEN
        WRITE (6, *) "Applied threshold to ", Nelth, ' elements'
      END IF

    END IF

  END SUBROUTINE HDG_applyThreshold

  SUBROUTINE HDG_FilterSolution()
    USE reference_element, ONLY: vandermonde_2d
    INTEGER*4             :: Ndim, Neq, Nel, Np, iel, i, Nu, ieq, Npm1, ind(refElPol%Nnodes2D)
    REAL*8, ALLOCATABLE    :: u(:, :)
    REAL*8                :: un(refElPol%Nnodes2D), um(refElPol%Nnodes2D), umho(refElPol%Nnodes2D)
    REAL*8                :: Vand(refElPol%Nnodes2D, refElPol%Nnodes2D), invVand(refElPol%Nnodes2D, refElPol%Nnodes2D)
    REAL*8, PARAMETER        :: tol = 1e-12
    REAL*8                :: se, s0
    REAL*8                :: soglia(phys%Neq)

    IF (.NOT. switch%filter) RETURN

    IF (MPIvar%glob_id .EQ. 0) THEN
      IF (utils%printint > 0) THEN
        WRITE (6, *) '*************************************************'
        WRITE (6, *) '*           Filtering solution                  *'
        WRITE (6, *) '*************************************************'
      END IF
    END IF

    Ndim = Mesh%ndim
    Neq = phys%Neq
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D
    Npm1 = refElPol%Ndeg*(refElPol%Ndeg + 1)/2

    soglia(1) = numer%thr
    soglia(2) = 0.
    soglia(3) = 1.5*phys%Mref*numer%thrpre*numer%thr
    soglia(4) = 1.5*phys%Mref*numer%thrpre*numer%thr

    s0 = LOG10(1./refElPol%Ndeg**4)

    nu = SIZE(sol%u)/neq
    ALLOCATE (u(nu, neq))
    u = TRANSPOSE(RESHAPE(sol%u, (/neq, nu/)))

    ! Vandermonde matrix
    CALL vandermonde_2d(Vand, refElPol)

    ! Invert Vandermonde matrix
    CALL invert_matrix(Vand, invVand)

    !*****************
    ! Loop in elements
    !*****************
    DO iel = 1, Nel
      ind = (iel - 1)*Np + (/(i, i=1, Np)/)

      DO ieq = 1, Neq

        ! Nodal solution for given element and equation
        un = u(ind, ieq)

        !  Convert solution into modal expansion
          um = MATMUL(invVand, un)

        ! Apply the filter only where there is no plasma:
        ! the solution if very small
          IF (SUM(um**2) > 1e-3) THEN
          CYCLE
        ENDIF
        ! Solution with only the ho mode
        umho = 0.
        umho(Npm1 + 1:Np) = um(Npm1 + 1:Np)

        ! Shock detector
          se = LOG10(tol + SUM(umho**2, 1)/(SUM(um**2, 1) + tol))

        ! coefficients
          s0 = LOG10(1./refElPol%Ndeg**4)

        IF (se > s0) THEN

          !           ! filtered solution for this element/equation
          !                                                                   uf = 0.
          !
          !                                                                   ! I take only the first node: constant solution
          !                                                                   uf(1) = um(1)
          !
          !                                                                   ! Convert back to nodal values
          !                                                                   u(ind,ieq) = matmul(Vand,uf)

          u(ind, ieq) = soglia(ieq)

        END IF
      END DO
    END DO

    ! Store filtered solution and deallcate
    sol%u = RESHAPE(TRANSPOSE(u), (/neq*nu/))
    DEALLOCATE (u)
  END SUBROUTINE HDG_FilterSolution

END MODULE HDG_LimitingTechniques

