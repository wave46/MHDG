!************************************************************
! project: MHDG
! file: preprocess.f90
! date: 16/12/2016
! Preprocess of the mesh: creates Tlin, F, N, flipFaces,
! innerfaces, outerfaces
!************************************************************
MODULE preprocess
  USE types
  USE globals
  USE printutils
  USE MPI_OMP

  IMPLICIT NONE
CONTAINS

  SUBROUTINE mesh_preprocess_serial(ierr)

    INTEGER, INTENT(OUT)            :: ierr
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          IF (MPIvar%glob_id .EQ. 0) THEN
             WRITE (6, *) '*************************************************'
             WRITE (6, *) '*                MESH PREPROCESS                *'
             WRITE (6, *) '*************************************************'
          ENDIF
       END IF
    ENDIF
    !**************************************
    ! Create the linear connectivity matrix
    !**************************************
    ALLOCATE (Mesh%Tlin(Mesh%Nelems, 1:refElPol%Nvertices))
    Mesh%Tlin = Mesh%T(:, 1:refElPol%Nvertices)

    !**********************************
    ! Create nodal connectivity matrix
    !**********************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Creating nodal connectivity "
       END IF
    ENDIF
    CALL createNodalConnectivity()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
       WRITE (6, *) "Nodal connectivity matrix N: "
       CALL displayMatrixInt(Mesh%N)
    END IF
    !**********************************
    ! Create outerfaces and innerfaces
    !**********************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Creating faces structures "
       END IF
    END IF
    CALL GetFaces(ierr)
    IF(ierr .EQ. 0) RETURN
    ! Print out stuff...
    IF(MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, '(A,I13)') " Number of interior faces:           ", Mesh%Nintfaces
          WRITE (6, '(A,I13)') " Number of faces:                    ", Mesh%Nfaces
       END IF
       IF (utils%printint > 1) THEN
          WRITE (6, *) "Exterior faces: "
          CALL displayMatrixInt(Mesh%extfaces)
          WRITE (6, *) "Interior faces: "
          CALL displayMatrixInt(Mesh%intfaces)
       END IF
    ENDIF

    !*******************************
    ! Boundary condition preprocess
    !*******************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Boundary conditions preprocess "
       END IF
    ENDIF
    CALL Bc_preprocess_serial()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
       WRITE (6, *) "Exterior faces after preprocess: "
       CALL displayMatrixInt(Mesh%extfaces)
       WRITE (6, *) "Boundary connectivity matrix Tb after preprocess: "
       CALL displayMatrixInt(Mesh%Tb)
    END IF

    !*****************************************
    ! Create face connectivity F and flipFaces
    !*****************************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Creating face connectivity "
       END IF
    ENDIF
    CALL CreateFaceConnectivity_serial()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
       WRITE (6, *) "Face connectivity matrix F: "
       CALL displayMatrixInt(Mesh%F)
       WRITE (6, *) "F-dir: "
       CALL displayMatrixLog(Mesh%Fdir)
    END IF

    !*****************************************
    ! Compute max and min element size
    !*****************************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Computing element size "
       END IF
    ENDIF
    CALL computeElementSize()
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Done! "
       END IF
    ENDIF

    !*****************************************
    ! Compute puff area
    !*****************************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Computing puff area"
       END IF
    ENDIF
    ! It does not work in parallel (TODO)
#ifndef PARALL
    CALL computePuffArea()
#endif
    ! It does not work in parallel (TODO)
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Puff area:  ", Mesh%puff_area*phys%lscale*phys%lscale, " m^2"
       END IF
    ENDIF
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Computing pump area"
       END IF
    ENDIF
#ifndef PARALL
    CALL computePumpArea()
#endif
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Pump area:  ", Mesh%pump_area*phys%lscale*phys%lscale, " m^2"
       END IF
    ENDIF
#ifndef PARALL
    CALL computeCoreArea()
#endif
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Core area:  ", Mesh%core_area*phys%lscale*phys%lscale, " m^2"
          WRITE (6, *) "Done! "
       END IF
    ENDIF
  END SUBROUTINE mesh_preprocess_serial

  SUBROUTINE mesh_preprocess(ierr)

    INTEGER, INTENT(OUT)            :: ierr

    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          IF (MPIvar%glob_id .EQ. 0) THEN
             WRITE (6, *) '*************************************************'
             WRITE (6, *) '*                MESH PREPROCESS                *'
             WRITE (6, *) '*************************************************'
          ENDIF
       END IF
    ENDIF
    !**************************************
    ! Create the linear connectivity matrix
    !**************************************
    ALLOCATE (Mesh%Tlin(Mesh%Nelems, 1:refElPol%Nvertices))
    Mesh%Tlin = Mesh%T(:, 1:refElPol%Nvertices)

    !**********************************
    ! Create nodal connectivity matrix
    !**********************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Creating nodal connectivity "
       END IF
    ENDIF
    CALL createNodalConnectivity()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
       WRITE (6, *) "Nodal connectivity matrix N: "
       CALL displayMatrixInt(Mesh%N)
    END IF
    !**********************************
    ! Create outerfaces and innerfaces
    !**********************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Creating faces structures "
       END IF
    END IF
    CALL GetFaces(ierr)
    IF(ierr .EQ. 0) RETURN
    ! Print out stuff...
#ifdef PARALL
    IF (utils%printint > 0) THEN
       WRITE (6, '(A,I13)') "Process: ", MPIvar%glob_id, " Number of interior faces:           ", Mesh%Nintfaces
       WRITE (6, '(A,I13)') "Process: ", MPIvar%glob_id, " Number of faces:                    ", Mesh%Nfaces
    END IF
#else
    IF (utils%printint > 0) THEN
       WRITE (6, '(A,I13)') " Number of interior faces:           ", Mesh%Nintfaces
       WRITE (6, '(A,I13)') " Number of faces:                    ", Mesh%Nfaces
    END IF
#endif
    IF (utils%printint > 1) THEN
       WRITE (6, *) "Exterior faces: "
       CALL displayMatrixInt(Mesh%extfaces)
       WRITE (6, *) "Interior faces: "
       CALL displayMatrixInt(Mesh%intfaces)
    END IF

    !*******************************
    ! Boundary condition preprocess
    !*******************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Boundary conditions preprocess "
       END IF
    ENDIF
    CALL Bc_preprocess()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
       WRITE (6, *) "Exterior faces after preprocess: "
       CALL displayMatrixInt(Mesh%extfaces)
       WRITE (6, *) "Boundary connectivity matrix Tb after preprocess: "
       CALL displayMatrixInt(Mesh%Tb)
    END IF

    !*****************************************
    ! Create face connectivity F and flipFaces
    !*****************************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Creating face connectivity "
       END IF
    ENDIF
    CALL CreateFaceConnectivity()
    ! Print out stuff...
    IF (utils%printint > 1) THEN
       WRITE (6, *) "Face connectivity matrix F: "
       CALL displayMatrixInt(Mesh%F)
       WRITE (6, *) "F-dir: "
       CALL displayMatrixLog(Mesh%Fdir)
    END IF

    !*****************************************
    ! Compute max and min element size
    !*****************************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Computing element size "
       END IF
    ENDIF
    CALL computeElementSize()
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Done! "
       END IF
    ENDIF

    !*****************************************
    ! Compute puff area
    !*****************************************
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Computing puff area"
       END IF
    ENDIF
    CALL computePuffArea()
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%puff_area, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Puff area:  ", Mesh%puff_area*phys%lscale*phys%lscale, " m^2"
       END IF
    ENDIF
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Computing pump area"
       END IF
    ENDIF
    CALL computePumpArea()
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%pump_area, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Pump area:  ", Mesh%pump_area*phys%lscale*phys%lscale, " m^2"
       END IF
    ENDIF
    CALL computeCoreArea()
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%core_area, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    IF (MPIvar%glob_id .EQ. 0) THEN
       IF (utils%printint > 0) THEN
          WRITE (6, *) "Core area:  ", Mesh%core_area*phys%lscale*phys%lscale, " m^2"
          WRITE (6, *) "Done! "
       END IF
    ENDIF
  END SUBROUTINE mesh_preprocess

  !********************
  ! Nodal Connectivity
  !********************
  SUBROUTINE createNodalConnectivity()

    INTEGER :: nc, k, iel, nn_Te(1:Mesh%Nnodesperelem)
    INTEGER :: nn(1:Mesh%Nnodes), Te(1:Mesh%Nnodesperelem)
    ! Find the number of colums of the N matrix
    nn = 1
    DO iel = 1, Mesh%Nelems
       Te = Mesh%T(iel, :)
       nn(Te) = nn(Te) + 1
    END DO
    nc = MAXVAL(nn)

    ALLOCATE (Mesh%N(Mesh%Nnodes, nc))
    Mesh%N = 0
    nn = 1
    nn_Te = 0
    DO iel = 1, Mesh%Nelems
       Te = Mesh%T(iel, :)
       nn_Te = nn(Te)
       DO k = 1, Mesh%Nnodesperelem
          Mesh%N(Te(k), nn_Te(k)) = iel
       END DO
       nn(Te) = nn(Te) + 1
    END DO
  END SUBROUTINE createNodalConnectivity

  !********************************************
  ! Create extFaces and intFaces: the exterior
  ! faces have the same numbering than in Tb
  !********************************************
  SUBROUTINE GetFaces(ierr)

    INTEGER :: iel, jel, ifa, jfa, ie, ii, node1
    INTEGER, INTENT(INOUT)          :: ierr
    INTEGER :: Efaces(refElPol%Nfaces, refElPol%Nfacenodeslin), nf(1:refElPol%Nfacenodeslin)
    INTEGER, ALLOCATABLE :: temp_intFaces(:, :)
    LOGICAL, ALLOCATABLE :: markE(:, :), markF(:)

    !Definition of the faces in the reference element
    SELECT CASE (refElPol%elemType)
    CASE (0) !triangle
       Efaces = RESHAPE((/1, 2, 3, 2, 3, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
    CASE (1) !quadrilaterals
       Efaces = RESHAPE((/1, 2, 3, 4, 2, 3, 4, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
    CASE (2) !thetrahedra
       !            Efaces = computeNodesFacesTetra(1);
    CASE (3) !hexahedra
       !        Efaces = computeNodesFacesHexa(1);
    END SELECT

    ALLOCATE (markE(1:Mesh%Nelems, 1:refElPol%Nfaces))
    ALLOCATE (markF(1:Mesh%Nextfaces))
    ALLOCATE (Mesh%extFaces(1:Mesh%Nextfaces, 2))
    ALLOCATE (temp_intFaces(2*Mesh%Nelems, 5))

    ie = 0
    ii = 0
    markE = .FALSE.
    markF = .FALSE.
    DO iel = 1, Mesh%Nelems
       DO ifa = 1, refElPol%Nfaces
          IF (.NOT. markE(iel, ifa)) THEN
             markE(iel, ifa) = .TRUE.
             nf = Mesh%Tlin(iel, Efaces(ifa, :))
             CALL FindElem(nf, iel, jel)
             IF (jel .NE. 0) THEN
                ! Inner element
                ii = ii + 1
                CALL FindFace(nf, Mesh%Tlin(jel, :), Efaces, jfa, node1);
                temp_intFaces(ii, 1) = iel
                temp_intFaces(ii, 2) = ifa
                temp_intFaces(ii, 3) = jel
                temp_intFaces(ii, 4) = jfa
                temp_intFaces(ii, 5) = node1
                markE(jel, jfa) = .TRUE.
             ELSE
                ! Boundary element
                CALL FindCorrespondingFace(iel, ifa, markF, ie, ierr)
                IF(ierr .EQ. 0) RETURN
                !                                                                                           ie = ie + 1
                Mesh%extFaces(ie, 1) = iel
                Mesh%extFaces(ie, 2) = ifa
             END IF
          END IF
       END DO
    END DO

    Mesh%Nintfaces = ii
    ALLOCATE (Mesh%intFaces(ii, 5))
    Mesh%intFaces = temp_intFaces(1:ii, :)
    Mesh%Nfaces = Mesh%Nextfaces + Mesh%Nintfaces
    DEALLOCATE (markE, markF, temp_intFaces)

  CONTAINS

    !*************************************************
    ! Find element: find the neighboring element to the
    ! element i connected by the nodes defined in nf
    !*************************************************
    SUBROUTINE FindElem(nf, iel, jel)
      INTEGER, INTENT(IN) :: nf(1:refElPol%Nfacenodeslin)
      INTEGER, INTENT(IN) :: iel
      INTEGER, INTENT(OUT):: jel
      INTEGER :: i, z, kel, check, k, zel
      INTEGER, DIMENSION(1:SIZE(Mesh%N, 2))::elems1, elemsi

      jel = 0
      elems1 = Mesh%N(nf(1), :)

      ! Loop in the element list associated to the 1st node
      DO k = 1, SIZE(Mesh%N, 2)
         kel = elems1(k)

         IF (kel .EQ. 0) EXIT
         IF (kel .EQ. iel) CYCLE ! I skip the element i

         ! Loop in the connecting node list
         DO i = 2, refElPol%Nfacenodeslin
            elemsi = Mesh%N(nf(i), :)
            check = 0

            ! Loop in the element list of the i-th node
            DO z = 1, SIZE(Mesh%N, 2)
               zel = elemsi(z)

               IF (zel .EQ. 0) EXIT
               IF (zel .EQ. iel) CYCLE ! I skip the element i
               IF (zel == kel) THEN
                  jel = kel
                  check = 1
                  EXIT
               END IF
            END DO

            IF (check == 0) jel = 0 ! if check=0 means I didn't find any equal element in this node
         END DO
         IF (jel .NE. 0) EXIT ! I found the element
      END DO
    END SUBROUTINE FindElem

    !*************************************************
    ! Find face: find the local face number in the
    ! element and the matching node for the first node
    !*************************************************
    SUBROUTINE FindFace(nfi, nodesE, Efaces, jfa, node1)
      INTEGER, INTENT(IN) :: nfi(:)
      INTEGER, INTENT(IN) :: nodesE(:)
      INTEGER, INTENT(IN) :: Efaces(:, :)
      INTEGER, INTENT(OUT):: jfa
      INTEGER, INTENT(OUT):: node1
      INTEGER :: check, fn, i, j
      INTEGER :: nfj(1:refElPol%Nfacenodeslin)

      fn = 0

      ! Find the corresponding face
      DO jfa = 1, refElPol%Nfaces
         nfj = nodesE(Efaces(jfa, :))
         ! check if the set of nodes nfj is equal to the set of nodes nfi
         fn = 0
         DO i = 1, refElPol%Nfacenodeslin
            check = 0
            DO j = 1, refElPol%Nfacenodeslin
               IF (nfi(i) .EQ. nfj(j)) THEN
                  check = 1
                  fn = fn + 1
                  EXIT
               END IF
            END DO
            IF (check .EQ. 0) THEN ! it is not the right face
               EXIT
            END IF
         END DO
         IF (fn .EQ. refElPol%Nfacenodeslin) THEN
            EXIT
         END IF
      END DO

      IF (fn .NE. refElPol%Nfacenodeslin) THEN
         WRITE (6, *) "Error! Corresponding face not found"
         STOP
      END IF

      ! Find the corresponding node
      DO i = 1, refElPol%Nfacenodeslin
         IF (nfj(i) .EQ. nfi(1)) EXIT
      END DO
      IF (i .EQ. refElPol%Nfacenodeslin + 1) THEN
         WRITE (6, *) "Error! Corresponding node not found"
         STOP
      END IF
      node1 = i

    END SUBROUTINE FindFace

    !*************************************************
    ! Find where to place the current exterior face in
    ! order to match the order of Tb
    !*************************************************
    SUBROUTINE FindCorrespondingFace(iel, ifa, markF, iextf, ierr)
      INTEGER, INTENT(IN)    :: iel, ifa
      INTEGER, INTENT(OUT)   :: iextf
      LOGICAL, INTENT(INOUT) :: markF(:)
      INTEGER, INTENT(INOUT) :: ierr
      INTEGER                :: i
      INTEGER                :: nodes(1:Mesh%Nnodesperface), nodes_check(1:Mesh%Nnodesperface)
      LOGICAL                :: isok

      nodes = Mesh%T(iel, refElPol%face_nodes(ifa, :))
      DO iextf = 1, Mesh%Nextfaces

         IF (markF(iextf)) CYCLE
         nodes_check = Mesh%Tb(iextf, :)! <--- Here I use Tb!!
         isok = .TRUE.

         ! Check if "nodes" and "nodes_check" have the same nodes
         DO i = 1, Mesh%Nnodesperface
            IF (nodes(i) .NE. nodes_check(i)) THEN
               isok = .FALSE.
               EXIT
            END IF
         END DO
         IF (isok) THEN
            markF(iextf) = .TRUE.
            EXIT
         END IF
      END DO

      IF (iextf == Mesh%Nextfaces + 1) THEN
         WRITE (6, *) 'Error! Corresponding face in Tb not found. Trying flipping faces.'
         ierr = 0
         RETURN
      END IF
      ierr = 1
    END SUBROUTINE FindCorrespondingFace

  END SUBROUTINE GetFaces

  SUBROUTINE GetFaces_mod(T, int_faces, ext_faces)

    INTEGER, INTENT(IN)                          :: T(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT)            :: int_faces(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT)            :: ext_faces(:,:)
    INTEGER                                      :: Efaces(refElPol%Nfaces, refElPol%Nfacenodeslin), nf(1:refElPol%Nfacenodeslin)
    INTEGER, ALLOCATABLE                         :: temp_intFaces(:, :), temp_extFaces(:, :)
    LOGICAL, ALLOCATABLE                         :: markE(:, :), markF(:)
    INTEGER                                      :: iel, jel, ifa, jfa, ie, ii, jj, node1

    !Definition of the faces in the reference element
    SELECT CASE (refElPol%elemType)
    CASE (0) !triangle
       Efaces = RESHAPE((/1, 2, 3, 2, 3, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
    CASE (1) !quadrilaterals
       Efaces = RESHAPE((/1, 2, 3, 4, 2, 3, 4, 1/), (/refElPol%Nfaces, refElPol%Nfacenodeslin/))
    CASE (2) !thetrahedra
       !            Efaces = computeNodesFacesTetra(1);
    CASE (3) !hexahedra
       !        Efaces = computeNodesFacesHexa(1);
    END SELECT

    ALLOCATE (markE(1:SIZE(T,1), 1:refElPol%Nfaces))
    ALLOCATE (markF(1:Mesh%Nextfaces))
    ALLOCATE (temp_intFaces(2*Mesh%Nelems, 5))
    ALLOCATE (temp_extFaces(2*Mesh%Nelems, 2))

    ie = 0
    ii = 0
    jj = 0
    markE = .FALSE.
    markF = .FALSE.
    DO iel = 1, SIZE(T,1)
       DO ifa = 1, refElPol%Nfaces
          IF (.NOT. markE(iel, ifa)) THEN
             markE(iel, ifa) = .TRUE.
             nf = T(iel, Efaces(ifa, :))
             CALL FindElem(T, nf, iel, jel)


             IF (jel .NE. 0) THEN
                ! Inner element
                ii = ii + 1

                CALL FindFace(nf, T(jel, :), Efaces, jfa, node1);
                temp_intFaces(ii, 1) = iel
                temp_intFaces(ii, 2) = ifa
                temp_intFaces(ii, 3) = jel
                temp_intFaces(ii, 4) = jfa
                temp_intFaces(ii, 5) = node1
                markE(jel, jfa) = .TRUE.

             ELSE
                ! Boundary element
                jj = jj + 1
                temp_extFaces(jj, 1) = iel
                temp_extFaces(jj, 2) = ifa

             END IF

          END IF
       END DO
    END DO

    ALLOCATE(int_faces(ii,5))
    ALLOCATE(ext_faces(jj,2))
    int_faces = temp_intFaces(1:ii, :)
    ext_faces = temp_extFaces(1:jj, :)
    DEALLOCATE (markE, markF, temp_intFaces, temp_extFaces)

  CONTAINS

    !*************************************************
    ! Find element: find the neighboring element to the
    ! element i connected by the nodes defined in nf
    !*************************************************
    SUBROUTINE FindElem(T, nf, iel, jel)
      INTEGER, INTENT(IN) :: nf(:)
      INTEGER, INTENT(IN) :: T(:,:)
      INTEGER, INTENT(IN) :: iel
      INTEGER, INTENT(OUT):: jel
      INTEGER, ALLOCATABLE :: T_temp(:,:), indices(:)
      INTEGER :: i, j, counter

      counter = 0
      DO i = 1, SIZE(T,1)
         IF(i .EQ. iel) CYCLE
         IF(ANY(T(i,:) .EQ. nf(1))) counter = counter + 1
      ENDDO

      ALLOCATE(T_temp(counter, SIZE(T,2)))
      ALLOCATE(indices(counter))

      counter = 1
      DO i = 1, SIZE(T,1)
         IF(i .EQ. iel) CYCLE
         IF(ANY(T(i,:) .EQ. nf(1))) THEN
            T_temp(counter,:) = T(i, :)
            indices(counter) = i
            counter = counter + 1
         ENDIF
      ENDDO

      jel = 0

      DO i = 2, SIZE(nf)
         DO j = 1, SIZE(T_temp,1)
            IF(ANY(T_temp(j,:) .EQ. nf(i))) THEN
               jel = indices(j)
               EXIT
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(T_temp)
      DEALLOCATE(indices)
    END SUBROUTINE FindElem

    !*************************************************
    ! Find face: find the local face number in the
    ! element and the matching node for the first node
    !*************************************************
    SUBROUTINE FindFace(nfi, nodesE, Efaces, jfa, node1)
      INTEGER, INTENT(IN) :: nfi(:)
      INTEGER, INTENT(IN) :: nodesE(:)
      INTEGER, INTENT(IN) :: Efaces(:, :)
      INTEGER, INTENT(OUT):: jfa
      INTEGER, INTENT(OUT):: node1
      LOGICAL, ALLOCATABLE :: check(:,:)
      INTEGER             :: i, j, ii

      jfa = -1
      node1 = -1

      ALLOCATE(check(SIZE(Efaces,2), SIZE(nfi)))

      DO i = 1, SIZE(Efaces,1)
         check = .FALSE.
         DO ii = 1, SIZE(Efaces,2)
            DO j = 1, SIZE(nfi)
               IF(nodesE(Efaces(i,ii)) .EQ. nfi(j)) THEN
                  check(ii,j) = .TRUE.
               ENDIF
            ENDDO
         ENDDO

         IF(ALL(ANY(check,2))) THEN
            jfa = i
            DO ii = 1, SIZE(Efaces,2)
               IF(nodesE(Efaces(i, ii)) .EQ. nfi(1)) THEN
                  node1 = ii
               ENDIF
            ENDDO
            EXIT
         ENDIF
      ENDDO

      IF ((jfa .EQ. -1) .OR. (node1 .EQ. -1)) THEN
         WRITE(*,*) "Corresponding face or node to found. STOP."
         STOP
      ENDIF

      DEALLOCATE(check)

    END SUBROUTINE FindFace
  END SUBROUTINE GetFaces_mod

  !********************************************
  ! Identify the number of boundaries and place
  ! Dirichlet boundaries at the bottom of the
  ! list
  !********************************************
  SUBROUTINE Bc_preprocess()
    INTEGER :: fl, nb, ndir, ifa, bt, idir
    INTEGER :: df(max_num_diff_bc)
    INTEGER, ALLOCATABLE :: aux_Tb(:, :), aux_extFaces(:, :), aux_boundaryflag(:)

    nb = 0
    df = 0
    ndir = 0

    ! Find how many different boundary types are present (different
    ! boundary conditions) and how many Dirichlet faces are used
    DO ifa = 1, Mesh%Nextfaces
       fl = Mesh%boundaryFlag(ifa) ! boundary flag of the current face


#ifdef PARALL
       IF (fl .EQ. 0) CYCLE
#endif
       IF (df(phys%bcflags(fl)) == 0) THEN
          df(phys%bcflags(fl)) = 1
          nb = nb + 1
       END IF
       IF (phys%bcflags(fl) == bc_dirichlet) THEN
          ndir = ndir + 1
       END IF
    END DO


    ! Here I place Dirichlet boundary faces at the bottom of Tb: I change
    ! the order of the faces only in case that at least 1 Dirichlet
    ! face exists and there are more then one type of boundary conditions.
    ! This guarantees compatibility with the Matlab version of the code
    IF (nb > 1 .AND. ndir > 0) THEN
       ALLOCATE (aux_Tb(SIZE(Mesh%Tb, 1), SIZE(Mesh%Tb, 2)))
       ALLOCATE (aux_extFaces(SIZE(Mesh%extFaces, 1), SIZE(Mesh%extFaces, 2)))
       ALLOCATE (aux_boundaryflag(Mesh%NextFaces))
       aux_boundaryflag = 0
       aux_Tb = 0
       aux_extFaces = 0
       idir = 0
       DO ifa = 1, Mesh%Nextfaces
          fl = Mesh%boundaryFlag(ifa) ! boundary flag of the current face
#ifdef PARALL
          IF (fl .EQ. 0) CYCLE ! Add by Benjamin
#endif
          IF (phys%bcflags(fl) == bc_dirichlet) THEN
             aux_Tb(ifa + Mesh%Nextfaces - ndir, :) = Mesh%Tb(ifa, :)
             aux_extFaces(ifa + Mesh%Nextfaces - ndir, :) = Mesh%extFaces(ifa, :)
             aux_boundaryflag(ifa + Mesh%Nextfaces - ndir) = Mesh%boundaryFlag(ifa)
             idir = idir + 1
          ELSE
             aux_Tb(ifa - idir, :) = Mesh%Tb(ifa, :)
             aux_extFaces(ifa - idir, :) = Mesh%extFaces(ifa, :)
             aux_boundaryflag(ifa - idir) = Mesh%boundaryFlag(ifa)
          END IF
       END DO
       Mesh%Tb = aux_Tb
       Mesh%extFaces = aux_extFaces
       Mesh%boundaryFlag = aux_boundaryflag
       DEALLOCATE (aux_Tb, aux_extFaces, aux_boundaryflag)
    END IF

    CALL set_boundary_flag_names()

    ! Periodic faces preprocess
    ALLOCATE(Mesh%periodic_faces(Mesh%Nextfaces))
    Mesh%periodic_faces = 0
    DO ifa = 1, Mesh%Nextfaces
       fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
       IF (fl .EQ. 0) CYCLE ! Add by Benjamin
#endif
       IF (phys%bcflags(fl) == bc_periodic) THEN
          CALL periodic_faces_preprocess()
          EXIT
       END IF
    END DO

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) '*************************************************'
       WRITE (6, *) '*          BOUNDARY CONDITIONS                  *'
       WRITE (6, *) '*************************************************'
    ENDIF
    bt = 0
    DO ifa = 1, Mesh%Nextfaces
       IF (Mesh%boundaryFlag(ifa) .EQ. bt) CYCLE
       fl = Mesh%boundaryFlag(ifa)
#ifdef PARALL
       IF (fl .EQ. 0) CYCLE
#endif
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) 'Boundary type: ', ADJUSTL(TRIM(bc_flag_name(phys%bcflags(fl)))), &
               ' on boundary flagged as: ', ADJUSTL(TRIM(bc_flag_type(fl)))
       ENDIF
       bt = fl
    END DO
    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) " "
    ENDIF

    Mesh%ndir = ndir
    Mesh%ukf = Mesh%Nfaces - ndir

    IF (utils%printint > 0) THEN
       WRITE (6, '(A,I13)') ' Number of different boundaries:     ', nb
       WRITE (6, '(A,I13)') ' Number of Dirichlet faces:          ', ndir
#ifdef PARALL
       WRITE (6, '(A,I13)') ' Number of Ghost faces:              ', Mesh%nghostfaces
#endif
       WRITE (6, *) ' '
    END IF

  CONTAINS


    SUBROUTINE periodic_faces_preprocess()
      INTEGER*4 :: i,j,ifa,ifan,fl,iel,ieln,ni,ne,nin,nen
      REAL*8    :: xi,xe,yi,ye,xin,xen,yin,yen
      DO i = 1, Mesh%Nextfaces
         fl = Mesh%boundaryFlag(i)
         IF (phys%bcflags(fl) == bc_periodic) THEN
            ! Find corresponding face
            iel = Mesh%extFaces(i,1)
            ifa = Mesh%extFaces(i,2)
            ni = refElPol%face_nodes(ifa,1)
            ne = refElPol%face_nodes(ifa,refElPol%Nfacenodes)
            xi=Mesh%X(Mesh%T(iel,ni),1)
            xe=Mesh%X(Mesh%T(iel,ne),1)
            yi=Mesh%X(Mesh%T(iel,ni),2)
            ye=Mesh%X(Mesh%T(iel,ne),2)
            DO j = 1, Mesh%Nextfaces
               IF (i==j) CYCLE
               ieln = Mesh%extFaces(j,1)
               ifan = Mesh%extFaces(j,2)
               nin = refElPol%face_nodes(ifan,1)
               nen = refElPol%face_nodes(ifan,refElPol%Nfacenodes)
               xin=Mesh%X(Mesh%T(ieln,nin),1)
               xen=Mesh%X(Mesh%T(ieln,nen),1)
               yin=Mesh%X(Mesh%T(ieln,nin),2)
               yen=Mesh%X(Mesh%T(ieln,nen),2)
               IF (((ABS(xi-xe) .GT. 1.e-12) .AND. (ABS(xi-xen) .LT. 1e-12) .AND. (ABS(xe-xin) .LT. 1.e-12)) .OR. ((ABS(yi-ye) .GT. 1e-12) .AND. (ABS(yi-yen) .LT. 1.e-12) .AND. (ABS(ye-yin) .LT. 1e-12))) THEN
                  Mesh%periodic_faces(i)=j
                  EXIT
               ENDIF
            END DO
         ENDIF
      END DO

    END SUBROUTINE periodic_faces_preprocess

  END SUBROUTINE Bc_preprocess

  SUBROUTINE Bc_preprocess_serial()
    INTEGER :: fl, nb, ndir, ifa, bt, idir
    INTEGER :: df(max_num_diff_bc)
    INTEGER, ALLOCATABLE :: aux_Tb(:, :), aux_extFaces(:, :), aux_boundaryflag(:)

    nb = 0
    df = 0
    ndir = 0

    ! Find how many different boundary types are present (different
    ! boundary conditions) and how many Dirichlet faces are used
    DO ifa = 1, Mesh%Nextfaces
       fl = Mesh%boundaryFlag(ifa) ! boundary flag of the current face

       IF (df(phys%bcflags(fl)) == 0) THEN
          df(phys%bcflags(fl)) = 1
          nb = nb + 1
       END IF
       IF (phys%bcflags(fl) == bc_dirichlet) THEN
          ndir = ndir + 1
       END IF
    END DO


    ! Here I place Dirichlet boundary faces at the bottom of Tb: I change
    ! the order of the faces only in case that at least 1 Dirichlet
    ! face exists and there are more then one type of boundary conditions.
    ! This guarantees compatibility with the Matlab version of the code
    IF (nb > 1 .AND. ndir > 0) THEN
       ALLOCATE (aux_Tb(SIZE(Mesh%Tb, 1), SIZE(Mesh%Tb, 2)))
       ALLOCATE (aux_extFaces(SIZE(Mesh%extFaces, 1), SIZE(Mesh%extFaces, 2)))
       ALLOCATE (aux_boundaryflag(Mesh%NextFaces))
       aux_boundaryflag = 0
       aux_Tb = 0
       aux_extFaces = 0
       idir = 0
       DO ifa = 1, Mesh%Nextfaces
          fl = Mesh%boundaryFlag(ifa) ! boundary flag of the current face
          IF (phys%bcflags(fl) == bc_dirichlet) THEN
             aux_Tb(ifa + Mesh%Nextfaces - ndir, :) = Mesh%Tb(ifa, :)
             aux_extFaces(ifa + Mesh%Nextfaces - ndir, :) = Mesh%extFaces(ifa, :)
             aux_boundaryflag(ifa + Mesh%Nextfaces - ndir) = Mesh%boundaryFlag(ifa)
             idir = idir + 1
          ELSE
             aux_Tb(ifa - idir, :) = Mesh%Tb(ifa, :)
             aux_extFaces(ifa - idir, :) = Mesh%extFaces(ifa, :)
             aux_boundaryflag(ifa - idir) = Mesh%boundaryFlag(ifa)
          END IF
       END DO
       Mesh%Tb = aux_Tb
       Mesh%extFaces = aux_extFaces
       Mesh%boundaryFlag = aux_boundaryflag
       DEALLOCATE (aux_Tb, aux_extFaces, aux_boundaryflag)
    END IF

    CALL set_boundary_flag_names()

    ! Periodic faces preprocess
    ALLOCATE(Mesh%periodic_faces(Mesh%Nextfaces))
    Mesh%periodic_faces = 0
    DO ifa = 1, Mesh%Nextfaces
       fl = Mesh%boundaryFlag(ifa)
       IF (phys%bcflags(fl) == bc_periodic) THEN
          CALL periodic_faces_preprocess()
          EXIT
       END IF
    END DO

    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) '*************************************************'
       WRITE (6, *) '*          BOUNDARY CONDITIONS                  *'
       WRITE (6, *) '*************************************************'
    ENDIF
    bt = 0
    DO ifa = 1, Mesh%Nextfaces
       IF (Mesh%boundaryFlag(ifa) .EQ. bt) CYCLE
       fl = Mesh%boundaryFlag(ifa)
       IF (MPIvar%glob_id .EQ. 0) THEN
          WRITE (6, *) 'Boundary type: ', ADJUSTL(TRIM(bc_flag_name(phys%bcflags(fl)))), &
               ' on boundary flagged as: ', ADJUSTL(TRIM(bc_flag_type(fl)))
       ENDIF
       bt = fl
    END DO
    IF (MPIvar%glob_id .EQ. 0) THEN
       WRITE (6, *) " "
    ENDIF

    Mesh%ndir = ndir
    Mesh%ukf = Mesh%Nfaces - ndir

    IF (utils%printint > 0) THEN
       WRITE (6, '(A,I13)') ' Number of different boundaries:     ', nb
       WRITE (6, '(A,I13)') ' Number of Dirichlet faces:          ', ndir
       WRITE (6, *) ' '
    END IF

  CONTAINS


    SUBROUTINE periodic_faces_preprocess()
      INTEGER*4 :: i,j,ifa,ifan,fl,iel,ieln,ni,ne,nin,nen
      REAL*8    :: xi,xe,yi,ye,xin,xen,yin,yen
      DO i = 1, Mesh%Nextfaces
         fl = Mesh%boundaryFlag(i)
         IF (phys%bcflags(fl) == bc_periodic) THEN
            ! Find corresponding face
            iel = Mesh%extFaces(i,1)
            ifa = Mesh%extFaces(i,2)
            ni = refElPol%face_nodes(ifa,1)
            ne = refElPol%face_nodes(ifa,refElPol%Nfacenodes)
            xi=Mesh%X(Mesh%T(iel,ni),1)
            xe=Mesh%X(Mesh%T(iel,ne),1)
            yi=Mesh%X(Mesh%T(iel,ni),2)
            ye=Mesh%X(Mesh%T(iel,ne),2)
            DO j = 1, Mesh%Nextfaces
               IF (i==j) CYCLE
               ieln = Mesh%extFaces(j,1)
               ifan = Mesh%extFaces(j,2)
               nin = refElPol%face_nodes(ifan,1)
               nen = refElPol%face_nodes(ifan,refElPol%Nfacenodes)
               xin=Mesh%X(Mesh%T(ieln,nin),1)
               xen=Mesh%X(Mesh%T(ieln,nen),1)
               yin=Mesh%X(Mesh%T(ieln,nin),2)
               yen=Mesh%X(Mesh%T(ieln,nen),2)
               IF (((ABS(xi-xe) .GT. 1.e-12) .AND. (ABS(xi-xen) .LT. 1e-12) .AND. (ABS(xe-xin) .LT. 1.e-12)) .OR. ((ABS(yi-ye) .GT. 1e-12) .AND. (ABS(yi-yen) .LT. 1.e-12) .AND. (ABS(ye-yin) .LT. 1e-12))) THEN
                  Mesh%periodic_faces(i)=j
                  EXIT
               ENDIF
            END DO
         ENDIF
      END DO

    END SUBROUTINE periodic_faces_preprocess

  END SUBROUTINE Bc_preprocess_serial

  SUBROUTINE CreateFaceConnectivity_serial
    USE MPI_OMP
    INTEGER :: ifa, igh
    INTEGER :: infoFace(5), infoFace_ex(2)
    LOGICAL :: isdir
    INTEGER :: ifan,ieln,Fi,Fi_per

    igh = 0
    ALLOCATE (Mesh%F(Mesh%Nelems, refElPol%Nfaces))
    ALLOCATE (Mesh%Fdir(Mesh%Nelems, refElPol%Nfaces))
    ALLOCATE (Mesh%flipFace(Mesh%Nelems, refElPol%Nfaces))
    Mesh%F = 0
    Mesh%Fdir = .FALSE.
    Mesh%flipFace = .FALSE.
    isdir = .FALSE.
    DO ifa = 1, Mesh%Nintfaces
       infoFace = Mesh%intFaces(ifa, :)
       Mesh%F(infoFace(1), infoFace(2)) = ifa
       Mesh%F(infoFace(3), infoFace(4)) = ifa
       IF (infoFace(1) < infoFace(3)) THEN
          Mesh%flipFace(infoFace(3), infoFace(4)) = .TRUE.
       ELSE
          Mesh%flipFace(infoFace(1), infoFace(2)) = .TRUE.
       END IF
    END DO

    DO ifa = 1, Mesh%Nextfaces
       infoFace_ex = Mesh%extFaces(ifa, :)
       IF (Mesh%boundaryFlag(ifa) == 0) THEN
          isdir = .FALSE.
       ELSE
          isdir = (phys%bcflags(Mesh%boundaryFlag(ifa)) == bc_dirichlet)
       ENDIF
       Mesh%F(infoFace_ex(1), infoFace_ex(2)) = ifa + Mesh%Nintfaces
       Mesh%Fdir(infoFace_ex(1), infoFace_ex(2)) = isdir
    END DO
    ! Modify flipface for periodic faces
    DO ifa = 1, Mesh%Nextfaces
       IF (Mesh%periodic_faces(ifa).NE.0) THEN
          infoFace_ex = Mesh%extFaces(ifa, :)
          Fi=Mesh%F(infoFace_ex(1), infoFace_ex(2))
          ieln = Mesh%extfaces(Mesh%periodic_faces(ifa),1)
          ifan = Mesh%extfaces(Mesh%periodic_faces(ifa),2)
          Fi_per = Mesh%F(ieln,ifan)
          IF (Fi>Fi_per) THEN
             Mesh%flipFace(infoFace_ex(1), infoFace_ex(2))=.TRUE.
          ENDIF

       ENDIF
    END DO

  END SUBROUTINE CreateFaceConnectivity_serial


  SUBROUTINE CreateFaceConnectivity()
    USE MPI_OMP
    INTEGER :: ifa, igh
    INTEGER :: infoFace(5), infoFace_ex(2)
    LOGICAL :: isdir
    INTEGER :: ifan,ieln,Fi,Fi_per

    igh = 0
    ALLOCATE (Mesh%F(Mesh%Nelems, refElPol%Nfaces))
    ALLOCATE (Mesh%Fdir(Mesh%Nelems, refElPol%Nfaces))
    ALLOCATE (Mesh%flipFace(Mesh%Nelems, refElPol%Nfaces))
    Mesh%F = 0
    Mesh%Fdir = .FALSE.
    Mesh%flipFace = .FALSE.
    isdir = .FALSE.
    DO ifa = 1, Mesh%Nintfaces
       infoFace = Mesh%intFaces(ifa, :)
       Mesh%F(infoFace(1), infoFace(2)) = ifa
       Mesh%F(infoFace(3), infoFace(4)) = ifa
       IF (infoFace(1) < infoFace(3)) THEN
          Mesh%flipFace(infoFace(3), infoFace(4)) = .TRUE.
       ELSE
          Mesh%flipFace(infoFace(1), infoFace(2)) = .TRUE.
       END IF
#ifdef PARALL
       IF (Mesh%ghostFaces(ifa) .EQ. 1) THEN
          igh = igh + 1
       END IF
#endif
    END DO

    DO ifa = 1, Mesh%Nextfaces
       infoFace_ex = Mesh%extFaces(ifa, :)
       IF (Mesh%boundaryFlag(ifa) == 0) THEN
          isdir = .FALSE.
       ELSE
          isdir = (phys%bcflags(Mesh%boundaryFlag(ifa)) == bc_dirichlet)
       ENDIF
       Mesh%F(infoFace_ex(1), infoFace_ex(2)) = ifa + Mesh%Nintfaces
       Mesh%Fdir(infoFace_ex(1), infoFace_ex(2)) = isdir
#ifdef PARALL
       IF (Mesh%ghostFaces(ifa + Mesh%Nintfaces) .EQ. 1) THEN
          igh = igh + 1
          IF (Mesh%ghostFlp(igh) .EQ. 1) THEN
             Mesh%flipFace(infoFace_ex(1), infoFace_ex(2)) = .TRUE.
          END IF
       END IF
#endif
    END DO
    ! Modify flipface for periodic faces
    DO ifa = 1, Mesh%Nextfaces
       IF (Mesh%periodic_faces(ifa).NE.0) THEN
          infoFace_ex = Mesh%extFaces(ifa, :)
          Fi=Mesh%F(infoFace_ex(1), infoFace_ex(2))
          ieln = Mesh%extfaces(Mesh%periodic_faces(ifa),1)
          ifan = Mesh%extfaces(Mesh%periodic_faces(ifa),2)
          Fi_per = Mesh%F(ieln,ifan)
          IF (Fi>Fi_per) THEN
             Mesh%flipFace(infoFace_ex(1), infoFace_ex(2))=.TRUE.
          ENDIF

       ENDIF
    END DO

  END SUBROUTINE CreateFaceConnectivity

  !********************************************
  ! Compute the element size in the
  ! mesh (the length of the edge 2D-3D is
  ! considered)
  !********************************************
  SUBROUTINE computeElementSize()
    INTEGER                         :: i
    REAL*8                          :: h1, h2, h3, h4
    REAL*8, DIMENSION(Mesh%ndim)     :: p1, p2, p3, p4

    ALLOCATE (Mesh%elemSize(Mesh%Nelems))
    Mesh%elemSize = 0.

    ! Loop in elements
    IF (refElPol%elemType .EQ. 0) THEN
       !$OMP PARALLEL PRIVATE(i,p1,p2,p3,h1,h2,h3)
       !$OMP DO
       DO i = 1, Mesh%Nelems
          p1 = Mesh%X(Mesh%Tlin(i, 1), :)
          p2 = Mesh%X(Mesh%Tlin(i, 2), :)
          p3 = Mesh%X(Mesh%Tlin(i, 3), :)
          h1 = NORM2(p1 - p2)
          h2 = NORM2(p1 - p3)
          h3 = NORM2(p3 - p2)
          Mesh%elemSize(i) = MIN(h1, h2, h3)
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
    ELSEIF (refElPol%elemType .EQ. 1) THEN
       !$OMP PARALLEL PRIVATE(i,p1,p2,p3,p4,h1,h2,h3,h4)
       !$OMP DO
       DO i = 1, Mesh%Nelems
          p1 = Mesh%X(Mesh%Tlin(i, 1), :)
          p2 = Mesh%X(Mesh%Tlin(i, 2), :)
          p3 = Mesh%X(Mesh%Tlin(i, 3), :)
          p4 = Mesh%X(Mesh%Tlin(i, 4), :)
          h1 = NORM2(p1 - p2)
          h2 = NORM2(p1 - p3)
          h3 = NORM2(p3 - p2)
          h4 = NORM2(p4 - p3)
          Mesh%elemSize(i) = MIN(h1, h2, h3, h4)
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
    END IF

  END SUBROUTINE computeElementSize


  SUBROUTINE computePuffArea()
    REAL*8   :: Xf(refElPol%Nfacenodes,2),xyg(refElPol%NGauss1D,2),xyg_d(refElPol%NGauss1D,2),dline
    INTEGER  :: i,el,fa,fl,g
    REAL*8   :: xyDerNorm_g


    Mesh%puff_area = 0.

    DO i = 1, Mesh%Nextfaces

       fl = Mesh%boundaryFlag(i)

#ifdef PARALL
       IF (fl .EQ. 0) CYCLE
#endif

       IF (phys%bcflags(fl) .NE. bc_BohmPuff) THEN
          CYCLE
       END IF

       el = Mesh%extfaces(i,1)
       fa = Mesh%extfaces(i,2)
       Xf = Mesh%X(Mesh%T(el,refElPol%face_nodes(fa,:)),:)
       xyg = MATMUL(refElPol%N1D,Xf)
       xyg_d = MATMUL(refElPol%Nxi1D,Xf)
#ifdef PARALL
       IF (Mesh%ghostFaces(Mesh%Nextfaces+i) .EQ. 0) THEN
#endif
          DO g = 1, refElPol%NGauss1D
             xyDerNorm_g = NORM2(xyg_d(g,:))
             dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
             dline = dline*xyg(g,1)
             Mesh%puff_area = Mesh%puff_area + 2*pi*dline
          END DO
#ifdef PARALL
       END IF
#endif
    END DO

  END SUBROUTINE computePuffArea

  SUBROUTINE computePumpArea()
    REAL*8   :: Xf(refElPol%Nfacenodes,2),xyg(refElPol%NGauss1D,2),xyg_d(refElPol%NGauss1D,2),dline
    INTEGER  :: i,el,fa,fl,g
    REAL*8   :: xyDerNorm_g


    Mesh%pump_area = 0.

    DO i = 1, Mesh%Nextfaces

       fl = Mesh%boundaryFlag(i)

#ifdef PARALL
       IF (fl .EQ. 0) CYCLE
#endif

       IF (phys%bcflags(fl) .NE. bc_BohmPump) THEN
          CYCLE
       END IF

       el = Mesh%extfaces(i,1)
       fa = Mesh%extfaces(i,2)
       Xf = Mesh%X(Mesh%T(el,refElPol%face_nodes(fa,:)),:)
       xyg = MATMUL(refElPol%N1D,Xf)
       xyg_d = MATMUL(refElPol%Nxi1D,Xf)
#ifdef PARALL
       IF (Mesh%ghostFaces(Mesh%Nextfaces+i) .EQ. 0) THEN
#endif
          DO g = 1, refElPol%NGauss1D
             xyDerNorm_g = NORM2(xyg_d(g,:))
             dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
             dline = dline*xyg(g,1)
             Mesh%pump_area = Mesh%pump_area + 2*pi*dline
          END DO
#ifdef PARALL
       END IF
#endif
    END DO

  END SUBROUTINE computePumpArea



  SUBROUTINE computeCoreArea()
    REAL*8   :: Xf(refElPol%Nfacenodes,2),xyg(refElPol%NGauss1D,2),xyg_d(refElPol%NGauss1D,2),dline
    INTEGER  :: i,el,fa,fl,g
    REAL*8   :: xyDerNorm_g


    Mesh%core_area = 0.

    DO i = 1, Mesh%Nextfaces

       fl = Mesh%boundaryFlag(i)

#ifdef PARALL
       IF (fl .EQ. 0) CYCLE
#endif

       IF (phys%bcflags(fl) > 10) THEN
          CYCLE
       END IF

       el = Mesh%extfaces(i,1)
       fa = Mesh%extfaces(i,2)
       Xf = Mesh%X(Mesh%T(el,refElPol%face_nodes(fa,:)),:)
       xyg = MATMUL(refElPol%N1D,Xf)
       xyg_d = MATMUL(refElPol%Nxi1D,Xf)
#ifdef PARALL
       IF (Mesh%ghostElems(el) .EQ. 0) THEN
#endif
          DO g = 1, refElPol%NGauss1D
             xyDerNorm_g = NORM2(xyg_d(g,:))
             dline = refElPol%gauss_weights1D(g)*xyDerNorm_g
             dline = dline*xyg(g,1)
             Mesh%core_area = Mesh%core_area + 2*pi*dline
          END DO
#ifdef PARALL
       END IF
#endif
    END DO

  END SUBROUTINE computeCoreArea


END MODULE preprocess
