MODULE domain_decomposition_module
   USE MPI_OMP
   USE globals
   USE adaptivity_common_module, ONLY: unique_1D, unique_stable
   USE Communications, ONLY: reshape_transpose_permute, reshape_transpose_permute_4D, flatten_row_major, flatten_row_major_4D
   USE preprocess, ONLY: GetFaces_mod, mesh_preprocess

   IMPLICIT NONE

   TYPE DD_mesh
      INTEGER, ALLOCATABLE                 :: TT(:, :), intfaces_loc(:, :), extfaces_loc(:, :), TTb_CON(:, :), Tb(:, :)
      REAL*8, ALLOCATABLE                  :: XX(:, :)
      INTEGER, ALLOCATABLE                 :: loc2glob_fa(:), loc2glob_el(:), reptFaces(:), ghostFaces(:), ghostElems(:), loc2glob_no(:), ghostFlp(:),&
      & mapextfa_loc(:), Tb_sizes(:), GhelsPro(:), GhelsLoc(:), GhostPro(:), GhostLoc(:), boundaryFlag(:)
      INTEGER                              :: Nfaces, Nnodes, Nelems, Ndim, Nnodesperelem, Nnodesperface, Nintfaces, NextFaces, elemType
   END TYPE DD_mesh

CONTAINS


   SUBROUTINE reorder_mesh(perm_elems)
      USE preprocess, ONLY: mesh_preprocess_serial, mesh_preprocess
      INTEGER, ALLOCATABLE, INTENT(OUT)       :: perm_elems(:)
      INTEGER                                 :: Nelems, Nnodesperelem, Nnodes
      INTEGER                                 :: iel
      INTEGER, ALLOCATABLE                    :: T_r(:, :)
      REAL*8, ALLOCATABLE                     :: X_r(:, :)
      INTEGER                                 :: ierr

      Nelems = Mesh%Nelems
      Nnodesperelem = Mesh%Nnodesperelem
      Nnodes = Mesh%Nnodes
#ifdef PARALL
      CALL mesh_preprocess_serial(ierr)
#else
      CALL mesh_preprocess(ierr)
#endif

      ALLOCATE(perm_elems(Nelems))
      perm_elems = 0

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "*           CONNECTIVITY REORDERING             *"
         WRITE (*, *) "*************************************************"

         ALLOCATE (T_r(Nelems, Nnodesperelem))
         ALLOCATE (X_r(Nnodes, mesh%Ndim))
         T_r = 0
         X_r = 0.d0

         CALL element_adjacency()

         DO iel = 1, Nelems
            T_r(iel, :) = Mesh%T(perm_elems(iel), :)
            X_r(T_r(iel, :), :) = Mesh%X(Mesh%T(perm_elems(iel), :), :)
         END DO


         Mesh%T = T_r
         Mesh%X = X_r

         DEALLOCATE (T_r, X_r)

      END IF

      CALL MPI_BCAST(Mesh%T,Nelems*Nnodesperelem,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Mesh%X,Nnodes*Mesh%Ndim,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(perm_elems,Nelems,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

   CONTAINS

      SUBROUTINE element_adjacency()

         INTEGER, ALLOCATABLE   :: edge(:,:), elem_of_edge(:)
         INTEGER                :: iel, k, tmp, ivert, jvert, e1, e2
         INTEGER, ALLOCATABLE   :: A(:,:)

         ALLOCATE(A(Nelems,Nelems))
         A = 0

         ALLOCATE(edge(3*Nelems,2))
         ALLOCATE(elem_of_edge(3*Nelems))

         k = 0

         DO iel = 1, Nelems

            ! Edge 1
            k = k + 1
            edge(k,:) = [Mesh%T(iel,1), Mesh%T(iel,2)]
            elem_of_edge(k) = iel

            ! Edge 2
            k = k + 1
            edge(k,:) = [Mesh%T(iel,2), Mesh%T(iel,3)]
            elem_of_edge(k) = iel

            ! Edge 3
            k = k + 1
            edge(k,:) = [Mesh%T(iel,3), Mesh%T(iel,1)]
            elem_of_edge(k) = iel

         END DO


         DO ivert = 1, 3*Nelems
            IF (edge(ivert,1) > edge(ivert,2)) THEN
               tmp = edge(ivert,1)
               edge(ivert,1) = edge(ivert,2)
               edge(ivert,2) = tmp
            END IF
         END DO

         DO ivert = 1, 3*Nelems-1
            DO jvert = ivert+1, 3*Nelems

               IF (edge(ivert,1) == edge(jvert,1) .AND. &
                  edge(ivert,2) == edge(jvert,2)) THEN

                  e1 = elem_of_edge(ivert)
                  e2 = elem_of_edge(jvert)

                  IF (e1 .NE. e2) THEN
                     A(e1,e2) = 1
                     A(e2,e1) = 1
                  END IF

               END IF

            END DO
         END DO

         CALL rcm(A, Nelems, perm_elems)

      END SUBROUTINE element_adjacency

      SUBROUTINE rcm(A, n, perm_elems)

         IMPLICIT NONE

         INTEGER, INTENT(IN)       :: n
         INTEGER, INTENT(IN)       :: A(n,n)
         INTEGER, INTENT(OUT)      :: perm_elems(n)

         INTEGER                   :: degree(n)
         INTEGER                   :: visited(n)
         INTEGER                   :: queue(n)
         INTEGER                   :: neighbors(n)
         INTEGER                   :: i, j
         INTEGER                   :: start, current
         INTEGER                   :: nvisited
         INTEGER                   :: nneigh
         INTEGER                   :: tmp
         INTEGER                   :: head, tail

         degree = 0

         !------------------------------------------------------------
         ! Compute degrees
         !------------------------------------------------------------
         DO i = 1, n
            DO j = 1, n
               degree(i) = degree(i) + A(i,j)
            END DO
         END DO

         visited = 0
         perm_elems = 0
         nvisited = 0

         !------------------------------------------------------------
         ! Process each connected component
         !------------------------------------------------------------
         DO WHILE (nvisited < n)

            ! Find unvisited node with minimum degree
            start = -1
            DO i = 1, n
               IF (visited(i) == 0) THEN
                  IF (start == -1) THEN
                     start = i
                  ELSE IF (degree(i) < degree(start)) THEN
                     start = i
                  END IF
               END IF
            END DO

            !--------------------------------------------------------
            ! BFS
            !--------------------------------------------------------
            queue = 0
            queue(1) = start
            visited(start) = 1

            head = 1
            tail = 1

            DO WHILE (head <= tail)

               current = queue(head)
               head = head + 1

               nvisited = nvisited + 1
               perm_elems(nvisited) = current

               ! Find unvisited neighbors
               nneigh = 0

               DO j = 1, n

                  IF (A(current,j) /= 0 .AND. visited(j) == 0) THEN

                     nneigh = nneigh + 1
                     neighbors(nneigh) = j

                     visited(j) = 1

                  END IF

               END DO

               !----------------------------------------------------
               ! Sort neighbors by increasing degree
               !----------------------------------------------------
               DO i = 1, nneigh-1
                  DO j = i+1, nneigh

                     IF (degree(neighbors(j)) < &
                        degree(neighbors(i))) THEN

                        tmp = neighbors(i)
                        neighbors(i) = neighbors(j)
                        neighbors(j) = tmp

                     END IF

                  END DO
               END DO

               !----------------------------------------------------
               ! Add neighbors to BFS queue
               !----------------------------------------------------
               DO i = 1, nneigh

                  tail = tail + 1
                  queue(tail) = neighbors(i)

               END DO

            END DO

         END DO

         !------------------------------------------------------------
         ! Reverse Cuthill-McKee
         !------------------------------------------------------------
         DO i = 1, n/2
            tmp = perm_elems(i)
            perm_elems(i) = perm_elems(n-i+1)
            perm_elems(n-i+1) = tmp
         END DO

      END SUBROUTINE rcm

   END SUBROUTINE reorder_mesh

#ifndef TOR3D
   SUBROUTINE reorder_solution(perm_elems)
      USE initialization, ONLY: extractFaceSolution
      INTEGER, INTENT(IN)                :: perm_elems(:)
      !integer, intent(in)                :: perm_faces(:)
      REAL*8, ALLOCATABLE, DIMENSION(:)  :: u, q
      INTEGER                            :: Nelems, Nnodesperelem, Neq, iel, i, ip, Ndim
      INTEGER, DIMENSION(phys%Neq)       :: idx_ord_u, idx_noord_u, indices_u
      INTEGER, DIMENSION(phys%Neq*Mesh%Ndim)  :: idx_ord_q, idx_noord_q, indices_q

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "*              SOLUTION REORDERING              *"
         WRITE (*, *) "*************************************************"
      END IF

      Nelems = Mesh%Nelems
      Nnodesperelem = Mesh%Nnodesperelem
      Neq = phys%Neq
      Ndim = Mesh%Ndim

      indices_u = (/ (i, i = 1, Neq) /)
      indices_q = (/ (i, i = 1, Neq*Ndim) /)

      ALLOCATE(u(Nelems*Nnodesperelem*Neq))
      ALLOCATE(q(Nelems*Nnodesperelem*Neq*Ndim))

      u = 0.
      q = 0.

      ! Reorder the solution vectors according to the reordered mesh
      DO iel = 1, Nelems
         DO ip = 1, Nnodesperelem
            idx_ord_u = (iel - 1)*Nnodesperelem*Neq + (ip-1)*Neq + indices_u
            idx_noord_u = (perm_elems(iel) - 1)*Nnodesperelem*Neq + (ip-1)*Neq + indices_u
            u(idx_ord_u) = sol%u(idx_noord_u)
         END DO
      END DO

      DO iel = 1, Nelems
         DO ip = 1, Nnodesperelem
            idx_ord_q = (iel - 1)*Nnodesperelem*Neq*Ndim + (ip-1)*Neq*Ndim + indices_q
            idx_noord_q = (perm_elems(iel) - 1)*Nnodesperelem*Neq*Ndim + (ip-1)*Neq*Ndim + indices_q
            q(idx_ord_q) = sol%q(idx_noord_q)
         END DO
      END DO

      sol%u = u
      sol%q = q

      DEALLOCATE (sol%u_tilde)
      ALLOCATE (sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
      CALL extractFaceSolution()

      DEALLOCATE(u, q)

   END SUBROUTINE reorder_solution

#else
   SUBROUTINE reorder_solution(perm_elems)
      USE initialization, ONLY: extractFaceSolution
      INTEGER, INTENT(IN)                                            :: perm_elems(:)
      INTEGER                                                        :: Nelems, Nnodesperelem, Nfaces, Nnodesperface, Neq, iel, i, ip, Ndim, start_range, end_range, itor, Ntor, Nnodesperelem_tor1d
      INTEGER, ALLOCATABLE, DIMENSION(:)                             :: indices_u, indices_q,idx_ord_u, idx_noord_u, idx_ord_q, idx_noord_q
      REAL*8, ALLOCATABLE, DIMENSION(:)                              :: u, q

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "              SOLUTION REORDERING                "
         WRITE (*, *) "*************************************************"
      END IF


      Nelems = Mesh%Nelems
      Nfaces = Mesh%Nfaces
      Nnodesperface = Mesh%Nnodesperface
      Nnodesperelem = Mesh%Nnodesperelem
      Neq = phys%Neq
      Ndim = Mesh%Ndim + 1
      Ntor = numer%Ntor
      Nnodesperelem_tor1d = numer%ptor + 1

      ALLOCATE(u(Ntor*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq))
      ALLOCATE(q(Ntor*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq*Ndim))

      ALLOCATE(indices_u(Nnodesperelem*Neq))
      ALLOCATE(indices_q(Nnodesperelem*Neq*Ndim))
      ALLOCATE(idx_ord_u(Nelems*Nnodesperelem*Neq))
      ALLOCATE(idx_ord_q(Nelems*Nnodesperelem*Neq*Ndim))
      ALLOCATE(idx_noord_u(Nelems*Nnodesperelem*Neq))
      ALLOCATE(idx_noord_q(Nelems*Nnodesperelem*Neq*Ndim))

      indices_u = (/ (i, i = 1, Neq*Nnodesperelem) /)
      indices_q = (/ (i, i = 1, Neq*Ndim*Nnodesperelem) /)

      u = 0.
      q = 0.

      DO itor = 1,Ntor
         DO ip = 1,Nnodesperelem_tor1d
            start_range = 1
            end_range = 1
            DO iel = 1, Nelems
               end_range = start_range + Neq*Nnodesperelem-1

               indices_u = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ip-1)*Neq*Nnodesperelem + (perm_elems(iel)-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + indices_u
               idx_noord_u(start_range:end_range) = indices_u
               indices_u = (/ (i, i = 1, Neq*Nnodesperelem) /)

               indices_u = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ip-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + indices_u
               idx_ord_u(start_range:end_range) = indices_u
               indices_u = (/ (i, i = 1, Neq*Nnodesperelem) /)

               start_range = end_range + 1
            END DO

            ! Reorder the solution vectors according to the reordered mesh
            u(idx_ord_u) = sol%u(idx_noord_u)
         END DO
      END DO

      DO itor = 1,Ntor
         DO ip = 1,Nnodesperelem_tor1d
            start_range = 1
            end_range = 1
            DO iel = 1, Nelems
               end_range = start_range + Neq*Nnodesperelem*Ndim-1
               indices_q = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq*Ndim + (ip-1)*Neq*Nnodesperelem*Ndim + (perm_elems(iel)-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d*Ndim + indices_q
               idx_ord_q(start_range:end_range) = indices_q
               indices_q = (/ (i, i = 1, Nnodesperelem*Neq*Ndim) /)
               indices_q = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq*Ndim + (ip-1)*Neq*Nnodesperelem*Ndim + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d*Ndim + indices_q
               idx_noord_q(start_range:end_range) = indices_q
               indices_q = (/ (i, i = 1, Nnodesperelem*Neq*Ndim) /)
               start_range = end_range + 1
            END DO

            ! Reorder the solution vectors according to the reordered mesh
            q(idx_ord_q) = sol%q(idx_noord_q)
         END DO
      END DO

      sol%u = u
      sol%q = q

      DEALLOCATE (sol%u_tilde)
      ALLOCATE (sol%u_tilde(Ntor*(Nelems*Nnodesperelem*Neq + Nfaces*Nnodesperface*Nnodesperelem_tor1d*Neq)))
      CALL extractFaceSolution()

      DEALLOCATE(idx_ord_u, idx_noord_u, idx_ord_q, idx_noord_q)
      DEALLOCATE(u, q)

   END SUBROUTINE reorder_solution
#endif

#ifdef PARALL
   SUBROUTINE domain_decomposition(mesh_name)
      USE in_out, ONLY: HDF5_save_mesh_serial
      CHARACTER(LEN=1024), INTENT(IN) :: mesh_name
      CHARACTER(LEN=255)              :: h5_filename, nid
      INTEGER                         :: ierr, i

      IF (MPIvar%glob_size .GT. 1) THEN
         ! split the mesh
         IF ((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
            CALL split_mesh(MPIvar%glob_size, 1, verbose=.FALSE.)
         ELSE
            CALL split_mesh(MPIvar%glob_size, 2, verbose=.FALSE.)
         END IF
         CALL mesh_preprocess(ierr)

         IF (switch%gmsh2h5) THEN
            WRITE (nid, '(I3.3)') MPIvar%glob_id
            h5_filename = TRIM(ADJUSTL(mesh_name)) // TRIM(ADJUSTL(nid))
            WRITE (nid, '(I3.3)') MPIvar%glob_size
            h5_filename = TRIM(ADJUSTL(h5_filename)) // TRIM(ADJUSTL(nid))
            Mesh%X = Mesh%X*phys%lscale
            CALL HDF5_save_mesh_serial(Mesh, h5_filename)
            Mesh%X = Mesh%X/phys%lscale
         END IF
      ELSE

#ifdef TOR3D
         IF (ASSOCIATED(Mesh%toroidal)) THEN
            DEALLOCATE (Mesh%toroidal)
            NULLIFY (Mesh%toroidal)
         END IF
#endif
         ALLOCATE (Mesh%ghostelems(Mesh%Nelems))
         ALLOCATE (Mesh%ghostfaces(Mesh%Nfaces))
         ALLOCATE (Mesh%ghostpro(Mesh%Nfaces))
         ALLOCATE (Mesh%ghostloc(Mesh%Nfaces))
         ALLOCATE (Mesh%loc2glob_el(Mesh%Nelems))
         ALLOCATE (Mesh%loc2glob_fa(Mesh%Nfaces))
         ALLOCATE (Mesh%loc2glob_nodes(Mesh%Nnodes))

         Mesh%ghostelems = 0
         Mesh%ghostfaces = 0
         Mesh%ghostpro = 0
         Mesh%ghostLoc = 0
         Mesh%loc2glob_el = [(i, i=1, Mesh%Nelems)]
         Mesh%loc2glob_fa = [(i, i=1, Mesh%Nfaces)]
         Mesh%loc2glob_nodes = [(i, i=1, Mesh%Nnodes)]

         Mesh%nghostelems = 0
         Mesh%nghostfaces = 0
         Mesh%Nel_glob = Mesh%Nelems
         Mesh%Nfa_glob = Mesh%Nfaces
         Mesh%Nno_glob = Mesh%Nnodes
         Mesh%Nintfaces_glob = Mesh%Nintfaces
         Mesh%Nextfaces_glob = Mesh%Nextfaces
         Mesh%Ndir_glob = Mesh%Ndir
         Mesh%Ngho_glob = 0
         Mesh%Nintfaces_nogho = Mesh%Nintfaces
         Mesh%Nextfaces_nogho = Mesh%Nextfaces

      END IF

   END SUBROUTINE domain_decomposition

   SUBROUTINE Bfield_Jtor_decomposition

      REAL*8, ALLOCATABLE       :: B_glob(:, :), magnetic_flux_glob(:), magnetic_psi_glob(:), Bperturb_glob(:, :), Jtor_glob(:)
      INTEGER                   :: i, Nnodes_local, Nnodes_glob, Nnodes_tor
      INTEGER, ALLOCATABLE      :: indices(:), indices_loc(:), indices_glob(:)

      IF (MPIvar%glob_size .EQ. 1) RETURN

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "*           BFIELD JTOR DECOMPOSITION           *"
         WRITE (*, *) "*************************************************"
      END IF

#ifdef TOR3D
      Nnodes_local = Mesh%Nnodes*Mesh%nnodes_toroidal
      Nnodes_glob = Mesh%Nno_glob*Mesh%nnodes_toroidal
      Nnodes_tor = Mesh%Nnodes_toroidal
#else
      Nnodes_local = Mesh%Nnodes
      Nnodes_glob = Mesh%Nno_glob
      Nnodes_tor = 1
#endif

      ALLOCATE (indices(Nnodes_local/Nnodes_tor))
      ALLOCATE (indices_loc(Nnodes_local/Nnodes_tor))
      ALLOCATE (indices_glob(Nnodes_local/Nnodes_tor))
      indices = [(i, i=1, Nnodes_local/Nnodes_tor)]
      indices_loc = 0
      indices_glob = 0


      ALLOCATE (B_glob(Nnodes_glob, 3))
      ALLOCATE (magnetic_flux_glob(Nnodes_glob))
      ALLOCATE (magnetic_psi_glob(Nnodes_glob))
      B_glob = phys%B

      magnetic_flux_glob = phys%magnetic_flux
      magnetic_psi_glob = phys%magnetic_psi

      IF (switch%rmp .OR. switch%ripple) THEN
         ALLOCATE (Bperturb_glob(Nnodes_glob, 3))
         Bperturb_glob = phys%Bperturb
         DEALLOCATE (phys%Bperturb); NULLIFY (phys%Bperturb)
         ALLOCATE (phys%Bperturb(Nnodes_local, 3))
      END IF
      IF (switch%OhmicSrc) THEN
         ALLOCATE (Jtor_glob(Nnodes_glob))
         Jtor_glob = phys%Jtor
         DEALLOCATE (phys%Jtor); NULLIFY (phys%Jtor)
         ALLOCATE (phys%Jtor(Nnodes_local))
      END IF

      DEALLOCATE (phys%B); NULLIFY (phys%B)
      DEALLOCATE (phys%magnetic_flux); NULLIFY (phys%magnetic_flux)
      DEALLOCATE (phys%magnetic_psi); NULLIFY (phys%magnetic_psi)
      ALLOCATE (phys%B(Nnodes_local, 3))
      ALLOCATE (phys%magnetic_flux(Nnodes_local))
      ALLOCATE (phys%magnetic_psi(Nnodes_local))

      DO i = 1, Nnodes_tor
         indices_loc = (i - 1)*Mesh%Nnodes + indices
         indices_glob = (i - 1)*Mesh%Nno_glob + Mesh%loc2glob_nodes(indices)
         phys%B(indices_loc, :) = B_glob(indices_glob, :)
         phys%magnetic_flux(indices_loc) = magnetic_flux_glob(indices_glob)
         phys%magnetic_psi(indices_loc) = magnetic_psi_glob(indices_glob)
      END DO

      IF (switch%rmp .OR. switch%ripple) THEN
         DO i = 1, Nnodes_tor
            indices_loc = (i - 1)*Mesh%Nnodes + indices
            indices_glob = (i - 1)*Mesh%Nno_glob + Mesh%loc2glob_nodes(indices)
            phys%Bperturb(indices_loc, :) = Bperturb_glob(indices_glob, :)
         END DO
      END IF

      IF (switch%OhmicSrc) THEN
         DO i = 1, Nnodes_tor
            indices_loc = (i - 1)*Mesh%Nnodes + indices
            indices_glob = (i - 1)*Mesh%Nno_glob + Mesh%loc2glob_nodes(indices)
            phys%Jtor(indices_loc) = Jtor_glob(indices_glob)
         END DO
      END IF

      DEALLOCATE (indices)
      DEALLOCATE (indices_loc)
      DEALLOCATE (indices_glob)
      DEALLOCATE (B_glob)
      DEALLOCATE (magnetic_flux_glob)
      DEALLOCATE (magnetic_psi_glob)
      IF (ALLOCATED(Bperturb_glob)) THEN
         DEALLOCATE (Bperturb_glob)
      END IF
      IF (ALLOCATED(Jtor_glob)) THEN
         DEALLOCATE (Jtor_glob)
      END IF

   END SUBROUTINE Bfield_Jtor_decomposition

#ifdef TOR3D
   SUBROUTINE solution_decomposition()
      REAL*8, POINTER           :: u_glob(:), u_tilde_glob(:), q_glob(:)
      REAL*8, ALLOCATABLE         :: u_glob_2d(:,:), u_loc_2d(:,:), u_loc(:), q_glob_2d(:,:), q_loc_2d(:,:), q_loc(:), utilde_glob_2d_pol(:,:), utilde_glob_2d_tor(:,:), utilde_loc_2d_pol(:,:), utilde_loc_2d_tor(:,:), utilde_loc(:)
      INTEGER, ALLOCATABLE      :: indices_u(:, :), indices_q(:, :), indices_utilde_pol(:, :), indices_utilde_tor(:, :)
      INTEGER                     :: Nelems_global, Nfaces_global, Nelems_local, Nfaces_local, Ntor, Nnodesperelem, Nnodesperface, Nnodesperelem_tor1d, Neq, size_utilde
      INTEGER                     :: iel, ifa, itor, index, index_1_u, index_2_u, index_start_u, index_end_u, index_1_q, index_2_q, index_start_q, index_end_q, index_start_utilde, index_end_utilde,index_1_utilde, index_2_utilde

      NULLIFY (u_glob, u_tilde_glob, q_glob)

      IF (mpivar%glob_size .EQ. 1) RETURN

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "*            SOLUTION DECOMPOSITION             *"
         WRITE (*, *) "*************************************************"
      END IF

      u_glob => sol%u
      u_tilde_glob => sol%u_tilde
      q_glob => sol%q

      Nelems_global = Mesh%Nel_glob
      Nfaces_global = Mesh%Nfa_glob
      Nelems_local = Mesh%Nelems
      Nfaces_local = Mesh%Nfaces
      Ntor = numer%ntor
      Nnodesperelem_tor1d = numer%ptor + 1
      Nnodesperelem = Mesh%Nnodesperelem
      Nnodesperface = Mesh%Nnodesperface
      Neq = phys%Neq
      ndim = 3

      ! size of utilde per equation local
      IF (MPIvar%ntor .GT. 1) THEN
         size_utilde = (Ntor*(Nfaces_local*Nnodesperface*Nnodesperelem_tor1d + Nelems_local*Nnodesperelem) + Nelems_local*Nnodesperelem)*neq
      ELSE
         size_utilde = Ntor*(Nfaces_local*Nnodesperface*Nnodesperelem_tor1d + Nelems_local*Nnodesperelem)*neq
      END IF

      ! allocate the indices array
      ALLOCATE (indices_u(Nelems_global, Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq))
      ALLOCATE (indices_q(Nelems_global, Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim))
      ALLOCATE (indices_utilde_pol(Nelems_global, Ntor*Nnodesperelem*Neq))
      ALLOCATE (indices_utilde_tor(Nfaces_global, Ntor*Nnodesperelem_tor1d*Nnodesperface*Neq))
      indices_u = 0
      indices_q = 0
      indices_utilde_pol = 0
      indices_utilde_tor = 0
      index_start_utilde = 0
      index_end_utilde = 0
      index_1_utilde = 0
      index_2_utilde = 0

      ! Loop over elements and toroidal slices
      DO iel = 1, Nelems_global
         DO itor = 1, ntor
            index_1_u = (itor - 1)*Nelems_global*Nnodesperelem_tor1d*Nnodesperelem*Neq + (iel - 1)*Nnodesperelem_tor1d*Nnodesperelem*Neq + 1
            index_2_u = index_1_u + Nnodesperelem_tor1d*Nnodesperelem*Neq - 1
            index_1_q = (itor - 1) * Nelems_global * Nnodesperelem_tor1d * Nnodesperelem * Neq * ndim + (iel - 1) * Nnodesperelem_tor1d * Nnodesperelem * Neq * ndim + 1
            index_2_q = index_1_q + Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim - 1
            index_1_utilde = (itor - 1) * (Nelems_global * Nnodesperelem + Nfaces_global * Nnodesperface * Nnodesperelem_tor1d) * Neq + (iel - 1)*Nnodesperelem*Neq + 1
            index_2_utilde = index_1_utilde + Nnodesperelem*Neq - 1

            index_start_u = (itor - 1)*Nnodesperelem_tor1d*Nnodesperelem*Neq + 1
            index_end_u = index_start_u + Nnodesperelem_tor1d*Nnodesperelem*Neq - 1
            index_start_q = (itor - 1)*Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim + 1
            index_end_q = index_start_q + Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim - 1
            index_start_utilde = (itor - 1)*Nnodesperelem*Neq + 1
            index_end_utilde = index_start_utilde + Nnodesperelem*Neq - 1

            indices_u(iel, index_start_u:index_end_u) = [(i, i=index_1_u, index_2_u)]
            indices_q(iel, index_start_q:index_end_q) = [(i, i=index_1_q, index_2_q)]
            indices_utilde_pol(iel, index_start_utilde:index_end_utilde) = [(i, i=index_1_utilde, index_2_utilde)]
         END DO
      END DO

      DO ifa = 1, Nfaces_global
         DO itor = 1, ntor
            index_1_utilde = (itor - 1) * (Nelems_global * Nnodesperelem + Nfaces_global * Nnodesperface * Nnodesperelem_tor1d) * Neq + (Nelems_global * Nnodesperelem + (ifa - 1) * Nnodesperface * Nnodesperelem_tor1d) * Neq + 1
            index_2_utilde = index_1_utilde + Nnodesperface*Nnodesperelem_tor1d*Neq - 1

            index_start_utilde = (itor - 1)*Nnodesperface*Nnodesperelem_tor1d*Neq + 1
            index_end_utilde = index_start_utilde + Nnodesperface*Nnodesperelem_tor1d*Neq - 1
            indices_utilde_tor(ifa, index_start_utilde:index_end_utilde) = [(i, i=index_1_utilde, index_2_utilde)]
         END DO
      END DO

      ALLOCATE (u_glob_2d(Nelems_global, Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq))
      ALLOCATE (q_glob_2d(Nelems_global, Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim))
      ALLOCATE (utilde_glob_2d_pol(Nelems_global, Ntor*Nnodesperelem*Neq))
      ALLOCATE (utilde_glob_2d_tor(Nfaces_global, Ntor*Nnodesperelem_tor1d*Nnodesperface*Neq))
      u_glob_2d = 0.
      q_glob_2d = 0.
      utilde_glob_2d_pol = 0.
      utilde_glob_2d_tor = 0.

      DO iel = 1, Nelems_global
         u_glob_2d(iel, :) = u_glob(indices_u(iel, :))
         q_glob_2d(iel, :) = q_glob(indices_q(iel, :))
         utilde_glob_2d_pol(iel, :) = u_tilde_glob(indices_utilde_pol(iel, :))
      END DO

      DO ifa = 1, Nfaces_global
         utilde_glob_2d_tor(ifa, :) = u_tilde_glob(indices_utilde_tor(ifa, :))
      END DO

      ALLOCATE (u_loc_2d(Nelems_local, Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq))
      ALLOCATE (q_loc_2d(Nelems_local, Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim))
      ALLOCATE (utilde_loc_2d_pol(Nelems_local, Ntor*Nnodesperelem*Neq))
      ALLOCATE (utilde_loc_2d_tor(Nfaces_local, Ntor*Nnodesperelem_tor1d*Nnodesperface*Neq))
      u_loc_2d = 0.
      q_loc_2d = 0.
      utilde_loc_2d_pol = 0.
      utilde_loc_2d_tor = 0.

      DO i = 1, Nelems_local
         index = Mesh%loc2glob_el(i)
         u_loc_2d(i, :) = u_glob_2d(index, :)
         q_loc_2d(i, :) = q_glob_2d(index, :)
         utilde_loc_2d_pol(i, :) = utilde_glob_2d_pol(index, :)
      END DO
      DO i = 1, Nfaces_local
         index = Mesh%loc2glob_fa(i)
         utilde_loc_2d_tor(i, :) = utilde_glob_2d_tor(index, :)
      END DO

      ALLOCATE (u_loc(Nelems_local*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq))
      ALLOCATE (q_loc(Nelems_local*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim))
      ALLOCATE (utilde_loc(size_utilde))
      u_loc = 0.
      q_loc = 0.
      utilde_loc = 0.

      DO itor = 1, ntor
         index_start_u = (itor - 1)*Nnodesperelem_tor1d*Nnodesperelem*Neq + 1
         index_end_u = index_start_u + Nnodesperelem_tor1d*Nnodesperelem*Neq - 1
         index_start_q = (itor - 1)*Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim + 1
         index_end_q = index_start_q + Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim - 1
         index_start_utilde = (itor - 1)*Nnodesperelem*Neq + 1
         index_end_utilde = index_start_utilde + Nnodesperelem*Neq - 1

         DO iel = 1, Nelems_local
            index_1_u = (itor - 1)*Nelems_local*Nnodesperelem*Nnodesperelem_tor1d*Neq + (iel - 1)*Nnodesperelem_tor1d*Nnodesperelem*Neq + 1
            index_2_u = index_1_u + Nnodesperelem_tor1d*Nnodesperelem*Neq - 1
            index_1_q = (itor - 1)*Nelems_local*Nnodesperelem *Nnodesperelem_tor1d *  Neq * ndim + (iel - 1 )*Nnodesperelem_tor1d * Nnodesperelem * Neq * ndim + 1
            index_2_q = index_1_q + Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim - 1
            index_1_utilde = (itor - 1)*(Nelems_local*Nnodesperelem + Nfaces_local*Nnodesperface*Nnodesperelem_tor1d) * Neq + (iel - 1) * Nnodesperelem * Neq + 1
            index_2_utilde = index_1_utilde + Nnodesperelem*Neq - 1

            u_loc(index_1_u:index_2_u) = u_loc_2d(iel, index_start_u:index_end_u)
            q_loc(index_1_q:index_2_q) = q_loc_2d(iel, index_start_q:index_end_q)
            utilde_loc(index_1_utilde:index_2_utilde) = utilde_loc_2d_pol(iel, index_start_utilde:index_end_utilde)
         END DO
      END DO

      DO itor = 1, Ntor
         index_start_utilde = (itor - 1)*Nnodesperface*Nnodesperelem_tor1d*Neq + 1
         index_end_utilde = index_start_utilde + Nnodesperface*Nnodesperelem_tor1d*Neq - 1
         DO ifa = 1, Nfaces_local
            index_1_utilde = Nelems_local*Nnodesperelem*Neq + (itor - 1)*(Nelems_local*Nnodesperelem + Nfaces_local*Nnodesperface*Nnodesperelem_tor1d) * Neq + (ifa - 1) * Nnodesperface*Nnodesperelem_tor1d * Neq + 1
            index_2_utilde = index_1_utilde + Nnodesperface*Nnodesperelem_tor1d*Neq - 1
            utilde_loc(index_1_utilde:index_2_utilde) = utilde_loc_2d_tor(ifa, index_start_utilde:index_end_utilde)
         END DO
      END DO

      DEALLOCATE (sol%u)
      DEALLOCATE (sol%q)
      DEALLOCATE (sol%u_tilde)
      NULLIFY (sol%u)
      NULLIFY (sol%q)
      NULLIFY (sol%u_tilde)

      ALLOCATE (sol%u(Nelems_local*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq))
      ALLOCATE (sol%q(Nelems_local*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq*ndim))
      ALLOCATE (sol%u_tilde(size_utilde))
      sol%u = u_loc
      sol%q = q_loc
      sol%u_tilde = utilde_loc

      DEALLOCATE (indices_u, indices_q, indices_utilde_pol, indices_utilde_tor)
      DEALLOCATE (u_glob_2d, u_loc_2d, u_loc)
      DEALLOCATE (q_glob_2d, q_loc_2d, q_loc)
      DEALLOCATE (utilde_glob_2d_pol, utilde_glob_2d_tor, utilde_loc_2d_pol, utilde_loc_2d_tor, utilde_loc)
      NULLIFY (u_glob, u_tilde_glob, q_glob)

   END SUBROUTINE solution_decomposition
#else
   SUBROUTINE solution_decomposition()
      REAL*8, POINTER           :: u_glob(:), u_tilde_glob(:), q_glob(:)
      REAL*8, ALLOCATABLE       :: u_3d(:, :, :), u_tilde_3d(:, :, :), q_4D(:, :, :, :)
      REAL*8, ALLOCATABLE       :: u_3d_local(:, :, :), u_tilde_3d_local(:, :, :), q_4d_local(:, :, :, :)
      INTEGER                   :: i, index

      NULLIFY (u_glob, u_tilde_glob, q_glob)

      IF (mpivar%glob_size .EQ. 1) RETURN

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "           SOLUTION DECOMPOSITION                "
         WRITE (*, *) "*************************************************"
      END IF

      u_glob => sol%u
      u_tilde_glob => sol%u_tilde
      q_glob => sol%q

      ALLOCATE (u_3d(Mesh%Nel_glob, Mesh%Nnodesperelem, phys%neq))
      ALLOCATE (u_tilde_3d(Mesh%Nfa_glob, Mesh%Nnodesperface, phys%neq))
      ALLOCATE (q_4d(Mesh%Nel_glob, Mesh%Nnodesperelem, phys%neq, Mesh%Ndim))

      u_tilde_3d = 0.
      u_3d = 0.
      q_4D = 0.

      CALL reshape_transpose_permute(u_glob, u_3d, phys%neq, Mesh%Nel_glob, Mesh%Nnodesperelem)
      CALL reshape_transpose_permute(u_tilde_glob, u_tilde_3d, phys%neq, Mesh%Nfa_glob, Mesh%Nnodesperface)
      CALL reshape_transpose_permute_4D(q_glob, q_4D, Mesh%Ndim, phys%neq, Mesh%Nel_glob, Mesh%Nnodesperelem)

      ALLOCATE (u_3d_local(Mesh%Nelems, Mesh%Nnodesperelem, phys%neq))
      ALLOCATE (u_tilde_3d_local(Mesh%Nfaces, Mesh%Nnodesperface, phys%neq))
      ALLOCATE (q_4d_local(Mesh%Nelems, Mesh%Nnodesperelem, phys%neq, Mesh%Ndim))
      u_3d_local = 0.
      u_tilde_3d_local = 0.
      q_4d_local = 0.

      DO i = 1, Mesh%Nelems
         index = Mesh%loc2glob_el(i)
         u_3d_local(i, :, :) = u_3d(index, :, :)
         q_4d_local(i, :, :, :) = q_4D(index, :, :, :)
      END DO

      DO i = 1, Mesh%Nfaces
         u_tilde_3d_local(i, :, :) = u_tilde_3d(Mesh%loc2glob_fa(i), :, :)
      END DO

      DEALLOCATE (sol%u, sol%u_tilde, sol%q)

      ALLOCATE (sol%u(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq))
      ALLOCATE (sol%u_tilde(Mesh%Nfaces*Mesh%Nnodesperface*phys%neq))
      ALLOCATE (sol%q(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq*Mesh%Ndim))
      sol%u = 0
      sol%u_tilde = 0
      sol%q = 0

      CALL flatten_row_major(u_3d_local, sol%u, SIZE(u_3d_local, 1), SIZE(u_3d_local, 2), SIZE(u_3d_local, 3))
      CALL flatten_row_major(u_tilde_3d_local, sol%u_tilde, SIZE(u_tilde_3d_local,1), SIZE(u_tilde_3d_local,2), SIZE(u_tilde_3d_local,3))
      CALL flatten_row_major_4D(q_4d_local, sol%q, SIZE(q_4d_local, 1), SIZE(q_4d_local, 2), SIZE(q_4d_local, 3), SIZE(q_4d_local, 4))

      DEALLOCATE (u_3d, u_tilde_3d, q_4d)
      DEALLOCATE (u_3d_local, u_tilde_3d_local, q_4d_local)
      NULLIFY (u_glob, u_tilde_glob, q_glob)
   END SUBROUTINE solution_decomposition
#endif
! TOR3D
#endif
! PARALL


   SUBROUTINE split_mesh(ndiv, divtype, verbose)

      INTEGER, INTENT(IN)                     :: nDiv, divType
      LOGICAL, INTENT(IN)                     :: verbose

      TYPE(DD_mesh)                           :: M(nDiv)

      INTEGER, POINTER, DIMENSION(:, :)        :: intfaces_glob => NULL(), extfaces_glob => NULL(), T => NULL(), TTb => NULL()
      INTEGER, POINTER, DIMENSION(:)          :: boundaryFlag => NULL()
      REAL*8, POINTER, DIMENSION(:, :)         :: X => NULL()

      INTEGER                                 :: Nvert, n_ext_faces, n_int_faces, Nelems, Nfaces, Nnodes, Nnodes_mean, Nelems_mean, Nel_tot
      INTEGER                                 :: Nel_loc, n_boundaries
      INTEGER                                 :: iel, it, i, j, ii, jj, counter, clg, offset, ifa, iface, ib, n, ierr, b1, b2, n1, n2, s
      INTEGER, PARAMETER                      :: imax = 1000
      INTEGER                                 :: Efaces(refElPol%Nfaces, refElPol%Nfacenodeslin)
      INTEGER, ALLOCATABLE, DIMENSION(:)      :: indgh, innext, ghel, ghel_temp, invmap, loc2glob_fa, loc2glob_el, reptFaces, ghostFaces, faceglo, nodes, unique_boundaryFlag, Tb_prov_unique, mapextf_unique, checkFaces
      INTEGER, ALLOCATABLE, DIMENSION(:)      :: nodesfix, mapTfix, Tb_sizes, mapextf, flipFace, indices, indices2, auxext, auxflp, boundaryFlag_loc
      INTEGER, ALLOCATABLE, DIMENSION(:,:)    :: auxgh, Tprov, Tfix, Tp, Tbprov_temp, intfaces, extfaces, Tb_bound, Tb_prov, T_CON, Tb_loc, loc2glob_faces
      REAL*8, ALLOCATABLE, DIMENSION(:, :)     :: Xprov
      LOGICAL, ALLOCATABLE, DIMENSION(:)      :: ind
      REAL*8                                  :: xcenter, ycenter, step, theta0, theta1, angle, tol = 1e-6, pi_loc
      REAL*8, ALLOCATABLE                     :: xb(:), yb(:), thetab(:), theta(:)
      LOGICAL                                 :: redo, check
      REAL*8                                  :: RotMat(2, 2), rot_angle

      NULLIFY (intfaces_glob, extfaces_glob, T, TTb, boundaryFlag, X)

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "*            DOMAIN DECOMPOSITION               *"
         WRITE (*, *) "*************************************************"
      END IF

      pi_loc = 4.D0*DATAN(1.D0)
      rot_angle = 1.

      Nvert = refElPol%Nvertices
      n_ext_faces = Mesh%Nextfaces
      n_int_faces = Mesh%Nintfaces

      intfaces_glob => Mesh%intfaces
      extfaces_glob => Mesh%extfaces

      Nelems = Mesh%Nelems
      Nfaces = Mesh%Nfaces
      Nnodes = Mesh%Nnodes

      Nnodes_mean = FLOOR(REAL(Nnodes/ndiv))
      Nelems_mean = FLOOR(REAL(Nelems/ndiv))

      ! rescale back to physical values
      Mesh%X = Mesh%X*phys%lscale

      X => Mesh%X
      T => Mesh%T

      ! Apply shift back if axisymmetric case
      IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) .OR. (switch%axisym .AND. MINVAL(Mesh%X(:,1)) < tol)) THEN
         X(:, 1) = X(:, 1) - geom%R0
      END IF

      boundaryFlag => Mesh%boundaryFlag
      TTb => Mesh%Tb

      CALL UNIQUE_STABLE(boundaryFlag, unique_boundaryFlag)
      n_boundaries = SIZE(unique_boundaryFlag)

      ALLOCATE (Tb_sizes(n_boundaries))
      Tb_sizes = 0
      counter = 1

      DO i = 1, n_boundaries
         Tb_sizes(i) = COUNT(boundaryFlag .EQ. unique_boundaryFlag(i))
      END DO

      angle = 0
      redo = .TRUE.

      xcenter = 0.5*(MAXVAL(X(:, 1)) + MINVAL(X(:, 1))) ! 3.4
      ycenter = 0.5*(MAXVAL(X(:, 2)) + MINVAL(X(:, 2))) ! 0

      ALLOCATE (xb(Nelems))
      ALLOCATE (yb(Nelems))
      ALLOCATE (thetab(Nelems))
      ALLOCATE (theta(nDiv + 1))
      ALLOCATE (auxgh(Nfaces, 2))
      ALLOCATE (indgh(Nfaces))
      ALLOCATE (ind(Nelems))

      xb = 0.
      yb = 0.
      thetab = 0.
      theta = 0.
      auxgh = 0
      indgh = 0
      ind = .FALSE.

      ! keep on rotating the mesh if no good partition is found
      DO WHILE (redo)

         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

         IF (MPIvar%glob_id .EQ. 0) THEN
            WRITE (*, *) "Angle of rotation: ", angle
         END IF

         IF (angle .GT. 360.d0) THEN
            IF (MPIvar%glob_id .EQ. 0) THEN
               WRITE (*, *) "No good partition in domain decomposition found. Stopping."
            END IF
            CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr)
         END IF

         !*********************** Mesh's info ************************!

         xb = 0.
         yb = 0.
         thetab = 0.
         theta = 0.
         auxgh = 0
         indgh = 0
         ind = .FALSE.

         DO iel = 1, Nelems
            xb(iel) = SUM(X(T(iel, 1:Nvert), 1))/REAL(Nvert)
            yb(iel) = SUM(X(T(iel, 1:Nvert), 2))/REAL(Nvert)
         END DO

         SELECT CASE (divType)
          CASE (1)
            ! Radial division
            DO iel = 1, Nelems
               thetab(iel) = ATAN(yb(iel), xb(iel))
               IF (thetab(iel) .LT. 0) THEN
                  thetab(iel) = thetab(iel) + 2*pi_loc
               END IF
            END DO
            theta = [((i - 1)*(2*pi_loc/REAL(nDiv)), i=1, nDiv + 1)]
          CASE (2)
            ! X division
            thetab = xb
            step = (MAXVAL(X(:, 1)) - MINVAL(X(:, 1)))/REAL(nDiv)
            theta = [(MINVAL(X(:, 1)) + (i - 1)*step, i=1, nDiv + 1)]
          CASE (3)
            ! Y division
            thetab = yb
            step = (MAXVAL(X(:, 2)) - MINVAL(X(:, 2)))/REAL(nDiv)
            theta = [(MINVAL(X(:, 2)) + (i - 1)*step, i=1, nDiv + 1)]
          CASE DEFAULT
            WRITE (*, *) "Option not supported for divType, domain decompositon. stop."
            STOP
         END SELECT

         !******************** Balancing Division *********************!

         ! !***** Division Loop->ndiv-1 *****!
         !

         DO WHILE (redo)

            IF (verbose) THEN
               IF (MPIvar%glob_id .EQ. 0) THEN
                  WRITE (*, *) "Balancing divisions"
               END IF
            END IF

            redo = .FALSE.
            Nel_tot = 0
            Nel_loc = 0

            ! Balance equation
            DO it = 1, nDiv - 1
               theta0 = theta(it)
               theta1 = theta(it + 2)

               DO i = 1, imax
                  Nel_loc = 0
                  DO j = 1, SIZE(thetab)
                     IF ((thetab(j) .GE. theta(it)) .AND. (thetab(j) .LE. theta(it + 1))) THEN
                        Nel_loc = Nel_loc + 1
                     END IF
                  END DO

                  IF ((Nel_loc - Nelems_mean) .GT. 5) THEN
                     theta1 = theta(it + 1)
                     theta(it + 1) = 0.5*(theta0 + theta1)
                  ELSEIF ((Nel_loc - Nelems_mean) .LT. -5) THEN
                     theta0 = theta(it + 1)
                     theta(it + 1) = 0.5*(theta0 + theta1)
                  ELSE
                     ! partitions are well balanced
                     EXIT
                  END IF
               END DO

               IF (i .EQ. imax + 1) THEN
                  redo = .TRUE.
               END IF
               Nel_tot = Nel_tot + Nel_loc
            END DO
         END DO

         Nel_tot = 0
         DO it = 1, SIZE(theta) - 1
            Nel_loc = 0
            DO i = 1, SIZE(thetab)
               IF ((thetab(i) .GE. theta(it)) .AND. (thetab(i) .LE. theta(it + 1))) THEN
                  Nel_loc = Nel_loc + 1
               END IF
            END DO
            IF (verbose) THEN
               IF (MPIvar%glob_id .EQ. 0) THEN
                  WRITE (*, *) "Number of elements per division:", Nelems_mean, Nel_loc
               END IF
            END IF

            Nel_tot = Nel_tot + Nel_loc
         END DO

         IF (Nel_tot .NE. Nelems) THEN
            IF (verbose) THEN
               IF (MPIvar%glob_id .EQ. 0) THEN
                  WRITE (*, *) 'Error in numbering of elements: ', Nel_tot, Nelems
               END IF
            END IF
            redo = .TRUE.
         END IF

#ifdef PARALL
         ! use MPI_Allreduce to check if any process encountered the exit condition
         CALL MPI_Allreduce(MPI_IN_PLACE, redo, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif

         IF (redo) THEN
            EXIT
         END IF

         CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

#ifndef PARALL
         DO it = 1, nDiv
#else
         DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1
#endif
            IF (verbose) THEN
               WRITE (*, *) "Division: ", it
            END IF

            ALLOCATE (loc2glob_el(COUNT((thetab .GE. theta(it)) .AND. (thetab .LE. theta(it + 1)))))
            loc2glob_el = 0
            Tb_sizes = 0

            ! set the index of the global element in the elemental loc2glob
            ind = .FALSE.

            counter = 1
            DO i = 1, SIZE(thetab)
               IF ((thetab(i) .GE. theta(it)) .AND. (thetab(i) .LE. theta(it + 1))) THEN
                  loc2glob_el(counter) = i
                  counter = counter + 1
                  ind(i) = .TRUE.
               END IF
            END DO

            IF (verbose) THEN
               WRITE (*, *) "Add ghost element to the division"
            END IF

            IF (nDiv .GT. 1) THEN

               auxgh = 0
               DO i = 1, n_int_faces
                  IF (ANY(intFaces_Glob(i, 1) .EQ. loc2glob_el)) auxgh(i, 1) = 1
                  IF (ANY(intFaces_Glob(i, 3) .EQ. loc2glob_el)) auxgh(i, 2) = 1
               END DO

               indgh = 0
               DO i = 1, Nfaces
                  IF (SUM(auxgh(i, :)) .EQ. 1) indgh(i) = 1
               END DO

               ALLOCATE (ghel_temp(2*SIZE(intFaces_Glob)))
               ghel_temp = -1
               counter = 1
               DO i = 1, SIZE(indgh)
                  IF (indgh(i) .NE. 0) THEN
                     IF (intFaces_Glob(i, 1) .NE. intFaces_Glob(i, 3)) THEN
                        ghel_temp(counter) = intFaces_Glob(i, 1)
                        ghel_temp(counter + 1) = intFaces_Glob(i, 3)
                        counter = counter + 2
                     ELSE
                        ghel_temp(counter) = intFaces_Glob(i, 1)
                        counter = counter + 1
                     END IF
                  END IF
               END DO

               DO i = 1, SIZE(ghel_temp)
                  IF (ghel_temp(i) .EQ. -1) ghel_temp(i) = ghel_temp(1)
               END DO

               CALL UNIQUE_1D(ghel_temp, ghel)
               DEALLOCATE (ghel_temp)

               counter = 1
               DO i = 1, SIZE(ghel)
                  IF (thetab(ghel(i)) .GT. theta(it + 1)) THEN
                     counter = counter + 1
                  END IF
               END DO

               ALLOCATE (innext(counter - 1))
               innext = 0
               counter = 1
               DO i = 1, SIZE(ghel)
                  IF (thetab(ghel(i)) .GE. theta(it + 1)) THEN
                     innext(counter) = i
                     counter = counter + 1
                  END IF
               END DO

               counter = 1
               DO i = 1, SIZE(innext)
                  ind(ghel(innext(i))) = .TRUE.
               END DO

               DEALLOCATE (loc2glob_el)
               ALLOCATE (loc2glob_el(COUNT(ind .EQV. .TRUE.)))
               loc2glob_el = 0

               counter = 1
               DO i = 1, SIZE(ind)
                  IF (ind(i)) THEN
                     loc2glob_el(counter) = i
                     counter = counter + 1
                  END IF
               END DO

               ALLOCATE (M(it)%ghostElems(SIZE(loc2glob_el)))
               M(it)%ghostElems = 0  ! Initialize all elements to false

               DO i = 1, SIZE(loc2glob_el)
                  DO j = 1, SIZE(innext)
                     IF (loc2glob_el(i) .EQ. ghel(innext(j))) THEN
                        M(it)%GhostElems(i) = 1  ! Set to true if loc2glob_el(i) is in ghel(innext)
                        EXIT
                     END IF
                  END DO
               END DO

               DEALLOCATE (ghel)
            ELSE

               ALLOCATE (M(it)%ghostElems(Nelems))
               M(it)%ghostElems = 0

            END IF

            IF (verbose) THEN
               IF (MPIvar%glob_id .EQ. 0) THEN
                  WRITE (*, *) "Generate X and T for the current division"
               END IF
            END IF

            ALLOCATE (Tprov(COUNT(ind .EQV. .TRUE.), Mesh%Nnodesperelem))
            ALLOCATE (Tfix(SIZE(Tprov, 1), Mesh%Nnodesperelem))
            Tprov = 0
            Tfix = 0

            counter = 1
            DO i = 1, SIZE(ind)
               IF (ind(i)) THEN
                  Tprov(counter, :) = T(i, :)
                  counter = counter + 1
               END IF
            END DO

            CALL UNIQUE_1D(RESHAPE(Tprov, [SIZE(Tprov, 1)*SIZE(Tprov, 2)]), nodes)

            ALLOCATE (Xprov(SIZE(nodes), Mesh%ndim))
            Xprov = X(nodes, :)

            ALLOCATE (mapTfix(SIZE(nodes)))
            mapTfix = nodes

            ALLOCATE (invmap(mapTfix(SIZE(mapTfix))))
            invmap = 0
            DO i = 1, SIZE(mapTfix)
               invmap(mapTfix(i)) = i
            END DO

            DO i = 1, SIZE(Tprov, 1)
               DO j = 1, SIZE(Tprov, 2)
                  Tfix(i, j) = invmap(Tprov(i, j))
               END DO
            END DO

            ! nodesfix is allocated in here
            CALL UNIQUE_1D(RESHAPE(Tfix, [SIZE(Tfix, 1)*SIZE(Tfix, 2)]), nodesfix)

            ! check that everything is ok so far...
            IF (ALL(nodesfix .NE. (/(j, j=1, MAXVAL(Tfix))/))) THEN
               IF (verbose) THEN
                  WRITE (*, *) 'Problem in fixed connectivity. Retrying.'
               END IF
               redo = .TRUE.
               !exit
            END IF

#ifdef PARALL
            ! use MPI_Allreduce to check if any process encountered the exit condition
            CALL MPI_Allreduce(MPI_IN_PLACE, redo, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            IF (redo) THEN
               EXIT
            END IF

            IF (verbose) THEN
               IF (MPIvar%glob_id .EQ. 0) THEN
                  WRITE (*, *) "Create loc2glob_fa: interior faces"
               END IF
            END IF

            ! Get faces takes ages!!
            CALL GetFaces_mod(Tfix(:, 1:Nvert), intfaces, extFaces)

            ALLOCATE (M(it)%intfaces_loc(SIZE(intfaces, 1), SIZE(intfaces, 2)))
            ALLOCATE (M(it)%extfaces_loc(SIZE(extfaces, 1), SIZE(extfaces, 2)))

            M(it)%intfaces_loc = intfaces
            M(it)%extfaces_loc = extfaces
            M(it)%Nfaces = SIZE(intfaces, 1) + SIZE(extfaces, 1)
            M(it)%Nintfaces = SIZE(intfaces, 1)
            M(it)%Nextfaces = SIZE(extfaces, 1)
            M(it)%Nnodesperface = refElPol%nDeg + 1
            M(it)%elemType = Mesh%elemType

            ALLOCATE (loc2glob_fa(M(it)%Nfaces))
            ALLOCATE (reptFaces(M(it)%Nfaces))
            ALLOCATE (faceglo(5))

            loc2glob_fa = 0
            reptFaces = 0
            faceglo = 0

            DO j = 1, SIZE(intfaces, 1)

               faceglo(1) = loc2glob_el(intfaces(j, 1))
               faceglo(2) = intfaces(j, 2)
               faceglo(3) = loc2glob_el(intfaces(j, 3))
               faceglo(4) = intfaces(j, 4)
               faceglo(5) = intfaces(j, 5)

               DO jj = 1, SIZE(intFaces_Glob, 1)
                  IF (SUM(ABS(intFaces_Glob(jj, :) - faceglo)) .EQ. 0) THEN
                     loc2glob_fa(j) = jj
                  END IF
               END DO

               IF (loc2glob_fa(j) .EQ. 0) THEN
                  WRITE (*, *) "Problem in the creation of loc2glob_fa. stop"
                  STOP
               END IF
            END DO
            IF (verbose) THEN
               IF (MPIvar%glob_id .EQ. 0) THEN
                  WRITE (*, *) "Generate boundary connectivity for the current division."
               END IF
            END IF

            offset = 0
            ALLOCATE (Tb_loc(SIZE(TTb, 1), SIZE(TTb, 2)))
            ALLOCATE (boundaryFlag_loc(SIZE(TTb, 1)))
            boundaryFlag_loc = 0
            Tb_sizes = 0
            Tb_loc = 0
            clg = 1

            DO ib = 1, n_boundaries
               n = COUNT(boundaryFlag .EQ. unique_boundaryFlag(ib))
               ALLOCATE (Tb_bound(n, Mesh%Nnodesperface))
               Tb_bound = 0

               counter = 1
               DO i = 1, SIZE(boundaryFlag)
                  IF (boundaryFlag(i) .EQ. unique_boundaryFlag(ib)) THEN
                     Tb_bound(counter, :) = TTb(i, :)
                     counter = counter + 1
                  END IF
               END DO

               ALLOCATE (Tbprov_temp(n, Mesh%Nnodesperface))
               Tbprov_temp = 0

               counter = 1

               DO ifa = 1, n
                  IF (ANY(Tprov .EQ. Tb_bound(ifa, 1)) .AND. ANY(Tprov .EQ. Tb_bound(ifa, SIZE(Tb_bound, 2)))) THEN
                     Tbprov_temp(counter, :) = Tb_bound(ifa, :)
                     counter = counter + 1

                     DO j = 1, SIZE(TTb, 1)
                        IF (ALL(TTb(j, :) .EQ. Tb_bound(ifa, :))) THEN
                           loc2glob_fa(clg + M(it)%Nintfaces) = j + n_int_faces
                           EXIT
                        END IF

                        IF (j .EQ. SIZE(TTb, 1)) THEN
                           WRITE (*, *) "Problem in the creation of loc2glob_fa. stop."
                           STOP
                        END IF
                     END DO

                     clg = clg + 1
                  END IF
               END DO

               Tb_sizes(ib) = counter - 1

               IF (counter .NE. 1) THEN
                  ALLOCATE (Tb_prov(counter - 1, Mesh%Nnodesperface))
                  Tb_prov = Tbprov_temp(1:counter - 1, :)

                  ! Tb_prov_unique is allocated here
                  CALL UNIQUE_1D(RESHAPE(Tb_prov, [SIZE(Tb_prov, 1)*SIZE(Tb_prov, 2)]), Tb_prov_unique)

                  DO i = 1, SIZE(Tb_prov_unique)
                     DO j = 1, SIZE(nodes)
                        IF (nodes(j) .EQ. Tb_prov_unique(i)) THEN
                           DO ii = 1, SIZE(Tb_prov, 1)
                              DO jj = 1, SIZE(Tb_prov, 2)
                                 IF (Tb_prov(ii, jj) .EQ. Tb_prov_unique(i)) THEN
                                    Tb_prov(ii, jj) = j
                                 END IF
                              END DO
                           END DO
                           EXIT
                        END IF
                     END DO

                     IF (j .EQ. SIZE(nodes) + 1) THEN
                        IF (verbose) THEN
                           WRITE (*, *) "Line was cut too close to boundary. Retrying."
                        END IF
                        redo = .TRUE.
                        EXIT
                     END IF
                  END DO

                  DEALLOCATE (Tb_prov_unique)

                  IF (redo) THEN
                     EXIT
                  END IF

                  Tb_loc(offset + 1:offset + SIZE(Tb_prov, 1), :) = Tb_prov
                  boundaryFlag_loc(offset + 1:offset + SIZE(Tb_prov, 1)) = unique_boundaryFlag(ib)
                  offset = offset + SIZE(Tb_prov, 1)
                  Tb_sizes(ib) = SIZE(Tb_prov, 1)
                  DEALLOCATE (Tb_prov)

               END IF

               DEALLOCATE (Tb_bound)
               DEALLOCATE (Tbprov_temp)

            END DO

#ifdef PARALL
            ! use MPI_Allreduce to check if any process encountered the exit condition
            CALL MPI_Allreduce(MPI_IN_PLACE, redo, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            ! if redo rotate the mesh
            IF (redo) THEN
               X(:, 1) = X(:, 1) - xcenter
               X(:, 2) = X(:, 2) - ycenter
               RotMat(1, 1) = COS(rot_angle*pi_loc/180.)
               RotMat(1, 2) = -SIN(rot_angle*pi_loc/180.)
               RotMat(2, 1) = SIN(rot_angle*pi_loc/180.)
               RotMat(2, 2) = COS(rot_angle*pi_loc/180.)
               X = TRANSPOSE(MATMUL(RotMat, TRANSPOSE(X)))
               X(:, 1) = X(:, 1) + xcenter
               X(:, 2) = X(:, 2) + ycenter
               angle = angle + rot_angle
            END IF

            ! if redo rotate deallocate all allocated structures
            IF (redo) THEN

               DEALLOCATE (loc2glob_fa)
               DEALLOCATE (intfaces)
               DEALLOCATE (extfaces)
               DEALLOCATE (reptFaces)
               DEALLOCATE (Tfix)
               DEALLOCATE (Tprov)
               DEALLOCATE (Xprov)
               DEALLOCATE (nodes)
               DEALLOCATE (loc2glob_el)
               DEALLOCATE (innext)
               DEALLOCATE (mapTfix)
               DEALLOCATE (invmap)
               DEALLOCATE (nodesfix)
               DEALLOCATE (faceglo)

#ifndef PARALL
               DO i = 1, it
                  DEALLOCATE (M(i)%ghostElems)
                  DEALLOCATE (M(i)%intfaces_loc)
                  DEALLOCATE (M(i)%extfaces_loc)
               END DO
#else

               DEALLOCATE (M(it)%ghostElems)
               DEALLOCATE (M(it)%intfaces_loc)
               DEALLOCATE (M(it)%extfaces_loc)
#endif

#ifndef PARALL
               DO i = 1, it - 1
#else
               DO i = it, it
#endif
                  IF (ALLOCATED(M(i)%loc2glob_fa)) DEALLOCATE (M(i)%loc2glob_fa)
                  IF (ALLOCATED(M(i)%reptFaces)) DEALLOCATE (M(i)%reptFaces)
                  IF (ALLOCATED(M(i)%TT)) DEALLOCATE (M(i)%TT)
                  IF (ALLOCATED(M(i)%XX)) DEALLOCATE (M(i)%XX)
                  IF (ALLOCATED(M(i)%loc2glob_no)) DEALLOCATE (M(i)%loc2glob_no)
                  IF (ALLOCATED(M(i)%loc2glob_el)) DEALLOCATE (M(i)%loc2glob_el)
                  IF (ALLOCATED(M(i)%Tb)) DEALLOCATE (M(i)%Tb)
                  IF (ALLOCATED(M(i)%boundaryFlag)) DEALLOCATE (M(i)%boundaryFlag)
                  IF (ALLOCATED(M(i)%Tb_sizes)) DEALLOCATE (M(i)%Tb_sizes)
               END DO

               IF (ALLOCATED(Tb_bound)) DEALLOCATE (Tb_bound)
               IF (ALLOCATED(Tbprov_temp)) DEALLOCATE (Tbprov_temp)
               IF (ALLOCATED(Tb_prov)) DEALLOCATE (Tb_prov)
               IF (ALLOCATED(Tb_loc)) DEALLOCATE (Tb_loc)
               IF (ALLOCATED(boundaryFlag_loc)) DEALLOCATE (boundaryFlag_loc)
               EXIT
            END IF

            ! store subdivision
            ALLOCATE (M(it)%loc2glob_fa(SIZE(loc2glob_fa)))
            ALLOCATE (M(it)%reptFaces(SIZE(reptFaces)))
            ALLOCATE (M(it)%TT(SIZE(Tfix, 1), SIZE(Tfix, 2)))
            ALLOCATE (M(it)%XX(SIZE(Xprov, 1), SIZE(Xprov, 2)))
            ALLOCATE (M(it)%loc2glob_no(SIZE(nodes)))
            ALLOCATE (M(it)%loc2glob_el(SIZE(loc2glob_el)))
            ALLOCATE (M(it)%Tb(SUM(Tb_sizes), SIZE(Tb_loc, 2)))
            ALLOCATE (M(it)%boundaryFlag(SUM(Tb_sizes)))
            ALLOCATE (M(it)%Tb_sizes(SIZE(Tb_sizes)))

            M(it)%Loc2Glob_fa = loc2glob_fa
            M(it)%ReptFaces = reptFaces
            M(it)%TT = Tfix
            M(it)%XX = Xprov
            M(it)%loc2glob_no = nodes
            M(it)%loc2glob_el = loc2glob_el
            M(it)%intfaces_loc = intfaces
            M(it)%extfaces_loc = extfaces
            M(it)%Tb = Tb_loc(1:offset, :)
            M(it)%boundaryFlag = boundaryFlag_loc(1:offset)
            M(it)%Tb_sizes = Tb_sizes
            M(it)%Nfaces = SIZE(intfaces, 1) + SIZE(extfaces, 1)
            M(it)%Nintfaces = SIZE(intfaces, 1)
            M(it)%Ndim = SIZE(M(it)%XX, 2)
            M(it)%Nnodes = SIZE(M(it)%XX, 1)
            M(it)%Nelems = SIZE(M(it)%TT, 1)
            M(it)%Nnodesperelem = SIZE(M(it)%TT, 2)

            IF (ALLOCATED(loc2glob_fa)) DEALLOCATE (loc2glob_fa)
            IF (ALLOCATED(Tb_loc)) DEALLOCATE (Tb_loc)
            IF (ALLOCATED(intfaces)) DEALLOCATE (intfaces)
            IF (ALLOCATED(extfaces)) DEALLOCATE (extfaces)
            IF (ALLOCATED(reptFaces)) DEALLOCATE (reptFaces)
            IF (ALLOCATED(Tfix)) DEALLOCATE (Tfix)
            IF (ALLOCATED(Tprov)) DEALLOCATE (Tprov)
            IF (ALLOCATED(Xprov)) DEALLOCATE (Xprov)
            IF (ALLOCATED(nodes)) DEALLOCATE (nodes)
            IF (ALLOCATED(loc2glob_el)) DEALLOCATE (loc2glob_el)
            IF (ALLOCATED(innext)) DEALLOCATE (innext)
            IF (ALLOCATED(mapTfix)) DEALLOCATE (mapTfix)
            IF (ALLOCATED(invmap)) DEALLOCATE (invmap)
            IF (ALLOCATED(nodesfix)) DEALLOCATE (nodesfix)
            IF (ALLOCATED(faceglo)) DEALLOCATE (faceglo)
            IF (ALLOCATED(boundaryFlag_loc)) DEALLOCATE (boundaryFlag_loc)

         END DO
      END DO

      ! Boundary connectivity for the overlapping faces

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

#ifndef PARALL
      DO it = 1, nDiv
#else
      DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1
#endif
         IF (verbose) THEN
            WRITE (*, *) "Division: ", it
         END IF

         ALLOCATE (Tp(SIZE(M(it)%TT, 1), SIZE(M(it)%TT, 2)))
         ALLOCATE (loc2glob_fa(SIZE(M(it)%loc2glob_fa)))
         ALLOCATE (loc2glob_el(SIZE(M(it)%loc2glob_el)))
         ALLOCATE (GhostFaces(M(it)%Nfaces))
         ALLOCATE (extfaces(SIZE(M(it)%extfaces_loc, 1), SIZE(M(it)%extfaces_loc, 2)))

         Tp = M(it)%TT
         loc2glob_fa = M(it)%loc2glob_fa
         loc2glob_el = M(it)%loc2glob_el
         ghostFaces = 0
         extfaces = M(it)%extfaces_loc

         ALLOCATE (T_con(SIZE(Tp, 1), SIZE(refElPol%Face_nodes, 2)))
         ALLOCATE (flipFace(SIZE(Tp, 1)))
         ALLOCATE (mapextf(SIZE(extFaces, 1)))

         T_con = 0
         flipFace = 0
         mapextf = 0
         clg = 0

         DO j = 1, SIZE(loc2glob_fa)
            IF (loc2glob_fa(j) .EQ. 0) THEN
               clg = j
               EXIT
            END IF
         END DO

         counter = 1
         DO ifa = 1, SIZE(extFaces, 1)
            iel = extFaces(ifa, 1);
            iface = extFaces(ifa, 2);
            n1 = Tp(iel, Efaces(iface, 1));
            n2 = Tp(iel, Efaces(iface, 2));
            check = .TRUE.
            s = 0

            DO ib = 1, n_boundaries

               IF (M(it)%Tb_sizes(ib) .NE. 0) THEN
                  ALLOCATE (Tb_bound(M(it)%Tb_sizes(ib), Mesh%Nnodesperface))
                  Tb_bound = M(it)%Tb(s + 1:s + M(it)%Tb_sizes(ib), :)

                  IF (ANY(Tb_bound .EQ. n1) .AND. ANY(Tb_bound .EQ. n2)) THEN

                     DO j = 1, SIZE(Tb_bound, 1)
                        IF ((n1 .EQ. Tb_bound(j, 1)) .AND. (n2 .EQ. Tb_bound(j, SIZE(Tb_bound, 2)))) THEN
                           mapextf(ifa) = j + s
                           check = .FALSE.
                           EXIT
                        ELSEIF ((n2 .EQ. Tb_bound(j, 1)) .AND. (n1 .EQ. Tb_bound(j, SIZE(Tb_bound, 2)))) THEN
                           mapextf(ifa) = j + s
                           check = .FALSE.
                           EXIT
                        END IF

                        IF (j .EQ. (SIZE(Tb_bound, 1) + 1)) THEN
                           WRITE (*, *) "Strange."
                           STOP
                        END IF
                     END DO
                  END IF
                  s = s + M(it)%Tb_sizes(ib)
                  DEALLOCATE (Tb_bound)
               END IF
            END DO

            IF (check) THEN
               T_CON(counter, :) = Tp(iel, refElPol%Face_nodes(iface, :))
               mapextf(ifa) = counter + s
               ghostFaces(counter + s + M(it)%Nintfaces) = 1

               b1 = 0
               b2 = 0
               DO j = 1, SIZE(intFaces_Glob, 1)
                  IF ((loc2glob_el(iel) .EQ. intFaces_Glob(j, 1)) .AND. (iface .EQ. intFaces_Glob(j, 2))) THEN
                     b1 = j
                  END IF
                  IF ((loc2glob_el(iel) .EQ. intFaces_Glob(j, 3)) .AND. (iface .EQ. intFaces_Glob(j, 4))) THEN
                     b2 = j
                  END IF
               END DO

               IF (b1 .NE. 0) THEN
                  loc2glob_fa(clg) = b1
                  IF (b2 .NE. 0) THEN
                     WRITE (*, *) "That is strange. stop."
                     STOP
                  END IF
                  IF (loc2glob_el(iel) .GT. intFaces_Glob(b1, 3)) THEN
                     flipFace(counter) = 1
                  END IF
               ELSEIF (b2 .NE. 0) THEN
                  loc2glob_fa(clg) = b2
                  IF (loc2glob_el(iel) .GT. intFaces_Glob(b2, 1)) THEN
                     flipFace(counter) = 1
                  END IF
               ELSE
                  WRITE (*, *) "That is very strange. stop."
                  STOP
               END IF

               counter = counter + 1
               clg = clg + 1
            END IF
         END DO

         ! little check
         IF (ANY(mapextf .EQ. 0)) THEN
            WRITE (*, *) "Error in mapextf. stop."
            STOP
         END IF

         ! watch out
         CALL UNIQUE_1D(mapextf, mapextf_unique)
         IF (SIZE(mapextf_unique) .NE. MAXVAL(mapextf)) THEN
            WRITE (*, *) "Error in mapextf. stop."
            STOP
         END IF
         DO j = 1, SIZE(mapextf_unique)
            IF (mapextf_unique(j) .NE. j) THEN
               WRITE (*, *) "Error in mapextf. stop."
               STOP
            END IF
         END DO

         ALLOCATE (M(it)%TTb_CON(counter - 1, Mesh%Nnodesperface))
         ALLOCATE (M(it)%ghostFlp(counter - 1))
         ALLOCATE (M(it)%ghostFaces(M(it)%Nfaces))
         ALLOCATE (M(it)%mapextfa_loc(SIZE(mapextf)))

         M(it)%TTb_CON = T_CON(1:counter - 1, :)
         M(it)%ghostFlp = flipFace(1:counter - 1)
         M(it)%loc2glob_fa = loc2glob_fa
         M(it)%mapextfa_loc = mapextf
         M(it)%ghostFaces = ghostFaces

         DEALLOCATE (T_CON)
         DEALLOCATE (flipFace)
         DEALLOCATE (mapextf)
         DEALLOCATE (ghostFaces)
         DEALLOCATE (loc2glob_fa)
         DEALLOCATE (loc2glob_el)
         DEALLOCATE (mapextf_unique)
         DEALLOCATE (Tp)
         DEALLOCATE (extfaces)

      END DO

      !  Find multiple occurrence of the same face across various processes, this part needs to be run un serial
      IF (verbose) THEN
         IF (MPIvar%glob_id .EQ. 0) THEN
            WRITE (*, *) "Find multiple occurrence of the same face across various processes"
         END IF
      END IF

      ALLOCATE (checkFaces(Nfaces))
      checkFaces = 0

#ifndef PARALL
      DO it = 1, ndiv
#else

      DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1

         IF (MPIvar%glob_id .NE. 0) THEN
            CALL MPI_Recv(checkFaces, Nfaces, MPI_INTEGER, MPIvar%glob_id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         END IF

#endif
         IF (verbose) THEN
            WRITE (*, *) "Division: ", it
         END IF

         ALLOCATE (loc2glob_faces(M(it)%Nintfaces, 5))
         ALLOCATE (indices(n_int_faces))
         indices = 0
         loc2glob_faces = 0

         loc2glob_faces(:, 1) = M(it)%loc2glob_el(M(it)%intfaces_loc(:, 1))
         loc2glob_faces(:, 2) = M(it)%intfaces_loc(:, 2)
         loc2glob_faces(:, 3) = M(it)%loc2glob_el(M(it)%intfaces_loc(:, 3))
         loc2glob_faces(:, 4) = M(it)%intfaces_loc(:, 4)
         loc2glob_faces(:, 5) = M(it)%intfaces_loc(:, 5)

         ! interior faces
         DO i = 1, n_int_faces
            DO j = 1, M(it)%Nintfaces
               IF (ALL(intFaces_Glob(i, :) .EQ. loc2glob_faces(j, :))) THEN
                  indices(i) = j
                  EXIT
               END IF
            END DO
         END DO

         DO i = 1, SIZE(indices)
            IF (indices(i) .NE. 0) THEN
               M(it)%reptFaces(indices(i)) = checkFaces(i)
               checkFaces(i) = 1
            END IF
         END DO

         DEALLOCATE (loc2glob_faces)
         DEALLOCATE (indices)
         ALLOCATE (loc2glob_faces(M(it)%NextFaces, 2))
         ALLOCATE (indices(n_ext_faces))
         loc2glob_faces = 0
         indices = 0

         ! what is this 2??
         loc2glob_faces(:, 1) = M(it)%loc2glob_el(M(it)%extfaces_loc(:, 1))
         loc2glob_faces(:, 2) = M(it)%extfaces_loc(:, 2)

         ! exterior faces
         DO i = 1, n_ext_faces
            DO j = 1, M(it)%NextFaces
               IF (ALL(extFaces_Glob(i, :) .EQ. loc2glob_faces(j, :))) THEN
                  indices(i) = j
                  EXIT
               END IF
            END DO
         END DO

         DO i = 1, SIZE(indices)
            IF (indices(i) .NE. 0) THEN
               M(it)%reptFaces(M(it)%Nintfaces + indices(i)) = checkFaces(n_int_faces + i)
               checkFaces(n_int_faces + i) = 1
            END IF
         END DO

         DEALLOCATE (loc2glob_faces)
         DEALLOCATE (indices)

         ALLOCATE (auxext(M(it)%Nextfaces))
         auxext = 0
         ! Change order in exterior faces
         counter = 1
         DO i = M(it)%Nintfaces + 1, M(it)%Nfaces
            auxext(counter) = M(it)%reptFaces(i)
            counter = counter + 1
         END DO

         counter = 1
         DO i = M(it)%Nintfaces + 1, M(it)%Nfaces
            auxext(M(it)%mapextfa_loc(counter)) = M(it)%reptFaces(i)
            counter = counter + 1
         END DO

         counter = 1
         DO i = M(it)%Nintfaces + 1, M(it)%Nfaces
            M(it)%reptFaces(i) = auxext(counter)
            counter = counter + 1
         END DO

         DO i = 1, M(it)%Nfaces
            IF (M(it)%reptFaces(i) .EQ. 1) THEN
               M(it)%ghostFaces(i) = 1
            END IF
         END DO

         DEALLOCATE (auxext)

         ! count the number of ghost faces
         counter = 0
         DO i = 1, M(it)%Nfaces
            IF (M(it)%ghostFaces(i) .EQ. 1) THEN
               counter = counter + 1
            END IF
         END DO

         ALLOCATE (auxflp(counter))
         auxflp = 0
         ! fix flipFaces dimension

         s = 1
         DO i = counter - SIZE(M(it)%ghostFlp) + 1, counter
            auxflp(i) = M(it)%ghostFlp(s)
            s = s + 1
         END DO

         M(it)%ghostFlp = auxflp

         DEALLOCATE (auxflp)

#ifdef PARALL
         IF (MPIvar%glob_id .NE. nDiv - 1) THEN
            CALL MPI_Send(checkFaces, Nfaces, MPI_INTEGER, MPIvar%glob_id + 1, 0, MPI_COMM_WORLD, ierr)
         END IF
#endif
      END DO

#ifdef PARALL
      ! before sending make sure all processes have finished
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
      ! Share Loc2Glob_el and Loc2Glob_fa with the other processes.
      DO it = 1, nDiv
         IF (MPIvar%glob_id + 1 .EQ. it) THEN
            ! if we are the root process, send our data to everyone
            DO i = 0, MPIvar%glob_size - 1
               IF (i .NE. MPIvar%glob_id) THEN
                  n = SIZE(M(it)%Loc2Glob_el)
                  CALL MPI_Send(n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, ierr)
                  CALL MPI_Send(M(it)%Loc2Glob_el, SIZE(M(it)%Loc2Glob_el), MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
                  n = SIZE(M(it)%Loc2Glob_fa)
                  CALL MPI_Send(n, 1, MPI_INT, i, 0, MPI_COMM_WORLD, ierr)
                  CALL MPI_Send(M(it)%Loc2Glob_fa, SIZE(M(it)%Loc2Glob_fa), MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
               END IF
            END DO
         ELSE
            ! if we are a receiver process, receive the data from the root (it-th process)
            CALL MPI_Recv(n, 1, MPI_INT, it - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            ALLOCATE (M(it)%Loc2Glob_el(n))
            CALL MPI_Recv(M(it)%Loc2Glob_el, SIZE(M(it)%Loc2Glob_el), MPI_INTEGER, it - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            CALL MPI_Recv(n, 1, MPI_INT, it - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            ALLOCATE (M(it)%Loc2Glob_fa(n))
            CALL MPI_Recv(M(it)%Loc2Glob_fa, SIZE(M(it)%Loc2Glob_fa), MPI_INTEGER, it - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
         END IF
      END DO
      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

      ! Create structure for the communications
      IF (verbose) THEN
         IF (MPIvar%glob_id .EQ. 0) THEN
            WRITE (*, *) "Create structure for the communications"
         END IF
      END IF

#ifndef PARALL
      DO it = 1, nDiv
#else
      DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1
#endif
         IF (verbose) THEN
            WRITE (*, *) "Division: ", it
         END IF

         ! Count number of ghostfaces
         n = COUNT(M(it)%ghostFaces .EQ. 1)

         IF (n .EQ. 0) THEN
            ALLOCATE (M(it)%GhostPro(1))
            ALLOCATE (M(it)%GhostLoc(1))
            M(it)%GhostPro = -1
            M(it)%GhostLoc = -1
         ELSE
            ALLOCATE (indices(n))
            ALLOCATE (indices2(n))
            ALLOCATE (M(it)%GhostPro(n))
            ALLOCATE (M(it)%GhostLoc(n))

            M(it)%GhostPro = 0
            M(it)%GhostLoc = 0
            indices = 0
            indices2 = 0

            counter = 1
            DO i = 1, M(it)%Nfaces
               IF (M(it)%ghostFaces(i) .EQ. 1) THEN
                  indices(counter) = i
                  counter = counter + 1
               END IF
            END DO

            DO i = 1, nDiv
               IF (i .EQ. it) CYCLE

               indices2 = 0
               DO j = 1, n
                  DO jj = 1, SIZE(M(i)%Loc2Glob_fa)
                     IF (M(it)%Loc2Glob_fa(indices(j)) .EQ. M(i)%Loc2Glob_fa(jj)) THEN
                        indices2(j) = jj
                        EXIT
                     END IF
                  END DO
               END DO

               DO j = 1, n
                  IF (indices2(j) .NE. 0) THEN
                     M(it)%GhostPro(j) = i
                     M(it)%GhostLoc(j) = indices2(j)
                  END IF
               END DO

            END DO
            DEALLOCATE (indices2)
            DEALLOCATE (indices)
         END IF

         ! Elements
         n = COUNT(M(it)%ghostElems .EQ. 1)

         IF (n .EQ. 0) THEN
            ALLOCATE (M(it)%GhelsPro(1))
            ALLOCATE (M(it)%GhelsLoc(1))
            M(it)%GhelsPro = -1
            M(it)%GhelsLoc = -1
         ELSE
            ALLOCATE (indices(n))
            ALLOCATE (indices2(n))
            ALLOCATE (M(it)%GhelsPro(n))
            ALLOCATE (M(it)%GhelsLoc(n))
            M(it)%GhelsPro = 0
            M(it)%GhelsLoc = 0
            indices = 0
            indices2 = 0

            counter = 1
            DO i = 1, M(it)%Nelems
               IF (M(it)%ghostElems(i) .EQ. 1) THEN
                  indices(counter) = i
                  counter = counter + 1
               END IF
            END DO

            DO i = 1, nDiv
               IF (i .EQ. it) CYCLE

               indices2 = 0
               DO j = 1, n
                  IF (indices(j) .NE. 0) THEN
                     DO jj = 1, SIZE(M(i)%Loc2Glob_el)
                        IF (M(it)%Loc2Glob_el(indices(j)) .EQ. M(i)%Loc2Glob_el(jj)) THEN
                           indices2(j) = jj
                           EXIT
                        END IF
                     END DO
                  END IF
               END DO

               DO j = 1, n
                  IF (indices2(j) .NE. 0) THEN
                     M(it)%GhelsPro(j) = i
                     M(it)%GhelsLoc(j) = indices2(j)
                  END IF
               END DO

            END DO
            DEALLOCATE (indices2)
            DEALLOCATE (indices)
         END IF
      END DO

#ifndef PARALL
      DO it = 1, nDiv
#else
      DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1
#endif

         ! Merge Tb_CON with Tb
         ALLOCATE (Tb_prov(SIZE(M(it)%Tb, 1), SIZE(M(it)%Tb, 2)))
         ALLOCATE (boundaryFlag_loc(SIZE(M(it)%Tb, 1)))

         Tb_prov = M(it)%Tb
         boundaryFlag_loc = M(it)%boundaryFlag

         DEALLOCATE (M(it)%Tb)
         DEALLOCATE (M(it)%boundaryFlag)

         ALLOCATE (M(it)%Tb(SIZE(M(it)%TTb_CON, 1) + SIZE(Tb_prov, 1), SIZE(Tb_prov, 2)))
         ALLOCATE (M(it)%boundaryFlag(SIZE(M(it)%TTb_CON, 1) + SIZE(Tb_prov, 1)))

         M(it)%Tb(1:SIZE(Tb_prov, 1), :) = Tb_prov
         M(it)%Tb(SIZE(Tb_prov, 1) + 1:SIZE(M(it)%Tb, 1), :) = M(it)%TTb_CON
         M(it)%boundaryFlag = 0
         M(it)%boundaryFlag(1:SIZE(Tb_prov, 1)) = boundaryFlag_loc

         DEALLOCATE (Tb_prov)
         DEALLOCATE (boundaryFlag_loc)
      END DO

      ! Rotate back the mesh
#ifndef PARALL
      DO it = 1, nDiv
#else
      DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1
#endif

         IF (angle .GT. 1.e-12) THEN
            M(it)%XX(:, 1) = M(it)%XX(:, 1) - xcenter
            M(it)%XX(:, 2) = M(it)%XX(:, 2) - ycenter
            RotMat(1, 1) = COS((-angle)*pi_loc/180.)
            RotMat(1, 2) = -SIN((-angle)*pi_loc/180.)
            RotMat(2, 1) = SIN((-angle)*pi_loc/180.)
            RotMat(2, 2) = COS((-angle)*pi_loc/180.)
            M(it)%XX = TRANSPOSE(MATMUL(RotMat, TRANSPOSE(M(it)%XX)))
            M(it)%XX(:, 1) = M(it)%XX(:, 1) + xcenter
            M(it)%XX(:, 2) = M(it)%XX(:, 2) + ycenter
         END IF

      END DO

#ifdef PARALL
      CALL send_to_processes(M)
#endif

#ifndef PARALL
      ! free the structures
      DO it = 1, nDiv
#else
      DO it = MPIvar%glob_id + 1, MPIvar%glob_id + 1
#endif
         DEALLOCATE (M(it)%TT)
         DEALLOCATE (M(it)%XX)
         DEALLOCATE (M(it)%intfaces_loc)
         DEALLOCATE (M(it)%extfaces_loc)
         DEALLOCATE (M(it)%TTb_CON)
         DEALLOCATE (M(it)%Tb)
         DEALLOCATE (M(it)%reptFaces)
         DEALLOCATE (M(it)%ghostFaces)
         DEALLOCATE (M(it)%ghostElems)
         DEALLOCATE (M(it)%loc2glob_no)
         DEALLOCATE (M(it)%ghostFlp)
         DEALLOCATE (M(it)%mapextfa_loc)
         DEALLOCATE (M(it)%Tb_sizes)
         DEALLOCATE (M(it)%GhelsPro)
         DEALLOCATE (M(it)%GhelsLoc)
         DEALLOCATE (M(it)%GhostPro)
         DEALLOCATE (M(it)%GhostLoc)
         DEALLOCATE (M(it)%boundaryFlag)

      END DO

      ! loc2glob_el and loc2glob_fa are shared by all processes
      DO it = 1, nDiv
         DEALLOCATE (M(it)%loc2glob_fa)
         DEALLOCATE (M(it)%loc2glob_el)
      END DO

      NULLIFY (intfaces_glob, extfaces_glob, T, TTb, boundaryFlag, X)

   END SUBROUTINE split_mesh

#ifdef PARALL

   SUBROUTINE send_to_processes(M)
      USE printutils

      IMPLICIT NONE

      TYPE(DD_mesh), INTENT(IN)        :: M(:)
      INTEGER                          :: ierr
      REAL*8                           :: xmin
      CHARACTER(10)                    :: str
      REAL*8                           :: tol = 1e-6
#ifdef TOR3D
      INTEGER                          :: i
#endif

      CALL free_mesh

      Mesh%Nelems = M(MPIvar%glob_id + 1)%Nelems
      Mesh%Nnodesperelem = M(MPIvar%glob_id + 1)%Nnodesperelem
      Mesh%Nnodesperface = M(MPIvar%glob_id + 1)%Nnodesperface
      Mesh%Nnodes = M(MPIvar%glob_id + 1)%Nnodes
      Mesh%ndim = M(MPIvar%glob_id + 1)%ndim
      Mesh%Nintfaces = M(MPIvar%glob_id + 1)%Nintfaces
      Mesh%Nextfaces = M(MPIvar%glob_id + 1)%Nextfaces
      Mesh%elemType = M(MPIvar%glob_id + 1)%elemType
      Mesh%Nfaces = M(MPIvar%glob_id + 1)%Nfaces

      ALLOCATE (Mesh%T(Mesh%Nelems, Mesh%Nnodesperelem))
      ALLOCATE (Mesh%X(Mesh%Nnodes, Mesh%ndim))
      ALLOCATE (Mesh%Tb(Mesh%Nextfaces, Mesh%Nnodesperface))
      ALLOCATE (Mesh%boundaryFlag(Mesh%Nextfaces))
      ALLOCATE (Mesh%ghostFaces(Mesh%Nfaces))
      ALLOCATE (Mesh%loc2glob_fa(Mesh%Nfaces))
      ALLOCATE (Mesh%loc2glob_el(Mesh%Nelems))
      ALLOCATE (Mesh%loc2glob_nodes(Mesh%Nnodes))
      ALLOCATE (Mesh%ghostElems(Mesh%Nelems))
      !allocate (Mesh%intFaces(Mesh%Nintfaces,Mesh%Nnodesperface))
      !allocate (Mesh%extFaces(Mesh%Nextfaces,Mesh%Nnodesperface))

      Mesh%T = M(MPIvar%glob_id + 1)%TT
      Mesh%X = M(MPIvar%glob_id + 1)%XX
      Mesh%Tb = M(MPIvar%glob_id + 1)%Tb
      Mesh%boundaryFlag = M(MPIvar%glob_id + 1)%boundaryFlag
      Mesh%ghostFaces = M(MPIvar%glob_id + 1)%ghostFaces
      Mesh%ghostElems = M(MPIvar%glob_id + 1)%ghostElems
      Mesh%loc2glob_fa = M(MPIvar%glob_id + 1)%loc2glob_fa
      Mesh%loc2glob_el = M(MPIvar%glob_id + 1)%loc2glob_el
      Mesh%loc2glob_nodes = M(MPIvar%glob_id + 1)%loc2glob_no
      !Mesh%intFaces        = M(MPIvar%glob_id+1)%intfaces_loc
      !Mesh%intFaces        = M(MPIvar%glob_id+1)%extfaces_loc

      ! Find the number of ghost faces
      Mesh%nghostfaces = SUM(Mesh%ghostFaces)
      ! Find the number of ghost elements
      Mesh%nghostElems = SUM(Mesh%ghostElems)

      IF (Mesh%nghostfaces .EQ. 0) THEN
         ALLOCATE (Mesh%ghostflp(1))
         ALLOCATE (Mesh%ghostpro(1))
         ALLOCATE (Mesh%ghostloc(1))
      ELSE
         ALLOCATE (Mesh%ghostflp(Mesh%nghostfaces))
         ALLOCATE (Mesh%ghostpro(Mesh%nghostfaces))
         ALLOCATE (Mesh%ghostloc(Mesh%nghostfaces))
      END IF

      IF (Mesh%nghostElems .EQ. 0) THEN
         ALLOCATE (Mesh%ghelspro(1))
         ALLOCATE (Mesh%ghelsloc(1))
      ELSE
         ALLOCATE (Mesh%ghelspro(Mesh%nghostElems))
         ALLOCATE (Mesh%ghelsloc(Mesh%nghostElems))
      END IF

      Mesh%ghostflp = M(MPIvar%glob_id + 1)%ghostflp
      Mesh%ghostpro = M(MPIvar%glob_id + 1)%ghostpro
      Mesh%ghostloc = M(MPIvar%glob_id + 1)%ghostloc
      Mesh%ghelspro = M(MPIvar%glob_id + 1)%ghelspro
      Mesh%ghelsloc = M(MPIvar%glob_id + 1)%ghelsloc

#ifdef TOR3D
      IF (MPIvar%ntor .GT. 1) THEN
         DO i = 1, SIZE(Mesh%ghostPro)
            IF (Mesh%ghostPro(i) .GT. -1) THEN
               Mesh%ghostPro(i) = Mesh%ghostPro(i) + (MPIvar%itor - 1)*MPIvar%npol
            END IF
         END DO
      END IF

      IF (MPIvar%ntor .GT. 1) THEN
         DO i = 1, SIZE(Mesh%ghelspro)
            IF (Mesh%ghelspro(i) .GT. -1) THEN
               Mesh%ghelsPro(i) = Mesh%ghelsPro(i) + (MPIvar%itor - 1)*MPIvar%npol
            END IF
         END DO
      END IF

#endif

      Mesh%Nextfaces_nogho = COUNT((Mesh%boundaryFlag .NE. 0) .AND. (Mesh%ghostFaces(Mesh%Nintfaces + 1:SIZE(Mesh%ghostFaces)) .EQ. 0))
      Mesh%Nintfaces_nogho = COUNT(Mesh%ghostFaces(1:Mesh%Nintfaces) .EQ. 0)

      CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_el), Mesh%Nel_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_fa), Mesh%Nfa_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MAXVAL(Mesh%loc2glob_nodes), Mesh%Nno_glob, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(Mesh%Nextfaces_nogho, Mesh%Nextfaces_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(Mesh%Nintfaces_nogho, Mesh%Nintfaces_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      !call MPI_ALLREDUCE(Mesh%ndir, Mesh%Ndir_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(Mesh%nghostfaces, Mesh%Ngho_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      xmin = MINVAL(Mesh%X(:, 1))

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, xmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      ! Apply shift if axisymmetric case
      IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80) .OR. (switch%axisym .AND. MINVAL(Mesh%X(:,1)) < tol)) THEN
         IF (MPIvar%glob_id .EQ. 0) THEN
            WRITE (6, *) "*** Applying translation in axisymmetric case!"
         END IF
         Mesh%X(:, 1) = Mesh%X(:, 1) + geom%R0
      END IF

      ! Apply length scale
      Mesh%X = Mesh%X/phys%lscale

      Mesh%xmax = MAXVAL(Mesh%X(:, 1))
      Mesh%xmin = MINVAL(Mesh%X(:, 1))
      Mesh%ymax = MAXVAL(Mesh%X(:, 2))
      Mesh%ymin = MINVAL(Mesh%X(:, 2))

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%xmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, Mesh%ymin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      IF (utils%printint > 0) THEN
         IF (MPIvar%glob_id .EQ. 0) THEN
            IF (Mesh%elemType == 0) THEN
               WRITE (str, '(A)') 'triangles'
            ELSEIF (Mesh%elemType == 1) THEN
               WRITE (str, '(A)') 'quads'
            ELSEIF (Mesh%elemType == 2) THEN
               WRITE (str, '(A)') 'thetra'
            ELSEIF (Mesh%elemType == 3) THEN
               WRITE (str, '(A)') 'hexa'
            END IF
            WRITE (6, *) '*************************************************'
            WRITE (6, *) '*                    MESH                       *'
            WRITE (6, *) '*************************************************'
            WRITE (6, '(A,I18)') ' Number of dimensions:         ', Mesh%ndim
            WRITE (6, '(A,A34)') ' Element type: ', TRIM(str)
            WRITE (6, '(A,I18)') ' Number of elements:           ', Mesh%Nelems
            WRITE (6, '(A,I18)') ' Number of nodes:              ', Mesh%Nnodes
            WRITE (6, '(A,I18)') ' Number of nodes per element:  ', Mesh%Nnodesperelem
            WRITE (6, '(A,I18)') ' Number of nodes per face:     ', Mesh%Nnodesperface
            WRITE (6, '(A,I18)') ' Number of exterior faces:     ', Mesh%Nextfaces
            WRITE (6, *) ' '
            WRITE (6, *) ' '

            IF (utils%printint > 1) THEN
               WRITE (6, *) "Connectivity matrix T:"
               CALL displayMatrixInt(Mesh%T)
               WRITE (6, *) "Boundary connectivity matrix Tb:"
               CALL displayMatrixInt(Mesh%Tb)
            END IF
         END IF
      END IF

   END SUBROUTINE send_to_processes
#endif

END MODULE domain_decomposition_module
