!*****************************************
! project: MHDG
! file: debug.f90
! date: 01/02/2017
! Collection of debugging routines
!*****************************************
MODULE debug

   USE globals
   USE Printutils
   USE in_out

   IMPLICIT NONE

CONTAINS

   SUBROUTINE print_matrices_tex()

      write (6, *) "Saving elemental matrices"
      CALL saveArray(elMat%iAqq, 'iAqq')
      CALL saveArray(elMat%Aqu, 'Aqu')
      CALL saveArray(elMat%Aql, 'Aql')
      CALL saveArray(elMat%Auq, 'Auq')
      CALL saveArray(elMat%Auu, 'Auu')
      CALL saveArray(elMat%Aul, 'Aul')
      CALL saveArray(elMat%Alq, 'Alq')
      CALL saveArray(elMat%Alu, 'Alu')
      CALL saveArray(elMat%All, 'All')
      IF (Mesh%ndir .GT. 0) THEN
         CALL saveMatrix(elMat%Aql_dir, 'Aql_dir')
         CALL saveMatrix(elMat%Aul_dir, 'Aul_dir')
      END IF

      CALL saveArray(elMat%UU, 'UU')
      CALL saveMatrix(elMat%U0, 'U0')
      CALL saveArray(elMat%LL, 'LL')
      CALL saveMatrix(elMat%L0, 'L0')
   END SUBROUTINE print_matrices_tex

   SUBROUTINE print_matrices_hdf5()

      write (6, *) "Saving elemental matrices in hdf5"
      CALL HDF5_save_array(elMat%iAqq, 'iAqq')
      CALL HDF5_save_array(elMat%Aqu, 'Aqu')
      CALL HDF5_save_array(elMat%Aql, 'Aql')
      CALL HDF5_save_array(elMat%Auq, 'Auq')
      CALL HDF5_save_array(elMat%Auu, 'Auu')
      CALL HDF5_save_array(elMat%Aul, 'Aul')
      CALL HDF5_save_array(elMat%Alq, 'Alq')
      CALL HDF5_save_array(elMat%Alu, 'Alu')
      CALL HDF5_save_array(elMat%All, 'All')
      CALL HDF5_save_matrix(elMat%S, 'S')
      CALL HDF5_save_matrix(elMat%fh, 'fh')
      IF (Mesh%ndir .GT. 0) THEN
         CALL HDF5_save_matrix(elMat%Aql_dir, 'Aql_dir')
         CALL HDF5_save_matrix(elMat%Aul_dir, 'Aul_dir')
      END IF
      CALL HDF5_save_array(elMat%UU, 'UU')
      CALL HDF5_save_matrix(elMat%U0, 'U0')
      CALL HDF5_save_array(elMat%LL, 'LL')
      CALL HDF5_save_matrix(elMat%L0, 'L0')
   END SUBROUTINE print_matrices_hdf5

   SUBROUTINE print_matrices_singleElem_tex(iel)
      INTEGER, INTENT(IN)       :: iel

      write (6, *) "Saving elemental matrices"
      CALL saveMatrix(elMat%M(:, :, iel), 'M')
      CALL saveMatrix(elMat%Cv(:, :, iel), 'Cv')
      CALL saveMatrix(elMat%H(:, :, iel), 'H')
      IF (Mesh%ndir .GT. 0) THEN
         CALL saveVector(elMat%Hdir(:, iel), 'Hdir')
         CALL saveVector(elMat%Edir(:, iel), 'Edir')
      END IF
      CALL saveMatrix(elMat%D(:, :, iel), 'D')
      CALL saveMatrix(elMat%E(:, :, iel), 'E')
      CALL saveVector(elMat%S(:, iel), 'S')
      CALL saveMatrix(elMat%UU(:, :, iel), 'UU')
      CALL saveVector(elMat%U0(:, iel), 'U0')
      CALL saveMatrix(elMat%Hf(:, :, iel), 'Hf')
      CALL saveMatrix(elMat%Df(:, :, iel), 'Df')
      CALL saveMatrix(elMat%Ef(:, :, iel), 'Ef')
      CALL saveVector(elMat%fH(:, iel), 'fH')
      CALL saveMatrix(elMat%B(:, :, iel), 'B')
      CALL saveMatrix(elMat%C(:, :, iel), 'C')
      IF (Mesh%ndir .GT. 0) THEN
         CALL saveVector(elMat%Cdir(:, iel), 'Cdir')
      END IF
      CALL saveMatrix(elMat%P(:, :, iel), 'P')
      CALL saveMatrix(elMat%G(:, :, iel), 'G')
      CALL saveMatrix(elMat%iL(:, :, iel), 'iL')
      CALL saveMatrix(elMat%Lf(:, :, iel), 'Lf')
      CALL saveMatrix(elMat%Qf(:, :, iel), 'Qf')
      CALL saveMatrix(elMat%LL(:, :, iel), 'LL')
      CALL saveVector(elMat%L0(:, iel), 'L0')
#ifdef TEMPERATURE
      CALL saveMatrix(elMat%TQ(:, :, iel), 'TQ')
      CALL saveMatrix(elMat%TQhf(:, :, iel), 'TQhf')
      CALL saveVector(elMat%Tfhf(:, iel), 'Tfhf')
      CALL saveVector(elMat%Tf(:, iel), 'Tf')
      IF (Mesh%ndir .GT. 0) THEN
         CALL saveVector(elMat%Thdir(:, iel), 'Thdir')
      END IF
#endif
   END SUBROUTINE print_matrices_singleElem_tex

   SUBROUTINE print_elmat_tex(K, f, iel)
      REAL, INTENT(IN)          :: K(:, :), f(:)
      INTEGER, INTENT(IN)       :: iel
      CHARACTER(LEN=10)         :: num

      WRITE (num, "(I10)") iel
      CALL saveMatrix(K, "K_"//ADJUSTL(num))
      CALL saveVector(f, "f_"//ADJUSTL(num))

   END SUBROUTINE print_elmat_tex

   SUBROUTINE print_elmat_tex_par(K, f, iel)
      USE MPI_OMP
      REAL, INTENT(IN)          :: K(:, :), f(:)
      INTEGER, INTENT(IN)       :: iel
      CHARACTER(LEN=60)         :: num

      WRITE (num, "(A,I0,A,I0,A,I0)") "iel", iel, "_ipr", MPIvar%glob_id, "_npr", MPIvar%glob_size
      CALL saveMatrix(K, "K_"//ADJUSTL(num))
      CALL saveVector(f, "f_"//ADJUSTL(num))

   END SUBROUTINE print_elmat_tex_par

   FUNCTION maxmat(V) RESULT(m)
      REAL                      :: V(:, :)
      REAL                      :: m
      INTEGER                   :: i, j

      m = 0.
      DO i = 1, size(V, 1)
         DO j = 1, size(V, 2)
            m = max(V(i, j), m)
         END DO
      END DO

   END FUNCTION maxmat

   SUBROUTINE findNaNArr(A)
      REAL*8, INTENT(IN)        :: A(:, :, :)
      INTEGER                   :: i, j, k

      DO k = 1, size(A, 3)
         DO j = 1, size(A, 2)
            DO i = 1, size(A, 1)
               IF (ISNAN(A(i, j, k))) THEN
                  WRITE (6, *) "NaN detected in array at pos: i=", i, " j=", j, " k=", k
                  STOP
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE findNaNArr

   SUBROUTINE findNaNMat(M)
      REAL*8, INTENT(IN)        :: M(:, :)
      INTEGER                   :: i, j

      DO j = 1, size(M, 2)
         DO i = 1, size(M, 1)
            IF (ISNAN(M(i, j))) THEN
               WRITE (6, *) "NaN detected in matrix at pos: i=", i, " j=", j
               STOP
            END IF
         END DO
      END DO
   END SUBROUTINE findNaNMat

   SUBROUTINE findNaNVec(V)
      REAL*8, INTENT(IN)        :: V(:)
      INTEGER                   :: i

      DO i = 1, size(V, 1)
         IF (ISNAN(V(i))) THEN
            WRITE (6, *) "NaN detected in vector at pos: i=", i
            STOP
         END IF
      END DO
   END SUBROUTINE findNaNVec

END MODULE debug
