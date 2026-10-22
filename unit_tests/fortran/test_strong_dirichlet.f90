PROGRAM test_strong_dirichlet
  USE hdg_assembly_module, ONLY: HDG_assembly
  USE global_system_solver, ONLY: solve_global_system
  USE globals
  USE preprocess, ONLY: Bc_preprocess, Bc_preprocess_serial, CreateFaceConnectivity
  USE LinearAlgebra, ONLY: col, tensorProduct
  USE MPI_OMP, ONLY: init_MPI_OMP, MPIvar
  USE solver_status, ONLY: solver_status_t, solve_success
#ifdef PARALL
  USE Communications, ONLY: init_Com
#endif
#ifdef WITH_PASTIX
  USE solve_pastix, ONLY: terminate_mat_PASTIX
#endif
  USE mpi, ONLY: MPI_FINALIZE
  IMPLICIT NONE
  INTEGER :: ierr, case_id, entry, face

  INTEGER, PARAMETER :: Neq=7, Npfl=2, Ndim=2
  INTEGER :: iel, i, j, ind_if(Npfl), ind_jf(Npfl), ind_asf(Npfl), ind_ash(Npfl)
  INTEGER :: ind_fe(Neq*Npfl), ind_ff(Neq*Npfl), ind_fg(Neq*Npfl*Ndim), mode, prescribed
  REAL*8 :: tau(Neq,Neq), NNif(Npfl,Npfl), Nif(Npfl), Ni(Npfl), kmult(Npfl,Npfl)
  REAL*8 :: uf(Neq), nodal(Neq*Npfl), expected(Neq*Npfl), jump(Neq), normal(2), field(2), qg(2*Neq)
  LOGICAL :: isdir

  CALL init_MPI_OMP()
  utils%printint = 0
  utils%timing = .FALSE.
  phys%bcflags = bc_dirichlet_weak_form
  phys%bcflags(1) = bc_dirichlet

  DO entry = 1,2
#ifdef PARALL
     IF (entry == 2) CYCLE
#endif
     DO case_id = 1,4
        CALL check_ordering(case_id,entry)
     ENDDO
  ENDDO
  CALL check_stabilization()
#ifdef WITH_PASTIX
  CALL check_partitioned_solve(.TRUE.)
  CALL check_partitioned_solve(.FALSE.)
#endif
  IF (MPIvar%glob_id == 0) WRITE(*,*) 'Strong Dirichlet checks: PASS'
  CALL MPI_FINALIZE(ierr)

CONTAINS

  SUBROUTINE check_ordering(pattern,entry)
    INTEGER, INTENT(IN) :: pattern,entry
    INTEGER :: expected(4), flags(4)

    Mesh%Nextfaces = 4
    Mesh%Nintfaces = 0
    Mesh%Nfaces = 4
    ALLOCATE(Mesh%Tb(4,2),Mesh%extfaces(4,2),Mesh%boundaryFlag(4))
    Mesh%Tb(:,1) = [1,2,3,4]
    Mesh%Tb(:,2) = [11,12,13,14]
    Mesh%extfaces = Mesh%Tb
    SELECT CASE(pattern)
    CASE(1)
       flags = [2,1,2,2]
       expected = [1,3,4,2]
    CASE(2)
       flags = [2,1,2,1]
       expected = [1,3,2,4]
    CASE(3)
       flags = 2
       expected = [1,2,3,4]
    CASE(4)
       flags = 1
       expected = [1,2,3,4]
    END SELECT
    Mesh%boundaryFlag = flags
    IF (entry == 1) THEN
       CALL Bc_preprocess_serial()
    ELSE
       CALL Bc_preprocess()
    ENDIF
    IF (ANY(Mesh%Tb(:,1) /= expected) .OR. ANY(Mesh%Tb(:,2) /= expected+10)) ERROR STOP 'boundary connectivity permutation'
    IF (ANY(Mesh%extfaces /= Mesh%Tb) .OR. ANY(Mesh%boundaryFlag /= flags(expected))) &
         ERROR STOP 'boundary metadata permutation'
    IF (Mesh%ndir /= COUNT(flags==1) .OR. Mesh%ukf /= 4-COUNT(flags==1)) ERROR STOP 'prescribed face count'
    CALL free_mesh_loc(Mesh)
  END SUBROUTINE check_ordering

  SUBROUTINE check_stabilization()

    ALLOCATE(elMat%Auu(Neq*Npfl,Neq*Npfl,1),elMat%Aul(Neq*Npfl,Neq*Npfl,1),elMat%S(Neq*Npfl,1))
    ALLOCATE(elMat%Aul_dir(Neq*Npfl,1),elMat%Aql_dir(Ndim*Neq*Npfl,1))
    iel = 1
    ind_asf = [0,Neq]
    ind_ash = [0,Neq*Ndim]
    ind_fe = [(i,i=1,Neq*Npfl)]
    ind_ff = ind_fe
    ind_fg = [(i,i=1,Ndim*Neq*Npfl)]
    Ni = [0.3d0,0.7d0]
    Nif = 2.d0*Ni
    NNif = 2.d0*tensorProduct(Ni,Ni)
    uf = [(1.d0+0.2d0*i,i=1,Neq)]
    nodal = [uf,uf]
    normal = [1.d0,0.d0]
    field = [0.d0,1.d0]
    qg = 0.d0
    DO mode = 1,3
       numer%stab = 1
       IF (mode == 2) numer%stab = 5
       IF (mode == 3) numer%stab = 6
       tau = 0.d0
       numer%tau = 0.5d0
       DO i=1,Neq
          tau(i,i) = 0.5d0
          IF (mode > 1) tau(i,i) = 2.d0+i
          IF (mode == 3 .AND. i < Neq) tau(i,i+1) = 0.125d0*i
       ENDDO
       DO prescribed = 0,1
          isdir = prescribed == 1
          elMat%Auu = 0.d0
          elMat%Aul = 0.d0
          elMat%S = 0.d0
          elMat%Aul_dir = 0.d0
          elMat%Aql_dir = 0.d0
          CALL assemble_stabilization()
          IF (isdir) CALL assembly_dirichletstr_bc(iel,ind_asf,ind_ash,ind_fe,ind_fg,Nif,qg,uf,field,normal)
          expected = MATMUL(elMat%Auu(:,:,1),nodal)+MATMUL(elMat%Aul(:,:,1),nodal) &
               -elMat%S(:,1)+elMat%Aul_dir(:,1)
          IF (MAXVAL(ABS(expected)) > 2.d-13) ERROR STOP 'equal element/trace state has nonzero stabilization'
          jump = MATMUL(tau,SPREAD(1.d0,1,Neq))
          expected = [jump*Nif(1),jump*Nif(2)]
          IF (MAXVAL(ABS(MATMUL(elMat%Auu(:,:,1),nodal+1.d0)+MATMUL(elMat%Aul(:,:,1),nodal) &
               -elMat%S(:,1)+elMat%Aul_dir(:,1)-expected)) > 2.d-13) ERROR STOP 'prescribed stabilization jump'
       ENDDO
    ENDDO
    DEALLOCATE(elMat%Auu,elMat%Aul,elMat%S,elMat%Aul_dir,elMat%Aql_dir)
  END SUBROUTINE check_stabilization

    SUBROUTINE assemble_stabilization()
#include "dirichlet_stabilization.F"
    END SUBROUTINE assemble_stabilization
    ! Zero convection isolates stabilization without duplicating its production assembly.
    SUBROUTINE jacobianMatricesFace(state,bn,A)
      REAL*8, INTENT(IN) :: state(:),bn
      REAL*8, INTENT(OUT) :: A(:,:)
      A=0.d0
    END SUBROUTINE jacobianMatricesFace
#include "dirichlet_prescribed.F"

#ifdef WITH_PASTIX
  SUBROUTINE check_partitioned_solve(strong)
    LOGICAL, INTENT(IN) :: strong
    INTEGER :: ids(4), expected_ids(4), block, element, global_face, unknown_faces
    REAL*8 :: wanted(4)
    TYPE(solver_status_t) :: status

    IF (MPIvar%glob_size > 2) ERROR STOP 'strong Dirichlet fixture supports one or two ranks'
    phys%Neq = 2
    Mesh%Ndim = 2
    Mesh%Nelems = 4
    Mesh%Nnodesperelem = 2
    Mesh%Nnodesperface = 2
    Mesh%Nfaces = 4
    Mesh%Nextfaces = 4
    Mesh%Nintfaces = 0
    refElPol%Nnodes2D = 2
    refElPol%Nfaces = 1
    refElPol%Nfacenodes = 2
    ALLOCATE(Mesh%Tb(4,2),Mesh%extfaces(4,2),Mesh%intfaces(0,5),Mesh%boundaryFlag(4))
    ids = [1,3,2,4]
    IF (MPIvar%glob_id == 1) ids = [2,3,1,4]
    DO face = 1,4
       Mesh%Tb(face,:) = [ids(face),ids(face)+10]
       Mesh%extfaces(face,:) = [face,1]
    ENDDO
    Mesh%boundaryFlag = [2,1,2,1]
    IF (.NOT. strong) Mesh%boundaryFlag = 2
#ifdef PARALL
    ALLOCATE(Mesh%ghostfaces(4),Mesh%loc2glob_fa(4))
    Mesh%loc2glob_fa = ids
    Mesh%ghostfaces = 0
    Mesh%nghostfaces = 0
    IF (MPIvar%glob_size == 2) THEN
       Mesh%nghostfaces = 2
       Mesh%ghostfaces = [0,0,1,1]
       IF (MPIvar%glob_id == 1) Mesh%ghostfaces = [0,1,1,0]
       Mesh%boundaryFlag(3) = 0
       IF (strong) Mesh%boundaryFlag(1) = 0
    ENDIF
    ALLOCATE(Mesh%ghostpro(MAX(1,Mesh%nghostfaces)),Mesh%ghostloc(MAX(1,Mesh%nghostfaces)), &
         Mesh%ghostflp(MAX(1,Mesh%nghostfaces)))
    Mesh%ghostpro = 0
    Mesh%ghostloc = 0
    Mesh%ghostflp = 0
    IF (MPIvar%glob_size == 2) THEN
       Mesh%ghostpro = 2-MPIvar%glob_id
       Mesh%ghostloc = [1,4]
       Mesh%ghostflp = [1,0]
       IF (MPIvar%glob_id == 1) THEN
          Mesh%ghostloc = [2,1]
          Mesh%ghostflp = [0,1]
       ENDIF
    ENDIF
#endif
    CALL Bc_preprocess()
    expected_ids = [1,2,3,4]
    IF (MPIvar%glob_id == 1) expected_ids = [2,1,3,4]
    IF (.NOT. strong) expected_ids = ids
    IF (ANY(Mesh%Tb(:,1) /= expected_ids)) ERROR STOP 'partition face ordering'
#ifdef PARALL
    IF (ANY(Mesh%loc2glob_fa /= expected_ids)) ERROR STOP 'global face map detached from connectivity'
    IF (MPIvar%glob_size == 2 .AND. strong) THEN
       IF (ANY(Mesh%ghostflp /= [1,0])) ERROR STOP 'compressed ghost metadata ordering'
       IF (ANY(Mesh%ghostloc /= [1,4-MPIvar%glob_id])) ERROR STOP 'remote owner indices not refreshed'
    ENDIF
#endif
    CALL CreateFaceConnectivity()
#ifdef PARALL
    CALL init_Com()
#endif
    block = 4
    ALLOCATE(elMat%ALL(block,block,4),elMat%Alu(block,block,4),elMat%Alq(block,2*block,4), &
         elMat%UU(block,block,4),elMat%U0(block,4),elMat%LL(2*block,block,4),elMat%L0(2*block,4),elMat%fh(block,4))
    elMat%ALL = 0.d0
    elMat%Alu = 0.d0
    elMat%Alq = 0.d0
    elMat%UU = 0.d0
    elMat%U0 = 0.d0
    elMat%LL = 0.d0
    elMat%L0 = 0.d0
    elMat%fh = 0.d0
    DO face = 1,4
       element = Mesh%extfaces(face,1)
       global_face = expected_ids(face)
       DO j = 1,block
          elMat%ALL(j,j,element) = 2.d0
          elMat%fh(j,element) = 2.d0*(10.d0*global_face+j)
       ENDDO
    ENDDO
    matK%start = .TRUE.
    CALL HDG_assembly()
    unknown_faces = Mesh%ukf
#ifdef PARALL
    unknown_faces = COUNT(Mesh%ghostfaces(1:Mesh%ukf)==0)
#endif
    IF (matK%n /= block*unknown_faces) ERROR STOP 'owned unknowns double-subtracted prescribed ghosts'
    lssolver%sollib = 1
    lssolver%timing = .FALSE.
    switch%ckeramp = .FALSE.
    numer%dumpnr = 1.d0
    ALLOCATE(sol%u_tilde(16))
    sol%u_tilde = -17.d0
#ifdef WITH_PETSC
    CALL solve_global_system(status,1)
#else
    CALL solve_global_system(status)
#endif
    IF (status%outcome /= solve_success) ERROR STOP 'partitioned strong Dirichlet solve failed'
    DO face = 1,Mesh%ukf
       wanted = 10.d0*expected_ids(face)+[(j,j=1,block)]
       IF (MAXVAL(ABS(sol%u_tilde((face-1)*block+1:face*block)-wanted)) > 1.d-11) &
            ERROR STOP 'unknown trace packing or ghost exchange'
    ENDDO
    IF (ANY(sol%u_tilde(Mesh%ukf*block+1:) /= -17.d0)) ERROR STOP 'linear solve overwrote prescribed trace slots'
    CALL terminate_mat_PASTIX()
    CALL free_el_mat()
    DEALLOCATE(sol%u_tilde)
    NULLIFY(sol%u_tilde)
    CALL free_mesh_loc(Mesh)
  END SUBROUTINE check_partitioned_solve
#endif
END PROGRAM test_strong_dirichlet
