PROGRAM test_hdg_residual
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value,ieee_quiet_nan,ieee_positive_inf
  USE globals
  USE LinearAlgebra, ONLY: invert_matrix
  USE hdg_residual_module
  USE balance_diagnostics, ONLY: balance_diag,balance_mode_equations
  USE hdg_assembly_module, ONLY: HDG_assembly
  USE MPI_OMP, ONLY: init_MPI_OMP,MPIvar
  USE omp_lib, ONLY: omp_get_max_threads
  USE mpi, ONLY: MPI_FINALIZE
  IMPLICIT NONE
  INTEGER, PARAMETER :: nv=6,nq=12,block=4,nf=12,ne=4,nh=32,volume_size=72,nt=104
  REAL*8 :: system(nt,nt),load(nt),state(nt),expected(nt),worst=0.d0
  TYPE(hdg_residual_t) :: residual
  TYPE(hdg_residual_norms_t) :: norms
  INTEGER :: seed,ierr,expected_ranks
  CHARACTER(16) :: rank_argument

  CALL init_MPI_OMP()
  expected_ranks=1
  CALL get_command_argument(1,rank_argument)
  IF(LEN_TRIM(rank_argument)>0)READ(rank_argument,*)expected_ranks
  IF(MPIvar%glob_size/=expected_ranks)ERROR STOP 'Unexpected MPI communicator size: check the launcher'
  IF(MPIvar%glob_size>2)ERROR STOP 'Residual fixture supports one or two ranks'
  CALL initialize_fixture()
  DO seed=1,4
    CALL check_fixture(seed,.FALSE.)
    CALL check_fixture(seed,.TRUE.)
  ENDDO
  CALL check_constant_representers()
  CALL check_convergence_policy()
  CALL check_no_trace_unknowns()
  CALL check_reporting_modes()
  IF(MPIvar%glob_id==0)WRITE(*,*) &
       'HDG residual PASS: 8 dense systems, assembler parity; ranks/threads=', &
       MPIvar%glob_size,omp_get_max_threads(),'; worst=',worst
  CALL MPI_FINALIZE(ierr)

CONTAINS

  SUBROUTINE initialize_fixture()
    INTEGER :: e,i
    REAL*8 :: t
    phys%neq=2;Mesh%Ndim=2;Mesh%Nelems=ne;Mesh%Nnodesperelem=3
    Mesh%Nnodesperface=2;Mesh%Nfaces=9;Mesh%Nintfaces=3;Mesh%Nextfaces=6
    Mesh%ukf=8;Mesh%Ndir=1
    refElPol%Nfaces=3;refElPol%Nfacenodes=2;refElPol%Nnodes2D=3
    utils%timing=.FALSE.;utils%printint=0;switch%axisym=.TRUE.
    refElPol%Ngauss2d=3;refElPol%Ngauss1d=2
    ALLOCATE(Mesh%X(6,2),Mesh%T(ne,3),refElPol%face_nodes(3,2), &
         refElPol%N2D(3,3),refElPol%Nxi2D(3,3),refElPol%Neta2D(3,3),refElPol%gauss_weights2D(3), &
         refElPol%N1D(2,2),refElPol%Nxi1D(2,2),refElPol%gauss_weights1D(2))
    Mesh%X(:,1)=[1.d0,2.d0,1.d0,2.d0,1.d0,2.d0]
    Mesh%X(:,2)=[0.d0,0.d0,1.d0,1.d0,2.d0,2.d0]
    Mesh%T(1,:)=[1,2,3];Mesh%T(2,:)=[2,4,3];Mesh%T(3,:)=[3,4,5];Mesh%T(4,:)=[4,6,5]
    refElPol%face_nodes(1,:)=[1,2];refElPol%face_nodes(2,:)=[2,3];refElPol%face_nodes(3,:)=[3,1]
    refElPol%N2D=1.d0/6.d0
    DO i=1,3
      refElPol%N2D(i,i)=2.d0/3.d0
    ENDDO
    DO i=1,3
      refElPol%Nxi2D(i,:)=[-1.d0,1.d0,0.d0];refElPol%Neta2D(i,:)=[-1.d0,0.d0,1.d0]
    ENDDO
    refElPol%gauss_weights2D=1.d0/6.d0
    t=1.d0/SQRT(3.d0)
    refElPol%N1D(1,:)=[(1.d0+t)/2,(1.d0-t)/2]
    refElPol%N1D(2,:)=refElPol%N1D(1,2:1:-1)
    refElPol%Nxi1D(:,1)=-.5d0;refElPol%Nxi1D(:,2)=.5d0
    refElPol%gauss_weights1D=1.d0
    numer%residual_scales=1.d0
    numer%residual_scales(2)=3.d0
    ALLOCATE(phys%conVarNam(2));phys%conVarNam=['density','energy ']
    ALLOCATE(Mesh%F(ne,3),Mesh%Fdir(ne,3),Mesh%flipface(ne,3), &
         Mesh%intfaces(3,4),Mesh%extfaces(6,2),Mesh%periodic_faces(6))
    Mesh%F(1,:)=[4,1,5];Mesh%F(2,:)=[6,2,1]
    Mesh%F(3,:)=[2,3,7];Mesh%F(4,:)=[8,9,3]
    Mesh%Fdir=.FALSE.;Mesh%Fdir(ne,2)=.TRUE.
    Mesh%flipface=.FALSE.;Mesh%flipface(2,3)=.TRUE.;Mesh%flipface(3,1)=.TRUE.;Mesh%flipface(4,3)=.TRUE.
    Mesh%intfaces(1,:)=[1,2,2,3];Mesh%intfaces(2,:)=[2,2,3,1];Mesh%intfaces(3,:)=[3,2,4,3]
    Mesh%extfaces(1,:)=[1,1];Mesh%extfaces(2,:)=[1,3]
    Mesh%extfaces(3,:)=[2,1];Mesh%extfaces(4,:)=[3,3]
    Mesh%extfaces(5,:)=[4,1];Mesh%extfaces(6,:)=[4,2]
#ifdef PARALL
    ALLOCATE(Mesh%ghostElems(ne),Mesh%ghostFaces(9),Mesh%loc2glob_fa(9))
    Mesh%loc2glob_fa=[1,2,3,4,5,6,7,8,9]
#endif
    ALLOCATE(sol%u(nv*ne),sol%q(nq*ne),sol%u_tilde(9*block))
    ALLOCATE(elMat%Auu(nv,nv,ne),elMat%Auq(nv,nq,ne),elMat%Aul(nv,nf,ne), &
         elMat%Aqu(nq,nv,ne),elMat%iAqq(nq,nq,ne),elMat%Aql(nq,nf,ne), &
         elMat%Alu(nf,nv,ne),elMat%Alq(nf,nq,ne),elMat%All(nf,nf,ne), &
         elMat%S(nv,ne),elMat%Aul_dir(nv,ne),elMat%Aql_dir(nq,ne),elMat%fh(nf,ne), &
         elMat%UU(nv,nf,ne),elMat%U0(nv,ne),elMat%LL(nq,nf,ne),elMat%L0(nq,ne))
  END SUBROUTINE initialize_fixture

  SUBROUTINE set_partition(periodic)
    LOGICAL, INTENT(IN) :: periodic
    Mesh%periodic_faces=0
    IF(periodic)THEN
      Mesh%periodic_faces(1)=5;Mesh%periodic_faces(5)=1
    ENDIF
#ifdef PARALL
    Mesh%ghostElems=0;Mesh%ghostFaces=0
    IF(MPIvar%glob_size==1)RETURN
    IF(MPIvar%glob_id==0)THEN
      Mesh%ghostElems=[0,0,1,1];Mesh%ghostFaces=[0,0,1,0,0,0,1,1,1]
      IF(periodic)Mesh%ghostFaces(8)=0
    ELSE
      Mesh%ghostElems=[1,1,0,0];Mesh%ghostFaces=[1,1,0,1,1,1,0,0,0]
      IF(periodic)Mesh%ghostFaces(8)=1
    ENDIF
#endif
  END SUBROUTINE set_partition

  SUBROUTINE set_matrices(seed)
    INTEGER, INTENT(IN) :: seed
    INTEGER :: e,i,j,k,indices(3)
    REAL*8 :: inverse(3,3)
    DO e=1,ne
      CALL pattern(elMat%Auu(:,:,e),.3d0,seed+e)
      CALL pattern(elMat%Auq(:,:,e),.2d0,seed+e+1)
      CALL pattern(elMat%Aul(:,:,e),.4d0,seed+e+2)
      CALL pattern(elMat%Aqu(:,:,e),.1d0,seed+e+3)
      CALL pattern(elMat%Aql(:,:,e),.3d0,seed+e+4)
      CALL pattern(elMat%Alu(:,:,e),.2d0,seed+e+5)
      CALL pattern(elMat%Alq(:,:,e),.1d0,seed+e+6)
      CALL pattern(elMat%All(:,:,e),.2d0,seed+e+7)
      elMat%iAqq(:,:,e)=0.d0
      DO i=1,nv
        elMat%Auu(i,i,e)=elMat%Auu(i,i,e)+3.d0
        elMat%S(i,e)=.1d0*(i+e)
        elMat%Aul_dir(i,e)=.017d0*(i-2*e)
      ENDDO
      DO i=1,nq
        elMat%Aql_dir(i,e)=.023d0*(i+e)
      ENDDO
      CALL invert_matrix(reference_element_mass(e),inverse)
      DO j=1,2
        DO k=1,2
          indices=[(k+(j-1)*2+(i-1)*4,i=1,3)]
          elMat%iAqq(indices,indices,e)=inverse
        ENDDO
      ENDDO
      DO i=1,nf
        elMat%All(i,i,e)=elMat%All(i,i,e)+2.d0
        elMat%fh(i,e)=.07d0*(i-e)
      ENDDO
      IF(e==ne)THEN
        elMat%Aul(:,5:8,e)=0;elMat%Aql(:,5:8,e)=0;elMat%All(:,5:8,e)=0
      ENDIF
    ENDDO
  END SUBROUTINE set_matrices

  SUBROUTINE check_fixture(seed,periodic)
    INTEGER, INTENT(IN) :: seed
    LOGICAL, INTENT(IN) :: periodic
    REAL*8 :: saved_system(nt,nt),saved_load(nt),expected_rms(3,2),expected_peak(3,2)
    INTEGER :: e,j,k,offset

    CALL set_partition(periodic)
#ifdef PARALL
    IF(MPIvar%glob_size==2.AND.seed==2.AND..NOT.periodic)THEN
      Mesh%ghostElems=0
      IF(MPIvar%glob_id==1)Mesh%ghostElems=1
    ENDIF
#endif
    CALL hdg_residual%invalidate_mesh()
    utils%timing=seed==4.AND.periodic
    timing%runtass=0.d0;timing%cputass=0.d0
    CALL set_matrices(seed)
    state=[(.5d0*SIN(REAL(j+seed,8)),j=1,nt)]
    DO e=1,ne
      offset=(e-1)*(nv+nq)
      sol%u((e-1)*nv+1:e*nv)=state(offset+1:offset+nv)
      sol%q((e-1)*nq+1:e*nq)=state(offset+nv+1:offset+nv+nq)
    ENDDO
    sol%u_tilde(1:nh)=state(volume_size+1:nt)
    sol%u_tilde(nh+1:)=1.d10 ! A prescribed slot is not an unknown column.
    CALL assemble_dense_reference()
    saved_system=system;saved_load=load
    expected=MATMUL(system,state)-load
    CALL evaluate_hdg_residual(residual)
    ! Local arrays are also needed on ghost elements; evaluating them requires no MPI reduction.
    DO e=1,ne
      offset=(e-1)*(nv+nq)
      CALL close(residual%u(:,e),expected(offset+1:offset+nv),'conservation residual')
      CALL close(residual%q(:,e),expected(offset+nv+1:offset+nv+nq),'inverse-mass gradient residual')
    ENDDO
    CALL hdg_residual%measure(residual,norms)
    IF(norms%evaluation_seconds<0.d0)ERROR STOP 'invalid residual timing'
    CALL reference_mass_norms(expected_rms,expected_peak)
    CALL close(RESHAPE(norms%rms,[6]),RESHAPE(expected_rms,[6]),'mass per-equation norms')
    CALL close(RESHAPE(norms%local_max,[6]),RESHAPE(expected_peak,[6]),'mass local maxima')
    CALL close([norms%total_rms],[SQRT(SUM(expected_rms**2)/6.d0)],'total scaled mass norm')
    CALL assemble_dense_reference()
    CALL close(RESHAPE(system-saved_system,[nt*nt]),SPREAD(0.d0,1,nt*nt),'matrix preservation')
    CALL close(load,saved_load,'affine RHS preservation')
    CALL check_assembler_trace()
  END SUBROUTINE check_fixture

  FUNCTION reference_element_mass(element) RESULT(mass)
    INTEGER, INTENT(IN) :: element
    REAL*8 :: mass(3,3),weight
    INTEGER :: g,a,b
    mass=0.d0
    ! Each fixture triangle has area 1/2. Integrate test polynomials directly.
    DO g=1,3
      weight=DOT_PRODUCT(refElPol%N2D(g,:),Mesh%X(Mesh%T(element,:),1))/6.d0
      DO b=1,3
        DO a=1,3
          mass(a,b)=mass(a,b)+weight*refElPol%N2D(g,a)*refElPol%N2D(g,b)
        ENDDO
      ENDDO
    ENDDO
  END FUNCTION reference_element_mass

  FUNCTION reference_face_mass(face) RESULT(mass)
    INTEGER, INTENT(IN) :: face
    REAL*8 :: mass(2,2),node(2,2),weight,x(2)
    INTEGER :: e,k,element,local_face,g,a,b
    element=0
    DO e=1,ne
      DO k=1,3
        IF(Mesh%F(e,k)/=face)CYCLE
        element=e;local_face=k
        EXIT
      ENDDO
      IF(element>0)EXIT
    ENDDO
    node=Mesh%X(Mesh%T(element,refElPol%face_nodes(local_face,:)),:)
    IF(Mesh%flipface(element,local_face))node=node(2:1:-1,:)
    mass=0.d0
    DO g=1,2
      x=MATMUL(refElPol%N1D(g,:),node)
      weight=x(1)*SQRT(SUM((node(2,:)-node(1,:))**2))/2.d0
      DO b=1,2
        DO a=1,2
          mass(a,b)=mass(a,b)+weight*refElPol%N1D(g,a)*refElPol%N1D(g,b)
        ENDDO
      ENDDO
    ENDDO
  END FUNCTION reference_face_mass

  SUBROUTINE reference_mass_norms(rms,peak)
    REAL*8, INTENT(OUT) :: rms(3,2),peak(3,2)
    REAL*8 :: mass(3,3),inverse(3,3),face_inverse(2,2),measure
    REAL*8 :: sums(3,2),measures(3),local(3),q(3),r(3)
    INTEGER :: e,j,k,f,offset
    sums=0.d0;peak=0.d0;measures=0.d0
    DO e=1,ne
      mass=reference_element_mass(e)
      measure=SUM(mass);measures(1:2)=measures(1:2)+measure
      CALL invert_matrix(mass,inverse)
      offset=(e-1)*(nv+nq)
      DO j=1,2
        r=expected(offset+j:offset+nv:2)
        local(1)=DOT_PRODUCT(r,MATMUL(inverse,r))
        local(2)=0.d0
        DO k=1,2
          q=expected(offset+nv+k+(j-1)*2:offset+nv+nq:4)
          local(2)=local(2)+DOT_PRODUCT(q,MATMUL(mass,q))/2.d0
        ENDDO
        sums(1:2,j)=sums(1:2,j)+local(1:2)
        peak(1:2,j)=MAX(peak(1:2,j),local(1:2)/measure)
      ENDDO
    ENDDO
    DO f=1,8
      measure=SUM(reference_face_mass(f));measures(3)=measures(3)+measure
      CALL invert_matrix(reference_face_mass(f),face_inverse)
      offset=volume_size+(f-1)*block
      DO j=1,2
        local(3)=DOT_PRODUCT(expected(offset+j:offset+block:2), &
             MATMUL(face_inverse,expected(offset+j:offset+block:2)))
        sums(3,j)=sums(3,j)+local(3)
        peak(3,j)=MAX(peak(3,j),local(3)/measure)
      ENDDO
    ENDDO
    DO j=1,2
      rms(:,j)=SQRT(sums(:,j)/measures)/numer%residual_scales(j)
      peak(:,j)=SQRT(peak(:,j))/numer%residual_scales(j)
    ENDDO
  END SUBROUTINE reference_mass_norms

  SUBROUTINE check_constant_representers()
    INTEGER :: e,j,k,f,face,offset,neighbours
    REAL*8 :: mass(3,3),face_mass(2,2)
    CALL set_partition(.FALSE.)
    CALL hdg_residual%invalidate_mesh()
    residual%trace=0.d0
    DO e=1,ne
      mass=reference_element_mass(e)
      DO j=1,2
        residual%u(j:nv:2,e)=SUM(mass,DIM=2)*numer%residual_scales(j)
        DO k=1,2
          residual%q(k+(j-1)*2:nq:4,e)=numer%residual_scales(j)
        ENDDO
      ENDDO
      DO f=1,3
        IF(Mesh%Fdir(e,f))CYCLE
        face=Mesh%F(e,f);offset=(f-1)*block
        neighbours=COUNT(Mesh%F==face);face_mass=reference_face_mass(face)
        DO j=1,2
          residual%trace(offset+j:offset+block:2,e)=SUM(face_mass,DIM=2)*numer%residual_scales(j)/neighbours
        ENDDO
      ENDDO
    ENDDO
    CALL hdg_residual%measure(residual,norms)
    CALL close(RESHAPE(norms%rms,[6]),SPREAD(1.d0,1,6),'constant representer global RMS')
    CALL close(RESHAPE(norms%local_max,[6]),SPREAD(1.d0,1,6),'constant representer peaks')
    CALL close([norms%total_rms],[1.d0],'constant representer total')
  END SUBROUTINE check_constant_representers

  SUBROUTINE check_convergence_policy()
    TYPE(hdg_residual_t) :: constant
    TYPE(hdg_residual_monitor_t) :: independent
    constant=residual ! The preceding constant-representer test has norm one.
    CALL hdg_residual%begin_solve()
    numer%newton_convergence='residual';numer%residual_atol=1.d-6;numer%residual_rtol=1.d-4
    CALL set_policy_residual(2.d0,constant)
    IF(hdg_residual%converged(0.d0))ERROR STOP 'initial nonzero residual accepted'
    CALL independent%measure(residual)
    CALL set_policy_residual(1.d-4,constant)
    IF(.NOT.hdg_residual%converged(1.d0))ERROR STOP 'residual convergence missing'
    IF(independent%converged(0.d0))ERROR STOP 'monitor histories are shared'
    numer%newton_convergence='both'
    IF(hdg_residual%converged(1.d0))ERROR STOP 'both ignores large update'
    numer%tNR=2.d-4
    IF(.NOT.hdg_residual%converged(1.d-5))ERROR STOP 'both rejects valid state'
    CALL hdg_residual%invalidate_mesh()
    IF(hdg_residual%converged(0.d0))ERROR STOP 'mesh reset retained convergence'
    CALL set_policy_residual(0.d0,constant)
    IF(hdg_residual%converged(0.d0))ERROR STOP 'both accepted without same-mesh update'
    numer%newton_convergence='residual'
    IF(.NOT.hdg_residual%converged(1.d0))ERROR STOP 'zero initial residual rejected'
    CALL close([norms%relative],[0.d0],'zero residual relative baseline')
    CALL set_policy_residual(1.d0,constant)
    IF(norms%relative/=HUGE(1.d0).OR.hdg_residual%converged(0.d0))ERROR STOP 'zero baseline masked a nonzero residual'
    CALL set_policy_residual(ieee_value(0.d0,ieee_quiet_nan),constant)
    IF(hdg_residual%converged(0.d0))ERROR STOP 'NaN residual accepted'
    CALL hdg_residual%begin_solve()
    CALL set_policy_residual(ieee_value(0.d0,ieee_positive_inf),constant)
    IF(hdg_residual%converged(0.d0))ERROR STOP 'Infinite initial residual accepted'
    numer%newton_convergence='update'
    IF(hdg_residual%converged(0.d0))ERROR STOP 'update stopping changed'
    CALL hdg_residual%begin_solve()
  END SUBROUTINE check_convergence_policy

  SUBROUTINE set_policy_residual(scale,constant)
    REAL*8, INTENT(IN) :: scale
    TYPE(hdg_residual_t), INTENT(IN) :: constant
    residual%u=scale*constant%u;residual%q=scale*constant%q;residual%trace=scale*constant%trace
    CALL hdg_residual%measure(residual,norms)
  END SUBROUTINE set_policy_residual

  SUBROUTINE assemble_dense_reference()
    INTEGER :: e,i,j,f,offset,global(nf),local(nf),face,partner
    system=0;load=0
    DO e=1,ne
      CALL trace_indices(e,global,local)
      offset=(e-1)*(nv+nq)
      system(offset+1:offset+nv,offset+1:offset+nv)=elMat%Auu(:,:,e)
      system(offset+1:offset+nv,offset+nv+1:offset+nv+nq)=elMat%Auq(:,:,e)
      system(offset+nv+1:offset+nv+nq,offset+1:offset+nv)=MATMUL(elMat%iAqq(:,:,e),elMat%Aqu(:,:,e))
      DO i=1,nq
        system(offset+nv+i,offset+nv+i)=1.d0
      ENDDO
      load(offset+1:offset+nv)=elMat%S(:,e)-elMat%Aul_dir(:,e)
      load(offset+nv+1:offset+nv+nq)=-MATMUL(elMat%iAqq(:,:,e),elMat%Aql_dir(:,e))
      DO j=1,nf
        IF(local(j)>nh)CYCLE
        system(offset+1:offset+nv,volume_size+local(j))=elMat%Aul(:,j,e)
        system(offset+nv+1:offset+nv+nq,volume_size+local(j))=MATMUL(elMat%iAqq(:,:,e),elMat%Aql(:,j,e))
      ENDDO
      DO f=1,3
        IF(Mesh%Fdir(e,f))CYCLE
        face=Mesh%F(e,f)
        CALL insert_trace_row(e,f,face,global)
        IF(face<=Mesh%Nintfaces)CYCLE
        partner=Mesh%periodic_faces(face-Mesh%Nintfaces)
        IF(partner==0)CYCLE
        partner=Mesh%F(Mesh%extfaces(partner,1),Mesh%extfaces(partner,2))
        global((f-1)*block+1:f*block)=(partner-1)*block+[(i,i=1,block)]
        CALL insert_trace_row(e,f,partner,global)
        CALL trace_indices(e,global,local)
      ENDDO
    ENDDO
  END SUBROUTINE assemble_dense_reference

  SUBROUTINE insert_trace_row(element,face,target,columns)
    INTEGER, INTENT(IN) :: element,face,target,columns(nf)
    INTEGER :: i,j,row,source,offset
    offset=(element-1)*(nv+nq)
    DO i=1,block
      row=volume_size+(target-1)*block+i;source=(face-1)*block+i
      system(row,offset+1:offset+nv)=system(row,offset+1:offset+nv)+elMat%Alu(source,:,element)
      system(row,offset+nv+1:offset+nv+nq)=system(row,offset+nv+1:offset+nv+nq)+elMat%Alq(source,:,element)
      load(row)=load(row)+elMat%fh(source,element)
      DO j=1,nf
        IF(columns(j)>nh)CYCLE
        system(row,volume_size+columns(j))=system(row,volume_size+columns(j))+elMat%All(source,j,element)
      ENDDO
    ENDDO
  END SUBROUTINE insert_trace_row

  SUBROUTINE check_assembler_trace()
    REAL*8 :: assembled
    INTEGER :: e,row,entry,global_row
    ! Constant local maps make K*U_hat-rhs exactly the current uncondensed trace residual.
    elMat%UU=0;elMat%LL=0
    DO e=1,ne
      elMat%U0(:,e)=sol%u((e-1)*nv+1:e*nv)
      elMat%L0(:,e)=sol%q((e-1)*nq+1:e*nq)
    ENDDO
    CALL HDG_assembly()
    DO row=1,MatK%n
      assembled=-rhs%vals(row)
      DO entry=MatK%rowptr(row),MatK%rowptr(row+1)-1
        assembled=assembled+MatK%vals(entry)*sol%u_tilde(MatK%cols(entry))
      ENDDO
      global_row=MatK%loc2glob(row)
      CALL close([assembled],[expected(volume_size+global_row)],'actual assembler trace row')
    ENDDO
  END SUBROUTINE check_assembler_trace

  SUBROUTINE trace_indices(element,global,local)
    INTEGER, INTENT(IN) :: element
    INTEGER, INTENT(OUT) :: global(nf),local(nf)
    INTEGER :: f,offset,i
    DO f=1,3
      offset=(f-1)*block
      global(offset+1:offset+block)=(Mesh%F(element,f)-1)*block+[(i,i=1,block)]
      local(offset+1:offset+block)=global(offset+1:offset+block)
      IF(Mesh%flipface(element,f))local(offset+1:offset+block)=global(offset+[3,4,1,2])
    ENDDO
  END SUBROUTINE trace_indices

  SUBROUTINE check_no_trace_unknowns()
    Mesh%Fdir=.TRUE.;Mesh%ukf=0;Mesh%Ndir=9;Mesh%periodic_faces=0
    CALL hdg_residual%invalidate_mesh()
    CALL evaluate_hdg_residual(residual)
    CALL hdg_residual%measure(residual,norms)
    IF(norms%measures(3)/=0.OR.ANY(norms%rms(3,:)/=0))ERROR STOP 'empty trace norm'
    IF(norms%total_rms<0.d0)ERROR STOP 'empty trace produced invalid total'
  END SUBROUTINE check_no_trace_unknowns

  SUBROUTINE check_reporting_modes()
    INTEGER :: unit,selected
    CHARACTER(9), PARAMETER :: modes(3)=[CHARACTER(9) :: 'off','summary','equations']
    utils%timing=.FALSE.
    IF(TRIM(hdg_residual%get_mode())/='off')ERROR STOP 'residual output default is active'
    IF(TRIM(numer%newton_convergence)/='update')ERROR STOP 'residual stopping default is active'
    ! Residual output works with physical balances off, including independent detail.
    DO selected=1,3
      CALL hdg_residual%configure(modes(selected))
      CALL hdg_residual%measure(residual)
      SELECT CASE(selected)
      CASE(1)
        CALL check_report_lines(0)
      CASE(2)
        CALL check_report_lines(2)
      CASE(3)
        CALL check_report_lines(5)
      END SELECT
    ENDDO
    ! Beginning another solve must suppress stale output until evaluation.
    CALL hdg_residual%begin_solve()
    CALL check_report_lines(0)
    ! Balance verbosity cannot enable residual reporting.
    CALL balance_diag%configure(balance_mode_equations)
    CALL hdg_residual%configure('off')
    CALL hdg_residual%measure(residual)
    CALL check_report_lines(0)
    ! Reporting off skips diagnostics under update stopping, but not evaluation
    ! required for residual acceptance. The zero state has exactly zero residual.
    sol%u=0.d0;sol%q=0.d0;sol%u_tilde=0.d0
    elMat%S=0.d0;elMat%Aul_dir=0.d0;elMat%Aql_dir=0.d0;elMat%fh=0.d0
    CALL hdg_residual%begin_solve()
    CALL hdg_residual%evaluate()
    numer%newton_convergence='residual'
    IF(hdg_residual%converged(0.d0))ERROR STOP 'inactive diagnostics evaluated a residual'
    CALL hdg_residual%evaluate()
    IF(.NOT.hdg_residual%converged(0.d0))ERROR STOP 'silent residual stopping skipped evaluation'
    CALL check_report_lines(0)
    numer%newton_convergence='update'
  END SUBROUTINE check_reporting_modes

  SUBROUTINE check_report_lines(expected_lines)
    INTEGER, INTENT(IN) :: expected_lines
    CHARACTER(256) :: line
    INTEGER :: unit,lines,status
    OPEN(NEWUNIT=unit,STATUS='scratch',ACTION='readwrite')
    CALL hdg_residual%report(unit)
    REWIND(unit);lines=0
    DO
      READ(unit,'(A)',IOSTAT=status)line
      IF(status/=0)EXIT
      lines=lines+1
      IF(lines==1.AND.INDEX(line,'before update: norm=')==0)ERROR STOP 'summary label'
    ENDDO
    CLOSE(unit)
    IF(MPIvar%glob_id/=0)THEN
      IF(lines/=0)ERROR STOP 'non-writer residual output'
    ELSEIF(lines/=expected_lines)THEN
      WRITE(*,*)'report lines actual/expected: ',lines,expected_lines
      ERROR STOP 'residual reporting policy'
    ENDIF
  END SUBROUTINE check_report_lines

  SUBROUTINE pattern(matrix,scale,seed)
    REAL*8, INTENT(OUT) :: matrix(:,:)
    REAL*8, INTENT(IN) :: scale
    INTEGER, INTENT(IN) :: seed
    INTEGER :: i,j
    DO j=1,SIZE(matrix,2)
      DO i=1,SIZE(matrix,1)
        matrix(i,j)=scale*SIN(REAL(3*i+7*j+seed,8))
      ENDDO
    ENDDO
  END SUBROUTINE pattern

  SUBROUTINE close(actual,wanted,label)
    REAL*8, INTENT(IN) :: actual(:),wanted(:)
    CHARACTER(*), INTENT(IN) :: label
    REAL*8 :: error
    error=MAXVAL(ABS(actual-wanted))/MAX(1.d0,MAXVAL(ABS(wanted)))
    worst=MAX(worst,error)
    IF(.NOT.(error<2.d-12))THEN
      WRITE(*,*)label,': rank/error=',MPIvar%glob_id,error
      ERROR STOP 'HDG residual comparison'
    ENDIF
  END SUBROUTINE close
END PROGRAM test_hdg_residual
