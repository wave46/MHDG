MODULE hdg_residual_module
  USE globals, ONLY: Mesh,phys,refElPol,elMat,sol,utils,numer,switch
  USE types, ONLY: Elmat_type,Sol_type
  USE LinearAlgebra, ONLY: invert_matrix
  USE solver_status, ONLY: finite_solver_values
  USE MPI_OMP, ONLY: MPIvar
#ifdef PARALL
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hdg_residual_t,hdg_residual_norms_t,hdg_residual_monitor_t,hdg_residual
  PUBLIC :: evaluate_element_residual,evaluate_hdg_residual

  TYPE hdg_residual_t
    REAL*8, ALLOCATABLE :: u(:,:),q(:,:),trace(:,:)
    REAL*8 :: evaluation_seconds=0.d0
  END TYPE hdg_residual_t

  TYPE hdg_residual_norms_t
    REAL*8, ALLOCATABLE :: rms(:,:),local_max(:,:)
    REAL*8 :: total_rms=0.d0,relative=1.d0,initial=0.d0,target=0.d0
    REAL*8 :: evaluation_seconds=0.d0
    REAL*8 :: measures(3)=0.d0
  END TYPE hdg_residual_norms_t

  TYPE hdg_residual_monitor_t
    PRIVATE
    INTEGER :: report_mode=0,samples=0
    TYPE(hdg_residual_norms_t) :: norms
    REAL*8, ALLOCATABLE :: volume_weights(:,:),face_inverse_mass(:,:,:),face_measure(:)
  CONTAINS
    PROCEDURE :: configure => configure_residual
    PROCEDURE :: get_mode => residual_report_mode
    PROCEDURE :: begin_solve => begin_residual_solve
    PROCEDURE :: invalidate_mesh => invalidate_residual_mesh
    PROCEDURE :: evaluate => evaluate_residual_monitor
    PROCEDURE :: measure => measure_residual
    PROCEDURE :: converged => residual_converged
    PROCEDURE :: report => report_residual
    PROCEDURE, PRIVATE :: prepare_geometry => prepare_residual_geometry
  END TYPE hdg_residual_monitor_t

  TYPE(hdg_residual_monitor_t) :: hdg_residual

CONTAINS

  PURE SUBROUTINE evaluate_element_residual(matrices,state,element,trace_local,trace_global,ru,rq,rtrace)
    TYPE(Elmat_type), INTENT(IN) :: matrices
    TYPE(Sol_type), INTENT(IN) :: state
    INTEGER, INTENT(IN) :: element
    REAL*8, INTENT(IN) :: trace_local(:),trace_global(:)
    REAL*8, INTENT(OUT) :: ru(:),rq(:),rtrace(:)
    INTEGER :: nv,nq,offset_u,offset_q

    nv=SIZE(matrices%Auu,1);nq=SIZE(matrices%Aqu,1)
    offset_u=(element-1)*nv;offset_q=(element-1)*nq
    ASSOCIATE(u=>state%u(offset_u+1:offset_u+nv),q=>state%q(offset_q+1:offset_q+nq))
      ru=MATMUL(matrices%Auu(:,:,element),u)+MATMUL(matrices%Auq(:,:,element),q)+ &
           MATMUL(matrices%Aul(:,:,element),trace_local)+matrices%Aul_dir(:,element)-matrices%S(:,element)
      ! iAqq already converts the gradient load to a gradient-value residual.
      rq=q+MATMUL(matrices%iAqq(:,:,element),MATMUL(matrices%Aqu(:,:,element),u)+ &
           MATMUL(matrices%Aql(:,:,element),trace_local)+matrices%Aql_dir(:,element))
      rtrace=MATMUL(matrices%Alu(:,:,element),u)+MATMUL(matrices%Alq(:,:,element),q)+ &
           MATMUL(matrices%All(:,:,element),trace_global)-matrices%fh(:,element)
    END ASSOCIATE
  END SUBROUTINE evaluate_element_residual

  SUBROUTINE evaluate_hdg_residual(residual)
    TYPE(hdg_residual_t), INTENT(OUT) :: residual
    INTEGER :: nv,nq,block,nf,e,f,k,j,face,offset,neq
    INTEGER :: permutation(Mesh%Nnodesperface*phys%neq)
    REAL*8 :: trace_local(refElPol%Nfaces*Mesh%Nnodesperface*phys%neq)
    REAL*8 :: trace_global(SIZE(trace_local))
    INTEGER(KIND=8) :: clock_start,clock_end,clock_rate

    IF(utils%timing)CALL SYSTEM_CLOCK(clock_start,clock_rate)
    neq=phys%neq
    nv=neq*Mesh%Nnodesperelem;nq=Mesh%Ndim*nv
    block=neq*Mesh%Nnodesperface;nf=refElPol%Nfaces*block
    ALLOCATE(residual%u(nv,Mesh%Nelems),residual%q(nq,Mesh%Nelems),residual%trace(nf,Mesh%Nelems))
    DO k=1,Mesh%Nnodesperface
      DO j=1,neq
        permutation((k-1)*neq+j)=(Mesh%Nnodesperface-k)*neq+j
      ENDDO
    ENDDO
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(e,f,face,offset,trace_local,trace_global)
    DO e=1,Mesh%Nelems
      DO f=1,refElPol%Nfaces
        face=Mesh%F(e,f);offset=(f-1)*block
        ! Prescribed states are already carried by Aul_dir/Aql_dir, not unknown columns.
        trace_global(offset+1:offset+block)=0.d0
        IF(.NOT.Mesh%Fdir(e,f)) &
             trace_global(offset+1:offset+block)=sol%u_tilde((face-1)*block+1:face*block)
        trace_local(offset+1:offset+block)=trace_global(offset+1:offset+block)
        IF(Mesh%flipface(e,f))trace_local(offset+1:offset+block)=trace_global(offset+permutation)
      ENDDO
      CALL evaluate_element_residual(elMat,sol,e,trace_local,trace_global, &
           ru=residual%u(:,e),rq=residual%q(:,e),rtrace=residual%trace(:,e))
    ENDDO
    !$OMP END PARALLEL DO
    IF(utils%timing)THEN
      CALL SYSTEM_CLOCK(clock_end)
      residual%evaluation_seconds=REAL(clock_end-clock_start,8)/REAL(clock_rate,8)
    ENDIF
  END SUBROUTINE evaluate_hdg_residual

  SUBROUTINE configure_residual(this,mode)
    CLASS(hdg_residual_monitor_t), INTENT(INOUT) :: this
    CHARACTER(*), INTENT(IN) :: mode
    SELECT CASE(TRIM(ADJUSTL(mode)))
    CASE('off')
      this%report_mode=0
    CASE('summary')
      this%report_mode=1
    CASE('equations')
      this%report_mode=2
    CASE DEFAULT
      ERROR STOP 'residual_diagnostics_mode must be off, summary or equations'
    END SELECT
#ifdef TOR3D
    IF(this%report_mode/=0)ERROR STOP 'HDG residual reporting currently supports 2D builds'
#endif
    CALL this%begin_solve()
  END SUBROUTINE configure_residual

  FUNCTION residual_report_mode(this) RESULT(mode)
    CLASS(hdg_residual_monitor_t), INTENT(IN) :: this
    CHARACTER(16) :: mode
    SELECT CASE(this%report_mode)
    CASE(0)
      mode='off'
    CASE(1)
      mode='summary'
    CASE(2)
      mode='equations'
    END SELECT
  END FUNCTION residual_report_mode

  SUBROUTINE begin_residual_solve(this)
    CLASS(hdg_residual_monitor_t), INTENT(INOUT) :: this
    ! New nonlinear problem: retain masses, discard its preceding history.
    this%samples=0
    this%norms=hdg_residual_norms_t()
  END SUBROUTINE begin_residual_solve

  SUBROUTINE invalidate_residual_mesh(this)
    CLASS(hdg_residual_monitor_t), INTENT(INOUT) :: this
    ! Mesh matrices are being rebuilt: their quadrature cache is stale too.
    IF(ALLOCATED(this%volume_weights)) &
         DEALLOCATE(this%volume_weights,this%face_inverse_mass,this%face_measure)
    CALL this%begin_solve()
  END SUBROUTINE invalidate_residual_mesh

  SUBROUTINE evaluate_residual_monitor(this)
    CLASS(hdg_residual_monitor_t), INTENT(INOUT) :: this
    TYPE(hdg_residual_t) :: residual
    IF(this%report_mode==0.AND.TRIM(numer%newton_convergence)=='update')RETURN
    CALL evaluate_hdg_residual(residual)
    CALL this%measure(residual)
  END SUBROUTINE evaluate_residual_monitor

  SUBROUTINE measure_residual(this,residual,norms)
    CLASS(hdg_residual_monitor_t), INTENT(INOUT) :: this
    TYPE(hdg_residual_t), INTENT(IN) :: residual
    TYPE(hdg_residual_norms_t), OPTIONAL, INTENT(OUT) :: norms
    REAL*8 :: initial
    TYPE(hdg_residual_norms_t) :: measured
    initial=this%norms%initial
    CALL compute_residual_norms(this,residual,measured)
    this%norms=measured
    IF(this%samples==0)initial=this%norms%total_rms
    this%samples=this%samples+1
    this%norms%initial=initial
    this%norms%relative=0.d0
    IF(initial>0.d0)THEN
      this%norms%relative=this%norms%total_rms/initial
    ELSEIF(this%norms%total_rms>0.d0)THEN
      this%norms%relative=HUGE(1.d0)
    ENDIF
    this%norms%target=numer%residual_atol+numer%residual_rtol*initial
    IF(PRESENT(norms))norms=this%norms
  END SUBROUTINE measure_residual

  LOGICAL FUNCTION residual_converged(this,update_error)
    CLASS(hdg_residual_monitor_t), INTENT(IN) :: this
    REAL*8, INTENT(IN) :: update_error
    residual_converged=.FALSE.
    IF(TRIM(numer%newton_convergence)=='update'.OR.this%samples==0)RETURN
    IF(.NOT.finite_solver_values([this%norms%total_rms,this%norms%initial]))RETURN
    residual_converged=this%norms%total_rms<=this%norms%target
    IF(TRIM(numer%newton_convergence)=='both') &
         residual_converged=residual_converged.AND.this%samples>1.AND.update_error<numer%tNR
  END FUNCTION residual_converged

  SUBROUTINE report_residual(this,output_unit)
    CLASS(hdg_residual_monitor_t), INTENT(IN) :: this
    INTEGER, OPTIONAL, INTENT(IN) :: output_unit
    IF(this%report_mode==0.OR.this%samples==0)RETURN
    CALL print_residual_norms(this%norms,this%report_mode==2,output_unit)
  END SUBROUTINE report_residual

  SUBROUTINE prepare_residual_geometry(this)
    CLASS(hdg_residual_monitor_t), INTENT(INOUT) :: this
    INTEGER :: e,f,face,g,nfp
    REAL*8 :: x(Mesh%Nnodesperelem,2),dx(refElPol%Ngauss2d,2),dy(refElPol%Ngauss2d,2)
    REAL*8 :: xf(Mesh%Nnodesperface,2),tangent(refElPol%Ngauss1d,2)
    REAL*8 :: weights(refElPol%Ngauss1d),mass(Mesh%Nnodesperface,Mesh%Nnodesperface)
    LOGICAL :: seen(Mesh%Nfaces)

    IF(ALLOCATED(this%volume_weights))RETURN
    nfp=Mesh%Nnodesperface
    ALLOCATE(this%volume_weights(refElPol%Ngauss2d,Mesh%Nelems), &
         this%face_inverse_mass(nfp,nfp,Mesh%Nfaces),this%face_measure(Mesh%Nfaces))
    this%volume_weights=0.d0;this%face_inverse_mass=0.d0;this%face_measure=0.d0;seen=.FALSE.
    DO e=1,Mesh%Nelems
#ifdef PARALL
      IF(Mesh%ghostElems(e)==0)THEN
#endif
        x=Mesh%X(Mesh%T(e,:),:)
        dx=MATMUL(refElPol%Nxi2D,x);dy=MATMUL(refElPol%Neta2D,x)
        this%volume_weights(:,e)=refElPol%gauss_weights2D*(dx(:,1)*dy(:,2)-dx(:,2)*dy(:,1))
        IF(switch%axisym)this%volume_weights(:,e)=this%volume_weights(:,e)*MATMUL(refElPol%N2D,x(:,1))
        IF(SUM(this%volume_weights(:,e))<=0.d0)ERROR STOP 'Nonpositive residual volume measure'
#ifdef PARALL
      ENDIF
#endif
      DO f=1,refElPol%Nfaces
        IF(Mesh%Fdir(e,f))CYCLE
        face=Mesh%F(e,f)
#ifdef PARALL
        IF(Mesh%ghostFaces(face)/=0)CYCLE
#endif
        IF(seen(face))CYCLE
        seen(face)=.TRUE.
        xf=Mesh%X(Mesh%T(e,refElPol%face_nodes(f,:)),:)
        IF(Mesh%flipface(e,f))xf=xf(nfp:1:-1,:)
        tangent=MATMUL(refElPol%Nxi1D,xf)
        weights=refElPol%gauss_weights1D*SQRT(SUM(tangent**2,DIM=2))
        IF(switch%axisym)weights=weights*MATMUL(refElPol%N1D,xf(:,1))
        IF(SUM(weights)<=0.d0)ERROR STOP 'Nonpositive residual face measure'
        mass=0.d0
        DO g=1,refElPol%Ngauss1d
          mass=mass+weights(g)*SPREAD(refElPol%N1D(g,:),2,nfp)*SPREAD(refElPol%N1D(g,:),1,nfp)
        ENDDO
        CALL invert_matrix(mass,this%face_inverse_mass(:,:,face))
        this%face_measure(face)=SUM(weights)
      ENDDO
    ENDDO
  END SUBROUTINE prepare_residual_geometry

  SUBROUTINE compute_residual_norms(this,residual,norms)
    CLASS(hdg_residual_monitor_t), INTENT(INOUT) :: this
    TYPE(hdg_residual_t), INTENT(IN) :: residual
    TYPE(hdg_residual_norms_t), INTENT(OUT) :: norms
    REAL*8 :: sums(3,phys%neq),peak(3,phys%neq),local_squared(2),measure
    REAL*8 :: gradient(refElPol%Ngauss2d),r(Mesh%Nnodesperelem)
    REAL*8 :: inverse_mass(Mesh%Nnodesperelem,Mesh%Nnodesperelem)
    REAL*8, ALLOCATABLE :: face_sum(:,:)
    INTEGER :: e,j,k,nv,nq,neq,ndim,block,face,ierr,active_blocks
    INTEGER :: mass_indices(Mesh%Nnodesperelem)
    INTEGER(KIND=8) :: clock_start,clock_end,clock_rate
    LOGICAL :: owned

    IF(utils%timing)CALL SYSTEM_CLOCK(clock_start,clock_rate)
    CALL this%prepare_geometry()
    neq=phys%neq;ndim=Mesh%Ndim
    nv=SIZE(residual%u,1);nq=SIZE(residual%q,1);block=neq*Mesh%Nnodesperface
    mass_indices=[(1+(k-1)*neq*ndim,k=1,Mesh%Nnodesperelem)]
    ALLOCATE(norms%rms(3,neq),norms%local_max(3,neq),face_sum(block,Mesh%Nfaces))
    sums=0.d0;peak=0.d0;face_sum=0.d0
    DO e=1,Mesh%Nelems
      owned=.TRUE.
#ifdef PARALL
      owned=Mesh%ghostElems(e)==0
#endif
      IF(owned)THEN
        measure=SUM(this%volume_weights(:,e))
        norms%measures(1:2)=norms%measures(1:2)+measure
        ! iAqq repeats the same scalar inverse mass for each equation/direction.
        inverse_mass=elMat%iAqq(mass_indices,mass_indices,e)
        DO j=1,neq
          r=residual%u(j:nv:neq,e)
          local_squared(1)=DOT_PRODUCT(r,MATMUL(inverse_mass,r))
          local_squared(2)=0.d0
          DO k=1,ndim
            gradient=MATMUL(refElPol%N2D,residual%q(k+(j-1)*ndim:nq:neq*ndim,e))
            local_squared(2)=local_squared(2)+SUM(this%volume_weights(:,e)*gradient**2)/ndim
          ENDDO
          sums(1:2,j)=sums(1:2,j)+local_squared
          peak(1:2,j)=MAX(peak(1:2,j),local_squared/measure)
        ENDDO
      ENDIF
      CALL accumulate_trace_rows(e,residual%trace(:,e),face_sum)
    ENDDO
    DO face=1,Mesh%ukf
#ifdef PARALL
      IF(Mesh%ghostFaces(face)/=0)CYCLE
#endif
      norms%measures(3)=norms%measures(3)+this%face_measure(face)
      DO j=1,neq
        measure=DOT_PRODUCT(face_sum(j:block:neq,face), &
             MATMUL(this%face_inverse_mass(:,:,face),face_sum(j:block:neq,face)))
        sums(3,j)=sums(3,j)+measure
        peak(3,j)=MAX(peak(3,j),measure/this%face_measure(face))
      ENDDO
    ENDDO
#ifdef PARALL
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,sums,SIZE(sums),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,norms%measures,3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,peak,SIZE(peak),MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif
    norms%rms=0.d0
    DO k=1,3
      IF(norms%measures(k)>0.d0)norms%rms(k,:)=SQRT(sums(k,:)/norms%measures(k))
    ENDDO
    norms%local_max=SQRT(peak)
    DO j=1,neq
      norms%rms(:,j)=norms%rms(:,j)/numer%residual_scales(j)
      norms%local_max(:,j)=norms%local_max(:,j)/numer%residual_scales(j)
    ENDDO
    active_blocks=COUNT(norms%measures>0.d0)
    IF(active_blocks>0)norms%total_rms=SQRT(SUM(norms%rms**2)/(neq*active_blocks))
    IF(utils%timing)THEN
      CALL SYSTEM_CLOCK(clock_end)
      norms%evaluation_seconds=residual%evaluation_seconds+REAL(clock_end-clock_start,8)/REAL(clock_rate,8)
#ifdef PARALL
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,norms%evaluation_seconds,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif
    ENDIF
  END SUBROUTINE compute_residual_norms

  SUBROUTINE accumulate_trace_rows(element,trace_residual,face_sum)
    INTEGER, INTENT(IN) :: element
    REAL*8, INTENT(IN) :: trace_residual(:)
    REAL*8, INTENT(INOUT) :: face_sum(:,:)
    INTEGER :: f,face,offset,block,partner,other_element,other_face
    REAL*8 :: trace_difference(SIZE(face_sum,1))

    block=SIZE(face_sum,1)
    DO f=1,refElPol%Nfaces
      IF(Mesh%Fdir(element,f))CYCLE
      face=Mesh%F(element,f);offset=(f-1)*block
#ifdef PARALL
      IF(Mesh%ghostFaces(face)/=0)CYCLE
#endif
      face_sum(:,face)=face_sum(:,face)+trace_residual(offset+1:offset+block)
      IF(face<=Mesh%Nintfaces)CYCLE
      partner=Mesh%periodic_faces(face-Mesh%Nintfaces)
      IF(partner==0)CYCLE
      other_element=Mesh%extfaces(partner,1);other_face=Mesh%extfaces(partner,2)
      partner=Mesh%F(other_element,other_face)
      ! Periodic assembly inserts this row twice, replacing its own trace column
      ! by the paired column on the second insertion. Current traces can differ.
      trace_difference=sol%u_tilde((partner-1)*block+1:partner*block)- &
           sol%u_tilde((face-1)*block+1:face*block)
      face_sum(:,partner)=face_sum(:,partner)+trace_residual(offset+1:offset+block)+ &
           MATMUL(elMat%All(offset+1:offset+block,offset+1:offset+block,element),trace_difference)
    ENDDO
  END SUBROUTINE accumulate_trace_rows

  SUBROUTINE print_residual_norms(norms,detailed,output_unit)
    TYPE(hdg_residual_norms_t), INTENT(IN) :: norms
    LOGICAL, INTENT(IN) :: detailed
    INTEGER, OPTIONAL, INTENT(IN) :: output_unit
    INTEGER :: destination,j

    IF(MPIvar%glob_id/=0)RETURN
    destination=6
    IF(PRESENT(output_unit))destination=output_unit
    WRITE(destination,'(A,ES11.3,A,ES11.3,A,ES11.3)') &
         'HDG residual before update: norm=',norms%total_rms, &
         '; relative=',norms%relative,'; target=',norms%target
    WRITE(destination,'(A,3(1X,ES11.3))') &
         '  local peaks (conservation, gradient, trace):',MAXVAL(norms%local_max,DIM=2)
    IF(detailed)THEN
      WRITE(destination,'(A)')'  equation       U RMS          Q RMS        trace RMS         U peak          Q peak       trace peak'
      DO j=1,SIZE(norms%rms,2)
        WRITE(destination,'(2X,A10,6(1X,ES14.5))')TRIM(phys%conVarNam(j)),norms%rms(:,j),norms%local_max(:,j)
      ENDDO
    ENDIF
    IF(utils%timing)WRITE(destination,'(A,ES11.3,A)') &
         '  HDG residual calculation wall time (max rank): ',norms%evaluation_seconds,' s'
  END SUBROUTINE print_residual_norms
END MODULE hdg_residual_module
