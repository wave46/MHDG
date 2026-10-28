PROGRAM test_stabilization
  USE extracted_faces
  USE neutral_flux_limiter, ONLY: neutral_tn_source_ti
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE
  REAL*8, PARAMETER :: perturbation=2.d-7,fd_tolerance=3.d-6,algebra_tolerance=1.d-12
  REAL*8, PARAMETER :: reference_trace(5)=[.8d0,-.15d0,1.6d0,2.2d0,.1d0]
  REAL*8 :: element_state(5),trace_state(5),gradient(10),b(3),n(2),tau(5,5),phi,held_heat_phi(2)
  REAL*8 :: worst
  INTEGER :: checks,states
  CALL initialize_fixture()

  CALL fixture('ordinary', reference_trace, .FALSE.)
  CALL fixture('physical heat limiting', reference_trace, .TRUE.)
  phys%T_fluxlim_maxi=.3d0;phys%T_fluxlim_maxe=.4d0
  CALL fixture('temperature caps',reference_trace,.TRUE.)
  phys%T_fluxlim_maxi=100.d0;phys%T_fluxlim_maxe=100.d0
  CALL fixture('advective velocity branch',[.8d0,3.2d0,8.d0,2.2d0,.1d0],.TRUE.)
  CALL fixture('Mach release',[1.d0,1.43d0,2.2d0,1.6d0,.1d0],.TRUE.)
  b=[.05d0,.05d0,SQRT(.995d0)]
  CALL fixture('grazing angular target',[.8d0,-.15d0,1.6d0,2.2d0,.1d0],.TRUE.)
  b=-b
  CALL fixture('inward normal projection',reference_trace,.TRUE.)
  b=[0.d0,0.d0,1.d0]
  CALL fixture('zero magnetic normal projection',reference_trace,.TRUE.)
  b=[.2d0,.4d0,SQRT(.8d0)]
#ifdef CONSTANTNEUTRALDIFF
  phys%diff_nn=.5d0*phys%diff_nn_min
  CALL fixture('supplied constant neutral diffusion below minimum',reference_trace,.TRUE.)
#endif
  CALL check_coefficient_corners()
  CALL check_scatter_maps()
  WRITE(*,*) 'Stabilization PASS; fixtures, checks, worst scaled discrepancy:',states,checks,worst
CONTAINS
  SUBROUTINE initialize_fixture()
    call set_model_layout()
    simpar%ndim=2;simpar%neq=5
    phys%Mref=1.3d0;phys%a=1.d0
    phys%diff_n=.3d0;phys%diff_u=.4d0;phys%diff_e=.5d0;phys%diff_ee=.6d0
    phys%diff_nn=100.d0;phys%diff_nn_min=1.d-4;phys%neutralp_ti_supp=.01d0
    phys%epn=2.5d0;phys%diff_pari=.4d0;phys%diff_pare=.7d0
    phys%Gmbohm=2.5d0;phys%Gmbohme=4.5d0
    phys%bohmth=.3d0;phys%bohm_energy_thresh=.1d0
    phys%c_fli=.2d0;phys%c_fle=.2d0;phys%T_fluxlim_maxi=100.d0;phys%T_fluxlim_maxe=100.d0
    simpar%refval_temperature=50.d0;simpar%refval_density=1.d19;simpar%refval_time=1.d-6
    call initialize_neutral_rate_runtime_constants()
    phys%alpha_iz=0.d0;phys%alpha_iz(1,1)=log(4.d5);phys%alpha_iz(2,1)=.15d0;phys%alpha_iz(1,2)=.03d0
    phys%alpha_cx=0.d0;phys%alpha_cx(1)=log(5.d5);phys%alpha_cx(2)=.2d0
    phys%neutral_flux_limiter_tn_source_id=neutral_tn_source_ti
    phys%neutral_flux_limiter_tn=0.d0;phys%neutral_flux_limiter_eps=0.d0
    phys%neutral_flux_limiter_fs_fraction=.2d0;phys%neutral_flux_limiter_fs_flux_min=0.d0
    call initialize_neutral_flux_limiter_config()
    phys%Re=.8d0;phys%Re_n=.7d0;phys%apply_trim=.false.;phys%bcflags(1)=bc_Bohm
    phys%recycle_plasma_diffusion=.false.
    switch%neutral_wall_sources_in_elements=.false.
    switch%testcase=54;switch%ME=.false.;switch%difcor=0;switch%transport_1d=.false.
    switch%decoup=.false.
    coefi=phys%diff_pari*(2.d0/(3.d0*phys%Mref))**(1.d0+phys%epn)
    coefe=phys%diff_pare*(2.d0/(3.d0*phys%Mref))**(1.d0+phys%epn)
    allocate(elMat%Auu(5,5,1),elMat%Aul(5,5,1),elMat%Auq(5,10,1),elMat%S(5,1))
    allocate(elMat%Alu(5,5,1),elMat%All(5,5,1),elMat%Alq(5,10,1),elMat%fh(5,1))

    allocate(Mesh%elemSize(1));Mesh%elemSize=.7d0;refElPol%ndeg=3;numer%stab=5
    switch%neutral_perpendicular_diffusion=.false.;numer%neutralp_lambda=.2d0
    switch%neutral_flux_limiter=.false.
    b=[.2d0,.4d0,sqrt(.8d0)];n=[.6d0,.8d0];phi=1.d0;checks=0;states=0;worst=0
  END SUBROUTINE initialize_fixture

  SUBROUTINE fixture(label,state,limited)
    CHARACTER(*), INTENT(IN) :: label
    REAL*8, INTENT(IN) :: state(5)
    LOGICAL, INTENT(IN) :: limited
    trace_state=state
    element_state=[1.2d0,.2d0,2.1d0,3.6d0,.3d0]
    gradient=[.2d0,-.1d0,.05d0,.12d0,.3d0,-.25d0,.04d0,.08d0,.5d0,-.4d0]
    switch%flux_limiter=limited
    CALL freeze_heat_limiter()
    states=states+1
    CALL check_coefficients(label)
    CALL check_faces(label,.TRUE.)
    CALL check_faces(label,.FALSE.)
    CALL check_prescribed_trace(label)
  END SUBROUTINE fixture

  SUBROUTINE freeze_heat_limiter()
    REAL*8 :: shi,she,face_gradient(2,5)
    face_gradient=RESHAPE(gradient,[2,5])
    CALL compute_spitzer_harm_flux_ions(trace_state,face_gradient,b,shi)
    CALL compute_spitzer_harm_flux_electrons(trace_state,face_gradient,b,she)
    CALL compute_heat_flux_limiters(trace_state,shi,she,held_heat_phi(1),held_heat_phi(2))
  END SUBROUTINE freeze_heat_limiter

  SUBROUTINE check_coefficients(label)
    CHARACTER(*), INTENT(IN) :: label
    REAL*8 :: derivative(5,5),tp(5,5),tm(5,5),trial(5),scale,heat(2),di(5,5),da(5,5),adv
    INTEGER :: i,j
    CALL get_tau(trace_state,gradient,tau,derivative)
    CALL thermal_parts(trace_state,heat)
    CALL coefficients(trace_state,1.d0,di,da)
    adv=ABS(DOT_PRODUCT(b(1:2),n))*MAX( &
         ABS(trace_state(2)/trace_state(1))+SQRT(ABS( &
         10.d0/9.d0*(trace_state(3)+trace_state(4))/trace_state(1) &
         -5.d0/9.d0*(trace_state(2)/trace_state(1))**2)), &
         5.d0/3.d0*ABS(trace_state(2)/trace_state(1)))
    DO i=1,5
      CALL close(tau(i,i),adv+di(i,i)*refElPol%ndeg/Mesh%elemSize(1)+ &
           SUM(heat, MASK=[i==3,i==4]),label//' coefficient value')
    ENDDO
    DO j=1,5
      scale=MAX(ABS(trace_state(j)),1.d-4)
      trial=trace_state;trial(j)=trace_state(j)+perturbation*scale;CALL get_tau(trial,gradient,tp)
      trial(j)=trace_state(j)-perturbation*scale;CALL get_tau(trial,gradient,tm)
      DO i=1,5
        CALL close((tp(i,i)-tm(i,i))/(2*perturbation),derivative(i,j)*scale,label//' coefficient derivative',fd_tolerance)
      ENDDO
    ENDDO
    CALL close(DOT_PRODUCT(derivative(1,:),trace_state),0.d0,label//' advective degree zero')
    DO i=3,4
      CALL close(DOT_PRODUCT(derivative(i,:),trace_state),-heat(i-2),label//' thermal degree minus one')
    ENDDO
    ! Q and physical limiter switches do not enter the coefficient.
    CALL get_tau(trace_state,100*gradient,tp)
    IF(ANY(tp/=tau)) ERROR STOP 'gradient entered stabilization'
    switch%flux_limiter=.NOT.switch%flux_limiter
    CALL get_tau(trace_state,gradient,tp)
    IF(ANY(tp/=tau)) ERROR STOP 'heat limiter entered stabilization'
    switch%flux_limiter=.NOT.switch%flux_limiter
  END SUBROUTINE check_coefficients

  SUBROUTINE check_faces(label,retain_diffusion)
    CHARACTER(*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: retain_diffusion
    REAL*8 :: r(5),expected(5),ju(5,5),jh(5,5),jq(5,10),jacobian(5,20)
    REAL*8 :: x(20),trial(20),plus(5),minus(5),scale
    INTEGER :: kind,j,i
    phys%bohm_momentum_always_include_diffusion=retain_diffusion
    phi=.4d0
    DO kind=1,3
      CALL assemble(kind,element_state,trace_state,gradient,phi,r,ju,jh,jq)
      CALL evaluate_face_residual(kind,element_state,trace_state,gradient,phi,expected)
      DO i=1,5
        CALL close(r(i),expected(i),label//' physical residual reconstruction')
      ENDDO
      jacobian(:,1:5)=ju;jacobian(:,6:10)=jh;jacobian(:,11:20)=jq
      x=[element_state,trace_state,gradient]
      DO j=1,20
        scale=MAX(ABS(x(j)),1.d-3)
        trial=x;trial(j)=x(j)+perturbation*scale
        CALL evaluate_face_residual(kind,trial(1:5),trial(6:10),trial(11:20),phi,plus)
        trial(j)=x(j)-perturbation*scale
        CALL evaluate_face_residual(kind,trial(1:5),trial(6:10),trial(11:20),phi,minus)
        DO i=1,5
          CALL close((plus(i)-minus(i))/(2*perturbation),jacobian(i,j)*scale,label//' assembled column',fd_tolerance)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE check_faces

  SUBROUTINE check_prescribed_trace(label)
    CHARACTER(*), INTENT(IN) :: label
    REAL*8 :: dt(5,5),di(5,5),da(5,5),base(5,5),baseq(5,10),rhs(5)
    INTEGER :: offsets(1),map(5),gradmap(10),i
    TYPE(balance_accumulator_type) :: diagnostics
    offsets=0;map=[1,2,3,4,5];gradmap=[(i,i=1,10)]
    CALL get_tau(trace_state,gradient,tau,dt);CALL coefficients(trace_state,phi,di,da)
    elMat%Auu=0;elMat%Aul=0;elMat%Auq=0;elMat%S=0
    CALL assemblyExtFacesContribution(1,.TRUE.,offsets,offsets,map,map,gradmap,b,1.d0,n,di,da,phi, &
         RESHAPE([1.d0],[1,1]),[1.d0],[DOT_PRODUCT(b(1:2),n)],element_state,trace_state,gradient,tau,0,n,.FALSE.,diagnostics)
    base=elMat%Auu(:,:,1);baseq=elMat%Auq(:,:,1);rhs=elMat%S(:,1)
    elMat%Auu=0;elMat%Aul=0;elMat%Auq=0;elMat%S=0
    CALL assemblyExtFacesContribution(1,.TRUE.,offsets,offsets,map,map,gradmap,b,1.d0,n,di,da,phi, &
         RESHAPE([1.d0],[1,1]),[1.d0],[DOT_PRODUCT(b(1:2),n)],element_state,trace_state,gradient,tau,0,n,.FALSE.,diagnostics,dtau_du=dt)
    IF(ANY(elMat%Auu(:,:,1)/=base).OR.ANY(elMat%Auq(:,:,1)/=baseq).OR.ANY(elMat%S(:,1)/=rhs)) &
         ERROR STOP 'prescribed trace gained coefficient derivatives'
    IF(ANY(elMat%Aul/=0.d0)) ERROR STOP 'prescribed trace gained unknown columns'
    CALL close(SUM(elMat%S(:,1)),SUM(rhs),label//' prescribed trace')
  END SUBROUTINE check_prescribed_trace

  SUBROUTINE check_coefficient_corners()
    REAL*8 :: dt(5,5),trial(5),tp(5,5),tm(5,5),heat(2),step,slope_plus,slope_minus
    REAL*8 :: conserved(1,5),physical(1,phys%npv)
    INTEGER :: row,side
    b=[1.d0,0.d0,0.d0];n=[1.d0,0.d0]
    phys%T_fluxlim_maxi=.3d0;phys%T_fluxlim_maxe=.4d0
    DO row=3,4
      DO side=-1,1,2
        trace_state=[1.d0,0.d0,1.5d0*phys%Mref*.24d0,1.5d0*phys%Mref*.32d0,.1d0]
        trace_state(row)=1.5d0*phys%Mref*MERGE(.3d0,.4d0,row==3)*(1+side*1.d-4)
        CALL check_coefficients('either side of temperature cap')
      ENDDO
    ENDDO
    DO row=3,4
      trace_state=[1.d0,0.d0,.5d0,.6d0,.1d0]
      conserved(1,:)=trace_state;CALL cons2phys(conserved,physical)
      IF(row==3)phys%T_fluxlim_maxi=physical(1,7)
      IF(row==4)phys%T_fluxlim_maxe=physical(1,8)
      CALL get_tau(trace_state,gradient,tau,dt)
      step=1.d-7;trial=trace_state;trial(row)=trace_state(row)-step
      CALL get_tau(trial,gradient,tm)
      CALL close(dt(row,row),(tau(row,row)-tm(row,row))/step,'one-sided derivative at cap',fd_tolerance)
    ENDDO
    trace_state=[1.d0,3.d0,6.d0,2.1d0,.1d0]
    CALL get_tau(trace_state,gradient,tau,dt)
    step=1.d-7;trial=trace_state;trial(4)=trace_state(4)+step;CALL get_tau(trial,gradient,tp)
    trial(4)=trace_state(4)-step;CALL get_tau(trial,gradient,tm)
    slope_plus=(tp(1,1)-tau(1,1))/step;slope_minus=(tau(1,1)-tm(1,1))/step
    IF(ABS(slope_plus-slope_minus)<.1d0) ERROR STOP 'fixture missed spectral MAX tie'
    IF(MIN(ABS(dt(1,4)-slope_plus),ABS(dt(1,4)-slope_minus))>1.d-5) ERROR STOP 'spectral MAX tie derivative'
    ! Subtract the common advective derivative to isolate thermal floor behavior.
    DO row=3,4
      trace_state=[1.d0,0.d0,1.d0,1.d0,.1d0];trace_state(row)=-1.d-6
      CALL get_tau(trace_state,gradient,tau,dt)
      CALL thermal_parts(trace_state,heat)
      CALL close(DOT_PRODUCT(dt(row,:)-dt(1,:),trace_state),-heat(row-2),'floor density contribution')
      CALL close(dt(row,row)-dt(1,row),0.d0,'floored temperature derivative')
    ENDDO
    trace_state=[1.d0,.2d0,-1.d0,-1.d0,.1d0]
    CALL check_coefficients('negative spectral radicand')
    trace_state=[1.d-21,3.d-21,1.d-20,1.d-20,1.d-22]
    CALL get_tau(trace_state,gradient,tau,dt)
    IF(.NOT.ALL(ieee_is_finite(dt))) ERROR STOP 'density floor produced nonfinite derivative'
    CALL thermal_parts(trace_state,heat)
    DO row=3,4
      CALL close((dt(row,1)-dt(1,1))*trace_state(1),-heat(row-2),'density floor derivative')
    ENDDO
    trace_state=[1.d0,0.d0,1.d0,-1.d0,.1d0]
    CALL get_tau(trace_state,gradient,tau,dt)
    IF(ANY(dt(1,:)/=0.d0)) ERROR STOP 'chi zero convention'
  END SUBROUTINE check_coefficient_corners

  SUBROUTINE check_scatter_maps()
    REAL*8 :: matrix(12,13),rhs(12),derivative(5,5),trace(5),jump(5),shape(2),nn(2,2),expected
    INTEGER :: rows(10),cols(10),offsets(2),i,j,a,c
    matrix=0;rhs=0;trace=[1.d0,2.d0,3.d0,4.d0,5.d0];jump=trace*.3d0
    derivative=RESHAPE([(i*.01d0,i=1,25)],[5,5]);shape=[.3d0,.7d0]
    nn=SPREAD(shape,2,2)*SPREAD(shape,1,2)
    rows=[10,9,8,7,6,5,4,3,2,1];cols=[2,4,6,8,10,1,3,5,7,9];offsets=[0,5]
    CALL add_tau_derivative(matrix,rhs,rows,cols,offsets,nn,shape,trace,jump,derivative)
    DO i=1,5
      DO a=1,2
        CALL close(rhs(rows(offsets(a)+i)),jump(i)*DOT_PRODUCT(derivative(i,:),trace)*shape(a),'affine mapped scatter')
        DO j=1,5
          DO c=1,2
            expected=jump(i)*derivative(i,j)*nn(a,c)
            CALL close(matrix(rows(offsets(a)+i),cols(offsets(c)+j)),expected,'matrix mapped scatter')
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    IF(ANY(matrix(11:,:)/=0.d0).OR.ANY(matrix(:,11:)/=0.d0).OR.ANY(rhs(11:)/=0.d0)) &
         ERROR STOP 'scatter changed unmapped entries'
  END SUBROUTINE check_scatter_maps

  SUBROUTINE close(actual,expected,label,tolerance)
    REAL*8, INTENT(IN) :: actual,expected
    CHARACTER(*), INTENT(IN) :: label
    REAL*8, OPTIONAL, INTENT(IN) :: tolerance
    REAL*8 :: error,allowed
    allowed=algebra_tolerance
    IF(PRESENT(tolerance)) allowed=tolerance
    error=ABS(actual-expected)/MAX(1.d0,ABS(actual),ABS(expected))
    checks=checks+1;worst=MAX(worst,error)
    IF(.NOT.(error<allowed))THEN
      WRITE(*,*) label,': actual, expected, error:',actual,expected,error
      ERROR STOP 'stabilization check'
    ENDIF
  END SUBROUTINE close

  SUBROUTINE get_tau(state,gradient,t,derivative)
    REAL*8, INTENT(IN) :: state(5),gradient(10)
    REAL*8, INTENT(OUT) :: t(5,5)
    REAL*8, OPTIONAL, INTENT(OUT) :: derivative(:,:)
    REAL*8 :: conserved(1,5),physical(1,phys%npv),di(5,5),da(5,5)
    conserved(1,:)=state;CALL cons2phys(conserved,physical)
    CALL coefficients(state,1.d0,di,da)
    CALL computeTauGaussPoints(physical(1,:),state,gradient,b(1:2),n,1,1.d0,[1.d0,1.d0],t,di,derivative)
  END SUBROUTINE get_tau

  SUBROUTINE thermal_parts(state,heat)
    REAL*8, INTENT(IN) :: state(5)
    REAL*8, INTENT(OUT) :: heat(2)
    REAL*8 :: conserved(1,5),physical(1,phys%npv),bn
    conserved(1,:)=state;CALL cons2phys(conserved,physical);bn=DOT_PRODUCT(b(1:2),n)
    heat=ABS(bn)*[phys%diff_pari,phys%diff_pare]* &
         MIN([phys%T_fluxlim_maxi,phys%T_fluxlim_maxe],physical(1,7:8))**2.5d0* &
         NORM2(b(1:2))/state(1)*refElPol%ndeg/Mesh%elemSize(1)
  END SUBROUTINE thermal_parts
  SUBROUTINE coefficients(state,neutral_phi,iso,ani)
    REAL*8, INTENT(IN) :: state(5),neutral_phi
    REAL*8, INTENT(OUT) :: iso(5,5),ani(5,5)
    REAL*8 :: conserved(1,5),di(5,5,1),da(5,5,1),xy(1,2),qc(1)
    conserved(1,:)=state;xy=1.d0;qc=1.d0
    CALL setLocalDiff(xy,conserved,di,da,qc)
    di(5,5,:)=neutral_phi*di(5,5,:)
    da(5,5,:)=neutral_phi*da(5,5,:)
    iso=di(:,:,1);ani=da(:,:,1)
  END SUBROUTINE coefficients

  SUBROUTINE assemble(kind,element,trace,gradient,neutral_phi,residual,ju,jh,jq)
    INTEGER, INTENT(IN) :: kind
    REAL*8, INTENT(IN) :: element(5),trace(5),gradient(10),neutral_phi
    REAL*8, INTENT(OUT) :: residual(5),ju(5,5),jh(5,5),jq(5,10)
    INTEGER :: offsets(1),map(5),gradmap(10),i
    REAL*8 :: di(5,5),da(5,5),nn(1,1),shape(1),normal_flux(1),dt(5,5)
    REAL*8 :: speed,ds(5),physical(1,phys%npv),conserved(1,5)
    REAL*8 :: p0,ps,pf,dp0(5),dps(5),dpf(5),active,dactive(5),target,dtarget(5)
    TYPE(balance_accumulator_type) :: diagnostics

    offsets=0;map=[1,2,3,4,5];gradmap=[(i,i=1,10)]
    nn=1;shape=1;normal_flux=DOT_PRODUCT(b(1:2),n)
    CALL coefficients(trace,neutral_phi,di,da)
    CALL get_tau(trace,gradient,tau,dt)
    elMat%Auu=0;elMat%Aul=0;elMat%Auq=0;elMat%S=0
    elMat%Alu=0;elMat%All=0;elMat%Alq=0;elMat%fh=0
    SELECT CASE(kind)
    CASE(1)
      CALL assemblyIntFacesContribution(1,offsets,offsets,map,map,gradmap,b,1.d0,n,di,da, &
           neutral_phi,nn,shape,normal_flux,trace,gradient,tau,0,n,uef=element,dtau_du=dt)
    CASE(2)
      CALL assemblyExtFacesContribution(1,.FALSE.,offsets,offsets,map,map,gradmap,b,1.d0,n,di,da, &
           neutral_phi,nn,shape,normal_flux,element,trace,gradient,tau,0,n,.FALSE.,diagnostics,dtau_du=dt)
    CASE(3)
      conserved(1,:)=trace
      CALL cons2phys(conserved,physical)
      CALL compute_cs(trace,speed)
      CALL compute_dcs_du(trace,ds)
      CALL evaluate_bohm_blend(trace,speed,ds,normal_flux(1),phys%bohmth,phys%bohm_energy_thresh, &
           p0,ps,pf,dp0,dps,dpf)
      active=p0+ps;dactive=dp0+dps
      target=SIGN(1.d0,normal_flux(1))*speed*ps
      dtarget=SIGN(1.d0,normal_flux(1))*(ps*ds+speed*dps)
      CALL assembly_bohm_bc(1,offsets,offsets,map,map,gradmap,nn,shape,gradient,trace,physical(1,:), &
           element,b(1:2),1.d0,n,tau,target,dtarget,active,pf,dactive,di,da,neutral_phi,1.d0,.TRUE.,0,n,dtau_du=dt)
    END SELECT
    IF(kind==3)THEN
      ju=elMat%Alu(:,:,1);jh=elMat%All(:,:,1);jq=elMat%Alq(:,:,1)
      residual=MATMUL(ju,element)+MATMUL(jh,trace)+MATMUL(jq,gradient)-elMat%fh(:,1)
      RETURN
    ENDIF
    ju=elMat%Auu(:,:,1);jh=elMat%Aul(:,:,1);jq=elMat%Auq(:,:,1)
    residual=MATMUL(ju,element)+MATMUL(jh,trace)+MATMUL(jq,gradient)-elMat%S(:,1)
    IF(kind==1)THEN
      IF(MAXVAL(ABS(residual-(MATMUL(elMat%Alu(:,:,1),element)+MATMUL(elMat%All(:,:,1),trace)+ &
           MATMUL(elMat%Alq(:,:,1),gradient)-elMat%fh(:,1))))>1.d-12) ERROR STOP 'element/trace flux differs'
    ENDIF
  END SUBROUTINE assemble

  ! Direct flux values, independent of the matrix/RHS reconstruction.
  SUBROUTINE evaluate_face_residual(kind,element,trace,gradient,neutral_phi,residual)
    INTEGER, INTENT(IN) :: kind
    REAL*8, INTENT(IN) :: element(5),trace(5),gradient(10),neutral_phi
    REAL*8, INTENT(OUT) :: residual(5)
    REAL*8 :: di(5,5),da(5,5),grad(2,5),w(5),g(5),f(5),pnormal(2),pn(2),bn
    REAL*8 :: shi,she,cs,target(5),wall_flux,beta,angle,activation,z,ps,pf
    INTEGER :: i
    CALL get_tau(trace,gradient,tau)
    grad=RESHAPE(gradient,[2,5]);bn=DOT_PRODUCT(b(1:2),n);pnormal=n-bn*b(1:2)
    CALL coefficients(trace,neutral_phi,di,da)
    DO i=1,5
      g(i)=di(i,i)*DOT_PRODUCT(grad(:,i),n)-da(i,i)*bn*DOT_PRODUCT(grad(:,i),b(1:2))
    ENDDO
    CALL compute_W2(trace,w,di(1,1),di(2,2))
    g(2)=g(2)+DOT_PRODUCT(MATMUL(grad,w),pnormal)
    CALL compute_W3(trace,w,di(1,1),di(2,2),di(3,3))
    g(3)=g(3)+DOT_PRODUCT(MATMUL(grad,w),pnormal)
    CALL compute_W4(trace,w,di(1,1),di(4,4))
    g(4)=g(4)+DOT_PRODUCT(MATMUL(grad,w),pnormal)
    pn=n
    IF(switch%neutral_perpendicular_diffusion) pn=pnormal
    CALL compute_W5p(trace,w)
    g(5)=g(5)+neutral_phi*DOT_PRODUCT(MATMUL(grad,w),pn)
    CALL compute_spitzer_harm_flux_ions(trace,grad,b,shi)
    CALL compute_spitzer_harm_flux_electrons(trace,grad,b,she)
    g(3)=g(3)-bn*held_heat_phi(1)*shi
    g(4)=g(4)-bn*held_heat_phi(2)*she
    IF(kind<3)THEN
      f=0
      f(1)=trace(2)
      f(2)=2.d0/3*(trace(2)**2/trace(1)+trace(3)+trace(4))
      f(3)=5.d0/3*trace(2)*trace(3)/trace(1)-1.d0/3*trace(2)**3/trace(1)**2
      f(4)=5.d0/3*trace(2)*trace(4)/trace(1)
      residual=bn*f-g+MATMUL(tau,element-trace)
      RETURN
    ENDIF
    target=element
    CALL compute_cs(trace,cs)
    beta=MAX(trace(3),1.d-20)/(MAX(trace(3),1.d-20)+phys%bohm_energy_thresh)
    z=MAX(0.d0,MIN(1.d0,ABS(bn)/phys%bohmth));angle=z*z*(3.d0-2.d0*z)
    z=MAX(0.d0,MIN(1.d0,(SIGN(1.d0,bn)*trace(2)/(trace(1)*cs)-1.d0/beta)/.20d0))
    activation=1.d0-z*z*(3.d0-2.d0*z)
    ps=activation*beta*angle;pf=1.d0-activation
    target(2)=ps*SIGN(1.d0,bn)*cs*trace(1)+pf*element(2)
    IF(.NOT.phys%bohm_momentum_always_include_diffusion) g(2)=pf*g(2)
    residual=MATMUL(tau,target-trace)-g
    residual(3)=residual(3)+bn*(5.d0-2*phys%Gmbohm)/3*trace(2)/trace(1)* &
         (trace(3)-.5d0*trace(2)**2/trace(1))
    residual(4)=residual(4)+bn*(5.d0-2*phys%Gmbohme)/3*trace(2)/trace(1)*trace(4)
    CALL compute_neutral_wall_flux(trace,phys%Re_n,wall_flux)
    residual(5)=residual(5)+phys%Re*bn*trace(2)-wall_flux
  END SUBROUTINE evaluate_face_residual
END PROGRAM test_stabilization
