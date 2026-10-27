! Compiled by check_bohm_momentum.py with verbatim extracted production blocks.
MODULE extracted_bohm
  USE, INTRINSIC :: iso_c_binding
  USE, INTRINSIC :: ieee_arithmetic
  USE bohm_boundary, ONLY: evaluate_bohm_blend, evaluate_isothermal_bohm_target
  IMPLICIT NONE
#ifdef NEUTRALGAMMA
  INTEGER, PARAMETER :: Neq=6
#else
  INTEGER, PARAMETER :: Neq=5
#endif
  INTEGER, PARAMETER :: Ndim=2, Npfl=1
  TYPE physics_type
    REAL*8 :: bohmth,bohm_energy_thresh
    LOGICAL :: bohm_momentum_always_include_diffusion=.TRUE.
  END TYPE
  TYPE(physics_type) :: phys
  INTEGER, PARAMETER :: equation_nu=2,term_diffusion=1,term_split_diffusion=2,term_tau_stabilization_inward=3
  TYPE balance_accumulator_type
    REAL*8 :: values(3)=0.d0
  END TYPE
! PRODUCTION_DIAGNOSTIC_TYPE
CONTAINS
! PRODUCTION_HELPERS
! PRODUCTION_DIAGNOSTIC_HELPERS

  SUBROUTINE add_bc_value(this,equation,term,value)
    CLASS(balance_accumulator_type),INTENT(INOUT) :: this
    INTEGER,INTENT(IN) :: equation,term
    REAL*8,INTENT(IN) :: value
    this%values(term)=this%values(term)+value
  END SUBROUTINE

  SUBROUTINE select_blend(uhat,bn,delta,free_weight,ddelta_du,target_coefficient,dtarget_du,weights)
    REAL*8, INTENT(IN) :: uhat(Neq),bn
    REAL*8, INTENT(OUT) :: delta,free_weight,ddelta_du(Neq),target_coefficient,dtarget_du(Neq),weights(3)
    REAL*8 :: ufg(1,Neq),upg(1,9),SoundSpeed,inc,sn,bohm_suppresion,energy_zero_threshold,dcs_du(Neq)
    REAL*8 :: p0,ps,pf,dp0(Neq),dps(Neq),dpf(Neq)
    INTEGER :: g
    g=1;ufg(1,:)=uhat;upg=0;upg(1,2)=uhat(2)/uhat(1);inc=bn;energy_zero_threshold=1.d-20
#ifdef TEMPERATURE
    CALL compute_cs(uhat,SoundSpeed)
    CALL compute_dcs_du(uhat,dcs_du)
#else
    SoundSpeed=1.d0;dcs_du=0.d0
#endif
! PRODUCTION_SELECTOR
#ifdef TEMPERATURE
    weights=[p0,ps,pf]
#else
    weights=[0.d0,delta,free_weight]
#endif
  END SUBROUTINE

  SUBROUTINE probe(u,uhat,q,iso,ani,b,n,tau_input,angle,energy,lag_velocity,tangent,retain,ju,jh,jq,rhs,weights,target,diagnostic) BIND(C)
    REAL(c_double), INTENT(IN) :: u(Neq),uhat(Neq),q(2,Neq),iso(Neq,Neq),ani(Neq,Neq),b(2),n(2),tau_input(Neq,Neq)
    REAL(c_double), VALUE :: angle,energy,lag_velocity
    INTEGER(c_int), VALUE :: tangent,retain
    REAL(c_double), INTENT(OUT) :: ju(Neq,Neq),jh(Neq,Neq),jq(Neq,2*Neq),rhs(Neq),weights(3),target,diagnostic(3)
    TYPE matrix_type
      REAL*8 :: All(Neq,Neq,1),Alu(Neq,Neq,1),Alq(Neq,2*Neq,1),fh(Neq,1)
    END TYPE
    TYPE(matrix_type) :: elMat
    TYPE(balance_accumulator_type) :: accumulator
    TYPE(balance_plasma_bc_type) :: plasma
    INTEGER :: i,j,k,idm,iel,Neqstab,Neqgrad,inn,ign
    INTEGER :: ind(1),indi(1),indj(1),ind_jf(1),ind_kf(1),ind_asf(1),ind_ash(1),ind_ff(Neq),ind_fe(Neq),ind_fg(2*Neq)
    REAL*8 :: ufg(Neq),uefg(Neq),upfg(9),Qpr(2,Neq),diffiso(Neq,Neq),diffani(Neq,Neq),bg(2),ng(2),bn,tau(Neq,Neq)
    REAL*8 :: delta,free_weight,ddelta_du(Neq),target_coefficient,dtarget_du(Neq),NiNi(1,1),Ni(1),kmult(1,1)
    REAL*8 :: target_affine,momentum_diffusion_factor,momentum_diffusive_flux
    REAL*8 :: W2(Neq),dW2_dU(Neq,Neq),QdW2(2,Neq)
    REAL*8 :: W3(Neq),dW3_dU(Neq,Neq),QdW3(2,Neq),W4(Neq),dW4_dU(Neq,Neq),QdW4(2,Neq)
    LOGICAL :: ntang
    REAL*8 :: diagnostic_target(Neq)
    phys%bohmth=angle;phys%bohm_energy_thresh=energy
    phys%bohm_momentum_always_include_diffusion=retain/=0
    bn=DOT_PRODUCT(b,n)
    CALL select_blend(uhat,bn,delta,free_weight,ddelta_du,target_coefficient,dtarget_du,weights)
    ufg=uhat;uefg=u
#ifdef TEMPERATURE
! PRODUCTION_DIAGNOSTIC_TARGET
#endif
#ifndef TEMPERATURE
    target=target_coefficient*u(1)+free_weight*u(2)
#endif
    ufg=uhat;uefg=u;upfg=0;upfg(2)=lag_velocity;Qpr=q;diffiso=iso;diffani=ani;bg=b;ng=n;tau=tau_input
    elMat%All=0;elMat%Alu=0;elMat%Alq=0;elMat%fh=0
    iel=1;Neqstab=Neq;Neqgrad=Neq;inn=0;ign=0
#ifdef NEUTRALGAMMA
    inn=5;ign=6
#endif
    ind_asf=0;ind_ash=0;ind_ff=[(i,i=1,Neq)];ind_fe=ind_ff;ind_fg=[(i,i=1,2*Neq)]
    NiNi=1.d0;Ni=1.d0
    ntang=tangent==0
    CALL compute_W2(ufg,W2,diffiso(1,1),diffiso(2,2))
    CALL compute_dW2_dU(ufg,dW2_dU,diffiso(1,1),diffiso(2,2));QdW2=MATMUL(Qpr,dW2_dU)
    CALL compute_W3(ufg,W3,diffiso(1,1),diffiso(2,2),diffiso(3,3))
    CALL compute_dW3_dU(ufg,dW3_dU,diffiso(1,1),diffiso(2,2),diffiso(3,3));QdW3=MATMUL(Qpr,dW3_dU)
    CALL compute_W4(ufg,W4,diffiso(1,1),diffiso(4,4))
    CALL compute_dW4_dU(ufg,dW4_dU,diffiso(1,1),diffiso(4,4));QdW4=MATMUL(Qpr,dW4_dU)
! PRODUCTION_DIFFUSION_FACTOR
! PRODUCTION_STABILIZATION
! PRODUCTION_DIFFUSION
! PRODUCTION_WEIGHT_DERIVATIVE
! PRODUCTION_DIAGNOSTIC_FACTOR
    diagnostic_target=u;diagnostic_target(2)=target
    CALL accumulate_momentum_closure(accumulator,uhat,diagnostic_target,q,n,b,bn,tau_input,iso,ani, &
         DOT_PRODUCT(MATMUL(q,W2),n-bn*b),plasma%momentum_diffusion_factor,1.d0)
    diagnostic=accumulator%values
    ju=elMat%Alu(:,:,1);jh=elMat%All(:,:,1);jq=elMat%Alq(:,:,1);rhs=elMat%fh(:,1)
  END SUBROUTINE
END MODULE
